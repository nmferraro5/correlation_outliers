#!/usr/bin/env Rscript
## Author: Joe Davis, 2017
## Modified: Nicole Ferraro, 2019

require(data.table)
require(ggplot2)
require(doMC)
require(foreach)
require(tidyr)
require(dplyr)
require(stringr)
require(optparse)

### Helper functions

# Function to make sure matrix with held out data isn't too sparse
# don't let any row or column have less than 10% observed data
# returns true is the matrix fails the sparsity conditions
too.sparse <- function(mat, sparsity = 0.10) {
    row_observed = apply(mat,1, function(x) sum(!is.na(x)))
    col_observed = apply(mat,2, function(x) sum(!is.na(x)))
    return (min(row_observed) < sparsity*ncol(mat) || min(col_observed) < sparsity*nrow(mat))
}

# Function that picks the data to hold out (tries 100 iterations to meet sparsity constraints)
# Args:
# * data: matrix, possibly with missing values
# * prop_holdout: the proportion of observed values to hold out
# Returns:
#  a list with the following elements
#  indices: the indices that were heldout
#  heldout: the values from those indices
#  mat: the resultsing matrix with those values removed
holdout.data <- function(data, prop_holdout = 0.05) {
    entries = which(!is.na(data))
    ok = FALSE
    n = 1
    while (!ok && n <= 100) {
        indices = sample(entries, ceiling(prop_holdout*length(entries)))
        heldout = data[indices]
        extra_missing = data
        extra_missing[indices] = NA
        if (!too.sparse(extra_missing)){
            ok = TRUE
        } else {
            n = n + 1
        }
    }
    if (n > 100) {
        cat('WARNING: sparsity conditions are probably too strict. Proceeding with a matrix that violates sparsity conditions.\n')
    }
    return(list(indices = indices, heldout = heldout, mat = extra_missing))
}

# Args:
# * holdout: list returned by holdout.data
# * FUN: imputation function
# * ...: optional arguments to the imputation function
# Returns:
#  the average squared error of the imputed values
reconstruction.error <- function(holdout, FUN, ...) {
    FUN = match.fun(FUN)
    result = FUN(holdout$mat, ...) # run the imputation
    complete = result$x
    error = error.function(complete, holdout)
    return(error)
}

# Function to calculate average squared error on heldout data
# Args:
# * imputed_data: completed data matrix
# * holdout: list returned by holdout.data
# Returns:
#   The average squared error, across heldout instances, of the imputed values
error.function <- function(imputed_data, holdout) {
    error = sum((imputed_data[holdout$indices] - holdout$heldout)^2) / length(holdout$indices)
    return(error)
}

# Function to select parameter values that minimize recon. error
# Args:
# * stat:
# * prefix:
# * drop_first:
# Returns:
# list of minimum parameter values
pick.params = function(stat, prefix, drop_first = FALSE) {
    stat_subset = mean_errors[str_detect(names(mean_errors), prefix)]
    if (drop_first) {
        stat_subset = stat_subset[-1]
    }
    stat_subset = sort(stat_subset)
    mins = which(stat_subset == min(stat_subset))
    return(names(mins))
}

### Argument parsing
option_list = list(make_option(c('--DATA.FILE'), type = 'character', default = NULL, help = 'path to the Z-score data'),
                   make_option(c('--OUTPUT.DIR'), type='character', default = NULL, help = 'path to output directory'),
                   make_option(c('--NSAMPLES'), type = 'numeric', default = 1000, help = 'number of samples (i.e. genes) to use for parameter evaluation'), 
                   make_option(c('--PROP.HOLDOUT'), type = 'numeric', default = 0.05, help = 'proportion data to holdout for evaluation'),
                   make_option(c('--NCORES'), type = 'numeric', default = 1, help = 'number of cores for parallelization'),
                   make_option(c('--SEED'), type = 'numeric', default = 1, help = 'set seed for reproducability'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Register backend for parallelization and set seed
registerDoMC(cores = opt$NCORES)
set.seed(opt$SEED)

## Source functions for various imputation methods
source('imputation/mvn_model_em.R')
source('imputation/mean_imputation.R')
source('imputation/knn_imputation.R')
source('imputation/svd_imputation.R')

## Load data
data = fread(opt$DATA.FILE)
setkey(data, Gene)

selected_genes = sample(unique(data[, Gene]), opt$NSAMPLES)

## Run imputation
ranks = c(5,10,15,20)
soft_its = c(5,10,15,20,30)
soft_pens = c(0,5,10,15,20,25)
pmd_pens = c(0.5,0.7,0.9,1)
pmd_maxit = 30

soft_params = expand.grid(soft_its, ranks, soft_pens)
pmd_params = expand.grid(ranks, pmd_pens)

em_tol = 1e-6
em_maxit = 100
knn_ks = c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,90,100,200,300)

reconstruction = foreach(i = 1:opt$NSAMPLES, .combine = rbind) %dopar% { 
    subset = as.matrix(t(data[Gene == selected_genes[i], 3:ncol(data)])) 
    holdout = holdout.data(subset, opt$PROP.HOLDOUT) 

    ## EM
    em_args = c(em_tol, em_maxit)
    names(em_args) = c('EM.TOL', 'EM.MAX.IT')
    em_results = em.mvn.impute(holdout$mat, em_args, return_all = TRUE)
    em_errs = unlist(lapply(em_results$x, error.function, holdout))
    em_errs_padded = rep(NA, em_maxit) # pad with NA in case of early convergence
    em_errs_padded[1:length(em_errs)] = em_errs

    ## MEAN
    mean_err = reconstruction.error(holdout, mean.impute)

    ## KNN
    missing = is.na(holdout$mat)
    M = dim(missing)[1]
    N = dim(missing)[2]
    colmax = max(colSums(missing)) / M + .01
    rowmax = max(rowSums(missing)) / N + .01
    
    knn_errs = rep(NA, length(knn_ks))
    for (j in 1:length(knn_ks)) {
        knn_args = knn_ks[j]
        names(knn_args) = 'KNN.K'
        knn_errs[j] = reconstruction.error(holdout, knn.impute, knn_args, rmax = rowmax, cmax = colmax)
    }

    ## SOFTIMPUTE
    soft_errs = rep(NA, nrow(soft_params))
    for (k in 1:nrow(soft_params)) {
        pmd_args = c(as.numeric(soft_params[k,]))
        names(pmd_args) = c('SOFT.MAX.IT', 'SOFT.RANK', 'SOFT.LAMBDA')
        soft_errs[k] = reconstruction.error(holdout, SVD.soft.impute, pmd_args)
    }
    
    ## PMDIMPUTE
    pmd_errs_padded = rep(NA, nrow(pmd_params) * pmd_maxit)
    for (l in 1:nrow(pmd_params)) {
        pmd_args = c(as.numeric(pmd_params[l,]), pmd_maxit)
        names(pmd_args) = c('PMD.RANK', 'PMD.LAMBDA', 'PMD.MAX.IT')
        pmd_results = PMD.impute(holdout$mat, pmd_args, return_all = TRUE)
        pmd_errs = unlist(lapply(pmd_results$x, error.function, holdout))
        pmd_errs_padded[((l-1)*pmd_maxit + 1):((l-1)*pmd_maxit + length(pmd_errs))] = pmd_errs
    }
    
    ## Write out results
    c(em_errs_padded, mean_err, knn_errs, soft_errs, pmd_errs_padded)
}


## Gather column names and add them
cnames = c(paste('EM', 1:em_maxit, sep = '_'),
           'Mean',
           paste('Knn', knn_ks, sep = '_'),
           paste('Soft', apply(soft_params, 1, paste, collapse = '_'), sep = '_'),
           paste('PMD', apply(expand.grid(1:pmd_maxit, ranks, pmd_pens),
                              1, paste, collapse = '_'), sep = '_'))

colnames(reconstruction) = cnames

## Make data frame long and split into different data frames by method for plotting
reconstruction_long = gather(as.data.frame(reconstruction) %>% mutate(Gene = rownames(reconstruction)),
                             Method, Error, -Gene)

recon_em_mean_knn = reconstruction_long %>%
    filter(str_detect(Method, 'EM|Mean|Knn')) %>%
    separate(Method, c('Method', 'Parameter'), '_', convert = TRUE, fill = 'right') %>%
    mutate(Parameter = ifelse(is.na(Parameter), 0, Parameter))

recon_soft = reconstruction_long %>%
    filter(str_detect(Method, 'Soft')) %>%
    separate(Method, c('Method', 'Num.it', 'Rank', 'Penalty'), '_', convert = TRUE) %>%
    mutate(Rank = factor(Rank), Penalty = factor(Penalty), Num.it = factor(Num.it))

recon_pmd = reconstruction_long %>%
    filter(str_detect(Method, 'PMD')) %>%
    separate(Method, c('Method', 'Num.it', 'Rank', 'Penalty'), '_', convert = TRUE) %>%
    mutate(Rank = factor(Rank), Penalty = factor(Penalty))

## Actual plotting
mytheme = theme_bw() + theme(legend.position = 'bottom',
                             strip.text = element_text(face = 'bold'),
                             strip.background = element_blank())

pdf(paste0(opt$OUTPUT.DIR, 'impute.parameter.plots.pdf'), width = 13, height = 5.5)

plot_em_mean_knn = ggplot(recon_em_mean_knn, aes(x = factor(Parameter), y = Error)) +
    geom_boxplot(aes(fill = Method)) + facet_grid(.~Method, scales = 'free') +
    xlab('Parameter') +
    scale_y_continuous(limits = c(0, 2)) + mytheme
print(plot_em_mean_knn)

plot_soft = ggplot(recon_soft, aes(x = Num.it, y = Error)) +
    geom_boxplot(aes(fill = Rank, alpha = Penalty)) +
    guides(alpha = guide_legend(override.aes = list(fill = 'darkgrey'))) +
    xlab('Number of iterations') +
    mytheme + ggtitle('SoftImpute')
print(plot_soft)

plot_soft2 = ggplot(recon_soft, aes(x = Penalty, y = Error)) +
    geom_boxplot(aes(fill = Rank, alpha = Num.it)) +
    guides(alpha = guide_legend(override.aes = list(fill = 'darkgrey'))) +
    mytheme + ggtitle('SoftImpute')
print(plot_soft2)

plot_pmd = ggplot(recon_pmd %>% filter(Num.it <= 10), aes(x = factor(Num.it), y = Error)) +
    geom_boxplot(aes(fill = Rank)) +
    xlab('Number of iterations') +
    facet_wrap(~Penalty, ncol = 2) + mytheme + ggtitle('PMD imputation')
print(plot_pmd)

plot_pmd2 = ggplot(recon_pmd %>% filter(Num.it <= 10), aes(x = Rank, y = Error)) +
    geom_boxplot(aes(fill = Rank, alpha = factor(Num.it))) +
    guides(alpha = guide_legend(override.aes = list(fill = 'darkgrey'))) +
    facet_wrap(~Penalty, ncol = 2) + mytheme + ggtitle('PMD imputation')
print(plot_pmd2)

dev.off()

mean_errors = colMeans(reconstruction)

em_min = pick.params(mean_errors, 'EM', drop_first = TRUE)
knn_min = pick.params(mean_errors, 'Knn')
soft_min = pick.params(mean_errors, 'Soft')
pmd_min = pick.params(mean_errors, 'PMD')

cat('Number of iterations for EM:\n')
print(str_split_fixed(em_min, '_', 2)[,2])
cat('Number of neighbors for knn:\n')
print(str_split_fixed(knn_min, '_', 2)[,2])
cat('Number of iterations, rank, penalty for SoftImpute:\n')
print(str_split_fixed(soft_min, '_', 4)[,2:4])
cat('Number of iterations, rank, penalty for PMD:\n')
print(str_split_fixed(pmd_min, '_', 4)[,2:4])

rm(data) # drop data for smaller file output

## Save workspace image
save.image(paste0(opt$OUTPUT.DIR, 'impute.args.RData'))
