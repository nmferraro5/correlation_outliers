#!/usr/bin/env Rscript

## Script to call outliers usign different methods given a file with Z-scores.
## The two first columns have the gene and tissue or phenotype, and the gene column is named "Gene".
## Subsequent columns have Z-score data for each individual and the individual ID is the column name.
## Author: Emily K. Tsang and Joe Davis, 2017
## Modified: Nicole Ferraro, 2019

rm(list = ls())


## Source the imputation methods
source('imputation/mvn_model_em.R')
source('imputation/knn_imputation.R')
source('imputation/mean_imputation.R')
source('imputation/svd_imputation.R')

## Load required packages
require(data.table)
require(optparse)
require(reshape2)
require(dplyr)
require(foreach)
require(doMC)
require(robustbase)

#--------------- FUNCTIONS

#### Function to calculate Mahalanobis distance for observations with missing values
## Input: dataset with missing observations (centered and scaled) and estimated covariance matrix
## Output: Mahalanobis distances, degrees of freedom and p-values
## Test statistic is assumed to be Chi-squared distributed under the null
mahala <- function(data, data_cov){
    ds = c()
    dfs = c()
    logps = c()
    for(i in 1:nrow(data)){
        pres_inds = !is.na(data[i, ])
        x = as.matrix(data[i, pres_inds], ncol = 1)
        sigma = data_cov[pres_inds, pres_inds]
        ds[i] = t(x) %*% solve(sigma) %*% x
        dfs[i] = sum(pres_inds)
        logps[i] = pchisq(ds[i], dfs[i], log.p = T, lower.tail = F)
    }
    return(data.frame(Ind = rownames(data), Dist = ds, Df = dfs, Pval = exp(logps), LogP = logps))
}

#### Function to pick top outliers per gene
## Input: Data frame of outlier test statistics and the number of data types required per test
## Output: Assignments of outliers and controls for genes with outliers
## Performs BF correction on the gene levels and Bh across genes
pick.outliers <- function(test_stats, ntypes, pthresh){
    test_stats = test_stats %>% filter(Df >= ntypes) %>%
        group_by(Gene) %>%
        mutate(N = n())
    test_stats = test_stats %>% mutate(BF = Pval * N) %>%
      mutate(BF = ifelse(BF > 1, 1, BF))
    outliers = test_stats %>% filter(Pval < pthresh) %>%
      ungroup() %>%
      mutate(FDR = p.adjust(BF, method = 'BH'), Y = 'outlier')
    controls = test_stats %>% filter(Pval >= pthresh) %>%
      ungroup() %>%
      mutate(FDR = NA, Y = 'control')
    out = rbind(outliers, controls) %>% select(Ind, Gene, N, Df, Dist, LogP, Pval, BF, FDR, Y)
    return(out)
}

#### Function to compute precision matrices for each gene and call outliers
## Input: Normalized data iwith the final N columns representing the N samples
## Output: Assignments of outliers and controls for genes with outliers
call.outliers <- function(data, ntypes, pthresh, impute_method, impute_args){
  genes = unique(sort(data$Gene))
  test_stats = foreach(i = 1:length(genes), .combine = rbind) %dopar% {
    gene = genes[i]
    gene_data = data %>% filter(Gene == gene)## %>% arrange(Phenotype)
    input_data = t(gene_data[, 3:ncol(gene_data)])
    imputed = impute_method(input_data, impute_args)
    test = mahala(input_data, imputed$C) %>% mutate(Gene = gene)
    test
  }
  outliers = pick.outliers(test_stats, ntypes, pthresh)
  return(outliers)
}

#--------------- MAIN

##-- Read command line arguments
option_list = list(make_option(c('--DATA.FILE'), type = 'character', default = NULL, help = 'path to the Z-score data'), 
  make_option(c('--OUTPUT.DIR'), type = 'character', default = NULL, help = 'output directory'),
	make_option(c('--OUT.PREFIX'), type = 'character', default = 'outliers', help = 'unique prefix for the outlier files'),
	make_option(c('--NCORES'), type = 'numeric', default = 1, help = 'number of cores for parallelization'),
	make_option(c('--NTYPES'), type = 'numeric', default = 5, help = 'number of observed measurements required to test for outlier'),
	make_option(c('--PTHRESH'), type = 'numeric', default = 0.01, help = 'nominal outlier p-value threshold'),
	make_option(c('--IMPUTE.METHOD'), type = 'character', default = 'KNN', help = 'method for data imputation'),
	make_option(c('--EM.MAX.IT'), type = 'numeric', default = 2, help = 'max iterations for EM'), 
	make_option(c('--EM.TOL'), type = 'numeric', default = 1e-6, help = 'convergence tolerance for EM'),
	make_option(c('--KNN.K'), type = 'numeric', default = 200, help = 'number of neighbors to use in KNN'), 
	make_option(c('--SOFT.MAX.IT'), type = 'numeric', default = 20, help = 'max iterations for SOFT'), 
	make_option(c('--SOFT.LAMBDA'), type = 'numeric', default = 15, help = 'penalty for SOFT'), 
	make_option(c('--SOFT.RANK'), type = 'numeric', default = 20, help = 'rank for SOFT'), 
	make_option(c('--PMD.MAX.IT'), type = 'numeric', default = 1, help = 'max iterations for PMD'), 
	make_option(c('--PMD.RANK'), type = 'numeric', default = 5, help = 'rank for PMD'), 
	make_option(c('--PMD.LAMBDA'), type = 'numeric', default = 0.5, help = 'penalty for PMD'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## Register parallel backend
registerDoMC(cores = opt$NCORES)

ntypes = opt$NTYPES
pthresh = opt$PTHRESH
prefix = opt$OUT.PREFIX
impute_method = toupper(opt$IMPUTE.METHOD)

##-- Analysis

## Read in the normalized data
data = as.data.frame(fread(paste0('zcat ', opt$DATA.FILE)))
impute_args = c(opt$EM.MAX.IT, opt$EM.TOL, opt$KNN.K, opt$SOFT.MAX.IT, opt$SOFT.LAMBDA, opt$SOFT.RANK, opt$PMD.MAX.IT, opt$PMD.RANK, opt$PMD.LAMBDA)
names(impute_args) = c('EM.MAX.IT', 'EM.TOL', 'KNN.K', 'SOFT.MAX.IT', 'SOFT.LAMBDA', 'SOFT.RANK', 'PMD.MAX.IT', 'PMD.RANK', 'PMD.LAMBDA')

##-- Call outliers using stated method
if (impute_method == 'KNN') {
  outliers = call.outliers(data, ntypes, pthresh, knn.impute, impute_args)
} else if (impute_method == 'EM') {
  outliers = call.outliers(data, ntypes, pthresh, em.mvn.impute, impute_args)
} else if (impute_method == 'MEAN') {
  outliers = call.outliers(data, ntypes, pthresh, mean.impute, impute_args)
} else if (impute_method == 'PMD') {
  outliers = call.outliers(data, ntypes, pthresh, PMD.impute, impute_args)
} else if (impute_method == 'SOFT') {
  outliers = call.outliers(data, ntypes, pthresh, SVD.soft.impute, impute_args)
} else {
  print('Invalid imputation method')
  print('Please select from: KNN, EM, MEAN, PMD, SOFT')
  quit(status = 2)
}

write.table(outliers, file=paste0(opt$OUTPUT.DIR, prefix, impute_method, '.txt'), sep = '\t', row.names = F, quote = F)



