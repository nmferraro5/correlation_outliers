#!/usr/bin/env Rscript

require(PMA)
require(softImpute)

## PMD imputation
## if quiet set to false, let you know when it hits max.it before convergence
## returns: completed matrix (or matrices from all iteration if return_all is TRUE) and sample precision and covariance matrices
PMD.impute <- function(x, impute_args, thresh = 1e-6, quiet = TRUE, return_all = FALSE) {
    maxit = impute_args['PMD.MAX.IT']
    penalty = impute_args['PMD.LAMBDA']
    rank = impute_args['PMD.RANK']
    
    missing_ind = is.na(x)
    x = apply(x, 2, function(cl) {cl[is.na(cl)] = mean(cl, na.rm = T); cl})
    x_new_all = list()
    
    for (i in 1:maxit) {
        ## takes care of missing data imputation
        j_pmd = PMD(x, type = 'standard', sumabs = penalty,
                    sumabsu = NULL, sumabsv = NULL, K = rank, trace = F)
        ## replace missing values by fitted values
        fitted = j_pmd$u %*% diag(j_pmd$d[1:rank], ncol = rank, nrow = rank) %*% t(j_pmd$v) + j_pmd$meanx
        ## put fitted values
        x_new = x
        x_new[missing_ind] = fitted[missing_ind]
        ## check for convergence
        if(norm((x-x_new), type = "F") / norm(x, type = "F") < thresh) {
            break
        } else {
            x = x_new
            x_new_all[[i]] = x_new
        }
    }
    if (i == maxit && !quiet) {
        cat('did not converge \n')
    }
    samp_cov = cov(x_new)
    samp_prec = solve(samp_cov)
    if (return_all) {
        samp_return = x_new_all
    } else {
        samp_return = x_new
    }
    return(list(x = samp_return, S = samp_prec, C = samp_cov))
}

## soft imputation (soft thresholding of svd)
## default is no penalty (lambda=0)
## returns completed matrix and sample precision and covariance matrices
SVD.soft.impute = function(x, impute_args) {
    maxit = impute_args['SOFT.MAX.IT']
    lambda = impute_args['SOFT.LAMBDA'] 
    rank = impute_args['SOFT.RANK']
    fit = softImpute(x, rank.max = rank, trace = FALSE, type = 'svd', lambda = lambda, maxit = maxit)
    x_new = complete(x, fit)
    samp_cov = cov(x_new)
    samp_prec = solve(samp_cov)
    return(list(x = x_new, S = samp_prec, C = samp_cov))
}
