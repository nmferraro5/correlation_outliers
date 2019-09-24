#!/usr/bin/env Rscript


# impute missing values with draws from a N(0,1)
random.impute = function(mat) {
    imputed = rnorm(sum(is.na(mat)))
    mat[is.na(mat)] = imputed
    covariance = cov(mat)
    precision = solve(covariance)
    return(list(x = mat,C = covariance,S = precision))
}
