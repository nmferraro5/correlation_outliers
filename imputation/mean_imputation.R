#!/usr/bin/env Rscript

# input matrix with missing data
# impute missing values by column means
# return completed data and the sample precision matrix of the imputed data
mean.impute = function(data, impute_args=NULL) {
  data = apply(data, 2, function(cl) {cl[is.na(cl)] = mean(cl, na.rm = TRUE); cl})
  samp.cov = cov(data)
  samp.prec = solve(samp.cov)
  return (list(x = data, S = samp.prec, C = samp.cov))
}
