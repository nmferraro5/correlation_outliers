#!/usr/bin/env Rscript

# Load required packages
require(Matrix)
require(mvtnorm)
require(matrixStats)
require(mvnmle)
require(psych)

# EM algorithm for learning of multivariate normal parameters with missing data
# Assumes data has already been centered
# Does not include a regularization term
# Input is matrix with tissues as columns and rows as individuals
# Output is learned covariance and precision matrices along with the imputed data and the convergence measure
# If return_all is TRUE, then return the imputed data from every step, not merely the last one.
em.mvn.impute <- function(samp, impute_args, return_all = FALSE){
  n = nrow(samp)
  t = ncol(samp)
  tol = impute_args['EM.TOL']
  maxit = impute_args['EM.MAX.IT']
  cov_change = c()  

  # Impute into samp.new
  samp_new_all = list()
  samp_new = samp
  samp_new[is.na(samp)] = 0

  # Initialize precision matrix
  # Start with sample covariance matrix--need to ensure PSD
  sigma = .5 * var(samp, use = 'pair') + .5 * diag(t) # Initialization idea from Dahl et al. 2016 (PHENIX)
  S = solve(sigma)

  # Iterate E and M steps until convergence
  counter = 1
  while(counter <= maxit){
    S.old = S
    # E-step
    T_S = matrix(0, nrow = t, ncol = t)
    for(i in 1:n){
      miss_inds = is.na(samp[i, ])
      
      xobs = matrix(0, nrow = t, ncol = 1)
      xobs[!miss_inds] = samp[i, !miss_inds]
      xmis = matrix(0, nrow = t, ncol = 1)
      
      sigma11 = matrix(sigma[miss_inds, miss_inds], nrow = sum(miss_inds), ncol = sum(miss_inds))
      sigma12 = matrix(sigma[miss_inds, !miss_inds], nrow = nrow(sigma11), ncol = t - nrow(sigma11))
      sigma22inv = solve(sigma[!miss_inds, !miss_inds])
      x2 = matrix(samp[i, !miss_inds], ncol = 1)
      xmis[miss_inds] = sigma12 %*% sigma22inv %*% x2
      sigmaMO = matrix(0, ncol = t, nrow = t)
      sigmaMO[miss_inds, miss_inds] = sigma11 - sigma12 %*% sigma22inv %*% t(sigma12)
      
      T_S = T_S + xobs %*% t(xobs) + xmis %*% t(xobs) + xobs %*% t(xmis) + sigmaMO + xmis %*% t(xmis)
      samp_new[i, ] = t(xobs) + t(xmis)
      samp_new_all[[counter]] = samp_new
    }
    # M-step
    sigma = (1 / n) * T_S
    S = solve(sigma)
    # Calculate log-likelihood and check for convergence
    cov_change[counter] = norm(S - S.old, type = 'F')^2 / norm(S, type = 'F')^2 # Convergence idea from Dahl et al. 2016 (PHENIX)
    if(cov_change[counter] < tol){
      break
    }
    counter = counter + 1
  }
  if (return_all) {
    samp_return = samp_new_all
  } else {
    samp_return = samp_new
  }
  return(list(C = sigma, S = S, x = samp_return, cov_change = cov_change))
}
