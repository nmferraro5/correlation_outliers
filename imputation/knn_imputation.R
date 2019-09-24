#!/usr/bin/env Rscript

# Load required packages
#require(Matrix)
#require(matrixStats)
require(impute)

# Function to impute missing data using KNN and then estimate precision matrix 
# Input: data matrix (rows = observations), k (number nearest neighbors), 
# row max and col max missingness rates
# Output: estimated precision matrix after imputation
knn.impute <- function(data, impute_args, rmax = 0.999, cmax = 0.999){
  k = impute_args['KNN.K']
	new_data = impute.knn(data, k = k, rowmax = rmax, colmax = cmax)
	covariance = cov(new_data$data)
	precision = solve(covariance)
	return(list(x = new_data$data, C = covariance, S = precision))
}
