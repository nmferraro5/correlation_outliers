## correlation_outliers
Identifying outlier samples in situations with multiple data measurements per sample based on the Mahalanobis distance from the measurement covariance matrix per sample. Applied in gene expression outlier calling where each gene has measurements across many tissues.
The method is further described in this publication: coming soon
Data used for development is available as part of the Genotype Tissue Expression (GTEx) project v8: https://gtexportal.org/home/
We recommend filtering samples with high missingness prior to outlier calls to prevent sparsity issues.
This code was originally written by Emily Tsang and Joe Davis and has been updated/modified by Nicole Ferraro for use in the above publication and release here.

## required R packages
data.table
dplyr
ggplot2
argparse
optparse
doMC
foreach
stringr
reshape2
robustbase

## set up directories
```
DATADIR=[insert]
DATAFILE=[insert]
OUTDIR=[insert]
```

## select parameters across potential imputation approaches (EM, MEAN, KNN, SOFT, PMD)
# Requires path to data file, output directory, number of samples (i.e. genes) used to select parameters, proportion of data to hold out for testing, number of cores for parallelization (default = 1), and a seed
```
Rscript pick_parameter_values_and_compare_methods.R \
        --DATA.FILE=$DATADIR/$DATAFILE
        --OUTPUT.DIR=$OUTDIR \
        --NSAMPLES=1000 \
        --PROP.HOLDOUT=0.05 \
        --NCORES=3 \
        --SEED=100
```
Generates an RData file containing the parameters tested and the values that minimize reconstruction error (impute.args.RData)
Generates a set of plots showing error across a range of parameter choices in impute.parameter.plots.pdf

## generate outlier statistics with given imputation method
# Requires path to data file, prefix for output, selected imputation method (valid = KNN, EM, MEAN, PMD, SOFT), minimum number of measurements needed per sample, number of cores (default = 1), and nominal outlier p-value threshold
# Data file should be gzipped file in the format Gene, Tissue/Data type, Samp1, Samp2, ...
# Either set parameters from above results, or however you prefer, and set here
# Optionally, include parameter settings for chosen method (defaults chosen from GTEx use case):
#	EM - max iterations (EM.MAX.IT, default = 2) and tolerance (EM.TOL, default = 1e-6)
#	KNN - number of neighbors to consider (KNN.K, default = 200)
#	SOFT - max iterations (SOFT.MAX.IT, default = 20), penalty (SOFT.LAMBDA, default = 15), rank (SOFT.RANK, default = 20)
# 	PMD - max iterations (PMD.MAX.IT, default = 1), penalty (PMD.LAMBDA, default = 0.5), rank (PMD.RANK, default = 5)
```
Rscript call_outliers.R \
        --DATA.FILE=$DATADIR/$DATAFILE \
        --OUTPUT.DIR=$OUTDIR \
        --OUT.FILE=outliers_ \
        --IMPUTE.METHOD=KNN \
        --NTYPES=5 \
        --NCORES=5 \
        --PTHRESH=0.0027 \
	--KNN.K=200
```
Generates a file with outlier statistics, including Sample ID, Gene, Number of individuals tested, Number of measurements/tissues for that sample, Mahalanobis distance, Log P-value, Pvalue, Bonferroni-corrected p-value, and FDR corrected p-value. Also includes a column indicating if the sample is an outlier at the nominal threshold, but additional thresholds can be defined downstream.

## assign outliers to specific group (considers all groups as separate, may want to collapse related, i.e. different Brain regions -> Brain)
# Requiires outlier file, output file, original data file (same as above), a threshold for determining a single measurement change and an outlier FDR threshold
```
Rscript filter_specific_outliers.R \
        --OUTLIER.FILE=$DATADIR/outliers_KNN.txt \
        --OUT.FILE=$DATADIR/specific_outliers_KNN.txt \
        --DATA.FILE=$DATADIR/$DATAFILE \
        --THRESH=3 \
        --FDR=0.001
```
Generates a file similar to the outlier file above, but with additional columns indicating if one measurement/tissue is driving the signal and the maximum value across all measurements. These are NA for controls or if no single measurement is driving the outlier call.
