# correlation_outliers
Identifying outlier samples in situations with multiple data measurements per sample based on the Mahalanobis distance from the measurement covariance matrix per sample. Applied in gene expression outlier calling where each gene has measurements across many tissues.
The method is further described in this publication: Ferraro, N. M., Strober, B. J., Einson, J., Li, X., Aguet, F., Barbeira, A. N., ... & Park, Y. (2019). Diverse transcriptomic signatures across human tissues identify functional rare genetic variation. bioRxiv, 786053.

Data used for development is available as part of the Genotype Tissue Expression (GTEx) project v8: https://gtexportal.org/home/. 
We recommend filtering samples with high missingness prior to outlier calls to prevent sparsity issues.
This code was originally written by Emily Tsang (https://github.com/ektsang) and Joe Davis (https://github.com/joed3) and has been updated/modified by Nicole Ferraro for use in the above publication and release here.

## Required R packages
* data.table
* dplyr
* ggplot2
* argparse
* optparse
* doMC
* foreach
* stringr
* reshape2
* robustbase
* Matrix
* mvtnorm
* matrixStats
* mvnmle
* psych
* PMA
* softImpute
* impute

## Testing environment
* R 3.6.0 
* data.table (1.12.8)
* dplyr (0.8.5)
* ggplot2 (3.3.0)
* argparse (2.0.1)
* optparse (1.6.4)
* doMC (1.3.6)
* foreach (1.4.8)
* stringr (1.4.0)
* reshape2 (1.4.3)
* robustbase (0.93-5)
* Matrix (1.2-18)
* mvtnorm (1.1-0)
* matrixStats (0.56.0)
* mvnmle (0.1-11.1)
* psych (1.9.12.31)
* PMA (1.2.1)
* softImpute (1.4)
* impute (1.60.0)

## Set up directories
Data file should be a gzipped file in the format of Gene, Tissue/Data type, Samp1, Samp2, ... (example in toy_data folder)
Data directory leads to where the original data is stored and output directory where to write the output files

```
DATADIR=[insert]
DATAFILE=[insert]
OUTDIR=[insert]
```

## Select parameters
* Tests several potential imputation approaches (EM, MEAN, KNN, SOFT, PMD)
* Requires path to data file, output directory, number of samples (i.e. genes) used to select parameters, proportion of data to hold out for testing, number of cores for parallelization (default = 1), and a seed
```
Rscript pick_parameter_values_and_compare_methods.R \
        --DATA.FILE=$DATADIR/$DATAFILE
        --OUTPUT.DIR=$OUTDIR \
        --NSAMPLES=1000 \
        --PROP.HOLDOUT=0.05 \
        --NCORES=3 \
        --SEED=100
```
Generates an RData file containing the parameters tested and the values that minimize reconstruction error (impute.args.RData) and a set of plots showing error across a range of parameter choices in impute.parameter.plots.pdf

## Generate outlier statistics
* Requires path to data file, prefix for output, selected imputation method (valid = KNN, EM, MEAN, PMD, SOFT), minimum number of measurements needed per sample, number of cores (default = 1), and nominal outlier p-value threshold
* Either set parameters from above results, or however you prefer, and set here
* Optional parameter settings for chosen method (defaults chosen from GTEx use case):
	*	EM - max iterations (EM.MAX.IT, default = 2) and tolerance (EM.TOL, default = 1e-6)
	*	KNN - number of neighbors to consider (KNN.K, default = 200)
	*	SOFT - max iterations (SOFT.MAX.IT, default = 20), penalty (SOFT.LAMBDA, default = 15), rank (SOFT.RANK, default = 20)
	* 	PMD - max iterations (PMD.MAX.IT, default = 1), penalty (PMD.LAMBDA, default = 0.5), rank (PMD.RANK, default = 5)
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
Generates a file with outlier statistics, including Sample ID, gene, number of individuals tested, number of measurements/tissues for that sample, mahalanobis distance, log p-value, p-value, bonferroni-corrected p-value, and FDR corrected p-value. Also includes a column indicating if the sample is an outlier at the nominal threshold, but additional thresholds can be defined downstream. Example in toy_data folder for a few genes.

## Assign outliers to specific group
* Requires outlier file, output file, original data file (same as above), a threshold for determining a single measurement change (assumes absolute value) and an outlier FDR threshold
* Considers all groups as separate, may want to collapse related, i.e. different Brain regions -> Brain
```
Rscript filter_specific_outliers.R \
        --OUTLIER.FILE=$DATADIR/outliers_KNN.txt \
        --OUT.FILE=$DATADIR/specific_outliers_KNN.txt \
        --DATA.FILE=$DATADIR/$DATAFILE \
        --THRESH=3 \
        --FDR=0.001
```
Generates a file similar to the outlier file above, but with additional columns indicating if one measurement/tissue is driving the signal and the maximum value across all measurements. These are NA for controls or if no single measurement is driving the outlier call. Example in toy_data folder for a few genes.
