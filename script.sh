################################
##### To execute this code, open a linux terminal, go to the working directory (the one containg the input files) and copy-paste the code lines on the
##### terminal. Output files will be placed in the working directory as well.
##### These scripts always take as arguments abundance table and metadata file and only samples represented in metadata will be considered.
##### All input tables are supposed to be column-tab separated with optionally one DESCRIPTION column. 

##### doTestDESeq2
# Implementation of DESeq2 library from Bioconductor.
# Test for differential expression by use of the negative binonial distribution and a shrinkage estimator for the distribution's variance.
# (parameter 1) abundance table file
# (parameter 2) metadata file
# (parameter 3) metadata column containing the levels of the groups to be compared
# (parameter 4) levels of the groups to be compared
# (parameter 5) tag used to build output files
# (parameter 6) if TRUE, low signal filtering is performed
# (parameter 7) name of feature DESCRIPTION column in abundance matrix
# (parameter 8) if TRUE, samples in groups are considered paired samples
# (parameter 9) if TRUE, Mean values in percentage scale will be reported for the two groups under comparisson in the results describing the significant features identified by DeSeq2. 
# (parameter 10) Type of algorithm ("parametric" or "local")

./doTestDESeq2.R counts.tsv metadata.tsv "patient_group" "R,NR" significant.features TRUE DESCRIPTION FALSE TRUE "parametric"

##### doBoruta
# Implementation of Boruta library from R
# A brief description of the algorithm:
# 1.- Create duplicate copies of all independent variables. When the number of independent variables in the original data is less than 5,
#	create at least 5 copies using existing variables.
# 2.- Shuffle the values of added duplicate copies to remove their correlations with the target variable. It is called shadow features or permuted copies.
# 3.- Combine the original ones with shuffled copies
# 4.- Run a random forest classifier on the combined dataset and performs a variable importance measure (the default is Mean Decrease Accuracy) to evaluate
#	the importance of each variable where higher means more important.
# 5.- Then Z score is computed. It means mean of accuracy loss divided by standard deviation of accuracy loss.
# 6.- Find the maximum Z score among shadow attributes (MZSA)
# 7.- Tag the variables as 'unimportant'  when they have importance significantly lower than MZSA. Then we permanently remove them from the process.
# 8.- Tag the variables as 'important'  when they have importance significantly higher than MZSA.
# 9.- Repeat the above steps for predefined number of iterations (random forest runs), or until all attributes are either tagged 'unimportant' or 'important', whichever comes first.
# This script takes as arguments abundance data (with description column), metadata and p-values table resulting from applying doTestDESeq2.
# (parameter 1) abundance table file
# (parameter 2) metadata file
# (parameter 3) DeSeq2 results
# (parameter 4) metadata column containing the levels of the groups to be compared
# (parameter 5) levels of the groups to be compared
# (parameter 6) tag used to build output files
# (parameter 7) metadata column containing the levels of the groups to be compared
# (parameter 8) if TRUE, Total Sum Scaling normalisation (TSS) will be performed at the begining
# (parameter 9) p-value cutoff. Only features with DESeq2 p-value lower than cutoff will be considered to perform boruta algorithm.

./doBoruta.R counts.tsv metadata.tsv DESeq2.results.tsv "patient_group" "NR,R" boruta.confirmed DESCRIPTION TRUE 0.05

##### doPrediction
# Validation samples status will be predicted by using a random forest model fitted on training samples. Then predicted values will be used
# to compute true positive and false positive rates for different discrimination thresholds so that ROC curve can be painted and AUC computed.
# Also, percentage of true positive at a given threshold are computed. 
# (parameter 1) abundance table file
# (parameter 2) metadata file for the training samples
# (parameter 3) metadata file for the validation samples
# (parameter 4) boruta results
# (parameter 5) metadata column containing the levels of the groups to be compared
# (parameter 6) levels of the groups to be compared
# (parameter 7) tag used to build output files
# (parameter 8) if TRUE, Total Sum Scaling normalisation (TSS) will be performed at the begining
# (parameter 9) Samples with predicted probability of belonging to one group or the other group under this cutoff will be removed from computations.
#			If the cutoff is set to 0.5 then it has no effect.
# (parameter 10) if TRUE, ROC curves computation will be done from probability tables previously filtered as determined by parameter 9

./doPrediction.R counts.tsv metadata.tsv metadata.validation.tsv results.boruta.tsv "patient_group" "NR,R" prediction.results TRUE 0.5 TRUE



