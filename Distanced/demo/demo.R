#Install and load packages
install.packages("devtools")
devtools::install_github(repo="thackmann/distanced", subdir="Distanced")
library(Distanced)

#Set parameters
sample.filepath=file.path=system.file("extdata", "Mock1V34run130401_test_data.fq", package = "Distanced")
reference.filepath=file.path=system.file("extdata", "reference_sequences_V34_no_primers.fasta", package = "Distanced")
n_sample_max.scalar=100
replace_ambiguous_letters.character=FALSE
random_seed.character=FALSE

#Run Distanced
Distanced(sample.filepath, reference.filepath, n_sample_max.scalar, replace_ambiguous_letters.character, random_seed.character)

#If Distanced ran successfully, output will match that below
#                        Mean pairwise distance
# Uncorrected                         0.2351698
# Estimated by Distanced              0.2236499
# Actual                              0.2249940
