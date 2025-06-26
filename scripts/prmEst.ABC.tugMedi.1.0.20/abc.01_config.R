
## Your working directory with R/, Input/, and abc.run.sh. 
base_dir <- './'


## For parallele execution ======================================
replicates <- 5
work_dir <- "./work/first"
run_script <- './abc.01_run.R'


## For ABC ======================================================
## Output
ABC_dir <- './ABC1'

## Sim VAFs
VAF_files_glob <- 'work/first/*/Output/VAF/VAF.txt'

## Obs VAFs
rint_file    <- 'Input/Samples/samples.Rint.txt'
rother_file  <- 'Input/Samples/samples.Rother.txt'
rother_regex <- "^_.*"
samples      <- c('TCGA-AZ-6608-01A-11D-1835-10')
sample_cols  <- c(6, 1, 11) ## sample, gene, vaf

##
tumor_contents <- c(1.0, 0.8, 0.6, 0.4)
min_vaf        <- 0.1

##
abc.tol    <- 0.1
abc.method <- 'rejection'


