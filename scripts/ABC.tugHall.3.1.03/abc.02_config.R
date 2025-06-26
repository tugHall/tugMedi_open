
## Your working directory with R/, Input/, and abc.run.sh. 
base_dir <- './'


## For parallele execution ======================================
replicates <- 5
work_dir <- "./work/TCGA-AZ-6608-01A-11D-1835-10"
run_script <- './abc.02_run.R'


## For ABC ======================================================
## Output
ABC_dir <- './ABC2'

## Sim VAFs
VAF_files_glob <- 'work/TCGA-AZ-6608-01A-11D-1835-10/*/Output/VAF/VAF.txt'

## Obs VAFs
rint_file    <- 'Input/Samples/samples.Rint.txt'
rother_file  <- 'Input/Samples/samples.Rother.txt'
rother_regex <- "^_.*"
samples      <- c('TCGA-AZ-6608-01A-11D-1835-10')
sample_cols  <- c(6, 1, 11) ## sample, gene, vaf

##
tumor_contents <- c(1.00, 0.8, 0.6, 0.4)
min_vaf        <- 0.1

##
abc.tol    <- 0.1
abc.method <- 'rejection'


## Selection based on ABC
selected_file <- 'ABC1/TCGA-AZ-6608-01A-11D-1835-10/TOLSel_rejection_tumor_content=0.40_tol=0.100.txt'
generator_dir <- 'work/first/%06d'


