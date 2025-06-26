
## Your working directory with R/, Input/, and abc.run.sh. 
base_dir <- './'


## For parallele execution ======================================
#replicates <- 2e3
replicates <- 10
work_dir <- "./work/TCGA-55-7903-01A-11D-2167-08"
run_script <- './test2.2_postABC.R'


## For ABC ======================================================
## Output
ABC_dir <- './ABC2'

## Sim VAFs
VAF_files_glob <- 'work/TCGA-55-7903-01A-11D-2167-08/*/Output/VAF/VAF.txt'

## Obs VAFs
rint_file    <- 'Input/Samples/samples.Rint.txt'
rother_file  <- 'Input/Samples/samples.Rother.txt'
rother_regex <- "^_.*"
samples      <- c('TCGA-55-7903-01A-11D-2167-08')
sample_cols  <- c(6, 1, 11) ## sample, gene, vaf

##
tumor_contents <- c(0.8)
min_vaf        <- 0.1

##
#abc.tol    <- 10e-3
abc.tol    <- 0.2
abc.method <- 'rejection'


## Selection based on ABC
selected_file <- 'ABC1/TCGA-55-7903-01A-11D-2167-08/TOLSel_rejection_tumor_content=0.80_tol=0.200.txt'
generator_dir <- 'work/first/%06d'


