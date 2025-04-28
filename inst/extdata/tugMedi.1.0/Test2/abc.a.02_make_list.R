#!/usr/bin/Rscript --slave --vanilla
## options(error = dump.frames)

config  <- 'config.R'
args <- commandArgs(trailingOnly = T)
if (length(args) >= 1) {
    config <- args[1]
}
source(config)

## this function might need to be modified
make_init_data <- function(base_dir, dest_dir, selected_file) {
    df_file <- read.table( sprintf('%s/%s', base_dir, selected_file), header=TRUE )
    input_dir <- ''
    if (length(df_file$param) == 1) {
        input_dir <- sprintf( paste0('%s/', generator_dir), base_dir, df_file$param[1] )
    }
    else {
        input_dir <- sprintf( paste0('%s/', generator_dir), base_dir, sample(df_file$param, 1) )
    }

    system(sprintf('cp -a %s/Input %s', input_dir, dest_dir))
    system(sprintf('cp -a %s/Output %s', input_dir, dest_dir))
}

run_list <- './run_list.txt'
base_dir <- normalizePath(base_dir)
run_script <- normalizePath(run_script)

if (file.exists(run_list)) {
    file.remove(run_list)
}

if (file.exists(work_dir)) {
    work_dir <- normalizePath(work_dir)
    old_dir <- sprintf('%s/old', dirname(work_dir))
    dir.create(old_dir, recursive = TRUE, , showWarnings = FALSE)
    file.rename(work_dir,
                sprintf('%s/%s_%s',
                        old_dir,
                        basename(work_dir),
                        format(Sys.time(), "%Y%m%d_%H%M%S")))
}

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
work_dir <- normalizePath(work_dir)

for (i in 1:replicates) {
    dest_dir <- sprintf('%s/%06d', work_dir, i)
    run_shell <- sprintf('%s/abc.run.sh', dest_dir)

    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    make_init_data(base_dir, dest_dir, selected_file)
    system(sprintf('cp -a %s %s/run.R', run_script, dest_dir))
    system(sprintf('cp -a %s/abc.run.sh %s/', base_dir, dest_dir))
    system(sprintf('ln -s %s/R %s', base_dir, dest_dir))

    cat(sprintf('%06d\n', i), file = sprintf('%s/run_number.txt', dest_dir))
    cat(run_shell, '\n', file = run_list, append=T)
}
