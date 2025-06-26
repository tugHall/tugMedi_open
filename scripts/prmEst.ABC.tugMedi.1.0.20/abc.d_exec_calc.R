#!/usr/bin/Rscript --slave --vanilla
library(abc)
library(dplyr)

config  <- 'config.R'
args <- commandArgs(trailingOnly = T)
if (length(args) >= 1) {
    config <- args[1]
}
source(config)

if (file.exists(ABC_dir)) {
    dest <- sprintf('old/%s_calc_%s',
                    basename(ABC_dir),
                    format(Sys.time(), "%Y%m%d_%H%M%S"))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    file.copy(ABC_dir, dest, recursive=TRUE)
}

uniq_sumstat <- function(samstat) {
    sumstat <- sumstat %>% distinct()
    return(sumstat)
}

for (sample in samples) {
    dir.abc <- sprintf('%s/%s', ABC_dir, sample)
    ## param_mat
    param_mat <- read.table(sprintf('%s/param_mat.txt', dir.abc),
                            header=TRUE, sep='\t')

    for (tumor_content in tumor_contents) {
        ## target
        target <- read.table(sprintf('%s/target_tumor_content=%.02f.txt',
                                     dir.abc,
                                     tumor_content),
                             header=TRUE, sep='\t', stringsAsFactors=F)
        genes <- target[,1]
        target <- target[2]
        ## rownames(target) <- genes
        rownames(target) <- sprintf('%s.%s', genes, ave(genes, genes, FUN=seq_along))

        ## sumstat
        sumstat_file <- sprintf('%s/sumstat_tumor_content=%.02f.txt',
                                dir.abc,
                                tumor_content)
        sumstat <- read.table(sumstat_file, header=TRUE, sep='\t', colClasses='character')
        sumstat[, 2:ncol(sumstat)] <- lapply(sumstat[, 2:ncol(sumstat)], as.numeric)
        sumstat_params <- sumstat[,1]
        sumstat <- sumstat[, 2:ncol(sumstat)]
        rownames(sumstat) <- sumstat_params
        
        ##print(ncol(target))
        ##print(nrow(target))
        ##print(head(target))

        ##print(ncol(sumstat))
        ##print(nrow(sumstat))
        ##print(head(sumstat))

        ##print(ncol(sumstat))
        ##print(nrow(sumstat))
        ##print(head(param_mat))

        if (nrow(uniq_sumstat(sumstat)) <= 1) {
            print(sprintf('!!! %s: same summaries', sumstat_file))
            sel_params <- c(rownames(sumstat)[1])
        }
        else{
            ## ABC
            abc.result <- abc(target=target, param=param_mat, sumstat=sumstat, tol=abc.tol, method=abc.method)
            dist <- abc.result$dist
            write.table(abc.result$dist,
                        sprintf('%s/ABC_dist_tumor_content=%.02f.txt', dir.abc, tumor_content),
                        quote=F, col.names=F, ,sep='\t')

            sel_params <- abc.result$unadj.values
        }
        write.table(data.frame(param=sel_params),
                    sprintf('%s/TOLSel_%s_tumor_content=%.02f_tol=%.03f.txt', dir.abc, abc.method, tumor_content, abc.tol),
                    quote=F, row.names=F, sep='\t')
    }
}
