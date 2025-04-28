#!/usr/bin/Rscript --slave --vanilla
## options(error = dump.frames)

config  <- 'config.R'
args <- commandArgs(trailingOnly = T)
if (length(args) >= 1) {
    config <- args[1]
}
source(config)

get_param_id <- function(file) {
    work_dir <- sub('/Output/VAF/VAF.txt', '', file)
    df_info <- read.table(sprintf('%s/run_number.txt', work_dir), header=FALSE, colClasses="character")
    return(df_info$V1[1])
}

get_sample <- function(file, sample) {
    df <- read.table(file, header=T, sep='\t', stringsAsFactors=FALSE, quote='')
    df <- df[, sample_cols]
    colnames(df) <- c('sample', 'gene', 'vaf')
    df <- df[(df$vaf >= min_vaf)
            & (df$sample == sample),]
    return(df)
}

get_extend_obs <- function(genes, df_obs, df_sim) {
    for (gene in genes) {
        obs_len <- sum(df_obs$gene == gene)
        df <- df_sim[df_sim$gene == gene, ]
        max_len <- 0
        if (nrow(df) > 0) {
            sim_lens <- as.vector(tapply(df$VAF_primary,
                                         df$param,
                                         length))
            max_len <- (max(max(sim_lens), obs_len))
        }
        if (obs_len < max_len) {
            diff <- max_len - obs_len
            df_obs <- rbind(df_obs,
                            data.frame(gene=rep(gene, diff), vaf=rep(0, diff)))
        }
    }
    df_obs_extend <- df_obs[order(df_obs$gene,
                                  df_obs$vaf,
                                  method='radix',
                                  decreasing=c(FALSE, TRUE)), ]
    return(df_obs_extend)
}

get_sumstat <- function(genes, target, df_sim, params) {
    sumstat <- matrix(0, nrow=length(params), ncol=length(target))
    dimnames(sumstat) <- list(params, names(target))
    for (param in params) {
        for (gene in genes) {
            df <- df_sim[(df_sim$param == param) & (df_sim$gene == gene), ]
            if (nrow(df) > 0) {
                sumstat_pos <- grep(gene, colnames(sumstat))
                sumstat[param, sumstat_pos] <- c(sort(df$VAF_primary, decreasing=TRUE),
                                                 rep(0, length(sumstat_pos) - length(df$VAF_primary)))
            }
        }
    }
    return(sumstat)
}

rename_gene_with_file <- function(file, sample, genes) {
    df <- read.table(file, header=T, sep='\t', stringsAsFactors=FALSE, quote='')
    df <- df[(df$VAF >= min_vaf)
            & (df$Sample == sample)
          , sample_cols]
    colnames(df) <- c('gene', 'vaf')
    genes <- ifelse(genes %in% df$gene, genes, '_Rother')
    return(genes)
}

rename_gene <- function(rother_regex, sample, genes) {
    genes <- gsub(rother_regex, '_Rother',  genes)
    return(genes)
}

## check target
for (sample in samples) {
    df_obs <- rbind(get_sample(rint_file, sample),
                    get_sample(rother_file, sample))
    print(head(df_obs, n=5))
}

if (file.exists(ABC_dir)) {
    dir.create('old/', recursive = TRUE, showWarnings = FALSE)
    file.rename(ABC_dir,
                sprintf('old/%s_make_%s',
                        basename(ABC_dir),
                        format(Sys.time(), "%Y%m%d_%H%M%S")))
}

for (tumor_content in tumor_contents) {
    ## read VAF
    VAF_files <- Sys.glob(VAF_files_glob)

    for (sample in samples) {
        cat(sprintf('%s tumor_content=%f\n', sample, tumor_content))
        df_sim <- do.call(rbind, lapply(VAF_files,
                                        function(file) {
                                            df <- read.table(file, header=T)
                                            time_max <- max(df$Time)
                                            df <- df[(df$Time == time_max)
                                                     & (df$VAF_primary >= min_vaf)
                                                     & (df$tumor_content == tumor_content),]
                                            df <- df[order(df$gene,
                                                           df$VAF_primary,
                                                           method='radix',
                                                           decreasing=c(FALSE, TRUE)), ]
                                            df$param <- get_param_id(file)
                                            return(df)
                                        }))
        params <- unique(sort(df_sim$param))

        dir.abc <- sprintf('%s/%s', ABC_dir, sample)
        dir.create(dir.abc, recursive = TRUE, showWarnings = FALSE)

        ## param_mat
        write.table(data.frame('param'=params),
                    sprintf('%s/param_mat.txt', dir.abc), row.names=F, quote=F, sep='\t')

        ## target
        df_obs <- rbind(get_sample(rint_file, sample),
                        get_sample(rother_file, sample))
        df_obs <- df_obs[, c(2, 3)]

        ## limit genes and renam
        genes <- unique(sort(df_sim$gene))
                                        # print(head(df_obs))
        df_obs <- df_obs[df_obs$gene %in% genes,]
                                        # print(df_obs$gene %in% genes)
        
        df_obs$gene <- rename_gene(rother_regex,
                                   sample,
                                   df_obs$gene)
        df_obs <- df_obs[order(df_obs$gene,
                               df_obs$vaf,
                               method='radix',
                               decreasing=c(FALSE, TRUE)), ]
        
        df_sim$gene <- rename_gene(rother_regex,
                                   sample,
                                   df_sim$gene)
        df_sim <- df_sim[order(df_sim$param,
                               df_sim$gene,
                               df_sim$VAF_primary,
                               method='radix',
                               decreasing=c(FALSE, FALSE, TRUE)), ]
        

        renamed_genes <- unique(sort(df_sim$gene))

        df_obs_extend <- get_extend_obs(renamed_genes, df_obs, df_sim)
        
        extended_obs_genes <- sprintf('%s.%s', df_obs_extend$gene, ave(df_obs_extend$gene, df_obs_extend$gene, FUN=seq_along))
        target <- df_obs_extend$vaf
        names(target) <- extended_obs_genes
        write.table(data.frame(gene=names(target), vaf=target),
                    sprintf('%s/target_tumor_content=%.02f.txt', dir.abc, tumor_content),
                    quote=F, row.names=F, sep='\t')

        ## sumstat
        write.table(df_sim,
                    sprintf('%s/sim_tumor_content=%.02f.txt', dir.abc, tumor_content),
                    row.names=F, quote=F, sep='\t')

        sumstat <- get_sumstat(renamed_genes,
                               target,
                               df_sim,
                               params)

        write.table(data.frame('param'=rownames(sumstat), sumstat),
                    sprintf('%s/sumstat_tumor_content=%.02f.txt', dir.abc, tumor_content),
                    row.names=F, quote=F, sep='\t')
    }
}
