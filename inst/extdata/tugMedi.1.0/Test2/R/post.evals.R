

post.evals.R             <- new.env( parent = emptyenv() )
post.evals.R$gene.Rother <- '_Rother'


#
calc_MNE <- function( obs, sim ) {
   if ( length( obs ) != 1 ) { stop( "obs should be scalar." ) }
   
   stat <- obs - mean( sim )
   
   return( stat )
}


#
calc_RMSE <- function( obs, sim ) {
   if ( length( obs ) != 1 ) { stop( "obs should be scalar." ) }
   
   ttl  <- sum( ( obs - sim )^2 ) / length( sim ) 
   stat <- sqrt ( ttl )
   
   return( stat )
}


#
.format_common_obs <- function( obs, LOD ) {

   tbl <- read.table( obs$file, header=T, sep="\t", 
                      quote="\"", stringsAsFactors = F ) # Note quote
   names(tbl)[ c( obs$col.gene, obs$col.VAF, obs$col.sample ) ] <- 
               c(        'gene',       'VAF',       'sample' )
   
   tbl2 <- tbl[ tbl$sample == obs$sample, ]
   if ( nrow(tbl2) == 0 ) stop( 'No data selected from obs with sample: ',
                                 obs$sample )
   tbl2 <- tbl2[ tbl2$VAF >= LOD,                              ]
   tbl2 <- tbl2[                , c( 'gene', 'VAF', 'sample' ) ]

   new <- data.frame()
   if ( ! is.null(obs$glob) ) 
      tbl2$gene[ grep( glob2rx( obs$glob ), tbl2$gene ) ] <- post.evals.R$gene.Rother
   for ( gene in unique(tbl2$gene) ) {
      tbl3 <- tbl2[ tbl2$gene == gene, ]
      
      tbl3 <- tbl3[ order( tbl3$VAF, decreasing=T ), ]
      tbl3 <- cbind( 'srtd_id'= paste( gene, 1:nrow(tbl3), sep='.' ), 
                     tbl3, 
                     stringsAsFactors = F )
      
      new <- rbind( new, tbl3 )
   }

   return( new )
}


#
.reformat_obs <- function( obs, LOD ) {
  
   obs2        <- obs$Rint
   obs2$sample <- obs$sample
   tbl.Rint    <- .format_common_obs( obs2, LOD )
   
   obs2        <- obs$Rother
   obs2$sample <- obs$sample
   tbl.Rother  <- .format_common_obs( obs2, LOD )

   new <- rbind( tbl.Rint, tbl.Rother )

   return( new )
}


#
.reformat_sim <- function( file, type, tc, glob ) {

   tbl  <- read.table( file, header=T, sep="\t", 
                       stringsAsFactors = F )

   tbl2      <- tbl[ tbl$Time          ==   max( as.numeric(tbl$Time) ) & 
                     tbl$tumor_content %in% tc, ]
   tbl2$Time <- NULL

   new <- data.frame()
   tbl2$gene[ grep( glob2rx( glob ), tbl2$gene ) ] <- post.evals.R$gene.Rother
   for ( tc in unique(tbl2$tumor_content) ) {
      tbl3 <- tbl2[ tbl2$tumor_content == tc, ]
   
      for ( gene in unique(tbl3$gene) ) {
         tbl4 <- tbl3[ tbl3$gene == gene, ]
         
         tbl4 <- tbl4[ order( tbl4[, type], decreasing=T ), ]
         tbl4 <- cbind( 'srtd_id'= paste( gene, 1:nrow(tbl4), sep='.' ), 
                        tbl4, 
                        stringsAsFactors = F )
         
         new <- rbind( new, tbl4 )
      }
   }

   return( new )
}


#
.get_big.sim <- function( sim, LOD )  {

   files <- Sys.glob( sim$glob )
   flog.info( paste( length(files), 'VAF.txt files are globbed for eval.') )

   big <- data.frame()
   cnt <- 0
   for ( file in files ) {
      cnt <- cnt + 1

      sim.tbl <- .reformat_sim( file, sim$type, sim$tc, sim$Rother$glob )
      sim.tbl <- cbind( 'file_no' = cnt, sim.tbl )
      
      big <- rbind( big, 
                    sim.tbl[ sim.tbl[ ,sim$type ] >= LOD, 
                             c( 'file_no', 'srtd_id', 'tumor_content', sim$type ) ], 
                    stringsAsFactors = F )
   }

   return( big ) 
}


# 
.get_obs_sims <- function( obs, sim, LOD ) {

   one.obs <- .reformat_obs( obs, LOD )
   big.sim <- .get_big.sim(  sim, LOD )

   if ( nrow(one.obs) == 0 ) stop(    "No data in obs for evaluations after LOD cut! LOD: ", LOD )
   if ( nrow(big.sim) == 0 ) stop( c( "No data in sim for evaluations after LOD cut! LOD: ", LOD, 
                                      "\n  First, check the glob pattern: ", sim$glob) )
   
   obs.smpl <- unique( one.obs$sample )
   if ( length(obs.smpl) != 1 ) stop( 'length(sample) is not 1!' )

   mat <- data.frame()
   for ( id in unique( c( one.obs$srtd_id, big.sim$srtd_id ) ) ) {
   
      obs.vaf <- one.obs[ one.obs$srtd_id == id, 'VAF' ]
      if ( length(obs.vaf) == 0 ) obs.vaf <- 0

      for ( tc in unique( big.sim$tumor_content ) ) {
         vafs <- rep( 0, sim$n.rep ) # for num of VAF files < n.rep
         for ( no in unique( big.sim$file_no ) ) {
            vaf <- big.sim[ big.sim$srtd_id       == id & 
                            big.sim$tumor_content == tc & 
                            big.sim$file_no       == no, 
                               sim$type ]
            if ( length(vaf) == 0 ) vaf <- 0
            
            vafs[ no ] <- vaf
         }
         row <- cbind( data.frame( tc, id ), 
                       data.frame( obs.vaf, t(vafs) ) )
         mat <- rbind( mat, row, stringsAsFactors = F )
      }
   }
   mat <- cbind( obs.smpl, mat ) 
   mat <- .format_obs_sims( mat, sim )
   
   return( mat )
}


#
.format_obs_sims <- function( mat, sim ) {

   names(mat) <- c( 'sample', 'tumor_content', 'id', 'obs', 
                    paste( 'sim', 1:sim$n.rep, sep='.' ) )
                    
   if ( 4+sim$n.rep < ncol(mat) ) 
      names(mat)[ (4+sim$n.rep+1):ncol(mat) ] <- 
         paste( 'out', (sim$n.rep+1):(ncol(mat)-4), sep='.' ) 
   
   mat$tumor_content <- as.numeric( mat$tumor_content )
   mat <- mat[ order( mat$id ), ]
   mat <- mat[ order( mat$tumor_content, decreasing = T ), ]

   return( mat )
}


# 
write_obs_sims <- function( obs_sims, out ) {
   
   if ( file.exists( out ) ) {
      out2 <- paste( out, 'SAVED', sep='.' )
      file.copy( out, out2, overwrite = T )
   }
   
   write.table( obs_sims, out, quote=F, sep="\t", row.names=F )
}


#
.get_stats <- function( mat ) {

   stats <- data.frame()
   for ( ii in 1:nrow(mat) ) {
      obs <- mat[ ii, 4           ] 
      sim <- mat[ ii, 5:ncol(mat) ]

      obs <-                            as.numeric( as.character( obs ) )
      sim <- sapply( X=sim, function(X) as.numeric( as.character(  X  ) ) )

      men  <- mean(     sim )
      sdv  <- sd(       sim )
      qq1  <- quantile( sim, prob=0.25 )
      med  <- quantile( sim, prob=0.5  )
      qq3  <- quantile( sim, prob=0.75 )
      mne  <- calc_MNE(  obs, sim )
      rmse <- calc_RMSE( obs, sim )

      row <- cbind( mat[ ii, c('sample', 'tumor_content', 'id') ], 
                    'obs'=obs, 'sim.mean'=men, 'sim.med'=med, 
                    'sim.sd'=sdv, 'sim.Q1'=qq1, 'sim.Q3'=qq3,
                    'ME'=mne, 'RMSE'=rmse )
      stats <- rbind( stats, row )
   }

   return( stats )
}


# 
write_stats <- function( obs_sims, out, survivorOnly = F ) {
   
   if ( survivorOnly == F ) {
      ;
   } else {
      out <- sub( "\\.([^\\.]*)$", ".survivor.\\1", out, perl=T )
      
      flog.info(paste( 'n.sim all      is:', ncol(obs_sims) - 4 ))
      obs_sims <- obs_sims[ , ! c( rep( F, 4 ), 
                                   apply( obs_sims[,-(1:4)], 2, sum ) == 0 
                                  ) 
                           ] # cols 1-4 are irrelevant to sum()
      flog.info(paste( 'n.sim survivor is:', ncol(obs_sims) - 4 ))
   }
   
   stats <- .get_stats( obs_sims )
   
   if ( file.exists( out ) ) {
      out2 <- paste( out, 'SAVED', sep='.' )
      file.copy( out, out2, overwrite = T )
   }
   
   write.table( stats, out, quote=F, sep="\t", row.names=F )
}


#' Function to write evals data to the file
#'
#' @param input File's name for input parameters
#'
#' @return NULL 
#' 
#' @export
#'
#' @examples
#' NULL
write_evals <- function( input ) {

   input.table <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )
   
   obs <- list()
   sim <- list()
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]
      
      # add quotes
      if ( length( grep( "^out\\.", V1, perl=T ) ) > 0 || 
           V1 == 'obs$sample'                          || 
           V1 == 'obs$Rint$file'                       || 
           V1 == 'obs$Rother$file'                     || 
           V1 == 'obs$Rother$glob'                     || 
           V1 == 'sim$glob'                            || 
           V1 == 'sim$type'                            || 
           V1 == 'sim$Rother$glob'
          ) V2 <- paste("\'", V2, "\'", sep='')
   
      text <- paste( V1, ' <- ', V2, sep='' )
      eval( parse( text = text ) )
   }

   obs_sims <- .get_obs_sims( obs, sim, LOD )
   write_obs_sims( obs_sims, out.obs_sims )
   write_stats(    obs_sims, out.stats    )
   if ( exists( 'survivorToo' ) && survivorToo == T ) 
      write_stats( obs_sims, out.stats, survivorOnly = T )
   flog.info( "Wrote eval files." )
}




##### unit test #########################################
if ( 0 ) { # change 1 to 0 for silence

library(futile.logger)

print("start")

print("---test---")

####### obs #######

LOD <- 0.1 # Limit Of Detection

obs <- list( sample = 'TCGA-DM-A28E-01A-11D-A16V-10'
             #sample = 'TCGA-AZ-6608-01A-11D-1835-10'
            )
obs$Rint   <- list( file = 'Input/Samples/samples.Rint.txt', 
                    col.gene   = 1, 
                    col.sample = 6, 
                    col.VAF    = 11 
                   )
obs$Rother <- list( file = 'Input/Samples/samples.Rother.txt', 
                    col.gene   = 1, 
                    col.sample = 6, 
                    col.VAF    = 11 
                   )

obs.tbl <- .reformat_obs( obs, LOD )

out.obs  <- './evals.obs.txt'
write.table( obs.tbl, out.obs, quote=F, sep="\t", row.names=F )


####### sim #######

sim <- list( glob   = './Output*/VAF/VAF.txt', 
             type  = 'VAF_primary', 
             tc    = c( 1.0, 0.9, 0.8 ), # tumor content
             n.rep = 4
            )

big.sim <- .get_big.sim( sim, LOD )

out.sims <- './evals.sims.txt'
write.table( big.sim, out.sims, quote=F, sep="\t", row.names=F )


####### obs_sims #######

obs_sims <- .get_obs_sims( obs, sim, LOD )

out.oss <- './evals.obs_sims.txt'
write_obs_sims( obs_sims, out.oss )


####### evals #######
print("---write test---")


write_evals( input = 'Input/ForPosts/evals.parameters.txt' )


print("---Done---")
}


