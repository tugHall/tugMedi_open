

#
get_small_ccds <- function( rother.file, ccds.db ) {

   rother   <- read.table( rother.file, header=T, sep="\t", stringsAsFactors = FALSE )
   ccds.ids <- as.character( rother$CDS_ID )

   table <- select_cds( ls = ccds.ids, f_in = ccds.db )
   return( table )
}


#
generate_cloneinit <- function( skeleton = NULL, pom, rother_file, ccds_file ) {
   if ( is.null( skeleton ) ) stop( "Not yet coded for empty skeleton!" )
   bone <- read.table( file = skeleton, header = TRUE, sep = '\t', stringsAsFactors = FALSE )
   
   ccds    <- get_small_ccds( rother_file, ccds_file )
   lengths <- get_lengths(    rother_file, ccds_file )

   ret <- data.frame()
   for ( clone in unique( bone$cloneID)  ) {
      slct <- bone[ bone$cloneID == clone, ]
      ret  <- rbind( ret, slct )

      if ( is.null( pom[[ as.character(clone) ]] ) ) next
      N_pom <- .get_random_number( opt = pom[[ as.character(clone) ]]$dist )
      if ( N_pom == 0 ) next
      for ( ii in 1:N_pom ) {
         mode <- 'pom'
         gene <- sample( names( lengths[[ mode ]] ), 1, prob = lengths[[ mode ]] )
      
         ccds.gene <- ccds[ ccds$Gene == gene, ]
         chr  <- ccds.gene$Chr[1]
         irow <- sample( 1:nrow(ccds.gene), 1, prob = ccds.gene$Len)
         pos  <- sample( ccds.gene[ irow, ]$Start:ccds.gene[ irow, ]$End, 1 )
      
         ret <- rbind( ret, 
                       data.frame( 'cloneID' = slct$cloneID[1], 
                                   'Ncells'  = slct$Ncells[1], 
                                   'Type'    = mode, 
                                   'Gene'    = gene, 
                                   'Chr'     = chr, 
                                   'Chr_stt' = pos, 
                                   'Chr_end' = 'NA', 
                                   'Pchr'    = sample( 1:2, 1 ), 
                                   'DrvPss'  = 'pss' ) )
      }
   }
   
   return( ret )
}


#
write_cloneinit <- function( input, force.stdout = F ) {

   input.table   <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )

   pom <- list()
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]

      # add quotes
      if ( length( grep( "^input\\.", V1, perl=T ) ) > 0 ||
           V1 == 'skeleton' ||
           V1 == 'out'      ||
           length( grep( "\\$spec$", V1, perl=T ) ) > 0
          ) V2 <- paste("\'", V2, "\'", sep='')
      if ( length( grep( "^pom\\[.+\\]", V1, perl=T ) ) > 0
          ) V1 <- sub( "\\[\\s*(\\S+)\\s*\\]", "[ \\'\\1\\' ]", V1, perl=T )

      text <- paste( V1, ' <- ', V2, sep='' )
      eval( parse( text = text ) )
   }

   table <- generate_cloneinit( skeleton = skeleton,
                                pom      = pom,
                                rother_file = input.Rother,
                                ccds_file   = input.CCDSdatabase )

   if ( force.stdout == T ) { print(table); return() }

   if ( file.exists( out ) ) {
      out2 <- paste( out, 'SAVED', sep='.' )
      file.copy( out, out2, overwrite = T )
   }

   write.table( table, out, quote=F, sep="\t", row.names=F )
}



# test
if ( 0 ) { # change 1 to 0 for silence

print("start")

source("R/generator.common.R")

write_cloneinit( input  = './Input/ForGenerators/cloneinit.parameters.txt',
                 force.stdout = T )

write_cloneinit( input  = './Input/ForGenerators/cloneinit.parameters.txt' )

}


