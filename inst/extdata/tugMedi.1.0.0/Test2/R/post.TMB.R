



#
.get_Nmut <- function( vaf.file, vaf.opt, rother.file ) {
   type <- vaf.opt$type
   LOD  <- vaf.opt$LOD
   
   vaf_all <- read.table( vaf.file, header=T, sep="\t", stringsAsFactors = FALSE )
   tm      <- unique( vaf_all$Time )
   tc      <- unique( vaf_all$tumor_content )
   
   rother  =  read.table( file = rother.file, header = TRUE, sep = '\t', stringsAsFactors = FALSE )
   
   vaf_all <- vaf_all[ vaf_all[ ,  type ]  >=   LOD,         ]
   vaf_all <- vaf_all[ vaf_all[ , 'gene' ] %in% rother$Gene, ]
   
   DF  =  data.frame( Time          = rep( tm, each = length( tc ) ), 
                      tumor_content = rep( tc, length( tm ) ),
                      TMB           = 0 )
   for( i in 1:nrow( DF ) ){
       tc_1 = DF$tumor_content[ i ]
       tm_1 = DF$Time[ i ]
       
       vaf <- vaf_all[ vaf_all$tumor_content == tc_1 & 
                       vaf_all$Time          == tm_1,  ]
       DF$TMB[ i ] <- nrow( vaf )
   }

   return( DF )
}


#
.get_Rsize <- function( rother.file, ccds.db  ) {

   rother <- read.table( rother.file, header=T, sep="\t", stringsAsFactors = FALSE )
   ccds.ids <- as.character( rother$CDS_ID )

   region <- make_map( ls = ccds.ids, f_in = ccds.db )

   rsize <- sum( region$Len )

   return( rsize )
}


#
get_TMB <- function( Rother.file, VAF.file, VAF.opt, CCDS.db, Mb = 1e6 ) {

   DF    <- .get_Nmut( VAF.file, VAF.opt, Rother.file )
   Rsize <- .get_Rsize( Rother.file, CCDS.db )

   DF$TMB <- ( DF$TMB / Rsize ) * Mb

   colnames(DF)[ 3 ] = paste0( 'TMBvaf', VAF.opt$LOD * 100, '%' )

   return( DF )
}


#
#' Function to write TMB data to the file
#'
#' @param input File's name for input parameters
#' @param force.stdout Logical to print parameters
#'
#' @return NULL 
#' 
#' @export
#'
#' @examples
#' NULL
write_TMB <- function( input, force.stdout = F ) {

   input.table <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )

   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]

      # add quotes
      if ( length( grep( "^input\\.", V1, perl=T ) ) > 0 ||
           V1 == 'out'                                   ||
           V1 == 'VAF.type'
          ) V2 <- paste("\'", V2, "\'", sep='')

      text <- paste( V1, ' <- ', V2, sep='' )
      eval( parse( text = text ) )
   }

   table <- get_TMB( Rother.file = input.Rother,
                     VAF.file    = input.VAF,
                     VAF.opt     = list( type = VAF.type,
                                         LOD  = VAF.LOD        # Limit Of Detection
                                        ),
                     CCDS.db     = input.CCDSdatabase )

   if ( force.stdout == T ) { print(table); return() }
   
   dir.create( dirname( out ), showWarnings = F, recursive = T )
   write.table( table, out, quote=F, sep="\t", row.names=F )
}


##### test #########################################
if ( 0 ) { # change 1 to 0 for silence

print("start")


library(stringr)

files <- c("../../tugHall.3/R/utils.R")
for (file in files) {
   message("For development loading: ", file)
   source(file)
}


#
print("---test---")

Rother.file <- 'Input/Rother.txt'

VAF.file <- 'Output/VAF/VAF.txt'
VAF.opt <- list( type = 'VAF_metastatic',
                 LOD  = 0.05               # Limit Of Detection
                )

CCDS.db <- 'Input/DATA/CCDS.current.txt'


TMB <- get_TMB( Rother.file, VAF.file, VAF.opt, CCDS.db )
print(TMB)


#
print("---write test---")

write_TMB( input = 'Input/ForPosts/TMB.parameters.txt',
           force.stdout = F )


print("---Done---")
}


