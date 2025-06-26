

#
generate_Rint <- function( skeleton = NULL, dist ) {

   if ( is.null( skeleton ) ) stop( "Not yet coded for empty skeleton!" )

   bone <- read.table( file = skeleton, header = TRUE, sep = '\t', stringsAsFactors = FALSE )
   ret <- bone

   for ( ii in 1:nrow(ret) ) {
   
      g1 <- as.vector( ret[ ii, ]$Gene )
      if ( g1 %in% names( dist$Type ) ) {
         opt <- dist$Type[[ g1 ]]
         ret[ ii, ]$Type <- .get_random_number( opt, mut = NULL )
      }
      
   }

   return( ret )
}


# 
#' Function to write Rint data to the file
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
write_Rint <- function( input, force.stdout = F ) {

   input.table   <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )

   dist <- list()
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]

      # add quotes
      if ( V1 == 'skeleton' ||
           V1 == 'out'      ||
           length( grep( "\\$spec$", V1, perl=T ) ) > 0 
          ) V2 <- paste("\'", V2, "\'", sep='')
      if ( length( grep( "\\$x$",    V1, perl=T ) ) > 0
          ) V2 <- gsub( "(s|o|dn)", "\\'\\1\\'", V2, perl=T )

      if ( length( grep( "^dist\\$Type", V1, perl=T ) ) > 0
          ) V1 <- sub( "\\[\\s*(\\S+)\\s*\\]", "[ \\'\\1\\' ]", V1, perl=T )

      text <- paste( V1, ' <- ', V2, sep='' )
      eval( parse( text = text ) )
   }

   table <- generate_Rint( skeleton = skeleton,
                           dist     = dist )

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

write_Rint( input = './Input/ForGenerators/Rint.parameters.txt',
            force.stdout = T )

write_Rint( input = './Input/ForGenerators/Rint.parameters.txt' )

}


