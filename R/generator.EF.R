

#
generate_EF.Rint <- function( skeleton = NULL, dist ) {

   if ( is.null( skeleton ) ) stop( "Not yet coded for empty skeleton!" )

   bone <- read.table( file = skeleton, header = TRUE, sep = '\t', stringsAsFactors = FALSE )
   ret <- bone

   opt <- dist[[ 'Time_step' ]]
   if ( ! is.null( opt ) ) {
      ret$Time_step <-
         sapply( X = as.character(bone$Mutation), function(X) {
                 round( .get_random_number( opt, mut = X ) )
                                                 } )
   }

   opt <- dist[[ 'Pchr' ]]
   if ( ! is.null( opt ) ) {
      ret$Pchr <-
         sapply( X = as.character(bone$Mutation), function(X) {
                 round( .get_random_number( opt, mut = NULL ) )
                                                 } )
   }

   return( ret )
}


# 
#' Function to write EF.Rint data to the file
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
write_EF.Rint <- function( input, force.stdout = F ) {

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

      text <- paste( V1, ' <- ', V2, sep='' )
      eval( parse( text = text ) )
   }

   table <- generate_EF.Rint( skeleton = skeleton,
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

#
print("---test pois---")
skeleton <- "Input/ForGenerators/EF.Rint.skeleton.txt"

EF.Rint <- generate_EF.Rint( skeleton = skeleton,
                             dist = list(
                                Time_step = list(
                                   spec   = 'pois',
                                   a      = 0,
                                   b      = Inf,
                                   lambda = c(M2=5, M3=20) ),
                                Pchr = list(
                                   spec = 'sample',
                                   x    = c(1, 2),
                                   prob = c(1, 1) ) ) )

print(EF.Rint)


print("---test2 norm---")

EF.Rint <- generate_EF.Rint( skeleton = skeleton,
                             dist = list(
                                Time_step = list(
                                   spec = 'norm',
                                   a    = 0,
                                   b    = Inf,
                                   mean = c(M2=5, M3=20),
                                   sd   = c(M2=1, M3=10) ),
                                Pchr = list(
                                   spec = 'sample',
                                   x    = c(1, 2),
                                   prob = c(1, 1) ) ) )

print(EF.Rint)


print("---test3 unif---")

EF.Rint <- generate_EF.Rint( skeleton = skeleton,
                             dist = list(
                                Time_step = list(
                                   spec = 'unif',
                                   min  = c(M2=0, M3=0),
                                   max  = c(M2=10, M3=100) ),
                                Pchr = list(
                                   spec = 'sample',
                                   x    = c(1, 2),
                                   prob = c(1, 1) ) ) )

print(EF.Rint)


#
print("---write---")

write_EF.Rint( input  = './Input/ForGenerators/EF.Rint.parameters.txt',
               force.stdout = T )

write_EF.Rint( input  = './Input/ForGenerators/EF.Rint.parameters.txt' )

}


