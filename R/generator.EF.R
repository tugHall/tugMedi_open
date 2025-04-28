

#
generate_EF.Rint <- function( skeleton = NULL, dist ) {

   if ( is.null( skeleton ) ) stop( 'Not yet coded for empty skeleton!' )

   bone <- read.table( file = skeleton, header = TRUE, sep = '\t', stringsAsFactors = FALSE )
   ret <- bone

   # timing
   timing <- intersect( c( 'When', 'Time_step' ), names( dist ) )
   if ( length( timing ) == 0 )    stop( 'Neither "When" nor "Time_step" in EF.Rint.parameters!' )
   if ( ! timing %in% names(ret) ) stop( c('No column of ', timing, ' in EF.Rint.skeleton!') )

   if      ( timing == 'When' )      { opt.list <-       dist[[ timing ]]   }
   else if ( timing == 'Time_step' ) { opt.list <- list( dist[[ timing ]] ) } # backward compatibility

   if ( ! is.null( opt.list ) ) {
      for ( ii in seq_along(bone$Mutation) ) {
         mut <- bone$Mutation[ ii ]
         opt <- .which_opt_with_mut( mut, opt.list )
         
         if ( ! is.null( opt ) ) {
            ret[[ timing ]][ ii ] <- round( .get_random_number( opt, mut = mut ) )
         } else {
            next
         }
         
      }
   }
   
   
   # Pchr
   opt <- dist[[ 'Pchr' ]]
   if ( ! is.null( opt ) ) {
      ret$Pchr <-
         sapply( X = as.character(bone$Mutation), function(X) {
                 round( .get_random_number( opt, mut = NULL ) )
                                                 } )
   }

   return( ret )
}


.which_opt_with_mut <- function( mut, opt.list ) {

   ret <- NULL
   
   for ( opt in opt.list ) {
      # rtrunc
      spec.stat <- c( 'lambda', 'mean', 'min', 'shape1', 'rate' )
      
      stat <- intersect( spec.stat, names( opt ) )
      if ( length( stat ) == 0 ) stop( 'No stat in spec.stat!' )
      
      if ( mut %in% names( opt[[ stat ]] ) ) {
         ret <- opt
         break
      }
      else {
         next
      }
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

   input.table <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )

   prm <- .get_prms_from_input.table( input.table )

   dist <- list()
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]

      # add quotes
      if ( length( grep( "^input\\.", V1, perl=T ) ) > 0 ||
           V1 == 'skeleton' ||
           V1 == 'out'      ||
           length( grep( "\\$spec$", V1, perl=T ) ) > 0 
          ) V2 <- paste("\'", V2, "\'", sep='')
      
      if ( length( grep( "^dist\\$When", V1, perl=T ) ) > 0
          ) V1 <- sub( "\\[\\s*(\\S+)\\s*\\]", "[ \\'\\1\\' ]", V1, perl=T )

      # prm.* to prm[ '*' ]
      if ( !is.null(prm) & length( grep( "prm\\.", V2, perl=T ) ) > 0 ) {
         V2 <- gsub( "prm\\.([a-zA-Z0-9_.]+)", "as.numeric( prm[ '\\1' ] )", V2 )
      }

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


.get_prms_from_input.table <- function( input.table ) {

   prms <- NULL
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]

      if ( V1 == 'input.prms' ) {
         prm_file <- V2
         prms <- .read_prms4gnrt( prm_file )
         break
      }
   }

   return( prms )
}


# test =================================================================
if ( 0 ) { # change 1 to 0 for silence

print("start")

library(futile.logger)

source("./generator.common.R")


#
print("--- test: old ver ---")
skeleton <- './EF.Rint.skeleton.txt'


if (0) {
dist <- list( Time_step = list( spec   = 'pois',
                                a      = 0,
                                b      = Inf,
                                lambda = c(M2=5, M3=20) 
                               ), 
              Pchr      = list( spec   = 'sample',
                                x      = c(1, 2),
                                prob   = c(1, 1) 
                               ) 
             )


EF.Rint <- generate_EF.Rint( skeleton = skeleton,
                             dist     = dist )

print(EF.Rint)
}


print("--- test2: new ver ---")


if (0) {
dist <- list( When = list( list( spec   = 'pois',
                                 a      = 0,
                                 b      = Inf,
                                 lambda = c(M10=3, M99=6) 
                                ), 
                           list( spec   = 'norm',
                                 a      = 0,
                                 b      = Inf,
                                 mean = c(M3=30), 
                                 sd   = c(M3=5) 
                                )
                          ), 
              Pchr = list( spec   = 'sample',
                           x      = c(1, 2),
                           prob   = c(1, 1) 
                          ) 
             )



EF.Rint <- generate_EF.Rint( skeleton = skeleton,
                             dist     = dist )


print(EF.Rint)
}


print("--- test3: old ver ---")


if (0) {
dist <- list( Time_step = list( spec = 'unif',
                                min  = c(M2=1, M3=4),
                                max  = c(M2=2, M3=5) 
                               ),
              Pchr      = list( spec   = 'sample',
                                x      = c(1, 2),
                                prob   = c(1, 1) 
                               ) 
             )


EF.Rint <- generate_EF.Rint( skeleton = skeleton,
                             dist     = dist )


print(EF.Rint)
}


#
print("---write---")
input <- './EF.Rint.parameters.txt'

write_EF.Rint( input  = input,
               force.stdout = T )

#write_EF.Rint( input  = input )


}


