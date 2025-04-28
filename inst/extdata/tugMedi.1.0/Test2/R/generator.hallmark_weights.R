



#
.get_dirichlet.mode <- function( alphas ) {
   
   mode <- vector()
   for ( ii in seq_along( alphas ) ) {
   
      alpha <-      alphas[  ii ]
      beta  <- sum( alphas[ -ii ] )
   
      if ( alpha > 1 && beta > 1 ) {
                   mode[ ii ] <- ( alpha - 1 ) / ( alpha + beta - 2 )
      } else {
                if ( isTRUE( all.equal( alpha, beta ) ) ) { 
                   stop("No mode in calculation of dirichlet mode")
         } else if ( alpha >  beta) {
                   mode[ ii ] <- 1
         } else if ( alpha <  beta) {
                   mode[ ii ] <- 0
         } else {
                   stop("Strange case in calculation of dirichlet mode")
         }
      }
   }
   if ( sum( mode ) == 0 ) { mode[ which.max( alphas ) ] <- 1 }
   
   det <- all.equal( sum(mode), 1 )
   if ( ! isTRUE( det ) ) stop( 'sum(mode) != 1, ie, ', sum(mode) )

   return( mode )
}


#
.generate_dirichlet_numbers <- function( spec, alphas ) {

   if ( spec == 'dirichlet' ) {
      gnums <- rdirichlet( 1, alphas )[ 1, ] }
      
   if ( spec == 'dirichlet.mode' ) {
      gnums <- .get_dirichlet.mode( alphas ) }
      
   if ( spec == 'dirichlet.mean' ) {
      gnums <- alphas / sum( alphas ) }

   return( gnums )
}


# 
generate_weights <- function( skeleton = NULL, dist ) {
   if ( is.null( skeleton ) ) stop( "Not yet coded for empty skeleton!" )
   bone <- read.table( file = skeleton, header = TRUE, sep = '\t', stringsAsFactors = FALSE )

   ret <- data.frame()
   for ( hallmark in unique( bone$Hallmark ) ) {
      genes <- as.character( bone[ bone$Hallmark == hallmark, 'Gene' ] )
      opt <- dist[[ hallmark ]]
      
      if ( is.null(opt) ) { 
         gnums <- bone[ bone$Hallmark == hallmark, 'Weight' ] }
      else {
      
         if ( is.null( opt$alpha.column ) || opt$alpha.column == 'NULL' ) {
            alphas <- rep( 1, length(genes) ) } 
         else {
            col    <- opt$alpha.column
            alphas <- as.numeric( bone[ bone$Hallmark == hallmark, col ] ) }
         
         gnums <- .generate_dirichlet_numbers( opt$spec, alphas )
      } 

      df <- data.frame( Hallmark = rep( hallmark, length(genes) ), 
                        Gene     =      genes, 
                        Weight   =      gnums )
      ret <- rbind( ret, df )
   }
   
   return( ret )
}


#
#' Function to write the hallmarks weights to the file
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
write_weights <- function( input, force.stdout = F ) {

   input.table <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )
   
   dist <- list()
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]
      
      # add quotes
      if ( V1 == 'skeleton' || 
           V1 == 'out'      || 
           length( grep( "\\$spec$",         V1, perl=T ) ) > 0 || 
           length( grep( "\\$alpha.column$", V1, perl=T ) ) > 0
          ) V2 <- paste("\'", V2, "\'", sep='')
   
      text <- paste( V1, ' <- ', V2, sep='' )
      eval( parse( text = text ) )
   }

   table <- generate_weights( skeleton = skeleton, 
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

library(dirmult)

#
print("---test---")
skeleton <- "Input/ForGenerators/hallmark_weights.skeleton.txt"

weights <- generate_weights( skeleton = skeleton, 
                             dist = list( 
                                       apoptosis = list( 
                                          spec         = 'dirichlet',
                                          alpha.column = NULL ), 
                                       growth = list( 
                                          spec         = 'dirichlet', 
                                          alpha.column = 'Weight' ),
                                       angiogenesis = list( 
                                          spec         = 'dirichlet.mode', 
                                          alpha.column = 'Weight' ), 
                                       immortalization = list( 
                                          spec         = 'dirichlet.mean', 
                                          alpha.column = 'Weight' ) ) )

print(weights)


print("---test2---")
skeleton <- "Input/ForGenerators/hallmark_weights.skeleton.txt.test2"

weights <- generate_weights( skeleton = skeleton, 
                             dist = list( 
                                       apoptosis = list( 
                                          spec         = 'dirichlet',
                                          alpha.column = NULL ), 
                                       angiogenesis = list( 
                                          spec         = 'dirichlet.mode', 
                                          alpha.column = 'Weight' ), 
                                       immortalization = list( 
                                          spec         = 'dirichlet.mean', 
                                          alpha.column = 'Weight' ) ) )

print(weights)


#
print("---write---")

write_weights( input = 'Input/ForGenerators/hallmark_weights.parameters.txt', 
               force.stdout = T  )

print("---write test2---")
write_weights( input = 'Input/ForGenerators/hallmark_weights.parameters.txt.test2', 
               force.stdout = T  )

print("---write test3---")
write_weights( input = 'Input/ForGenerators/hallmark_weights.parameters.txt.test3', 
               force.stdout = T  )

write_weights( input = 'Input/ForGenerators/hallmark_weights.parameters.txt' )
                         
}


