

#
generate_CF <- function( skeleton = NULL, dist ) {

   if ( is.null( skeleton ) ) stop( "Not yet coded for empty skeleton!" )

   bone <- read.table( file = skeleton, header = TRUE, sep = '\t', stringsAsFactors = FALSE )
   vals          <- as.character( bone[, 2] )
   names(vals)   <- as.character( bone[, 1] )

   ret <- data.frame()   
   for ( prm in names(vals) ) {
      opt <- dist[[ prm ]]
      
      if ( is.null(opt) ) { 
         val <- vals[ prm ]
         names(val) <- NULL } 
      else {
         val <- .get_random_number( opt )
      } 

      df <- data.frame( Hallmark         = prm, 
                        CompactionFactor = val, 
                        stringsAsFactors = F )
      ret <- rbind( ret, df )
   }
   
   return( ret )
}


#
#' Function to write the compaction factors to the file
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
write_CF <- function( input, force.stdout = F ) {

   input.table <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )
   
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

   table <- generate_CF( skeleton = skeleton, 
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
print("---test---")
skeleton <- "Input/ForGenerators/CF.skeleton.txt"

CF <- generate_CF( skeleton = skeleton, 
                   dist = list( 
                             apoptosis = list( 
                                spec   = 'beta',
                                a      = 0, 
                                b      = 1, 
                                shape1 = 1, 
                                shape2 = 1 ), 
                             growth = list( 
                                spec = 'sample', 
                                x    = c(0.2, 0.1, 0.02), 
                                prob = rep(1, 3) ), 
                             immortalization = list( 
                                spec = 'norm', 
                                a    = 0,
                                b    = 1,
                                mean = 0.3, 
                                sd   = 0.2 ) ) )
                                  
print(CF)

#
print("---write---")

write_CF( input = 'Input/ForGenerators/CF.parameters.txt', force.stdout = T )

print("---write test2---")

write_CF( input = 'Input/ForGenerators/CF.parameters.txt.test2', force.stdout = T )


write_CF( input = 'Input/ForGenerators/CF.parameters.txt')
                         
}


