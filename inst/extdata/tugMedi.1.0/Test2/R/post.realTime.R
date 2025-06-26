

post.realTime.R            <- new.env( parent = emptyenv() )
post.realTime.R$Ncells.mm3 <- 1e6
post.realTime.R$Tsize.Serr <- 1 # tumor size scale error


#
.get_N1N2 <- function( realTumorSize, 
                       Ncells.mm3 = post.realTime.R$Ncells.mm3 ) {
                       
   LD1.mm <- realTumorSize$LD1.mm
   SD1.mm <- realTumorSize$SD1.mm
   LD2.mm <- realTumorSize$LD2.mm
   SD2.mm <- realTumorSize$SD2.mm
   V1.mm3 <- realTumorSize$V1.mm3
   V2.mm3 <- realTumorSize$V2.mm3

   if        ( length( c( LD1.mm, SD1.mm, LD2.mm, SD2.mm ) ) == 4 ) {
      V1.mm3 <- ( LD1.mm * (SD1.mm^2) ) / 2
      V2.mm3 <- ( LD2.mm * (SD2.mm^2) ) / 2
   } else if ( length( c( V1.mm3, V2.mm3 )                 ) == 2 ) {
      ;
   } else {
      stop("Measurements are not complete to calculate number of cells!")
   } 

   N1 <- V1.mm3 * Ncells.mm3
   N2 <- V2.mm3 * Ncells.mm3

   return( c('N1'=N1, 'N2'=N2) )
}


#
.get_correspond_simTime <- function( clone, realTumorSize, 
                                     Tsize.Serr = post.realTime.R$Tsize.Serr ) {

   NN <- .get_N1N2( realTumorSize )

   Time_Ncells <- clone[, c( 'Time', realTumorSize$Ntype ) ]
   Time_Ncells <- unique( Time_Ncells )

   Ncells <- Time_Ncells[, realTumorSize$Ntype ]
   scl1 <- min( abs( log10(Ncells) - log10(NN[ 'N1' ]) ) ) 
   scl2 <- min( abs( log10(Ncells) - log10(NN[ 'N2' ]) ) ) 
   if ( scl1 > Tsize.Serr || 
        scl2 > Tsize.Serr    ) 
      stop( paste( 'Real tumor sizes from simulation\'s are too different!', 
                   paste( 'Scale ratio:', scl1, scl2, 'in log10' ), 
            sep="\n" ) )

   spi <- spline( Time_Ncells[, realTumorSize$Ntype ], Time_Ncells[, 'Time' ], 
                  xout = c( NN[ 'N1' ], NN[ 'N2' ] ), ties = mean )
   if ( spi$y[2] - spi$y[1] < 1 ) 
      stop( paste( 'Difference of interpolated times is < 1!', 
                   paste( 'At time 2,', spi$y[2] ), 
                   paste( 'At time 1,', spi$y[1] ), 
                   paste( 'Dif:', spi$y[2] - spi$y[1] ), 
            sep="\n" ) )
   min_max.sim <- c( min( Time_Ncells[, realTumorSize$Ntype ] ), 
                     max( Time_Ncells[, realTumorSize$Ntype ] ) )
   if ( any( NN < min_max.sim[1] | min_max.sim[2] < NN ) ) 
      stop( paste( 'Extrapolation happened!', 
                   paste( c( 'Real cell number @ t1 t2:', NN ),          collapse=' ' ), 
                   paste( c( 'Siml cell number min max:', min_max.sim) , collapse=' ' ), 
            sep="\n" ) )

   tt <- c( spi$y[1], spi$y[2] )
   names( tt ) <- c( '1', '2' )

   return( tt )
}


#
get_conversion_rate <- function( clone, 
                                 VDT           = NULL, 
                                 realTumorSize = NULL, 
                                 Days.month  = 365/12, 
                                 Months.year = 12 ) {

   if        ( !is.null( VDT ) && is.null( realTumorSize ) ) {

      nume <- VDT$VDT.days * log2( VDT$n2 / VDT$n1 )
      deno <- VDT$t2 - VDT$t1

   } else if ( is.null( VDT ) && !is.null( realTumorSize ) ) {

      nume <- realTumorSize$deltaT.days
      tt   <- .get_correspond_simTime( clone, realTumorSize )
      deno <- tt[ '2' ] - tt[ '1' ]

   } else {
      stop("Either of VDT or realTumorSize is necessary!")
   }
   
   tau <- ( nume / deno ) / c( 'day'   = 1, 
                               'month' = Days.month, 
                               'year'  = Days.month * Months.year )

  return( tau )
}


#
get_realTime_clone <- function( clone, 
                                tau, outUnit = 'year', 
                                col.replace = F ) {

   realTimeStep <- clone$Time * tau[ as.character( outUnit ) ]
   realTimeStep <- round( realTimeStep, digit = 2 )

   clone.new <- cbind( clone[, 1], realTimeStep, clone[, -1] )
   colnames( clone.new )[1] <- colnames( clone )[1]
   colnames( clone.new )[2] <- paste( 'Time', outUnit, sep = '.' )
   
   if ( col.replace == T ) clone.new <- clone.new[, -1]
   
   return( clone.new )
}


#
#' Function to rewrite output file with real time info
#'
#' @param input File's name for input parameters
#' @param force.stdout Logical to print parameters
#' @param cnvRate Conversion rate
#'
#' @return NULL 
#' 
#' @export
#'
#' @examples
#' NULL
write_realTime_clone <- function( input, force.stdout = F, cnvRate = NULL ) {

   input.table <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )
   
   VDT           <- NULL
   realTumorSize <- NULL
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]
      
      # add quotes
      if ( length( grep( "^input\\.", V1, perl=T ) ) > 0 || 
           V1 == 'out'                                   || 
           V1 == 'realTumorSize$Ntype'                   || 
           V1 == 'outUnit'
          ) V2 <- paste("\'", V2, "\'", sep='')
   
      text <- paste( V1, ' <- ', V2, sep='' )
      eval( parse( text = text ) )
   }

   clone <- read.table( input.cloneout, header=T, sep="\t", stringsAsFactors = FALSE )
   if ( is.null(cnvRate) ) 
      cnvRate <- get_conversion_rate( clone, 
                                      VDT, realTumorSize
                                     ) 
   table <- get_realTime_clone( clone, 
                                cnvRate, outUnit, 
                                col.replace
                               )

   if ( force.stdout == T ) { print(table); return( cnvRate ) }
   
   write.table( table, out, quote=F, sep="\t", row.names=F )
   return( cnvRate )
}


# unit test ###############################################
if ( 0 ) { # change 1 to 0 for silence

print("start")

clonefile <- 'Output/cloneout.txt'
rows <- 1:15
cols <- 1:8

old <- read.table( clonefile, header=T, sep="\t", stringsAsFactors = FALSE )
print( old[ rows, cols ] )


#
print("---test---")


VDT <- list( VDT.days = 200, 
             n1 = 1e5, n2 = 1e6, # sim
             t1 = 20,  t2 = 30   # sim
            )
realTumorSize <- NULL

new <- get_realTime_clone( clonefile   = clonefile, 
                           VDT         = VDT, 
                           outUnit     = 'year', 
                           col.replace = F )
print( new[ rows, cols ] )


#
print("---write test---")

VDT           <- NULL
realTumorSize <- NULL

write_realTime_clone( input = 'Input/ForPosts/realTime.parameters.txt', force.stdout = F )

print("---write test2---")
write_realTime_clone( input = 'Input/ForPosts/realTime.parameters.txt.test2' )

print("---write test3---")
write_realTime_clone( input = 'Input/ForPosts/realTime.parameters.txt.test3' )


print("---Done---")
}


