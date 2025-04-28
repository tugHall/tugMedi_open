

#
.read_prms4gnrt <- function( input ) {

   input.table <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )

   prms <- list()
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]

      # add quotes
      if ( length( grep( "^input\\.", V1, perl=T ) ) > 0 ||
           V1 == 'out'                                   ||
           length( grep( "\\$spec$",  V1, perl=T ) ) > 0 ||
           V1 == 'Parameter'                            ||
           V1 == 'metastatic_model'
          ) V2 <- paste("\'", V2, "\'", sep='')

      text <- paste0( V1, ' <- ', V2 )
      text <- paste0( 'prms$', text )
      eval( parse( text = text ) )
   }

   return( prms )
}


# 
.insert_poms_at <- function( prms1, prms2, lengths ) {
   if ( is.null( prms1$N_pom$dist ) ) return( )

   N_pom <- .get_random_number( opt = prms1$N_pom$dist )

   tbl <- data.frame()
   for ( ii in 1:N_pom ) {
      ndiv <- prms1$N_pom$ndiv
      mode <- 'pom'
      gene <- sample( names( lengths[[ mode ]] ), 1, prob = lengths[[ mode ]] )
   
      tbl <- rbind( tbl, 
                    data.frame( 'num.d' = ndiv, 
                                'mode'  = mode, 
                                'gene'  = gene ) )
   }

   return( tbl )
}


# To keep the mean same
.scale_correction <- function( dist ) {

   if        ( dist$EEL$spec == 'gamma' )   { 
      cor <- 1 / dist$EEL$shape 
   } else if ( dist$EEL$spec == 'weibull' ) { 
      cor <- 1 / gamma( 1 + 1 / dist$EEL$shape ) 
   } else if ( dist$EEL$spec == 'exp')      {
      cor <- 1
   } else                                   { 
      stop( 'Cannot correct scale parameter!' ) 
   }

   return( cor )
}


#
.get_waiting_divisions <- function( prms1, prms2, lengths ) {

   rates <- list( 'pom' = .get_mut_mal_rates( 'pom', prms2 ), 
                  'del' = .get_mut_mal_rates( 'del', prms2 ), 
                  'dup' = .get_mut_mal_rates( 'dup', prms2 ) )
   if ( prms1$CNA_presence == F ) { 
      rates[[ 'del' ]][ 'm' ] <- 0
      rates[[ 'dup' ]][ 'm' ] <- 0
   }

   modes <- sample( c( 'pom', 'del', 'dup' ), prms1$N_mut, rep=T, 
                    prob = sapply( X=c( 'pom', 'del', 'dup' ), function(X) 
                       4 * rates[[ X ]][ 'm' ] * ( 1 - rates[[ X ]][ 'u' ] ) * 
                       sum( unlist( lengths[[ X ]] ) ) ) # sum across genes 
                   )

   wd <- data.frame()
   for ( ii in 1:prms1$N_mut ) {
      mode <- modes[[ ii ]]
      gene <- sample( names( lengths[[ mode ]] ), 1, prob = lengths[[ mode ]] )

      s.lam <- 4 * rates[[ mode ]][ 'm' ] * ( 1 - rates[[ mode ]][ 'u' ] ) * 
               lengths[[ mode ]][[ gene ]]
      opt.w <- list( spec  = prms1$dist$EEL$spec, a = -Inf, b = Inf,
                     rate  = s.lam,
                     shape = prms1$dist$EEL$shape,
                     scale = ( 1 / s.lam ) * 
                                 .scale_correction( prms1$dist ) )

      num.d <- .get_random_number( opt = opt.w ) 

      wd <- rbind( wd, data.frame( num.d, mode, gene ) )
   }
   wd$num.d <- ceiling( wd$num.d )

   return( wd )
}


#
.get_mut_mal_rates <- function( mode, prms2 ) {

   if        ( mode == 'pom' ) {

      mm <- prms2$m0
      uu <- prms2$uo + prms2$us

   } else if ( mode == 'dup' ) {

      mm <- prms2$m_dup
      uu <- prms2$uo_dup + prms2$us_dup

   } else if ( mode == 'del' ) {

      mm <- prms2$m_del
      uu <- prms2$uo_del + prms2$us_del

   } else {
      stop( 'Exceptional mode! mode: ', mode )
   }

   return( c( 'm'=mm, 'u'=uu ) )
}


#
.arrange_wds <- function( table ) {

   colnames(table) <- c( 'Waiting_division', 'Type', 'Gene' )

   table$Waiting_division <- as.numeric( table$Waiting_division )
   table$Type             <- as.character( table$Type )
   table$Gene             <- as.character( table$Gene )

   table <- table[ order( table$Waiting_division,
                          table$Gene, 
                          table$Type 
                         ), ]
                         
   table$Mutation <- paste0( 'A', seq( 1:nrow(table) ) )
   table <- table[, c(1, 4, 2, 3) ]
   
   return( table )
}


#
#' Function to write EF.Rother data to the file
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
write_EF.Rother <- function( input, force.stdout = F ) {

   prms1   <- .read_prms4gnrt( input )
   prms2   <- .read_prms4gnrt( prms1$input.prms )
   lengths <- get_lengths( prms1$input.Rother, prms1$input.CCDSdatabase )

   table1 <- .insert_poms_at( prms1, prms2, lengths )
   table2 <- .get_waiting_divisions( prms1, prms2, lengths )
   table  <- rbind( table1, table2 )
   if ( nrow(table) == 0 ) stop( 'No waiting division data!' )
   table <- .arrange_wds( table )

   if ( force.stdout == T ) { print(table); return() }

   out <- prms1$out
   if ( file.exists( out ) ) {
      out2 <- paste( out, 'SAVED', sep='.' )
      file.copy( out, out2, overwrite = T )
   }

   write.table( table, out, quote=F, sep="\t", row.names=F )
}


##### test #########################################
if ( 0 ) { # change 1 to 0 for silence

print("start")


library("truncdist")
source( "R/generator.common.R" )


#
print("---test---")


#sink("tmp")
write_EF.Rother( input = 'Input/ForGenerators/EF.Rother.parameters.txt',
                 force.stdout = F )
#sink()


print("---Done---")
}


