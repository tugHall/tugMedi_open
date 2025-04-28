

load_lib_functions4generators <- function( ){

   pkgs = list( 
      truncdist = c( 'rtrunc' ),
      dirmult   = c( 'rdirichlet' ),
      MASS      = c( 'mvrnorm' )
               )

   for( pck in names( pkgs ) ){
      require( package = pck, character.only = TRUE, include.only = pkgs[[ pck ]])
   }

}

load_lib_functions4generators()


#
#
.get_random_number <- function( opt, mut = NULL ) {

   if ( opt$spec == 'sample' ) {
      if ( length( opt$x ) != 1 ) {
         rnum <- sample( opt$x, 1, prob = opt$prob, rep=T ) 
      } else {
         flog.info( paste0( 'x length of sample(x, ...) is 1. ', 
                            'Just used x: ', opt$x ) )
         rnum <- opt$x
      }
   }
   if ( opt$spec == 'pois' ) {
      if ( ! is.null( mut ) ) { opt$lambda <- opt$lambda[ mut ] }

      rnum <- rtrunc( 1, spec = opt$spec, a = opt$a, b = opt$b,
                         lambda = opt$lambda ) 
   }
   if ( opt$spec == 'norm' ) {
      if ( ! is.null( mut ) ) { opt$mean <- opt$mean[ mut ];
                                opt$sd   <- opt$sd[ mut ] }

      rnum <- rtrunc( 1, spec = opt$spec, a = opt$a, b = opt$b,
                         mean = opt$mean, sd = opt$sd ) 
   }
   if ( opt$spec == 'unif' ) {
      if ( ! is.null( mut ) ) { opt$min <- opt$min[ mut ];
                                opt$max <- opt$max[ mut ] }

      rnum <- rtrunc( 1, spec = opt$spec, a = -Inf, b = Inf,
                         min = opt$min, max = opt$max ) 
   }
   if ( opt$spec == 'beta' ) {
      if ( ! is.null( mut ) ) { opt$shape1 <- opt$shape1[ mut ];
                                opt$shape2 <- opt$shape2[ mut ] }

      rnum <- rtrunc( 1, spec = opt$spec, a = opt$a, b = opt$b,
                         shape1 = opt$shape1, shape2 = opt$shape2 ) 
   }


   if ( opt$spec == 'exp' ) {
      if ( ! is.null( mut ) ) { opt$rate <- opt$rate[ mut ] }
   
      rnum <- rtrunc( 1, spec = opt$spec, a = opt$a, b = opt$b,
                         rate = opt$rate )
   }
   if ( opt$spec == 'gamma' ) {
      rnum <- rtrunc( 1, spec = opt$spec, a = opt$a, b = opt$b,
                         shape = opt$shape, scale = opt$scale )
   }
   if ( opt$spec == 'weibull' ) {
      rnum <- rtrunc( 1, spec = opt$spec, a = opt$a, b = opt$b,
                         shape = opt$shape, scale = opt$scale )
   }


   if ( opt$spec == 'NULL' || is.null( opt$spec ) ) {
      rnum <- rtrunc( 1, spec = 'beta', a = 0, b = 1,
                         shape1 = 1, shape2 = 1 ) 
   }


   return( rnum )
}


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
get_lengths <- function( rother.file, ccds.db  ) {

   rother <- read.table( rother.file, header=T, sep="\t", stringsAsFactors = FALSE )
   ccds.ids <- as.character( rother$CDS_ID )

   table <- select_cds( ls = ccds.ids, f_in = ccds.db )

   lengths <- list()
   for ( id in unique( table$Gene ) ) {
      table.1 <- table[ table$Gene == id, ]

      # size exons
      lengths[[ 'pom' ]][[ id ]] <- sum( table.1$Len )
      # size cds
      lengths[[ 'dup' ]][[ id ]] <- max( table.1$End ) - min( table.1$Start ) + 1
      lengths[[ 'del' ]][[ id ]] <- lengths[[ 'dup' ]][[ id ]]
   }

   return( lengths )
}


#
select_cds <- function( ls, f_in ) {

    df <- read.delim( file = f_in, sep = '\t', header = TRUE, stringsAsFactors = FALSE )

    if ( any( (ls %in% df$ccds_id) == F ) )
       stop( "Not all of given CCDS IDs are found in the database!" )

    ###  Get the only needed data
    dt <- df[ df$ccds_id %in% ls, ]
    rm( df )

    ### get a vector of strings:
    dfALL  <-  NULL
    for ( ii in 1:nrow(dt) ) {
        lct <- as.character( dt[ ii, 'cds_locations' ] )
        lct <- gsub( '[\\[\\]\\s]', '', lct, perl=T )
        lct <- strsplit( lct, split = ',' )
        lct <- lct[[1]]   # returned as [[1]]
        lct <- strsplit( lct, split = '-' )

        df <- data.frame( matrix( unlist(lct), nrow=length(lct), byrow=TRUE ),
                          stringsAsFactors=FALSE )
        df[ , 1 ] <- as.numeric( df[ , 1 ] )
        df[ , 2 ] <- as.numeric( df[ , 2 ] )
        names( df ) <- c( 'Start', 'End' )

        dfALL  <-  rbind( dfALL,
                          data.frame( Chr     = as.character( dt$X.chromosome[ii] ),
                                      CDS_ID = as.character( dt$ccds_id[ii] ),
                                      Gene    = as.character( dt$gene[ii] ),
                                      df,
                                      stringsAsFactors = FALSE )
                         )
    }

    dfALL$Len  <-  dfALL$End  -  dfALL$Start + 1

    return( dfALL )
}


### test ###

if ( 0 ) { # change 1 to 0 for silence

opt <- list( spec = 'exp', 
             a    = 1, 
             b    = Inf, 
             rate = 10
            )

rn <- .get_random_number( opt )

print(rn)

}


