

#' Function to make a gene_map data.frame with information of genes' locations
#' @keywords internal
make_map  <-  function( ls   =  c( 'CCDS4107.1',  'CCDS8702.1',
                                   'CCDS43171.1', 'CCDS11118.1' ),
                        f_in =  'Input/CCDS.current.txt'
                       ) {

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


#' Function to order info in gene_map data.frame with information of genes' locations
#' @keywords internal
order_gene_map  <-  function( gene_map ){
    gm = NULL
    for( t in   unique(gene_map$Chr) ) {

        w = which( gene_map$Chr == t )
        gm2  =  gene_map[ w, ]
        gmt  =  gm2[ order( gm2$Start ), ]
        gm   =  rbind( gm, gmt )

    }
    rownames(gm) <- 1:length(gm$Chr)
    return( gm )
}

#' Function to get length of CDS and whole gene from gene_map data.frame
#' @keywords internal
get_len_cds_rna  <-  function( gene_map ){

    name0  =  unique( gene_map$Gene )
    cds0   =  NULL
    rna0   =  NULL
    # Get length of CDS and RNA from gene_map:
    for ( i in 1:length(name0) ) {
        w     =  which( gene_map$Gene == name0[ i ] )
        cds0  =  c( cds0, sum( gene_map$End[w]  -  gene_map$Start[w] + 1 ) )
        rna0  =  c( rna0, max( gene_map$End[w] ) - min( gene_map$Start[w]) + 1  )
    }

    if ( !is.null( pck.env$onco ) ){
        w = unlist( sapply( name0, FUN = function(x) which( pck.env$onco$name  ==  x ) ) )
    
        name0  =  name0[w]
        rna0   =  rna0[w]
        cds0   =  cds0[w]
    }

    return( list( Name = name0, CDS = cds0, LEN_Genes = rna0 ) )
}


#' Function to read a file
#' @keywords internal
read_file  <-  function( file_name = '', stringsAsFactors = FALSE, header = TRUE ){
    if ( !file.exists( file_name ) ) {
        warning( paste0('The file ', file_name, ' does not exist. ' ) )
        return( NULL )
    }
    if ( file.size( file_name )  < 10 ) return( NULL )
    return( read.table( file = file_name, stringsAsFactors  =  stringsAsFactors ,
                        sep="\t", header = header ))
}


#' Check the installation of a package for some functions
#' @keywords internal
.check_pkg  <-  function( pkg ){
    msg  =  paste0( 'Package ', pkg, ' must be installed to use this function. \n ' )
    if ( !requireNamespace( pkg , quietly = TRUE ) )    stop( msg, call. = FALSE )
}


#' Check the installation of packages and attach them with corresponding functions
#' @keywords internal
check_packages  <-  function( pkgs = NULL ){

   if ( is.null( pkgs ) ) {
      pkgs = list( dplyr  = '%>%',
                   methods = c( 'new', 'setRefClass', 'initialize', 'is' ), 
                   stats = c( 'rbinom', 'rexp', 'rnorm', 'runif' ),
                   stringr = c( 'str_split', 'str_replace_all' ),
                   utils = c( 'read.delim', 'read.table', 'write.table' ),
                   withr = c( 'local_environment' ),
                   truncdist = c( 'rtrunc' ),
                   dirmult = c( 'rdirichlet' ),
                   MASS = c( 'mvrnorm' ), 
                   futile.logger = c( 'flog.debug',
                                      'flog.info',
                                      'flog.warn',
                                      'flog.error',
                                      'flog.appender',
                                      'flog.threshold'
                                     )
                  )
   }

    ### Attach the packages
    for( pck in names( pkgs ) ){
        .check_pkg( pkg = pck )
        require( package = pck, character.only = TRUE, include.only = pkgs[[ pck ]])
    }
}


#' Foolproof function allows to checking the consistency of some input parameters
#' @keywords internal
foolproof  <-  function(){

    # 1. List of genes in different objects: onco, gene_map
    genes_onco  =  sort( unique( pck.env$onco$name) )
    genes_map   =  sort( unique( pck.env$gene_map$Gene ) )

    if ( length( genes_onco ) == length( genes_map ) ){ 
        if ( !all.equal( sort( genes_onco ), sort( genes_map ) ) ){
            warning( paste( c( 'Names of genes defined for onco: ', genes_onco ), collapse = '  ' ) )
            warning( paste( c( 'Names of genes defined for genes location/gene_map: ', genes_map ), collapse = '  ' ) )
            stop( 'Names of genes defined in onco and gene_map are different.' )
        }

    } else {
        message("Please, be sure that all definitions of hallmarks and genes names are correct!")
        message( paste( c( 'Names of genes defined for onco object: ', genes_onco ), collapse = '  ' ) )
        message( paste( c( 'Names of genes defined for genes location/gene_map: ', genes_map ), collapse = '  ' ) )

    }

    # 2. NA and NULL check in onco and hall:

    ### Check onco object:
    flds  =  c( "name", "cds_1", "cds_2", "rna_1", "rna_2", "p0_1", "p0_2", "prob_1", "prob_2",
                "sum_prob_1", "sum_prob_2", "onsp", "len" )

    for( v in flds ){
        if ( length(  pck.env$onco[[ v ]] ) == 0 |
             any( is.na(   pck.env$onco[[ v ]] ) )   |
             is.null( pck.env$onco[[ v ]] ) ){
            stop( paste0( 'The parameter ', v, ' in the onco object is not defined.'))
        }
    }

    if ( !all( pck.env$onco$onsp %in% c("dn", "o",  "s" ) ) ){
        stop( 'Not all the genes are onco (o) or suppressor (s) or dominant negative (dn).' )
    }

    # Check the definition of compaction factors if they are used:
    if ( pck.env$compaction_factor ){
        if ( is.null( pck.env$CF ) ) stop( 'Compaction factors are NOT defined yet. They are needed to define hallmarks weights.' )
        if ( is.null( pck.env$CF$Ha  ) ) stop( 'Compaction factor for apoptosis is not defined. ' )
        if ( is.null( pck.env$CF$Hd  ) ) stop( 'Compaction factor for division is not defined. ' )
        if ( is.null( pck.env$CF$Hi  ) ) stop( 'Compaction factor for immortalization is not defined. ' )
        if ( is.null( pck.env$CF$Him ) ) stop( 'Compaction factor for invasion/metastasis transformation is not defined. ' )
    }

    ### Check hall object:
    flds  =  c( "Hd", "Hd_w" )
    if ( pck.env$tumbler_for_apoptosis_trial )       flds  =  c( flds, "Ha",  "Ha_w"  )
    if ( pck.env$tumbler_for_angiogenesis_trial )    flds  =  c( flds, "Hb",  "Hb_w"  )
    if ( pck.env$tumbler_for_immortalization_trial ) flds  =  c( flds, "Hi",  "Hi_w"  )
    if ( pck.env$tumbler_for_metastasis_trial )      flds  =  c( flds, "Him", "Him_w" )

    for( v in flds ){
        if ( length(  pck.env$hall[[ v ]] ) == 0 |
             any( is.na(   pck.env$hall[[ v ]] ) )    |
             is.null( pck.env$hall[[ v ]] ) ){
            message( paste( c ( 'In hallmarks definition you have parameter ', v, ' = ', pck.env$hall[[ v ]]  ), collapse = ' ' ) )
            stop( paste0( 'The parameter ', v, ' in the hall object is not defined.'))
        }
    }

    if ( pck.env$Fb == 0 ) {
        stop( 'The factor Fb for angiogenesis hallmark should be non-zero.' )
    }

    return( NULL )
}


#' Function to calculate binomial distribution including BIG NUMBERS like 10^12 and more using approximation with normal distribution
#'
#' @param tr Length of vector with successes trials
#' @param n Number of independent Bernoulli trials
#' @param p Probability to get successes in trials
#'
#' @return Vector of integer numbers of successes trials
#'
#' @keywords internal
calc_binom <- function(tr,n,p){
    if (n*p < 10^8){
        ou <- rbinom(tr,n,p)
    } else {
        m <- n * p
        s <- sqrt(  n*p* (1-p)  )
        ou <- rnorm( tr, mean = m, sd = s )
    }

    return(  round( ou )  )
}




