
#' Copy the test data of this package under the current directory. 
#'
#' @param subdir Directory where the test data are stored. By default, 'extdata'
#'
#' @export
#'
#' @examples
#' copy_test_data()
copy_test_data <- function( subdir = 'extdata' ){

   path.dir  <- system.file( subdir, package = 'tugHall.3', mustWork = T )
   path.glob <- paste0( path.dir, '/*' )
   dirs      <- Sys.glob( path.glob )
   
   for ( dir in dirs ) {
      file.copy( dir, '.', recursive = T, overwrite = F )
   }

}


#' Environment of the package 'tugHall.3' to store all the objects of a simulation
#'
#' @description \code{pck.env } is environment of the package 'tugHall.3'
#' where all the objects of a simulation are stored and used
#'
#' @export
#'
pck.env = new.env( parent = emptyenv() )
attr( pck.env, "name" ) = "tugHall.Environment"


#' Remove all the objects from the tugHall environment pck.env
#'
#' @keywords internal
clear_tugHall.Environment  <-  function( ) {

    remove( list = ls( envir = pck.env ), envir = pck.env, inherits = FALSE )

    return( NULL )
}


#' Get results of simulation stored in the tugHall environment pck.env
#'
#' @keywords internal
get_tugHall.Environment  <-  function(){

    results  =  list()
    l = ls( pck.env )
    if ( length( l ) > 0 ){
        for( i in 1:length( l ) ){
            results[[ l[i] ]]  =  pck.env[[ l[i] ]]
        }
    }

    return( results )
}


#' Load previous results of simulation to the environment pck.env
#'
#' @param results List of results of a simulation
#'
#' @keywords internal
load_tugHall.Environment  <-  function( results ){

    l = ls( results )
    if ( length( l ) > 0 ){
        for( i in 1:length( l ) ){
            pck.env[[ l[i] ]]  =  results[[ l[i] ]]
        }
    }

    return( NULL )
}


#########################################################
#
read_files4sim <- function( in.file, check_output = T ) {

   input.table <- read.table( in.file, header=F, sep="\t", stringsAsFactors = FALSE )

   files <- list()
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]

      # add quotes
      V2 <- paste("\'", V2, "\'", sep='')

      text <- paste( V1, ' <- ', V2, sep='' )
      text <- paste0( 'files$', text )
      eval( parse( text = text ) )
   }
   
   if ( check_output ) {
        .check_previous_output( files$Output_dir )
        .make_dirs( files )
   }
   
   return( files )
}


# 
.check_previous_output  <-  function( out ){

    stmp  =  str_replace_all( as.character( Sys.time() ), ' ', '_')
    stmp  =  str_replace_all( stmp, '-', '_' )
    stmp  =  str_replace_all( stmp, ':', '_' )

    if ( dir.exists( out ) ){

        check_files = list.files( out )
        if ( length( check_files ) > 0 ){
            file.rename( from = out, to = paste0( out, '_', stmp ) )
            dir.create( out )
            flog.info( paste0( 'Output folder was renamed with time stamp.' ) )
        }

    }

    return( NULL )
}


#
.make_dirs  <-  function( files ){

    sbdrs  =  c( files$Input_dir,
                 files$Output_dir,

                 files$OutInfo_dir,
                 files$OutVAF_dir,
                 files$OutMut_dir )
    for ( sbdr in sbdrs ){
        if ( ! file.exists( sbdr ) ){
            dir.create( sbdr )
        }
    }

}


### Define the PARAMETERS ------------------------------------------------

#' Define all the parameters for a simulation
#'
#' @return Values of all the parameters
#' 
#' @keywords internal
read_parameters  <-  function( file_prm, file_fixed ){

   prms  <- .read_prms( file_fixed )
   prms2 <- .read_prms( file_prm )
   # if present in both, overridden by file_prm
   prms[ names(prms2) ] <- prms2[ names(prms2) ] 

   if ( is.na( prms$kN ) ) {
      sgmd  =  1 / ( 1 + exp( -prms$sN * ( 0 - 0.5 ) ) )
      if ( prms$tumbler_for_apoptosis_trial ){
          if (prms$dN > sgmd ) {
              prms$kN  =  prms$dN - sgmd
          } else {
              prms$kN  =  prms$dN
          }
      } else {
          prms$kN  =  prms$dN
      }
      if ( prms$kN < 0 ) prms$kN = 0
   } 

   if ( prms$uo_del != 0 ) warning( 'uo_del is not zero.' )
   if ( prms$us_dup != 0 ) warning( 'us_dup is not zero.' )

   return( prms )
}


#
.read_prms <- function( file ) {

   input.table <- read.table( file, header=T, sep="\t", stringsAsFactors=F )
   if (  ! all( colnames(input.table) ==  c('Parameter', 'Value') ) )
      stop( 'Column names in the parameter file should be Parameter and Value')
   
   prms <- list()
   for( ii in 1:nrow( input.table ) ) {
      V1 <- input.table[ ii, 1 ]
      V2 <- input.table[ ii, 2 ]
      
      # add quotes
      if ( V1 == 'metastatic_model' || 
           V1 == 'growth'           || 
           length( grep( "\\$spec$", V1, perl=T ) ) > 0
          ) V2 <- paste("\'", V2, "\'", sep='')
   
      text <- paste( 'prms$', V1, ' <- ', V2, sep='' )
      eval( parse( text = text ) )
   }
   
   return( prms )
}


# 
define_parameters  <-  function( prms ){

   # The list of parameters:
   vrs.number  <- c( 'm0', 'uo', 'us', 
                     'dN', 'kN','sN', 'K_N', 'Fb', 'ctmax', 

                     'm_dup',            'm_del', 
                     'uo_dup', 'us_dup', 'uo_del', 'us_del', 
                     'ave_len_dup',      'ave_len_del', 
                    
                     'Zim', 'kappa', 
                    
                     'control.censor_cell_number', 
                     'control.censor_time_step', 
                     'control.censor_real_time' )
   vrs.logical <- c( 'compaction_factor', 
                     'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
                     'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial',
                     'tumbler_for_event_enforcement', 
                     'control.monitor' )
   vrs.char    <- c( ) 

   vrs = c( vrs.number, vrs.logical, vrs.char )
   for( vr in vrs ){
      if ( is.null( prms[[ vr ]] ) || 
           is.na(   prms[[ vr ]] ) || 
           length(  prms[[ vr ]] ) == 0 ) {
          stop( paste0( 'The parameter ', vr, ' is not defined. '))
      }
   }
   
   sapply( X=vrs.number,  function(X) pck.env[[ X ]] <- as.numeric(  prms[ X ]) )
   sapply( X=vrs.logical, function(X) pck.env[[ X ]] <- as.logical(  prms[ X ]) )
   sapply( X=vrs.char,    function(X) pck.env[[ X ]] <- as.character(prms[ X ]) )
}


#' Define compaction factor
#'
#' @return Data frame with with compaction factors for all the hallmarks
#'
#' @keywords internal
define_compaction_factor  <-  function( file_name ){

    data_log  =  read.table( file = file_name, sep = '\t', header = TRUE, stringsAsFactors = FALSE )
    if ( ! all( names( data_log ) == c( 'Hallmark', 'CompactionFactor' ) ) )
       stop( 'Column names in compaction factor file should be Hallmark and CompactionFactor' )

    cf = data.frame( Ha = 1, Hd = 1, Hi = 1, Him = 1 )
    cf$Ha   =  data_log$CompactionFactor[ data_log$Hallmark == 'apoptosis' ]
    cf$Hd   =  data_log$CompactionFactor[ data_log$Hallmark == 'growth' ]
    cf$Hi   =  data_log$CompactionFactor[ data_log$Hallmark == 'immortalization' ]
    cf$Him  =  data_log$CompactionFactor[ data_log$Hallmark == 'invasion' ]

    pck.env$CF  =  cf
}


#
define_regions <- function( file_Rint, file_Rother ) {

   data_int = read.table( file = file_Rint, sep = '\t', stringsAsFactors = FALSE, header = TRUE )
   if (  ! all( names( data_int ) ==  c( 'Gene', 'CDS_ID', 'Type' ) ) )
      stop( 'Column names in Rint file should be Gene, CDS_ID, Type. ')

   data_other  =  read.table( file = file_Rother, sep = '\t', stringsAsFactors = FALSE, header = TRUE )
   if (  ! all( names( data_other ) ==  c( 'Gene', 'CDS_ID' ) ) )
      stop( 'Column names in Rother file should be Gene, CDS_ID. ')

   if ( any( data_int$Gene %in% data_other$Gene ) )
      stop( 'Genes in Rint and Rother files should be different.' )

   pck.env$ls_genes  =  c( data_int$CDS_ID, data_other$CDS_ID )

}


#' Define genes' location in chromosome
#' 
#' @param genes_list is a list of genes' names like CCDS4107.1 in the CCDS database.
#' @param file_database is a name of file to input where the information about genes location is defined. It can be loaded from CCDS database https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/
#' @param file_selected is a name of file to output where the information about genes location is saved.
#' @param save.out is a logical parameter to save the output to the file
#'
#' @return table of genes' locations in DNA
#' @keywords internal
define_gene_location  <-  function( genes_list, file_database, file_selected, save.out = F ) {

    gene_map = make_map( ls   = genes_list,
                         f_in = file_database )

    ### We have to be sure in the sorting of start position for each chromosome
    gene_map  = order_gene_map( gene_map )
    
    pck.env$gene_map = gene_map 
    if ( save.out ){
        write.table( pck.env$gene_map, file = file_selected, sep = "\t",
                     col.names = TRUE, row.names = FALSE )
    }
}


#' Function to define information about all the genes regarding probabilities of \code{u} for deletion, duplication and point mutation
#'
#' @keywords internal
define_dominant_recessive  <-  function( gfile, name, onsp ){
    if ( is.null(pck.env$onco) ) 
       stop( 'For define_dominant_recessive function onco object should be define first.' )
    if ( is.null(pck.env$hall) ) 
       stop( 'For define_dominant_recessive function hall object should be define first.' )
    
    gtype <- read.table( gfile, header = T, sep = '\t', 
                         stringsAsFactors = F)
    
    DR  =  data.frame( gene  =  name,
                       dr    =  onsp,
                       pom   =  0,
                       del   =  0,
                       dup   =  0,
                       role_pom   =  0,
                       role_del   =  0,
                       role_dup   =  0,
                       stringsAsFactors = FALSE )

    for ( os in c( 'o', 's', 'dn' ) ) {
       w  =  which( DR$dr == os )
       DR$pom[ w ]      = pck.env[[ gtype[ gtype$type==os, 'u.pom' ] ]]
       DR$del[ w ]      = pck.env[[ gtype[ gtype$type==os, 'u.del' ] ]]
       DR$dup[ w ]      = pck.env[[ gtype[ gtype$type==os, 'u.dup' ] ]]
       DR$role_pom[ w ] =           gtype[ gtype$type==os, 'mode.pom' ]
       DR$role_del[ w ] =           gtype[ gtype$type==os, 'mode.del' ]
       DR$role_dup[ w ] =           gtype[ gtype$type==os, 'mode.dup' ]
    }

    return( DR )
}
