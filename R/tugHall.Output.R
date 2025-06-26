# Functions related to output of a simulation -----------------------------------------

#' Function to write parameters
#'
#' @keywords internal
write_parameters <- function( file, prms, prms2 = NULL ) {
   
   data <- data.frame( 'Parameter' = 'Working_folder', 'Value' = getwd(), 
                       stringsAsFactors = F )
   for ( name in names(prms) ) {
      if ( length(    prms[[ name ]] ) != 1 ) next
      if ( is.atomic( prms[[ name ]] ) != T ) next

      data <- rbind( data, 
                     data.frame( 'Parameter' = name, 'Value' = prms[[ name ]],
                                 stringsAsFactors = F  ) 
                    )
   }
   for ( name in setdiff( names(prms2), names(prms) ) ) {
      if ( ! is.atomic( prms2[[ name ]] ) ) next
      data <- rbind( data, 
                     data.frame( 'Parameter' = name, 'Value' = prms2[[ name ]],
                                 stringsAsFactors = F  ) 
                    )
   }
   
   data <- data[ order( nchar( data$Parameter ), 
                               data$Parameter ), ]

   dir.create( dirname( file ), showWarnings = F, recursive = T )
   write.table( data, file, sep = "\t", col.names = F, quote = F, row.names = F)
}


#' Function to write info about Hallmark data
#'
#' @param outfile File name for output info
#' @param hall Object of class "Hallmark"
#' @param compaction_factor Compaction factor, logical type only. True means 'to use', False means 'do not use'.
#' @param CF Vector with values of compaction factor for each hallmark
#'
#' @keywords internal
write_geneout <- function( outfile, hall, compaction_factor, CF ) {

    hallmark <- NULL
    gene <- NULL
    weight_used <- NULL
    weight_woCF <- NULL
    
    if ( pck.env$tumbler_for_apoptosis_trial )       {
        hallmark <- c( hallmark, rep("apoptosis",       length(pck.env$onco$name[ hall$Ha  ])) )
        gene     <- c( gene,     pck.env$onco$name[ hall$Ha   ])
        weight_used <- c( weight_used, hall$Ha_w )
        if ( compaction_factor ) weight_woCF <- c( weight_woCF, hall$Ha_w  / CF$Ha )
    }
    if ( pck.env$tumbler_for_immortalization_trial ) {
        hallmark <- c( hallmark, rep("immortalization", length(pck.env$onco$name[ hall$Hi  ])) )
        gene     <- c( gene,     pck.env$onco$name[ hall$Hi   ])
        weight_used <- c( weight_used, hall$Hi_w )
        if ( compaction_factor ) weight_woCF <- c( weight_woCF, hall$Hi_w  / CF$Hi )
    }
    
    hallmark <- c( hallmark,  rep("growth",          length(pck.env$onco$name[ hall$Hd  ])) )
    gene <- c( gene, pck.env$onco$name[ hall$Hd  ] )
    weight_used <- c( weight_used, hall$Hd_w )
    if ( compaction_factor ) weight_woCF <- c( weight_woCF, hall$Hd_w  / CF$Hd )
    
    if ( pck.env$tumbler_for_angiogenesis_trial )   {
        hallmark <- c( hallmark, rep("angiogenesis",    length(pck.env$onco$name[ hall$Hb  ])) )
        gene     <- c( gene,     pck.env$onco$name[ hall$Hb   ])
        weight_used <- c( weight_used, hall$Hb_w )
        if ( compaction_factor ) weight_woCF <- c( weight_woCF, hall$Hb_w )
    }
    if ( pck.env$tumbler_for_metastasis_trial )   {
        hallmark <- c( hallmark, rep("invasion",        length(pck.env$onco$name[ hall$Him ])) )
        gene     <- c( gene,     pck.env$onco$name[ hall$Him  ])
        weight_used <- c( weight_used, hall$Him_w )
        if ( compaction_factor ) weight_woCF <- c( weight_woCF, hall$Him_w / CF$Him )
    } 
    
    if ( !compaction_factor ) weight_woCF = weight_used 

    
    df <- data.frame( hallmark = hallmark, gene=gene, 
                      weight_woCF = weight_woCF,
                      weight_used = weight_used )
    dir.create( dirname( outfile ), showWarnings = F, recursive = T )
    write.table( df, outfile, quote = F, row.names = FALSE, sep = '\t' )
}


#' Function to write data to cloneout file at a time step
#'
#' @param outfile File name for output info
#' @param env Object of class 'Environ'
#' @param clones List of objects of class 'Clone'
#' @param onco_clones List of objects of class 'OncoGene'
#'
#' @keywords internal
write_cloneout <- function( outfile, env, clones, onco_clones ) {
    debugging <- FALSE

    dir.create( dirname( outfile ), showWarnings = F, recursive = T )
    if ( !file.exists( outfile )) .write_header( outfile, pck.env$env, onco_clones )

    intact_normal  =  sum( unlist( sapply( clones, FUN = function( cl ) 
                                              ifelse( cl$CNA_ID[ 1 ] == 0 & cl$PointMut_ID[ 1 ] == 0, cl$N_cells, 0 ) ) ) )

    common1 <- c( env$T, 'avg', '-',  '-', '-', '-', env$c, env$d, env$i, env$im, env$a, env$k, env$E,
                  intact_normal, env$N - intact_normal,  
                  env$Nmax, env$P, env$M, env$Ha, env$Him, env$Hi, env$Hd, env$Hb, '-', env$mutden,
                  env$n_divisions 
                 )
    if ( debugging == F ){
        data  =  c( common1, rep( '-', 4                                ) )
    } else {
        data  =  c( common1, rep( '-', 10 + 4*length(pck.env$onco$name) ) )
    }
    columns_structure  =  c( 1:3, 5, 6, 24, 4, 17, 18, 14, 15, 7, 8, 12, 11, 9, 10, 13, 16, 22, 19, 21, 20, 23, 25:length( data ) )
    write( data[ columns_structure ], outfile, append=TRUE, ncolumns = length(data), sep="\t")

    if (length(clones) > 0 ) {
        for (i in 1:length( clones ) ) {
            clone1  =  clones[[i]]
            onco1   =  onco_clones[[i]]
            type    =  get_type( clone1 = clone1 )
            common2 <- c( env$T, i, clone1$id, clone1$N_cells, clone1$parent, clone1$birthday, clone1$c, clone1$d,
                          clone1$i, clone1$im, clone1$a, clone1$k, clone1$E,
                          intact_normal, env$N - intact_normal, 
                          clone1$Nmax, env$P, env$M,
                          clone1$Ha, clone1$Him, clone1$Hi, clone1$Hd, clone1$Hb, type, 
                          clone1$mutden,
                          env$n_divisions,
                          paste( pck.env$onco$name[ which( clone1$gene    == 1 ) ], collapse = ', ' ), 
                          paste( pck.env$onco$name[ which( clone1$pasgene == 1 ) ], collapse = ', ' ), 
                          paste(clone1$PointMut_ID, collapse = ', '),
                          paste(clone1$CNA_ID,      collapse = ', ')
                         ) 
            if ( debugging == F ){
                data    =  c( common2 )
            } else {
                data    =  c( common2, onco1$cds_1, onco1$rna_1, onco1$prob_1,
                                       onco1$cds_2, onco1$rna_2, onco1$prob_2 )
            }
            write( data[ columns_structure ], outfile, append=TRUE, ncolumns = length(data), sep="\t")
        }
    }
}


#' Function to write the header to a file
#'
#' @param outfile File name for output info
#' @param env Object of class 'Environ'
#' @param onco Object of class "OncoGene"
#'
#' @keywords internal
.write_header <- function(outfile, env, onco) {
    debugging <- F
    
    common <- c('Time', 'AvgOrIndx', 'ID', 'N_cells', 'ParentID', 'Birth_time', 'ct', 'd', 'i', 'im', 'a',
                'k', 'K', 'N_normal_intact',	'N_normal_speckled', 'Nmax', 'N_primary', 'N_metastatic',
                'Ha', 'Him', 'Hi', 'Hd', 'Hb', 'Type', 'mutden',
                'total_divIdx',
                'genesDysfunc', 'genesWmutsNoDys',
                'PointMut_ID', 'CNA_ID' )
    if ( debugging == F ){
        header <- c( common )
    } else {
        header <- c( common,
                     paste('Pchr1_CDS_', onco$name, sep=''), paste('Pchr1_RNA_', onco$name, sep=''),
                     'Pchr1_prob_point_mut', 'Pchr1_prob_del', 'Pchr1_prob_dup' ,
                     paste('Pchr2_CDS_', onco$name, sep=''), paste('Pchr2_RNA_', onco$name, sep=''),
                     'Pchr2_prob_point_mut', 'Pchr2_prob_del', 'Pchr2_prob_dup' )
    }
    columns_structure  =  c( 1:3, 5, 6, 24, 4, 17, 18, 14, 15, 7, 8, 12, 11, 9, 10, 13, 16, 22, 19, 21, 20, 23, 25:length( header ) )
    write( header[ columns_structure ], outfile, append=FALSE, ncolumns = length(header), sep="\t")
}



#' Function to write the the reason of stop of the simulations into a file
#' @keywords internal
.write.log  <- function( fl, l.clones, time.step, cells_number, df.time , censors ){
    
    stop.reason = 'Undefined'
    if ( l.clones < 1 )                                       stop.reason = 'Extinction'
    if ( time.step >= censors$control.censor_time_step )      stop.reason = 'censor_time_step'
    if ( cells_number > censors$control.censor_cell_number )  stop.reason = 'censor_cell_number'
    if ( df.time > censors$control.censor_real_time )         stop.reason = 'censor_real_time'
    
    if ( !file.exists( fl ) ) {
        write.table( paste('exit_reason', stop.reason, sep = "\t" ), fl, sep = "\t", col.names = F, row.names = F, quote = F )
    } else {
        write.table( paste('exit_reason', stop.reason, sep = "\t" ), fl, sep = "\t", col.names = F, row.names = F, quote = F, append = T )
    }
}

#
write_genetic_type <- function( file, dr ) {
   
   dir.create( dirname( file ), showWarnings = F, recursive = T )
   write.table( dr, file, sep = "\t", 
                col.names = c( 'gene', 'type', 'pom', 'del', 'dup', 
                               'mode.pom', 'mode.del', 'mode.dup' ),
                row.names = F, quote = F )

}


#' Function to write the point mutation info for all clones for all time steps, used at the last time step or after simulation
#'
#' @param pnt_clones List of objects of class 'Point_Mutations'
#' @param file_out  File name to write pnt_clones for allele B
#' @param file_outA File name to write pnt_clones for allele A
#' @param safe_to_file Logical parameter to safe data frame to the file or do not safe. By default it is TRUEile
#'
#' @keywords internal
write_pnt_clones <- function( pnt_clones, file_out, file_outA, safe_to_file = TRUE ){

    dir.create( dirname( file_out ), showWarnings = F, recursive = T )
    pn  <-  NULL
    if ( !is.null( pnt_clones ) ){
        for (i in 1:length(pnt_clones)) {
            pnt1 <-  unlist( pnt_clones[[i]] )
            pn1  <-  safe_pnt_mut( pnt1 )     ### pnt1$safe()
            pn   <-  rbind( pn, pn1)
        }
    }

    if ( safe_to_file ){
        # Save allele B
        w  =  which( !is.na( pn$MalfunctionedByPointMut) )
        w_name  =  which( !( 'mut_order' == names( pn ) ) )
        write.table( pn[ w, w_name ], file = file_out, append=FALSE, sep="\t", row.names = FALSE)

        # Save allele A
        w  =  which( is.na( pn$MalfunctionedByPointMut) )
        w_name  =  which( !( 'mut_order' == names( pn ) ) )
        write.table( pn[ w, w_name ], file = file_outA, append=FALSE, sep="\t", row.names = FALSE)
        
        
    }
    return( pn )
}


#' Function to write the CNA mutation info for all clones for all time steps, used at the last time step or after simulation
#'
#' @param cna_clones List of objects of class 'CNA_Mutations'
#' @param file_out File name to write
#'
#' @keywords internal
write_cna_clones <- function( cna_clones, file_out ){

    dir.create( dirname( file_out ), showWarnings = F, recursive = T )
    cn  <-  NULL
    if ( !is.null( cna_clones ) ){
        for (i in 1:length(cna_clones)) {
            cna1 <-  unlist( cna_clones[[ i ]] )
            cn1  <-  safe_pnt_mut( cna1 )     ### pnt1$safe()
            cn   <-  rbind( cn, cn1)
        }
    }

    w_name  =  which( !( 'mut_order' == names( cn ) ) )
    write.table( cn[ , w_name ], file = file_out, append=FALSE, sep="\t", row.names = FALSE)
}


#' Function to get type of the clone: normal, primary or metastatic
#'
#' @return One of characters 'normal', 'primary' or 'metastatic'
#'
#' @keywords internal
get_type  <-  function( clone1 ){
    if ( clone1$invasion ) return( 'metastatic' )
    if ( !clone1$primary ) return( 'normal' )
    return( 'primary' )
}


##### monitors #########################

#' Function to write a simulation monitoring data
#'
#' @keywords internal
monitor_Ncells.Nmuts <- function( outfile, start = FALSE, env, clones, cna_clones, tumber_EEL ){

    if ( start ) {
        header <- c( 'Time', 'N_clones', 'N_normal_intact', 
                     'N_normal_speckled', 'N_primary', 'N_metastatic',
                     'N_point_mutations', 'N_duplications', 'N_deletions' )
        if ( tumber_EEL ) header  =  c( header, 'N_divisions' )
        
        write( header, outfile, append = FALSE, ncolumns = length( header ), sep="\t" )
    } else {
        if ( length( clones )  >  0 ) {
            point_mut_list  =  sort( unique( as.numeric( unlist( 
               sapply( X = 1:length( clones ), FUN = function( x ) 
                  clones[[ x ]]$PointMut_ID ) ) ) ) )
            
            if ( point_mut_list[ 1 ] == 0 ) l_pm  =  length( point_mut_list ) - 1 
            else                            l_pm  =  length( point_mut_list )
            
            cna_list  =  sort( unique( as.numeric( unlist( 
               sapply( X = 1:length( clones ), FUN = function( x ) 
                  clones[[ x ]]$CNA_ID ) ) ) ) )
            if ( cna_list[ 1 ] == 0 ) cna_list  =  cna_list[ -1 ]

            dupdel  =  unlist( sapply( X = cna_list, FUN = function( x ) { 
                                           cna_clones[[ x ]]$dupOrdel } ) )

            l_dup   =  length( which( dupdel  ==  'dup' ) )
            l_del   =  length( which( dupdel  ==  'del' ) )

            # Get intact and speckled normal cells:
            i_n  =  which( unlist( sapply( 1:length(clones), FUN = function(x) 
               get_type( clones[[ x ]] ) ) ) == 'normal')
               
            if ( length( i_n ) > 0 ){
                int = unlist( sapply( i_n, FUN = function( x ) sum( clones[[ x ]]$pasgene ) ) )
                if ( length( which( int == 0 ) ) > 0 ){
                    N_intact  =  unlist( sum( sapply( which( int == 0 ), FUN = function( x ) 
                                         clones[[ i_n[ x ] ]]$N_cells ) ) )
                } else { 
                    N_intact  =  0
                }
                N_speckled  =  env$N - N_intact
            } else {
                N_intact    =  0
                N_speckled  =  0
            }

            data <- c( env$T, length( clones ), N_intact, N_speckled, 
                       env$P, env$M, l_pm, l_dup, l_del )
            if ( tumber_EEL ) data  =  c( data, env$n_divisions )

            write(data, outfile, append=TRUE, ncolumns = length(data), sep="\t")
        }


    }

}


#
monitor_pck.env <- function( file, first, time, prev.obj = NULL, obj ) {

   obj  <- sort( obj )
   obj2 <- obj[ ! ( obj %in% prev.obj ) ]

   if ( length( obj2 ) == 0 ) { return( obj ) }
   tbl <- data.frame( 'time'    = time,
                      'new_pck.env' = obj2 )

   col.flg <- ifelse( length( obj ) != length( obj2 ) || first == F,
                      F, T )
   write.table( tbl, file, append = !col.flg, col.names = col.flg,
                quote = F, sep = "\t", row.names = F,  )

   return( obj )
}


