


# FUNCTIONS ---------------------------------------------------------------



#' Function to get data about simulation from cloneoutfile
#'
#' @return list of data.frames
#' 
#' @keywords internal
get_flow_data <- function( cloneoutfile, pom.file, pomA.file, cna.file, time = NA ) {

    # Get data about onco and hallmarks

    ## Get data from output file
    data_out  <-  read.table( cloneoutfile, sep="\t", header = TRUE, stringsAsFactors = FALSE )
    data_out[is.na(data_out)]  <-  ""

    # average data
    data_avg  =  data_out[ which( data_out$AvgOrIndx  ==  "avg"), ]
    data_avg$N_normal  =  data_avg$N_normal_intact  +  data_avg$N_normal_speckled

    # data without averaging - flow data
    data_flow  =  data_out[ which( !data_out$AvgOrIndx == "avg" ), ]

    nms  =  names( data_flow )
    nms_exclude  =  c( 'Type', 'genesMalfunc', 'genesWmutsNoExp', 'PointMut_ID', 'CNA_ID' )
    w_nms  =  as.integer( sapply( X = nms_exclude, FUN = function( x ) which( nms == x ) ) )

    clmns  =  ( 1:ncol( data_flow ) )[-w_nms]  #  numeric columns
    data_flow[ , clmns ]  =  sapply( X = clmns,  FUN = function( x ) {
                                                    as.numeric( data_flow[ , x ] )
                                                  } )
    data_flow$N_normal  =  data_flow$N_normal_intact  +  data_flow$N_normal_speckled

    # the data at the defined time step
    if ( is.na(time) ) {
        time = max( data_flow$Time )
    }
    else if ( time >= max( data_flow$Time ) ) {
        return( NA )
    }
    else if ( time < min( data_flow$Time ) ) {
        return( NA )
    }
    data_last  =  data_flow[ which( data_flow$Time  ==  time ), ]
    rm( data_out )
    if ( file.size( cna.file ) < 5 ){
        cna_mut  =  NULL
    } else {
        cna_mut  =  read.table( file = cna.file, sep = '\t', header = TRUE, stringsAsFactors = FALSE )
        }
    if ( file.size( pom.file ) < 5 ){
        pnt_mut  =  NULL
    } else {
        pnt_mut  =  read.table( file = pom.file, sep = '\t', header = TRUE, stringsAsFactors = FALSE )
    }
    
    if ( file.size( pomA.file ) < 5 ){
        pnt_mut_A  =  NULL
    } else {
        pnt_mut_A  =  read.table( file = pomA.file, sep = '\t', header = TRUE, stringsAsFactors = FALSE )
    }

    return( list( data_avg = data_avg,
                  data_flow = data_flow, time_max = time, data_last = data_last,
                  cna_mut = cna_mut, 
                  pnt_mut = pnt_mut, 
                  pnt_mut_A = pnt_mut_A
                  ) )
}


#' Function to calculate the Variant allele frequencies (VAF)
#'
#' @param pnt_mut data.frame with point mutation info for allele B
#' @param pnt_mut_A data.frame with point mutation info for allele A
#' @param data_last data.frame with data of simulation at the defined time step
#'
#' @return data.frame with info about Variant allele frequencies
#' 
#' @keywords internal
get_VAF  <-  function( pnt_mut, pnt_mut_A, data_last ){

    pnt_mut_B  =  pnt_mut

    if ( nrow(data_last) == 0 ) return( NULL )

    ids  =  str_split( data_last$PointMut_ID, pattern = ',' )
    ids  =  lapply( X = 1:nrow( data_last), FUN = function(x) as.numeric( ids[[ x ]] ) )

    nqu  =  sort( unique( unlist( ids ) ) )  # unique IDs of point mutations
    if ( nqu[1] == 0 ) nqu = nqu[ -1 ]       # exclude intact normal cells

    VAF  =  NULL

    if ( length( nqu ) == 0 ) return( NULL )

    for( j in 1:length( nqu ) ){
        # wc - which clones have an ID of point mutation
        wc  =  unlist( sapply( X = 1:length( ids ), FUN = function( x ) is.element( nqu[j] , ids[[ x ]] ) ) )

        VAF_1  =  pnt_mut_B[ which( pnt_mut_B$PointMut_ID == nqu[ j ] ) , ]

        ### number of cells: speckled normal cells,
        ###     primary tumor cells and metastatic cells:
        if ( any( data_last[ which( wc ), 'Type' ] == 'normal' ) ){
            sm  =  sum( as.numeric( data_last[ which( wc & data_last$Type ==  'normal' ),  'N_cells'  ] ) )
            if ( length( sm ) == 0 ) sm = 0
            VAF_1$N_speckled_normal =  sm
        } else VAF_1$N_speckled_normal =  0

        if ( any( data_last[ which( wc ), 'Type' ] == 'primary' ) ){
            VAF_1$N_primary  =  sum( as.numeric( data_last[ which( wc & data_last$Type ==  'primary' ),  'N_cells'  ] ) )
        } else VAF_1$N_primary  =  0

        if ( any( data_last[ which( wc ), 'Type' ] == 'metastatic' ) ){
            VAF_1$N_metastatic =  sum( as.numeric( data_last[ which( wc & data_last$Type ==  'metastatic' ),  'N_cells'  ] ) )
        } else VAF_1$N_metastatic =  0

        # add copy number of original allele A:
        VAF_A  =  pnt_mut_A[ which( pnt_mut_A$PointMut_ID == nqu[ j ] ) , ]

        VAF_1$Copy_number_A  =  VAF_A$Copy_number

        VAF  =  rbind( VAF, VAF_1)
    }

    # Total number of speckled normal cells (no driver mutations but at least one passenger )
    if ( any( data_last[ , 'Type' ] == 'normal' & ( data_last[ , 'CNA_ID'] != '0' | data_last[ , 'PointMut_ID'] != '0' ) ) ) {
      w = which( data_last$Type == 'normal' & ( data_last[ , 'CNA_ID'] != '0' | data_last[ , 'PointMut_ID'] != '0' ) )
      VAF$N_speckled_normal_total = sum( as.numeric( data_last[ w, 'N_cells' ] ) )
    } else  VAF$N_speckled_normal_total  =  0

    # Total number of primary tumor cells (at least one driver)
    if ( any( data_last[ , 'Type' ] == 'primary' ) ) {
        VAF$N_primary_total = sum( as.numeric( data_last[ which( data_last$Type == 'primary' ), 'N_cells' ] ) )
    } else  VAF$N_primary_total  =  0

    # Total number of metastatic cells
    if ( any( data_last[ , 'Type' ] == 'metastatic' ) ) {
        VAF$N_metastatic_total = sum( as.numeric( data_last[ which( data_last$Type == 'metastatic' ), 'N_cells' ] ) )
    } else  VAF$N_metastatic_total  =  0

    VAF  = VAF[ , -which(names(VAF) == 'rst.ratio')]
    ### Save VAF to the file:
    nms  =  names( VAF)[ c( 1:6, 14, 7:8, 11:13, 15:17 ) ]
    VAF_out  =  VAF[ , nms ]
    names( VAF_out )  =  c( nms[ c( 1:7) ], 'Copy_number_B', nms[ c( 9:15 ) ] )
    VAF_out$Chr  =  as.character( VAF_out$Chr )

    return( VAF_out )
}


#' Function to get Variant allele frequencies (VAF) based on rho input parameters
#'
#' @param vf data.frame getting from get_VAF() function
#' @param rho Vector of rho parameter in the range (0,1)
#'
#' @return VAF for different rho with separation for metastatic cells and (primary tumor + speckled normal) cells
#' 
#' @keywords internal
get_rho_VAF  <-  function( vf = NULL, rho = c( 0.0, 0.1, 0.5 ) ){

    if ( min(rho) < 0 | max(rho) >1 ) stop( 'rho values should be in the range [0,1]' )
    nq_i  =  unique( vf$Ref_pos )
    if ( length(nq_i) < 1 ) return( NULL )

    N_speckled_normal_total  =  vf$N_speckled_normal_total[1]
    N_primary_total          =  vf$N_primary_total[1]
    N_metastatic_total       =  vf$N_metastatic_total[1]

    VAF  =  NULL
    for( k in 1:length( rho ) ){
        # Scale for admixture rate of intact normal cells rho[ k ]:
        k_scale  =  ( 1 - rho[ k ] )   # rho[ k ] * ( N_primary_total + N_speckled_normal_total ) / N_primary_total

        for( i in nq_i ){
            w  =  which( vf$Ref_pos == i )
                                        # for primary tumor and speckled normal cells
            if ( ( N_primary_total + N_speckled_normal_total ) > 0 ){
                numenator_N    =  sum( vf[ w , c( 'N_speckled_normal', 'N_primary' ) ] * vf[ w, 'Copy_number_B'] ) / ( N_primary_total + N_speckled_normal_total )
                denominator_N  =  2 * ( 1 - ( sum( vf[ w , c( 'N_speckled_normal', 'N_primary' ) ] ) / ( N_primary_total + N_speckled_normal_total ) ) )  +
                                  sum( vf[ w , c( 'N_speckled_normal', 'N_primary' ) ] * ( vf[ w, 'Copy_number_B'] + vf[ w, 'Copy_number_A'] ) ) / ( N_primary_total + N_speckled_normal_total )
            } else {
                numenator_N    =  0
                denominator_N  =  0
            }
                                        # for metastatic cells
            if ( N_metastatic_total > 0 ){
                numenator_M    =  sum( vf[ w , 'N_metastatic'] * vf[ w, 'Copy_number_B'] ) / N_metastatic_total
                denominator_M  =  2 * ( 1 - ( sum( vf[ w , 'N_metastatic']  ) / N_metastatic_total ) )  +
                                  sum( vf[ w , 'N_metastatic'] * ( vf[ w, 'Copy_number_B'] + vf[ w, 'Copy_number_A'] ) ) / N_metastatic_total
            } else {
                numenator_M    =  0
                denominator_M  =  0
            }
                                        # VAF calculations:
            if ( numenator_N == 0 ){
                VAF_N_rho  =  0
            } else {
                VAF_N_rho  =  k_scale * numenator_N / ( 2*( 1 - k_scale ) + k_scale * denominator_N )
            }

            if ( numenator_M == 0 ){
                VAF_M_rho  =  0
            } else {
                VAF_M_rho  =  k_scale * numenator_M / ( 2*( 1 - k_scale ) + k_scale * denominator_M )
            }
                                        # save to data.frame:
            VAF_1 = data.frame( Chr = vf[ w[1], 'Chr' ] ,
                                site = i,
                                gene = vf[ w[1], 'Gene_name' ] ,
                                tumor_content = 1 - rho[ k ],
                                VAF_primary = VAF_N_rho,
                                VAF_metastatic = VAF_M_rho,

                                PointMut_ID    =  paste( vf[ w, 'PointMut_ID' ], collapse = ', ' )
                    )
            VAF_1[ is.na.data.frame( VAF_1 ) ]  =  0  #  division by 0 if rho = 1
            VAF  =  rbind( VAF, VAF_1 )
        }
    }

    return( VAF )
}


#
#' Function to write VAF data to the file
#'
#' @param input File's name for input parameters
#'
#' @return NULL 
#' 
#' @export
#'
#' @examples
#' NULL
write_VAF <- function( input ) {
    
    input.table <- read.table( input, header=F, sep="\t", stringsAsFactors = FALSE )

    for( ii in 1:nrow( input.table ) ) {
        V1 <- input.table[ ii, 1 ]
        V2 <- input.table[ ii, 2 ]

        # add quotes
        if ( length( grep( "^input\\.",  V1, perl=T ) ) > 0 ||
             length( grep( "^output\\.", V1, perl=T ) ) > 0
            ) V2 <- paste("\'", V2, "\'", sep='')

        text <- paste( V1, ' <- ', V2, sep='' )
        eval( parse( text = text ) )
    }
    timeVAF  =  time

    df_VAF_list <- list()
    df_rho_VAF_list <- list()
    for (time in c(timeVAF, NA)) {
        dtst = get_flow_data( cloneoutfile = input.cloneout,
                              pom.file     = input.pom,
                              pomA.file    = input.pomA, 
                              cna.file     = input.cna,
                              time         = time
                             )
        if (is.list(dtst)) {
            vf = get_VAF( pnt_mut = dtst$pnt_mut, pnt_mut_A = dtst$pnt_mut_A, data_last = dtst$data_last )
            if ( length( vf ) != 0 ){
                vf = cbind(Time = dtst$time_max, vf)
                df_VAF_list <- c(df_VAF_list, list(vf))
           
                rho_vf = get_rho_VAF( vf = vf, rho = ( 1 - tumor_content ) )
                rho_vf = cbind(Time = dtst$time_max, rho_vf)
                df_rho_VAF_list <- c(df_rho_VAF_list, list(rho_vf))
            }
        }
    }
    dfForVAF  =  do.call(rbind, df_VAF_list)
    dfVAF     =  do.call(rbind, df_rho_VAF_list)
    dfVAF     =  dfVAF[ , c(1, 5, 8, 2, 3, 4, 6, 7 ) ]
    dir.create( dirname( output.ForVAF ), showWarnings = FALSE, recursive = TRUE )
    dir.create( dirname( output.VAF    ), showWarnings = FALSE, recursive = TRUE )
    write.table( dfForVAF, file = output.ForVAF, append = FALSE, sep = "\t", row.names = FALSE )
    write.table( dfVAF,    file = output.VAF,    append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
    flog.info( paste0( 'VAF is saved in the file ', output.VAF ))
}


# CODE TO RUN -------------------------------------------------------------

if ( 0 ) { # set 0 for silence


library(stringr)

write_VAF( input = './Input/ForPosts/VAF.parameters.txt' )


}

