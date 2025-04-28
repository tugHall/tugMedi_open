# Functions related to Point Mutations:  --------------------------------------

#' Function to save 1 point mutation in a data frame
#'
#' @keywords internal
safe_pnt_mut  <-  function(pnt){
    return( pnt$safe() )
}

#' Function to copy of pnt1 without mutation info for allele A
#'
#' @return Object of class 'Point_Mutations' for another chromosome
#'
#' @keywords internal
copy_pnt_no_mutation  <-  function( pnt1 ){
    pnt2  =  Point_Mutations$new()
    pnt2$PointMut_ID  =  pnt1$PointMut_ID
    pnt2$Allele = 'A'
    pnt2$Parental_1or2 = ifelse( pnt1$Parental_1or2 == 2, as.integer(1), as.integer(2) )
    pnt2$Chr  =  pnt1$Chr
    pnt2$Ref_pos  =  pnt1$Ref_pos
    pnt2$Phys_pos  =  NA
    pnt2$Delta     =  NA
    pnt2$Copy_number  =  1
    pnt2$Gene_name    =  pnt1$Gene_name
    pnt2$MalfunctionedByPointMut  =  NA
    pnt2$mut_order    =  pnt1$mut_order
    pnt2$rst.ratio    =  0

    return( pnt2 )
}

#' Function to copy of point mutation info
#'
#' @return The same object of class 'Point_Mutations' with the same ID
#'
#' @keywords internal
copy_pnt  <-  function( pnt1 ){
    ### Function to copy of pnt1
    pnt2  =  Point_Mutations$new()
    pnt2$PointMut_ID  =  pnt1$PointMut_ID
    pnt2$Allele       = pnt1$Allele
    pnt2$Parental_1or2  =  pnt1$Parental_1or2
    pnt2$Chr  =  pnt1$Chr
    pnt2$Ref_pos  =  pnt1$Ref_pos
    pnt2$Phys_pos  =  pnt1$Phys_pos
    pnt2$Delta     =  pnt1$Delta
    pnt2$Copy_number  =  pnt1$Copy_number
    pnt2$Gene_name    =  pnt1$Gene_name
    pnt2$MalfunctionedByPointMut  =  pnt1$MalfunctionedByPointMut
    pnt2$mut_order    =  pnt1$mut_order
    pnt2$Ovlp_CNA_ID  =  pnt1$Ovlp_CNA_ID
    pnt2$rst.ratio    =  pnt1$rst.ratio

    return( pnt2 )
}

#' Function to copy CNA info
#'
#' @param CNA1 Object of class 'CNA_Mutations'
#'
#' @return The same object of class 'CNA_Mutations'
#'
#' @keywords internal
copy_CNA  <-  function( CNA1 ){
    CNA2  =  CNA_Mutations$new()
    CNA2$CNA_ID  =  CNA1$CNA_ID
    CNA2$Parental_1or2  =  CNA1$Parental_1or2
    CNA2$dupOrdel = CNA1$dupOrdel
    CNA2$Chr  =  CNA1$Chr
    CNA2$Ref_start  =  CNA1$Ref_start
    CNA2$Ref_end    =  CNA1$Ref_end
    CNA2$Gene_names =  CNA1$Gene_names
    CNA2$MalfunctionedByCNA  =  CNA1$MalfunctionedByCNA

    CNA2$mut_order  =  CNA1$mut_order
    CNA2$rst.ratio  =  CNA1$rst.ratio

    return( CNA2 )
}

#' Generation point mutation info
#'
#' @param onco1 Object of class 'OncoGene'
#' @param gm_1_2 List of two data frames (for 1st and 2nd parental chromosomes) with genes' location information
#' 
#' @keywords internal
get_point_mutation <- function( onco1, gm_1_2 ){

    prntl  =   sample( c(1,2), size = 1, replace = TRUE, prob = c( sum(onco1$cds_1), sum(onco1$cds_2) ))
    gene   =   ifelse( prntl == 1,
                       sample( onco1$name, size = 1, replace = TRUE, prob = onco1$cds_1/sum(onco1$cds_1) ),
                       sample( onco1$name, size = 1, replace = TRUE, prob = onco1$cds_2/sum(onco1$cds_2) ))

    gm     =   gm_1_2[[ prntl ]]

    sp  =  which( gm$Gene == gene)   
    p   =  sample( sp, size = 1, replace = TRUE,
                   prob = gm[sp,'Len'] / sum( gm[sp,'Len'] ) )   
    pos =  sample( gm$Start[p]:gm$End[p], 1, replace=FALSE)
    Chr =  gm$Chr[p]

    return( list( prntl, gene, pos, Chr ) )
}

#' Generation point mutation info for the particular gene
#'
#' @param onco1 Object of class 'OncoGene'
#' @param gm_1_2 List of two data frames (for 1st and 2nd parental chromosomes) with genes' location information
#' @param gene Gene's name where point mutation should be occurred
#' 
#' @keywords internal
get_point_mutation_for_gene <- function( onco1, gm_1_2, gene ){
    prntl  =   sample( c(1,2), size = 1, replace = TRUE, prob = c( sum(onco1$cds_1), sum(onco1$cds_2) ))

    gm     =   gm_1_2[[ prntl ]]

    sp  =  which( gm$Gene == gene)   
    p   =  sample( sp, size = 1, replace = TRUE,
                   prob = gm[sp,'Len'] / sum( gm[sp,'Len'] ) )   
    pos =  sample( gm$Start[p]:gm$End[p], 1, replace=FALSE)
    Chr =  gm$Chr[p]

    return( list( prntl, gene, pos, Chr ) )
}


#' Function to generate an object of class 'Point_Mutations'
#'
#' @param prntl Parental chromosome, could be 1 or 2
#' @param gene Gene name
#' @param pos Position of point mutation
#' @param onco1 Object of class 'OncoGene'
#' @param Chr Chromosome name
#' @param mutation If mutation is NOT NA then MalfunctionedByPointMut = TRUE, else it is defined by corresponding probabilities
#' @param rst.ratio Resist ratio against drug intervention
#'
#' @return Object of class 'Point_Mutations'
#'
#' @keywords internal
generate_pnt  <-  function( onco1, prntl, gene, pos, Chr, malfunc = NA, rst.ratio = 0 ) {
    
    pnt0 = Point_Mutations$new()
    pnt0$Gene_name = gene
    pnt0$PointMut_ID =  ifelse( length( pck.env$pnt_clones ) == 0, 1,
                                pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID + 2 )
    pnt0$Allele = 'B'  # Mutation occurs on allele B
    pnt0$Parental_1or2  =  as.integer(prntl)
    pnt0$Chr = Chr
    pnt0$Ref_pos  = pos
    pnt0$Phys_pos = pos
    pnt0$Delta = 0
    pnt0$Copy_number = 1
    u = ifelse( onco1$onsp[ which(onco1$name == gene) ] == 'o', uo, us)
    if ( is.na( malfunc ) ) {
        pnt0$MalfunctionedByPointMut  =  ifelse( (runif(1) < u), TRUE, FALSE )
    } else {
        pnt0$MalfunctionedByPointMut  =  malfunc
    }
    pck.env$mut_order  =  pck.env$mut_order  +  1
    pnt0$mut_order  =  pck.env$mut_order
    if ( is.null( rst.ratio ) ) {
        pnt0$rst.ratio  =  0
    } else {
        if ( is.na( rst.ratio ) ) {
            pnt0$rst.ratio  =  0
        } else {
            pnt0$rst.ratio  =  rst.ratio
           }
    }
    pck.env$pnt_clones = c( pck.env$pnt_clones, pnt0 )

    pnt1  =  copy_pnt_no_mutation( pnt0 )
    pck.env$pnt_clones = c( pck.env$pnt_clones, pnt1 )

    return( pnt0 )
}

#' Function to generate the same object of class 'Point_Mutations' with coping all information from input object
#'
#' @param pnt Object of class 'Point_Mutations'
#'
#' @return The same object of class 'Point_Mutations' with different ID
#'
#' @keywords internal
generate_to_copy_pnt  <-  function( pnt ) {

    pnt0 = Point_Mutations$new()
    pnt0$Gene_name = pnt$Gene_name
    pnt0$PointMut_ID =  ifelse( length( pck.env$pnt_clones ) == 0, 1,
                                pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID + 2 )
    pnt0$Allele = pnt$Allele
    pnt0$Parental_1or2  =  pnt$Parental_1or2
    pnt0$Chr = pnt$Chr
    pnt0$Ref_pos  = pnt$Ref_pos
    pnt0$Phys_pos = pnt$Phys_pos
    pnt0$Delta = pnt$Delta
    pnt0$Copy_number = pnt$Copy_number
    pnt0$MalfunctionedByPointMut  =  pnt$MalfunctionedByPointMut

    pnt0$mut_order  =  pnt$mut_order
    pnt0$Ovlp_CNA_ID  =  pnt$Ovlp_CNA_ID
    pnt0$rst.ratio  =  pnt$rst.ratio

    pck.env$pnt_clones  =  c( pck.env$pnt_clones, pnt0 )

    pnt1  =  copy_pnt_no_mutation( pnt0 )
    pck.env$pnt_clones  =  c( pck.env$pnt_clones, pnt1 )

    return( pnt0 )
}



# Functions related to CNA ------------------------------------------------


#' Function to choose probability of CNA mutation for several genes
#'
#' @param genes Names of genes, vector of names
#' @param dupOrdel It could be 'dup' or 'del' to denote duplication or deletion
#'
#' @return Single value of maximal probability from probabilities for several genes
#'
#' @keywords internal
get_u_cna <- function( genes, dupOrdel ){
    # input: genes, u_del or u_dup for oncogenes and suppressors
    u = NULL
    for (gene in genes ) {
        os = pck.env$onco$onsp[ which( pck.env$onco$name == gene) ]
        u1 = ifelse( dupOrdel == 'dup',
                     ifelse( os == 'o', uo_dup, us_dup ),   # for duplication
                     ifelse( os == 'o', uo_del, us_del ) )  # for deletion
        u = c( u, u1 )
    }

    return( max( u ) )
}


#' Generation CNA mutation info
#'
#' @param onco1 Object of class 'OncoGene'
#' @param dupOrdel It could be 'dup' or 'del' to denote duplication or deletion
#' @param gm_1_2 List of two data frames (for 1st and 2nd parental chromosomes) with genes' location information
#'
#' @return List of ( prntl - 1 or 2 parental chromosome, Chr - name of chromosome, genes - genes names, start_end - vector with start and end positions of CNA, w_cna - rows of CNA in gene_map data frame )
#'
#' @keywords internal
get_cna_mutation <- function( onco1, dupOrdel, gm_1_2 ){

    prntl  =  sample( c(1,2), size = 1, replace = TRUE, prob = c( sum(onco1$cds_1), sum(onco1$cds_2) ))
    lambda =  ifelse( dupOrdel == 'dup', ave_len_dup, ave_len_del )
    l_cna  =  round( rexp(1, 1/lambda) ) +1 # rpois( n = 1, lambda ) + 1   # length of CNA

    gm     =  gm_1_2[[ prntl ]]

    w1  =  sample( 1:length( gm$Start ), size = 1, prob = gm$Len )  # choose the row in gene_map
    Chr = gm$Chr[ w1 ]

    pos_st =  sample( gm$Start[w1]:gm$End[w1], 1, replace=FALSE)   ### position of start CNA
    pos_end  = pos_st  +  l_cna
    if ( pos_end > max( gm[ which( gm$Chr == Chr ), 'End'   ] ) ) pos_end  = max( gm[ which( gm$Chr == Chr ), 'End'   ] )
    
    w3  =  which( gm$Start < pos_end & gm$Chr == Chr )
    w_cna  = w3[ which(w3 >= w1)]    ### rows of CNA in gene_map
    start_end  =  c(0,0)
    start_end[1]  =  pos_st
    start_end[2]  =  ifelse( pos_end < gm$End[ max(w_cna) ], pos_end, gm$End[ max(w_cna) ] )
    genes   =   unique( gm$Gene[w_cna] )   ### Genes of CNA at same Chr

    return( list( prntl, Chr, genes, start_end, w_cna ) )
}


#' Generation CNA mutation info (new algorithm)
#'
#' @param onco1 Object of class 'OncoGene'
#' @param dupOrdel It could be 'dup' or 'del' to denote duplication or deletion
#' @param gm_1_2 List of two data frames (for 1st and 2nd parental chromosomes) with genes' location information
#' @param gene Name of a gene where CNA will start
#'
#' @return List of ( prntl - 1 or 2 parental chromosome, Chr - name of chromosome, genes - genes names, start_end - vector with start and end positions of CNA, w_cna - rows of CNA in gene_map data frame )
#'
#' @keywords internal
get_cna_mutation_for_gene <- function( onco1, dupOrdel, gm_1_2, gene ){

    w  =  which( onco1$name  == gene )
    CDS1  =  onco1$cds_1[ w ]
    CDS2  =  onco1$cds_2[ w ]
    if ( CDS1 <= 0 & CDS2 <= 0 ) stop( 'That is impossible to generate CNA in a gene when the gene is deleted.' )
    if ( CDS1 <= 0 ) prntl  =  2
    if ( CDS2 <= 0 ) prntl  =  1
    if ( CDS1 > 0 & CDS2 > 0 ){
        prntl  =  sample( c(1,2), size = 1, replace = TRUE, prob = c( CDS1, CDS2 ) / sum( CDS1 + CDS2 ) )
    }

    lambda =  ifelse( dupOrdel == 'dup', ave_len_dup, ave_len_del )
    l_cna  =  round( rexp(1, 1/lambda) ) + 1 # rpois( n = 1, lambda ) + 1   # length of CNA

    gm     =  gm_1_2[[ prntl ]]
    # Define Chromosome based on gene of interest:
    Chr    =  gm$Chr[ which( gm$Gene == gene ) ][ 1 ]

    gm_Chr      =  gm[ which(gm$Chr == Chr ), ]
    if ( length(gm_Chr$Start ) > 1 ){
        for( i in 2:( length(gm_Chr$Start )  ) ){
            gm_Chr[ i, 'gap_before' ]  =  gm_Chr$Start[ i ]  -  gm_Chr$End[ i - 1 ]
        }
    }
    gm_Chr[ 1, 'gap_before' ]  =  - l_cna

    # Get the ranges for a chromosome:
    rngs  =  data.frame( Start = 0, End = 0 )
    i_r   =  c( unique( c( 1, which( l_cna < gm_Chr$gap ) ) ), length( gm_Chr$End ) )
    if ( length( i_r ) == 0 ) stop( 'There is no range for CNA generation.' )
    for( i in  1:( length( i_r ) - 1 ) ){
        rngs[ i, ]  =  c( gm_Chr$Start[ i_r[ i ] ] - l_cna, gm_Chr$End[ i_r[ i + 1 ] - 1 ] )
    }
    rngs$End[ i ]  =  gm_Chr$End[ i_r[ i + 1 ] ]
    rngs$length  =  rngs$End - rngs$Start

    w1   =  sample( x = 1:length( rngs$End ), size = 1, replace = TRUE, prob = rngs$length )    ### Area for generation of CNA
    pos_st =  sample( rngs$Start[ w1 ] : rngs$End[ w1 ], 1, replace=FALSE)   ### position of start CNA
    pos_end  = pos_st  +  l_cna
    if ( pos_st  < min( gm[ which( gm$Chr == Chr ), 'Start' ] ) ) pos_st   = min( gm[ which( gm$Chr == Chr ), 'Start' ] )
    if ( pos_end > max( gm[ which( gm$Chr == Chr ), 'End'   ] ) ) pos_end  = max( gm[ which( gm$Chr == Chr ), 'End'   ] )
    
    w3  =  min( which( gm$End   >= pos_st  & gm$Chr == Chr ) )
    w4  =  max( which( gm$Start <= pos_end & gm$Chr == Chr ) )
    w_cna  = w3 : w4  ### rows of CNA in gene_map

    genes   =   unique( gm$Gene[ w_cna ] )   ### Genes of CNA at the same chromosome

    start_end  =  c(0,0)
    start_end[ 1 ]  =  pos_st
    start_end[ 2 ]  =  pos_end

    return( list( prntl, Chr, genes, start_end, w_cna ) )
}


#' Function to get order of mutation for all possible types
#'
#' @return data.frame with fields order, type, ID
#'
#' @keywords internal
mixed_mut_order   <-   function( clone1 ) {

    order_tp_num  <-  data.frame( order = NULL, type = NULL, ID = NULL )
    i  =  0
    if ( length( clone1$PointMut_ID ) > 0  & clone1$PointMut_ID[ 1 ]  != 0 ){
        for (i in 1:length( clone1$PointMut_ID )) {
            order_tp_num[i,'type']   =  'pnt'
            order_tp_num[i,'ID']     =  as.numeric( clone1$PointMut_ID[i] )
            order_tp_num[i,'order']  =  pck.env$pnt_clones[[ order_tp_num[i,'ID'] ]]$mut_order 
        }
    }

    if ( length( clone1$CNA_ID ) > 0 & clone1$CNA_ID[ 1 ]  != 0 ){
        for (j in 1:length( clone1$CNA_ID ) ) {
            order_tp_num[j+i,'ID']      =   clone1$CNA_ID[ j ]
            cn1   =   pck.env$cna_clones[[ order_tp_num[ j+i, 'ID' ] ]]
            order_tp_num[j+i,'order']   =   cn1$mut_order 
            order_tp_num[j+i,'type']    =   cn1$dupOrdel  
        }
    }

    if ( length(order_tp_num)  !=  0 ){
        order_tp_num  <-  order_tp_num[ order( order_tp_num$order), ]
        row.names(order_tp_num) <- 1:length( order_tp_num[,'type'])
    } else order_tp_num = NULL

    return( order_tp_num )
}

#'  Function to generate object of CNA mutation
#'
#' @return Object of class 'CNA_Mutations'
#'
#' @keywords internal
generate_cna  <-  function( onco1, prntl, genes, start_end, Chr, malfunc = NA, dupOrdel, rst.ratio = 0) {

    cna0 = CNA_Mutations$new()
    cna0$CNA_ID =  ifelse( length( pck.env$cna_clones ) == 0, 1,
                           pck.env$cna_clones[[ length( pck.env$cna_clones ) ]]$CNA_ID + 1 )
    cna0$Parental_1or2  =  prntl
    cna0$dupOrdel = dupOrdel
    cna0$Chr = Chr
    cna0$Ref_start  = start_end[1]
    cna0$Ref_end    = start_end[2]
    cna0$Gene_names = genes
    u = get_u_cna( genes, dupOrdel )
    if ( is.na( malfunc ) ) {
        cna0$MalfunctionedByCNA  =  ifelse( ( runif(1) < u ), TRUE, FALSE )
    } else {
        cna0$MalfunctionedByCNA  =  malfunc
    }
    pck.env$mut_order  =  pck.env$mut_order  +  1
    cna0$mut_order  =  pck.env$mut_order
    if ( is.null( rst.ratio ) ) {
        cna0$rst.ratio  =  0
    } else {
        if ( is.na( rst.ratio ) ) {
            cna0$rst.ratio  =  0
        } else {
            cna0$rst.ratio  =  rst.ratio
        }
    }
    pck.env$cna_clones = c( pck.env$cna_clones, cna0 )

    return( cna0 )
}

#' Function to add the mutations to the data.frame gene_map
#'
#' @return list( gm1, gm2 ), where gm1 and gm2 are data.frames gene_maps with mutation information
#'
#' @keywords internal
modify_gene_map  <-  function( clone1 , onco1 ){
    gm1    =  gm2  =  gene_map 
    gm1$pnts  =  gm2$pnts  = ''
    mixed_order  =  mixed_mut_order( clone1 = clone1 )
    if ( is.null(mixed_order) ) return( list( gm1, gm2 ) )

    for (l in 1:nrow( mixed_order ) ) {
        cn1  =  NULL
        if ( mixed_order[l, 'type'] == 'del' ) {

            cn1        =  pck.env$cna_clones[[  mixed_order$ID[l] ]]   # cn[ mixed_order$ID[l], ]
            Ref_start  =  cn1$Ref_start
            Ref_end    =  cn1$Ref_end
            Chr        =  cn1$Chr

            # change gene map only one of two chromosomes: 1 or 2
            ifelse( cn1$Parental_1or2 == 1, gm <- gm1, gm <- gm2 )

            gm  =  add_deletion( gm = gm, Ref_start = Ref_start,
                                 Ref_end = Ref_end, Chr = Chr )

            ### come back to gene map for certain chromosome 1 or 2
            ifelse( cn1$Parental_1or2 == 1, gm1 <- gm, gm2 <- gm )
            rm( gm )
        }

        if ( mixed_order[l, 'type'] == 'dup' ) {

            cn1        =  pck.env$cna_clones[[  mixed_order$ID[l] ]]
            Ref_start  =  cn1$Ref_start
            Ref_end    =  cn1$Ref_end
            Chr        =  cn1$Chr

            # change gene map only one of two chromosomes: 1 or 2
            ifelse( cn1$Parental_1or2 == 1, gm <- gm1, gm <- gm2 )

            gm  =  add_duplication( gm = gm, Ref_start = Ref_start,
                                    Ref_end = Ref_end, Chr = Chr )

            ### come back to gene map for the certain chromosome 1 or 2
            ifelse( cn1$Parental_1or2 == 1, gm1 <- gm, gm2 <- gm )
            rm( gm )
        }

        if ( mixed_order[l, 'type'] == 'pnt' ) {

            pn1        =  pck.env$pnt_clones[[ mixed_order$ID[l] ]] 
            pos_pnt1   =  pn1$Ref_pos
            Chr        =  pn1$Chr

            # change gene map only one of two chromosomes: 1 or 2
            ifelse( pn1$Parental_1or2 == 1, gm <- gm1, gm <- gm2 )

            gm  =  add_pnt_mutation( gm = gm, pos_pnt = pos_pnt1, Chr = Chr )

            ### come back to gene map for the certain chromosome 1 or 2
            ifelse( pn1$Parental_1or2 == 1, gm1 <- gm, gm2 <- gm )
            rm( gm )

        }

    }

    return( list( gm1, gm2 ) )
}

#' Function to add point mutation to data.frame gene_map (chromosomal location data frame)
#'
#' @param gm Chromosomal location data frame
#' @param pos_pnt Position of point mutation
#' @param Chr Chromosome name
#'
#' @return Chromosomal location data frame with additional point mutation info
#'
#' @keywords internal
add_pnt_mutation   <-  function( gm = gm, pos_pnt, Chr = Chr ){

    w = which( gm$Chr == Chr & gm$Start <= pos_pnt  & gm$End >= pos_pnt )

    if ( length( w ) > 1 ) w = w[ 1 ]
    if ( length(w) == 0 ) {
        pos_pnt_original = pos_pnt
        pos_pnt  <-  .correct_pos_pnt( gm, pos_pnt, Chr )
        if ( pos_pnt != pos_pnt_original ){
            flog.info(msg = paste0('Position of point mutation was corrected from ', pos_pnt_original, 
                                   ' to new one ', pos_pnt,   ' due to CNA presence.') )
        }
        w = which( gm$Chr == Chr & gm$Start <= pos_pnt  & gm$End >= pos_pnt )
    }
    
    if ( length(w) == 0 ) {
        write( paste( 'Chr:', Chr, '| pos:', pos_pnt ), 
               '__tmp.addPnt__', sep="\t" )
        write( paste( 'gm Chr:', gm$Chr, '| gm Pos:', gm$Start, '-', gm$End ), 
               '__tmp.addPnt__', sep="\t", append=T )
        write( paste( 'w:',       w       ), '__tmp.addPnt__', sep="\t", append=T )
        write( paste( 'pos_pnt:', pos_pnt ), '__tmp.addPnt__', sep="\t", append=T )
        write.table( gm,                     '__tmp.addPnt__', sep="\t", append=T, quote=F, row.names=F, )
        stop( cat( ' Point mutation does not meet into parental chr ranges, ', 
                   ' which may differ between clones.', 
                   ' Wrote info in __tmp.addPnt__ file.', sep="\n" )  )
    } else {
        gm[w, 'pnts']  =  ifelse( gm[w, 'pnts'] ==  '',
                                  gm[w, 'pnts']  <-  as.character(pos_pnt),
                                  gm[w, 'pnts']  <-  paste(gm[w, 'pnts'], pos_pnt, sep = ',') )
    }

    return( gm )
}


.correct_pos_pnt  <-  function( gm, pos_pnt, Chr ){
    
    gm_original  =  pck.env$gene_map[ which( pck.env$gene_map$Chr == Chr ), ]
    gm           =  gm[ which( gm$Chr == Chr ), ]
    
    w = which( gm_original$Chr == Chr & gm_original$Start <= pos_pnt  & gm_original$End >= pos_pnt )
    if ( length( w ) == 0 ){
        flog.info( paste0('Position  ', pos_pnt, ' can not be corrected. ' ) )
        return( pos_pnt )
    }
    w_origin  =  w 
    
    
    nr  =  min( nrow( gm ), nrow( gm_original ) )
    
    nq  =  unique( gm_original$Start[ 1:nr ]  -  gm$Start[ 1:nr ] )
    nq  =  nq[ ! (nq %in% 0 ) ]
    
    for ( pos in ( pos_pnt - nq ) ){
        w = which( gm$Chr == Chr & gm$Start <= pos  & gm$End >= pos )
        if ( length( w ) != 0 ) return( pos )
    }
    
    # if it is could not be assigned exactly, let's to assign it as close as possible to the original position:
    st_original  =  min( gm_original$Start[ w_origin ] )
    w  =  which.min( abs( gm$Start -  st_original ) )[ 1 ]
    
    pos_pnt_corrected  =  round( runif( n = 1, min = gm$Start[ w ], max = gm$End[ w ] ) )
    
    return( pos_pnt_corrected )
}

#' Function to add deletion to gene map (chromosomal location data frame)
#' 
#' @return Chromosomal location data frame with additional deletion info
#'
#' @keywords internal
add_deletion  <-  function( gm, Ref_start, Ref_end, Chr ){

    if (Ref_end < Ref_start)  stop( 'End should be larger or equal then Start for CNA' )
    w1  =  max( which( gm$Start <= Ref_start  &  gm$Chr == Chr ) )
    if ( length( w1) == 0 ) w1 =  min( which( gm$Chr == Chr ) )
    w2  =  max( which( gm$Start <= Ref_end    &  gm$Chr == Chr ) )

    if ( w1 == w2 ) {
        gm  =  rbind( gm[1:w1,], gm[ w1:nrow(gm), ] )  # duplicate the row w1=w2
        w2  =  w1 + 1
    }
    ### Correction of rows w1 and w2:
    if ( gm[ w1, 'Start']   <   Ref_start ){  # to delete part of row or leave as before
        gm[ w1, 'End']    =   ifelse( Ref_start <=  gm[ w1, 'End'], Ref_start - 1, gm[ w1, 'End']  )
        gm[ w1, 'Len']    =   gm[ w1, 'End']  -  gm[ w1, 'Start']    +  1
        gm[ w1, 'pnts']   =   check_pnts( gm[ w1,  ] )
    } else {
        w1  =  w1  -  1   # to delete whole row
    }

    if ( gm[ w2, 'End']   >   Ref_end ){   # to delete part of row or leave as before
        gm[ w2, 'Start']    =   Ref_end + 1
        gm[ w2, 'Len']      =   gm[ w2, 'End']  -  gm[ w2, 'Start']  +  1
        gm[ w2, 'pnts']     =   check_pnts( gm[ w2,  ] )
    } else {
        w2  =  w2  +  1    # to delete whole row
    }

    ### delete part of gm related to deletion
    if ( w2 - w1 > 1 )      gm   =    gm[ -c( (w1+1):(w2-1) ), ]

    # To subtract the delta from positions
    w    =  which( gm$Start > Ref_end & gm$Chr == Chr )
    dlt  =  Ref_end  -  Ref_start  +  1
    if ( length(w) > 0 ){
        gm[w, 'Start']    =  gm[w, 'Start']  -  dlt
        gm[w, 'End']      =  gm[w, 'End']    -  dlt
        gm[w, 'pnts']     =  unlist( sapply( w, FUN = function(x)  pnts_add_dlt( gm[ x, ], dlt  =  -dlt)  ) )
    }

    return( gm )
}

#' Function to add duplication to gene map (chromosomal location data frame)
#'
#' @param gm Chromosomal location data frame
#' @param Ref_start Starting position of duplication
#' @param Ref_end Final position of duplication
#' @param Chr Chromosome name
#'
#' @return Chromosomal location data frame with additional duplication info
#'
#' @keywords internal
add_duplication  <-  function( gm, Ref_start, Ref_end, Chr ){

    if (Ref_end < Ref_start)  stop( 'End should be larger or equal then Start for CNA' )
    dlt        =  Ref_end  -  Ref_start  +  1  # delta for all next chromosomal positions

    ### Change gene_map with CNA duplication:
    w1  =  max( which( gm$Start <= Ref_start  &  gm$Chr == Chr ) )
    if ( length( w1) == 0 ) w1 =  min( which( gm$Chr == Chr ) )
    w2  =  max( which( gm$Start <= Ref_end    &  gm$Chr == Chr ) )

    # ADD rows for duplication
    if ( (w2+1) <= nrow(gm) ){
        gm  =  rbind( gm[1:w2, ], gm[w1:w2, ] , gm[(w2+1):nrow(gm), ])
    } else {
        gm  =  rbind( gm[1:w2, ], gm[w1:w2, ] )
    }

    w3  =  w2 + 1

    w5  =  max( which( gm$Chr  ==  Chr ) )
    # change the final position BEFORE duplication
    gm[w2, 'End']    =  ifelse( gm[w2, 'End'] > Ref_end, Ref_end, gm[w2, 'End'] )
    gm[w2, 'Len']    =  gm[w2, 'End'] - gm[w2, 'Start'] + 1
    gm[w2, 'pnts']   =  check_pnts( gm[ w2,  ] )

    # Add delta to all positions after duplication
    gm[w3:w5, 'Start']  =  gm[w3:w5, 'Start'] + dlt
    gm[w3:w5, 'End']    =  gm[w3:w5, 'End']   + dlt
    gm[w3:w5, 'pnts']   =  unlist( sapply( w3:w5, FUN = function(x) pnts_add_dlt( gm[ x, ], dlt )  ) )

    # Correct Start of w3 row if it's necessary
    if ( gm[w3, 'End']  >= (Ref_start + dlt) ){
        gm[w3, 'Start']  =  Ref_start + dlt
        gm[w3, 'Len']    =  gm[w3, 'End'] - gm[w3, 'Start'] + 1
        gm[w3, 'pnts']   =  check_pnts( gm[ w3,  ] )
    } else { 
        # delete the row if duplication started from outside of gene positions
        gm  =  gm[-w3,]
    }

    return( gm )
}

#' Function to subtract delta from position of point mutations
#' 
#' @param gm_w1 A row from data.frame gene_map
#' @param dlt Delta to subtract from positions of point mutations
#'
#' @return Return the pnts - dlt for one row of data.frame gene_map
#'
#' @keywords internal
pnts_add_dlt  <-  function( gm_w1 , dlt ){

    if ( is.null(gm_w1) )  stop( 'The input is null' )
    pnts = gm_w1$pnts
    if ( length(pnts)  ==  0 )  stop( 'The input is empty' )
    if ( pnts == '' ) return( '' )
    pnts  =  as.numeric( unlist( strsplit( pnts, split = ',') ) )
    if ( !is.numeric( pnts ) ) stop( 'Incorrect format of points mutation: should be numeric')
    pnts  =  pnts  +  dlt
    pnts  =  paste( as.character( pnts ) , collapse = ',')

    return( pnts )
}

#' Function to check what pnts do fall into the range
#' @param gm_w1 A row from data.frame gene_map
#'
#' @return Return the point mutations which fall into the range
#'
#' @keywords internal
check_pnts  <-  function( gm_w1 ){

    if ( is.null(gm_w1) )  stop( 'The input is null' )
    pnts = gm_w1$pnts
    if ( length(pnts)  ==  0 )  stop( 'The input is empty' )
    if ( pnts == '' ) return( '' )
    pnts  =  as.numeric( unlist( strsplit( pnts, split = ',') ) )
    if ( !is.numeric( pnts ) ) stop( 'Incorrect format of points mutation: should be numeric')

    check_log  =  unlist( sapply( pnts, FUN = function( x ) ( x <= gm_w1$End & x >= gm_w1$Start ) ) )
    pnts  =  paste( as.character( pnts[ check_log ] ) , collapse = ',')
    if ( length(pnts) == 0 ) pnts = ''

    return( pnts )
}


#' Function to get length of CDS and of genes from data.frame gene_map and related probabilities
#' @param gm data.frame gene_map with info about genes' location
#'
#' @keywords internal
get_cds_rna  <-  function( gm ){

    name0  =  unique( pck.env$gene_map$Gene )
    rna0  =  cds0  =  NULL
    for ( i in 1:length(name0) ) {
        w     =  which( gm$Gene == name0[ i ] )
        if ( length( w ) != 0 ){
            cds0  =  c( cds0, sum( gm$End[w]  -  gm$Start[w] + 1 ) )
            rna0  =  c( rna0, max( gm$End[w] ) - min( gm$Start[w]) + 1  )
        } else {
            cds0  =  c( cds0, 0 )
            rna0  =  c( rna0, 0 )
        }
    }

    sum_prob   =  sum( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) )
    prob       =  c(   m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) ) / sum_prob
    p0         =   (1-m0)**sum(cds0) * (1-m_dup)**sum(rna0) * (1-m_del)**sum(rna0)

    return( list( names = name0, CDS = cds0, RNA = rna0, PROB = prob, SUM = sum_prob, P0 = p0 ) )

}
