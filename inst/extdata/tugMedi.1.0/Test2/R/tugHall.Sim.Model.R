

# Functions related to trial  ---------------------------------------------



#' Aggregate data of a clone for environment object
#'
#' @param env Object of class 'Environ'
#' @param clones List of all the objects of class 'Clone'
#'
#' @keywords internal
sum_cell <- function(env, clones) {
    if (length(clones) > 0) {
        avg = apply(matrix(unlist(lapply(clones, sum_mutation)),ncol=length(clones)),1,sum)
        sNMP  =  pck.env$env$N + pck.env$env$M + pck.env$env$P
        pck.env$env$c = avg[1] / sNMP
        pck.env$env$d = avg[2] / sNMP
        pck.env$env$i = avg[3] / sNMP
        pck.env$env$a = avg[4] / sNMP
        pck.env$env$k = avg[5] / sNMP
        pck.env$env$E = avg[6] / sNMP
        pck.env$env$Nmax = avg[7] / sNMP
        pck.env$env$im = avg[8] / sNMP
        pck.env$env$Ha = avg[9] / sNMP
        pck.env$env$Him = avg[10] / sNMP
        pck.env$env$Hi = avg[11] / sNMP
        pck.env$env$Hb = avg[12] / sNMP
        pck.env$env$Hd = avg[13] / sNMP
        pck.env$env$type = avg[14] / sNMP
        pck.env$env$mutden = avg[15] / sNMP
    } else {
        pck.env$env$M = 0
        pck.env$env$N = 0
        pck.env$env$P = 0
        pck.env$env$c = 0
        pck.env$env$d = 0
        pck.env$env$i = 0
        pck.env$env$a = 0
        pck.env$env$k = 0
        pck.env$env$E = 0
        pck.env$env$Nmax = 0
        pck.env$env$im = 0
        pck.env$env$Ha = 0
        pck.env$env$Him = 0
        pck.env$env$Hi = 0
        pck.env$env$Hb = 0
        pck.env$env$Hd = 0
        pck.env$env$type = 0
    }

}

#' Serve function for sum_cell() function
#'
#' @keywords internal
sum_mutation <- function(clone1) {

    return(c(clone1$c*clone1$N_cells,      clone1$d*clone1$N_cells,    clone1$i*clone1$N_cells,
             clone1$a*clone1$N_cells,      clone1$k*clone1$N_cells,    clone1$E*clone1$N_cells,
             clone1$Nmax*clone1$N_cells,   clone1$im*clone1$N_cells,   clone1$Ha*clone1$N_cells,
             clone1$Him*clone1$N_cells,    clone1$Hi*clone1$N_cells,   clone1$Hb*clone1$N_cells,
             clone1$Hd*clone1$N_cells,     ifelse(clone1$invasion,1,0),
             clone1$mutden*clone1$N_cells))
           
}

#' Function to calculate number of cells = normal + primary tumor + metastatic
#'
#' @param env Object of class 'Environ'
#' @param clones List of all the objects of class 'Clone'
#'
#' @return Number of all the cells in a simulation (normal + primary tumor + metastatic)
#'
#' @keywords internal
sum_N_P_M <- function(env, clones) {
    if (length(clones) > 0) {
        avg = apply(matrix(unlist(lapply(clones, number_N_P_M)),ncol=length(clones)),1,sum)  #  /length(clones)
        pck.env$env$N = avg[ 1 ]
        pck.env$env$P = avg[ 2 ]
        pck.env$env$M = avg[ 3 ]
        sumNPM = pck.env$env$N + pck.env$env$P + pck.env$env$M
        if ( length( sumNPM ) == 0 ) sumNPM = 0
        return( sumNPM )
    }
}

#' Function to get number of cells of a clone with indicator (normal, primary tumor or metastatic)
#'
#' @param clone1 Object of class 'Clone'
#'
#' @return Vector c( N_normal, N_primary, N_metastatic )
#'
#' @keywords internal
number_N_P_M <- function(clone1) {
    indicator = c(  ifelse(clone1$invasion,0,1) * ifelse(clone1$primary,0,1),
                    ifelse(clone1$invasion,0,1) * ifelse(clone1$primary,1,0),
                    ifelse(clone1$invasion,1,0) )
    return( clone1$N_cells * indicator )
}


##### trial #############################################################
#' 
#' Function to conduct trials
#'
#' @return Number of new clones originated by clone1
#'
#' @keywords internal
conduct_trials <- function( clone1, onco1, model_param, drug_int_param=NULL ) {

   # Constant cell death trial
   N_die <- .get_rtrunc4trial( spec = model_param$trial$spec, max = clone1$N_cells,
                               prob = clone1$k )

   # Apoptosis trial
   if ( pck.env$tumbler_for_apoptosis_trial ){
      N_die <- N_die + .get_rtrunc4trial( spec = model_param$trial$spec, max = clone1$N_cells,
                                          prob = clone1$a )
   }

   # Drug intervention trial
   if ( ( ! is.null( drug_int_param ) ) && clone1$drug.on ){
      N_die = N_die + .get_rtrunc4trial( spec = model_param$trial$spec, max = clone1$N_cells,
                                         prob = drug_int_param[[ 'kill_prob' ]] )
   }

   # Invasion / metastasis trial
   if ( pck.env$tumbler_for_metastasis_trial ) {
      N_im <- .im_trial( model_param, clone1 )
      N_die <- N_die + N_im$die
      if ( N_im$suc > 0 ) { clone1$invasion.suc <- N_im$suc
                            clone1$invasion     <- TRUE     }
   }
   
   # Cell division trial
   N_divisions = .get_rtrunc4trial( spec = model_param$trial$spec, max = clone1$N_cells,
                                    prob = clone1$d * 
                                           ( 1 - pck.env$m0    
                                               - pck.env$m_del 
                                               - pck.env$m_dup ) )
                                           # Mutated clones are later addressed outside.

   # Replication limit trial
   if ( pck.env$tumbler_for_immortalization_trial ) {
      frac <- 1 - pnorm( pck.env$ctmax, mean = clone1$c, sd = sqrt(clone1$c.var) )
      N_divisions <-                round( N_divisions * (1 - frac) ) + 
                     calc_binom( 1, round( N_divisions *      frac  ), 1 - clone1$i )
   }

   # Per-cell division counter. Distribution specified by the mean and variance, 
   # where the binomials with different probs approximated to the normals
   if ( clone1$N_cells > 0 ) {
      pp <- 1 / clone1$N_cells
      clone1$c     = clone1$c     + N_divisions * pp
      clone1$c.var = clone1$c.var + N_divisions * pp * ( 1 - pp )
   }


   # Updates ----

   # Total division counter
   pck.env$env$n_divisions  =  pck.env$env$n_divisions + N_divisions

   # Cell number, be updated later altogether for <= 0
   clone1$N_cells = clone1$N_cells + N_divisions - N_die

   names( N_divisions ) <- clone1$id
   return( N_divisions )
}


#
.get_rtrunc4trial <- function( spec, max, prob ) {

   if ( spec == 'pois' )  {
      ret <- rtrunc( 1, spec="pois",  a=-Inf, b=max,
                        lambda = max * prob * pck.env$kappa )
   }

   if ( spec == 'binom' ) {
      ret <- calc_binom( 1, max, prob )
   }

   return( ret )
}


#
.im_trial <- function( model_param, clone1 ) {

   N_suc <- 0; N_die <- 0
   if ( ! is.nan(clone1$im) && clone1$im > 0 && 
        ! clone1$invasion ) {

      if ( model_param$metastatic_model == 'proportional_metastatic') {
         # cells trying im
         N_im <- .get_rtrunc4trial( spec = model_param$trial$spec, max = clone1$N_cells,
                                       prob = clone1$im )

         # success and fail
         N_suc <- rbinom( 1, N_im, prob = pck.env$Zim )
         N_die <- N_im - N_suc
      }

      # Under test
      if ( model_param$metastatic_model == 'threshold_metastatic') {
         # im prob. decreases death rate
         N_die <- .get_rtrunc4trial( spec = model_param$trial$spec, max = clone1$N_cells,
                                     prob = 1 - clone1$im )

         # If over the threshold, all of the survivor will transform.
         N_suc = ifelse( clone1$im >= pck.env$Zim,
                         clone1$N_cells - N_die, 0 )
      }
   }

   return( list( suc = N_suc, die = N_die ) )
}


# Mutagenesis trial -------------------------------------------------------

#' Function to insert mutation
#' 
#' @return Changed object clone1, add related mutations to the lists of point mutations and/or CNA mutations
#'
#' @keywords internal
insert_mutations <- function( clone1, onco1, mut_info ) {
    gm   <-  modify_gene_map( clone1 , onco1 )
    ncol <- list( base   = 1, 
                  Rother = 3, 
                  Rint   = 9 )

    for ( ii in 1:nrow( mut_info ) ) {
        mut_row1 <- mut_info[ ii, ]
        type <- mut_row1[ 1, 'Type' ]

        if ( type  ==  'pom' ) {
            pom_info <- .get_pom_info( onco1, gm, ncol, mut_row1 )
            prntl    <- pom_info$prntl
            gene     <- pom_info$gene
            pos      <- pom_info$pos
            Chr      <- pom_info$Chr
            pnt0     <- pom_info$pnt0

            if ( clone1$PointMut_ID[1] == 0 ) {
                id =     pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID
            } else  { 
                id =  c( clone1$PointMut_ID, 
                         pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID ) 
            }
            clone1$PointMut_ID = id 

            if ( pnt0$MalfunctionedByPointMut ){
                ck_dr  =  check_driver( clone1 = clone1, gene = gene )
                if ( ck_dr ){
                    clone1$gene[ which( pck.env$onco$name == gene ) ] = 1
                } else {
                    clone1$pasgene[ which( pck.env$onco$name == gene ) ] = 1
                }
            } else {
                clone1$pasgene[ which( pck.env$onco$name == gene ) ] = 1
            }

            gm[[ prntl ]]  =  add_pnt_mutation( gm = gm[[ prntl ]], pos_pnt = pos, Chr = Chr )
        }


        if ( type == 'dup' | type == 'del' ) {
            cna_info  <- .get_cna_info( onco1, gm, ncol, mut_row1 )
            prntl     <- cna_info$prntl
            genes     <- cna_info$genes
            start_end <- cna_info$start_end
            Chr       <- cna_info$Chr
            cna0      <- cna_info$cna0

            if ( clone1$CNA_ID[1] == 0 ) {
                id  =  pck.env$cna_clones[[ length( pck.env$cna_clones ) ]]$CNA_ID
            } else  id  =  c( clone1$CNA_ID, pck.env$cna_clones[[ length( pck.env$cna_clones ) ]]$CNA_ID )
            clone1$CNA_ID  =  id

            if ( cna0$MalfunctionedByCNA ){
                for( gene in genes ){
                    ck_dr  =  check_driver( clone1 = clone1, gene = gene )
                    if ( ck_dr ){
                        clone1$gene[ which( pck.env$onco$name == gene ) ] = 1
                    } else {
                        clone1$pasgene[ which( pck.env$onco$name == gene ) ] = 1
                    }

                }
            } else {
                clone1$pasgene[ unlist( sapply(genes, FUN = function(x) which(pck.env$onco$name == x) ) ) ] = 1
            }

            ### Change the gene_map for the corresponding chromosome
            ifelse( type == 'dup',
                    gm[[ prntl ]] <- add_duplication( gm = gm[[ prntl ]], 
                                        Ref_start = start_end[1], Ref_end = start_end[2], 
                                        Chr = Chr ),
                    gm[[ prntl ]] <- add_deletion( gm = gm[[ prntl ]], 
                                        Ref_start = start_end[1], Ref_end = start_end[2], 
                                        Chr = Chr )  )

            ### Check which point mutations match into the CNA
            sp   = FALSE
            sp_A = FALSE
            if ( clone1$PointMut_ID[ 1 ] != 0 ){
                sp = unlist( sapply( clone1$PointMut_ID , FUN = function( x )  {
                    chk_pnt_mut( pnt1  =  pck.env$pnt_clones[[ x ]], Ref_start = start_end[1],
                                 Ref_end = start_end[2], Chr = Chr, prntl  =  prntl )
                }) )
                ### to check original allele A do/don't match into CNA:
                prntl_inv  =  ifelse( prntl  ==  1, 2, 1 )
                sp_A = unlist( sapply( clone1$PointMut_ID , FUN = function( x )  {
                    chk_pnt_mut( pnt1  =  pck.env$pnt_clones[[ x ]], Ref_start = start_end[1],
                                 Ref_end = start_end[2], Chr = Chr, prntl  =  prntl_inv )
                } ) )


            }

            if ( any( sp ) | any( sp_A ) ){
                ### Before changing point mutation for NEW clone
                ### we have to copy it and avoid any changes in OTHER (parental) clones

                CNA_ID_last  =  clone1$CNA_ID[ length( clone1$CNA_ID ) ]

                ### copy pnt for a NEW clone
                for( q in clone1$PointMut_ID[ sp | sp_A ] ){
                    pnt0 = generate_to_copy_pnt( pnt = pck.env$pnt_clones[[ q ]] )
                    x  =  which( clone1$PointMut_ID == q )
                    clone1$PointMut_ID[ x ]  =  pnt0$PointMut_ID
                }

                sapply( clone1$PointMut_ID[ sp ],
                        FUN = function( x ) change_pnt_by_cna( pnt1 = pck.env$pnt_clones[[ x ]],
                                                               start_end, type, CNA_ID_last ) )
                sapply( clone1$PointMut_ID[ sp_A ],
                        FUN = function( x ) change_allele_A_by_cna( pnt1 = pck.env$pnt_clones[[ x+1 ]],
                                                                    start_end, type, CNA_ID_last ) )
            }

        }
    }

    onco_update( onco1, gm )
    if ( sum( clone1$gene ) > 0 ) clone1$primary = TRUE
}


#' Function to check that malfunction of a gene is driver or not
#'
#' @param clone1  Object of class 'Clone'
#' @param gene Gene name to check the malfunction of this gene cases driver state or does not
#'
#' @return logical value about malfunction of a gene is driver or not
#'
#' @keywords internal
check_driver  <-  function( clone1, gene ){

    driver   =  FALSE   # Factor of driver for a gene
    driver1  =  FALSE   # Factor of malfunction for chromosome 1 of the gene
    driver2  =  FALSE   # Factor of malfunction for chromosome 2 of the gene

    il  =  which( pck.env$DR$gene == gene )

    # Check the point mutations (can be dominant or recessive):
    ids  =  clone1$PointMut_ID
    if ( ( ids != 0 )[ 1 ] ){
        for( id in ids ){
            gn  =  pck.env$pnt_clones[[ id ]]$Gene_name
            mf  =  pck.env$pnt_clones[[ id ]]$MalfunctionedByPointMut
            if ( mf & gene == gn & ( pck.env$DR[ il, 'role_pom' ] != 'N' ) ){
                if ( pck.env$pnt_clones[[ id ]]$Parental_1or2 == 1 ){
                    driver1  =  TRUE
                } else {
                    driver2  =  TRUE
                }

                if ( pck.env$DR[ il, 'role_pom' ] == 'D' ) driver = TRUE

            }
        }
    }

    # Check the CNA mutations (can be dominant, recessive or nothing):
    ids  =  clone1$CNA_ID
    if ( ( ids != 0 )[ 1 ] ){
        for( id in ids ){
            gns  =  pck.env$cna_clones[[ id ]]$Gene_names
            mf  =  pck.env$cna_clones[[ id ]]$MalfunctionedByCNA
            nm  =  pck.env$cna_clones[[ id ]]$dupOrdel
            nm  =  paste0( 'role_', nm, collapse = '' )

            if ( mf & gene %in% gns & ( pck.env$DR[ il, nm ] != 'N' ) ){

                if ( pck.env$cna_clones[[ id ]]$Parental_1or2 == 1 ){
                    driver1  =  TRUE
                } else {
                    driver2  =  TRUE
                }

                if ( pck.env$DR[ il, nm ] == 'D' ) driver = TRUE

            }

        }
    }

    # Check combination of drivers 1 and 2:
    if ( driver1 & driver2 ) driver = TRUE

    return( driver )
}


.get_pom_info <- function( onco1, gm, ncol, mut_row1 ) {

   pom_info <- list()
   if        ( ncol( mut_row1 ) == ncol$base ) { # the base algorithm
      pm <- get_point_mutation( onco1, gm )
      
      pom_info$prntl = unlist( pm[[1]] )
      pom_info$gene  = unlist( pm[[2]] )
      pom_info$pos   = unlist( pm[[3]] )
      pom_info$Chr   = unlist( pm[[4]] )
      MalfunctionedByPointMut <- NA   # if NA, randomly generated later in generate_pnt
      
   } else if ( ncol( mut_row1 ) == ncol$Rother ) {
      pm <- get_point_mutation_for_gene( onco1, gm, mut_row1[ 1, 'Gene' ] )
     
      pom_info$prntl = unlist( pm[[1]] )
      pom_info$gene  = unlist( pm[[2]] )
      pom_info$pos   = unlist( pm[[3]] )
      pom_info$Chr   = unlist( pm[[4]] )
      MalfunctionedByPointMut <- F
      
   } else if ( ncol( mut_row1 ) == ncol$Rint ) {
      # type casting is needed for safe
      pom_info$prntl <- as.numeric(   mut_row1[ 1, 'Pchr' ] )
      pom_info$gene  <- as.character( mut_row1[ 1, 'Gene' ] )
      pom_info$pos   <- as.numeric(   mut_row1[ 1, 'Chr_stt' ] )
      pom_info$Chr   <- as.character( mut_row1[ 1, 'Chr' ] )
      MalfunctionedByPointMut <- ifelse( mut_row1[ 1, 'DrvPssRst' ] != 'pss', TRUE, FALSE )
      
   } else { 
      stop( 'Strange case!' )
   }

   pom_info$pnt0 <- generate_pnt( onco1, pom_info$prntl, pom_info$gene, 
                                  pom_info$pos, pom_info$Chr,
                                  malfunc = MalfunctionedByPointMut, 
                                  rst.ratio = mut_row1$rst.ratio )

   return( pom_info )
}


#
.get_cna_info <- function( onco1, gm, ncol, mut_row1 ) {
   type <- mut_row1[ 1, 'Type' ]
   
   cna_info <- list()
   if        ( ncol( mut_row1 ) == ncol$base ) { # the base algorithm
      cna_mut <- get_cna_mutation( onco1, dupOrdel = type, gm )
      
      cna_info$prntl     = unlist( cna_mut[[1]] )
      cna_info$Chr       = unlist( cna_mut[[2]] )
      cna_info$genes     = unlist( cna_mut[[3]] )
      cna_info$start_end = unlist( cna_mut[[4]] )
      MalfunctionedByCNA <- NA   # if NA, randomly generated later in generate_cna
      
   } else if ( ncol( mut_row1 ) == ncol$Rother ) {
      cna_mut <- get_cna_mutation_for_gene( onco1, dupOrdel = type, gm, mut_row1[ 1, 'Gene' ] )

      cna_info$prntl     = unlist( cna_mut[[1]] )
      cna_info$Chr       = unlist( cna_mut[[2]] )
      cna_info$genes     = unlist( cna_mut[[3]] )
      cna_info$start_end = unlist( cna_mut[[4]] )
      MalfunctionedByCNA <- F
      
   } else if ( ncol( mut_row1 ) == ncol$Rint ) {
      # type casting is needed for safe
      cna_info$prntl     <-       as.numeric( mut_row1[ 1, 'Pchr' ] )
      cna_info$Chr       <-     as.character( mut_row1[ 1, 'Chr' ] )
               genes     <-     as.character( mut_row1[ 1, 'Gene' ] )
      cna_info$genes     <- unlist( strsplit( genes, '[, ]+', perl=T ) )
      cna_info$start_end <-    as.numeric( c( mut_row1[ 1, 'Chr_stt' ], 
                                              mut_row1[ 1, 'Chr_end' ] ) )
      MalfunctionedByCNA <- ifelse( mut_row1[ 1, 'DrvPssRst' ]  != 'pss', TRUE, FALSE )
      
   } else { 
      stop( 'Strange case!' )
   }

   cna_info$cna0 = generate_cna( onco1, cna_info$prntl, cna_info$genes, 
                                 cna_info$start_end, cna_info$Chr, 
                                 malfunc = MalfunctionedByCNA, type, mut_row1$rst.ratio )

   return( cna_info )
}


# Functions related to mutations ------------------------------------
#' Function to change the point mutation due to CNA
#'
#' @param pnt1 Object of class 'Point_Mutations'
#' @param start_end Vector with initial and final positions of CNA
#' @param t 'dup' or 'del' for duplication or deletion respectively
#' @param CNA_ID_last The last ID of clone's CNA_ID which should change the point mutation
#'
#' @keywords internal
change_pnt_by_cna  <-  function( pnt1, start_end, t, CNA_ID_last ) {

    ntrs  =  intersect( which( pnt1$Phys_pos >= start_end[1] ),
                        which( pnt1$Phys_pos <= start_end[2] )  )
    cf    =  length( ntrs )  ### coefficient for CN

    if ( pnt1$Copy_number != 0 )  {
        pn_cn  =  pnt1$Copy_number + cf * ifelse(t == 'dup', 1, -1 )
        pnt1$Copy_number  =  ifelse( pn_cn >= 0, pn_cn, 0 )
    }

    pos_pnt       =   pnt1$Phys_pos[ ntrs ]

    dlt  =  ifelse( t == 'dup', start_end[2]  -  start_end[1], start_end[1]  -  start_end[2] )
    pnt1$Phys_pos   =  c( pnt1$Phys_pos, pos_pnt + dlt )
    pnt1$Delta      =    pnt1$Phys_pos  -  pnt1$Ref_pos

    if ( pnt1$Ovlp_CNA_ID == 0 ){
        pnt1$Ovlp_CNA_ID  =  CNA_ID_last
    } else {
        pnt1$Ovlp_CNA_ID  =  c( pnt1$Ovlp_CNA_ID, CNA_ID_last )
    }

}

#' Function to change copy number of the allele A of the point mutation at the allele B due to CNA
#'
#' @param pnt1 Object of class 'Point_Mutations'
#' @param start_end Vector with initial and final positions of CNA
#' @param t 'dup' or 'del' for duplication or deletion respectively
#' @param CNA_ID_last The last ID of clone's CNA_ID which should change the point mutation
#'
#' @keywords internal
change_allele_A_by_cna  <-  function( pnt1, start_end, t, CNA_ID_last ) {
    if ( pnt1$Copy_number != 0 )  {
        pn_cn  =  pnt1$Copy_number + ifelse(t == 'dup', 1, -1 )
        pnt1$Copy_number  =  ifelse( pn_cn >= 0, pn_cn, 0 )
    }
    if ( pnt1$Ovlp_CNA_ID == 0 ){
        pnt1$Ovlp_CNA_ID  =  CNA_ID_last
    } else {
        pnt1$Ovlp_CNA_ID  =  c( pnt1$Ovlp_CNA_ID, CNA_ID_last )
    }
}

#' Function to update onco1 after mutation (for usage in trial_mutagenesis() function)
#'
#' @param onco1 Object of class 'OncoGene'
#' @param gm data.frame gene_map
#'
#' @return onco1 with updated info
#'
#' @keywords internal
onco_update  <-  function( onco1, gm ){

    lst1  =  get_cds_rna( gm[[1]] )
    rd1   =  as.integer( unlist( sapply( onco1$name, FUN = function(x) which( x  ==  lst1[[1]]) ) ) )

    lst2  =  get_cds_rna( gm[[2]] )
    rd2   =  as.integer( unlist( sapply( onco1$name, FUN = function(x) which( x  ==  lst2[[1]]) ) ) )

    # change the onco1 related to new gene_map:
    onco1$cds_1  =  lst1$CDS[ rd1 ]
    onco1$cds_2  =  lst2$CDS[ rd2 ]

    onco1$rna_1  =  lst1$RNA[ rd1 ]
    onco1$rna_2  =  lst2$RNA[ rd2 ]

    onco1$prob_1 =  lst1$PROB
    onco1$prob_2 =  lst2$PROB

    onco1$sum_prob_1  =  lst1$SUM
    onco1$sum_prob_2  =  lst2$SUM

    onco1$p0_1   =  lst1$P0
    onco1$p0_2   =  lst2$P0

    return( onco1 )
}

#' Function to check point mutations match or don't match into duplication or deletion
#'
#' @param pnt1 Object of class 'Point_Mutations'
#' @param Ref_start Initial position of CNA
#' @param Ref_end Final position of CNA
#' @param Chr Chromosome name
#' @param prntl Parental chromosome 1 or 2
#'
#' @return Logical: TRUE if point mutation matches CNA, FALSE if it doesn't match
#'
#' @keywords internal
chk_pnt_mut  <-  function( pnt1 , Ref_start, Ref_end, Chr, prntl ){

    for( X in pnt1$Phys_pos ){
        if ( pnt1$Chr == Chr &  X <= Ref_end & X >= Ref_start & pnt1$Parental_1or2 == prntl ) {
            return( TRUE )}
    }
    return( FALSE )
}


#' Function to read file with initial clones
#'
#' @param clonefile File to read
#' @param n_genes Number of genes
#'
#' @return List of objects of class 'Clone
#'
#' @keywords internal
init_clones <- function( clonefile, n_genes ) {

   data = read.table( clonefile, sep="\t", header=T, stringsAsFactors = F )
   if ( nrow(data) == 0 ) stop( 'There is no data for the initial clones.' )
   data <- .check.init.clones( data )
   
   clone1 <- clone$new( gene_size = n_genes ) # empty clone
   clones <- c()
   for ( id in unique( data$cloneID ) ) {
      clone2 <- clone_copy( clone1 )
      n_cells <- unique( data[ data$cloneID == id, 'Ncells' ] )
      if ( length(n_cells) > 1 ) stop( 'Not unique Ncells in the initial clone file!' )
      clone2$N_cells <- as.numeric( n_cells )
      clone2$tmp     <- data[ data$cloneID == id, ]
      clones         <- c( clones, clone2 )
   }

   cnt <- 0
   for ( clone in clones ) {
      cnt <- cnt + 1
      clone$id = cnt
      clone$parent = 0
      clone$birthday = 0

      df <- clone$tmp; 
      # Driver mutations
      for ( gene1 in unique( df[ df$DrvPssRst != 'pss', ]$Gene ) ) { 
         drv1    <- df[ df$DrvPssRst != 'pss' & df$Gene == gene1 & 
                        as.character(df$Pchr) == '1', ]$Type
         drv2    <- df[ df$DrvPssRst != 'pss' & df$Gene == gene1 & 
                        as.character(df$Pchr) == '2', ]$Type
      
         ww <- which( pck.env$DR$gene == gene1 )
         r1 <- pck.env$DR[ ww, paste0('role_', drv1[1] ) ]
         r2 <- pck.env$DR[ ww, paste0('role_', drv2[1] ) ]
         driver <- F
         if        ( length( drv1 ) != 0 && r1 == 'R' && 
                     length( drv2 ) != 0 && r2 == 'R' ) {
            driver <- T
            
         } else if ( length( drv1 ) != 0 && r1 == 'D' ) {
            driver <- T
            
         } else if ( length( drv2 ) != 0 && r2 == 'D' ) {
            driver <- T
            
         }

         ifelse( driver, clone$gene[    which( pck.env$onco$name == gene1 ) ] <- 1, 
                         clone$pasgene[ which( pck.env$onco$name == gene1 ) ] <- 1 )
      }
      if ( sum( clone$gene ) > 0 ) clone$primary = TRUE
                
      # Passenger mutations
      for ( gene1 in unique( df[ df$DrvPssRst == 'pss', ]$Gene) ) { 
         clone$pasgene[ which( pck.env$onco$name == gene1 ) ] <- 1
      }
   }

   pck.env$env$last_id = length( clones )
   return( as.list( clones ) )
}

.check.init.clones <- function( data ){
    if ( ncol(data) == 10 ) return( data )
    if ( ncol(data) == 9 ) {
        data.new = data
        names(data.new)[ which( names(data.new) == 'DrvPss')] = 'DrvPssRst'
        data.new$Ratio = NA
        return( data.new )
    } else {
        stop( 'The number of columns in the initial clone file is not correct.' )
    }
    
}

#' Function to make list of onco_clones (object of class 'OncoGene') and generate initial onco settings for all clones
#'
#' @param onco1 Object of class 'OncoGene'
#' @param clones List of objects of class 'Clone'
#'
#' @keywords internal
init_onco_clones <- function( onco1, clones ) {
    onco_clones = NULL
    for ( i in 1:length( clones ) ) {
        clone1 = clones[[i]]
        onco_clone2 = onco_copy( onco1 )
        onco_clone2$id = clone1$id
        
        onco_clones = c( onco_clones, onco_clone2 )
    }

    return( as.list( onco_clones ) )
}


#' Function to generate point mutations and CNAs for initial clones
#'
#' @param clones List of objects of class 'Clone'
#' @param onco_clones List of objects of class 'OncoGene'
#'
#' @keywords internal
init_pnt_cna_clones   <- function( clones, onco_clones ) {
   if ( is.null( clones ) )  return( NULL )

   for( i in 1:length( clones ) ){
      clone1  =  clones[[ i ]]
      onco1   =  onco_clones[[ i ]]

      df <- clone1$tmp; clone1$tmp <- data.frame()
      if ( all( is.na( df[ df$DrvPssRst != 'pss' | 
                           df$DrvPssRst == 'pss', ] ) ) ) {
         warning( 'No mutations in the initial cells.' )
         return()
      }

      df$Mutation <- paste0( 'I', 1:nrow(df) )
      mut_info <- df[ , c( 'Mutation', 'Type', 'Gene', 
                           'Chr', 'Chr_stt', 'Chr_end', 
                           'Pchr', 'DrvPssRst', 'Ratio' ) ] # same as EF table
      insert_mutations( clone1, onco1, mut_info )
   }

}


