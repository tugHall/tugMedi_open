#' Simulation for start with parameters from Input folder
#'
#' @description \code{start_simulation()} makes a simulation with parameters from Input folder
#' and results save in pck.env as well in some RDS file in \code{work_dir} folder
#'
#' @param first Logical for the first simulation or continue simulation
#' @param input_files The path to the file with input files' names, for example 'Input/DATA/FILES.txt'
#' @param saveTo The RDS file name to save results of simulation, for example 'Output/Results.sim.02.RDS'
#' @param seed Numeric type to set seed for a simulation, if seed = NA then it will be skipped
#' @param loadFrom The RDS file name to load results of the previous simulation if first = F
#' @param change_parameters The list of parameters to change for a simulation
#' @param change_files The list of files to change for a simulation
#' @param drug_int_param Parameters for drug intervention trials
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' NULL
start_simulation <- function( first = T, 
                              input_files, 
                              saveTo, 
                              seed = NA,
                              # only work if first = FALSE
                              loadFrom = NULL,
                              change_parameters = NULL,
                              change_files      = NULL,
                              drug_int_param = NULL
                             ) {
    check_packages()
    files <- read_files4sim( input_files, check_output = first )
    if ( file.exists( files$out.clone ) ) flog.info( "The cloneout.txt file already exists, \n make sure that this is according to a simulation target!" )
    if ( !is.na( seed ) ) set.seed( seed )
    model_param <- read_parameters( files$in.prm, files$in.fixed_prm )

    if ( first == T ) {
        clear_tugHall.Environment()
        
        define_parameters( model_param ) # for pck.env
        define_compaction_factor( files$in.cf )
        define_regions( files$in.Rint, files$in.Rother )
        define_gene_location( pck.env$ls_genes, files$in.CDS, files$out.CDS, model_param$save.CDS.info )
    } else {
        local_environment( env = pck.env )
        load_tugHall.Environment( results = readRDS( file = loadFrom ) )
        
        lapply( X = names(change_parameters),
           function(X) pck.env[[ X ]] <- change_parameters[[ X ]] )
        files[ names(change_files) ] <- change_files[ names(change_files) ]
    }

    forNext <- 
       do_simulation( first, files, model_param, drug_int_param )
    pck.env$clones      <- forNext$clones
    pck.env$onco_clones <- forNext$oncos
    pck.env$ttl_divs    <- forNext$ttl_divs

    write_pnt_clones( pck.env$pnt_clones, file_out = files$out.poms, file_outA = files$out.pomsA )
    write_cna_clones( pck.env$cna_clones, file_out = files$out.CNAs )

    res  =  get_tugHall.Environment()
    saveRDS( object = res, file = saveTo )
}


#' Main function to simulate clones' evolution
#' 
#' @param first  Logical parameter if it is the preliminary simulation or continuing one
#' @param model_param List of model's parameters to use or change
#' @param drug_int_param List of parameter for drug intervention trial
#'
#' @keywords internal
do_simulation <- function( first, files, model_param, drug_int_param ){

    local_environment( env = pck.env )

    if ( first ) {
       pck.env$onco = oncogene$new()        # make the vector onco about the hallmarks
       pck.env$onco$read( files$in.Rint )
       pck.env$hall = hallmark$new()        # make a vector hall with hallmarks parameters
       tumblers = list( apoptosis       = pck.env$tumbler_for_apoptosis_trial, 
                        angiogenesis    = pck.env$tumbler_for_angiogenesis_trial, 
                        immortalization = pck.env$tumbler_for_immortalization_trial, 
                        metastasis      = pck.env$tumbler_for_metastasis_trial )
       
       pck.env$hall$read( file          = files$in.weights, 
                          normalization = TRUE, 
                          tumblers      = tumblers, 
                          compaction_factor = pck.env$compaction_factor )
       pck.env$env = environ$new(Fb)        # new vector for average values of cells
       pck.env$pnt_clones = NULL
       pck.env$cna_clones = NULL
       pck.env$mut_order  =  0
       # Names of the driver genes involved into hallmarks
       w_drivers  =  unique( c( pck.env$hall$Ha, pck.env$hall$Hb, pck.env$hall$Hd, pck.env$hall$Hi, pck.env$hall$Him ) )
       w_drivers  =  w_drivers[ !is.na( w_drivers ) ]
       pck.env$DR = define_dominant_recessive( files$in.genetic_type, 
                                               pck.env$onco$name[ w_drivers ], 
                                               pck.env$onco$onsp[ w_drivers ] )
       # Define the names of genes in the Rother region:
       pck.env$Rother_genes  =  setdiff(x = pck.env$onco$name, y = pck.env$DR$gene)  

       # initial settings
       clones      = init_clones( files$in.clone, length( pck.env$onco$cds_1 ) )
       onco_clones = init_onco_clones( pck.env$onco, clones )
                     init_pnt_cna_clones( clones, onco_clones ) # pom and CNA
    }
    foolproof()
    write_parameters( files$out.prm, pck.env, model_param ) 
    write_geneout( files$out.weights, pck.env$hall, compaction_factor, pck.env$CF )
    write_genetic_type( files$out.genetic_type, pck.env$DR )

    if ( model_param$growth == 'exponential' )
       lapply( X=clones, function( X ) { X$invasion = T; X$im = NaN } )

    cells_number <- sum_N_P_M( pck.env$env, clones )  # to calculate cells numbers for initial cells
    # update Hallmarks and probabilities
    if ( length( clones ) > 0 ){
        lapply( 1:length( clones ), FUN = function( id ) {
            update_Hallmarks( clone1 = clones[[ id ]],
                              onco1  = onco_clones[[ id ]] )  
        } )
    }
    pck.env$hall$updateEnviron( pck.env$env, clones ) # to average probabilities and hallmarks

    if ( first ) write_cloneout( files$out.clone, pck.env$env, clones, onco_clones )  #  write initial clones
    if ( control.monitor ) {
       monitor_Ncells.Nmuts( outfile = files$out.monitor.N, start = first,
                             env = env, clones = clones, 
                             cna_clones = pck.env$cna_clones,
                             tumbler_for_event_enforcement )
       prev.pck <- monitor_pck.env( files$out.monitor.pck, first,
                                    pck.env$env$T,
                                    prev.obj = NULL, obj = ls( envir=pck.env ) )
    }
    
    # -------------------------------------------------

    time_start  =  Sys.time()

    # Counter for waiting divisions
    if ( first ) ttl_divs <- 0

    flog.info( paste( "Start of simulation time step:", pck.env$env$T ) )
    EF.Rint.data   <- .get_EF_data( files$in.EF.Rint,   tumbler_for_event_enforcement )
    EF.Rother.data <- .get_EF_data( files$in.EF.Rother, tumbler_for_event_enforcement )
    while( 0                < length( clones )
           && cells_number  < control.censor_cell_number
           && pck.env$env$T < control.censor_time_step
           && ( as.numeric( difftime( Sys.time(), time_start, units = 'secs') ) < control.censor_real_time ) ) {
        pck.env$env$T = pck.env$env$T + 1

        # Trial
        if ( ! first ) {
            drug_c_o    <- .druggalize_clones_oncos( clones, onco_clones, drug_int_param )
            clones      <- drug_c_o$clones
            onco_clones <- drug_c_o$oncos
        }
        if ( length( clones ) != length( onco_clones ) ) stop( "Lengths of clones and oncos differ!" )
        N_divisions = unlist( mapply( conduct_trials,
                                      clone1 = clones, onco1 = onco_clones,
                                      MoreArgs = list( model_param = model_param,
                                                       drug_int_param = drug_int_param ) ) )
        if ( ! first ) {
            restored_c_o <- .restore_druggalized_clones_oncos( clones, onco_clones, drug_int_param )
            clones       <- restored_c_o$clones
            onco_clones  <- restored_c_o$oncos
        }

        # new clones by invasion / metastasis
        newIMc_o    <- .make_IM_clones_oncos( clones = clones, oncos = onco_clones,
                                              model_param = model_param )
        IMc_o       <- .initialize_IM_states( clones, onco_clones )
        clones      <- c( IMc_o$clones, newIMc_o$clones )
        onco_clones <- c( IMc_o$oncos,  newIMc_o$oncos )

        # new clones by mutations
        if ( ! tumbler_for_event_enforcement ) {
           # mutation insertion by the base algorithm
           warning( "tumbler_for_event_enforcement == F, ", 
                    "the base algorithm is not kept validated!" )
           newBSc_o    <- .make_BS_clones_oncos( N_divisions, model_param$trial$spec,
                                                 clones = clones, oncos = onco_clones )
           clones      <- c( clones,      newBSc_o$clones )
           onco_clones <- c( onco_clones, newBSc_o$oncos )
        } else {
           # mutation insertion by EF
           newRIc_o    <- .make_EF_clones_oncos( N_divisions, pck.env$env$T, EF.Rint.data,
                                                 clones = clones, oncos = onco_clones )
           clones      <- c( clones,      newRIc_o$clones )
           onco_clones <- c( onco_clones, newRIc_o$oncos )

           # by "EEL"
           newROc_o    <- .make_EF_clones_oncos( N_divisions, ttl_divs,      EF.Rother.data,
                                                 clones = clones, oncos = onco_clones )
           clones      <- c( clones,      newROc_o$clones )
           onco_clones <- c( onco_clones, newROc_o$oncos )
        }

        # Updates ----
        ttl_divs    <- ttl_divs + sum( N_divisions )
        TorF        <- sapply( X=clones, function(X) X$N_cells > 0 )
        clones      <- clones[ TorF ]
        onco_clones <- onco_clones[ TorF ]

        cells_number <- sum_N_P_M( pck.env$env, clones )  # to calculate cells numbers for next step
        # update Hallmarks and probabilities
        if ( length( clones ) > 0 ){
            lapply( 1:length( clones ), FUN = function( id ) {
                update_Hallmarks( clone1 = clones[[ id ]],
                                  onco1  = onco_clones[[ id ]] )  
            } )
        }
        pck.env$hall$updateEnviron( pck.env$env, clones ) # to average probabilities and hallmarks

        write_cloneout( files$out.clone, pck.env$env, clones, onco_clones )
        if ( control.monitor ) {
           monitor_Ncells.Nmuts( outfile = files$out.monitor.N, start = FALSE,
                                 env = pck.env$env, clones = clones,
                                 cna_clones = pck.env$cna_clones,
                                 tumbler_for_event_enforcement )
           prev.pck <- monitor_pck.env( files$out.monitor.pck, first,
                                        pck.env$env$T,
                                        prev.pck, obj = ls( envir=pck.env ) )
        }

        if ( pck.env$env$T %% 20 == 0 ) 
           flog.info( paste( "Updated simulation time step:", pck.env$env$T ) )
    }

    .write.log( fl = files$out.log, l.clones = length( clones ), time.step = pck.env$env$T, 
               cells_number = ifelse( length( cells_number) > 0 , cells_number, 0), 
               df.time = as.numeric( difftime( Sys.time(), time_start, units = 'secs') ), 
               censors = as.list( c( control.censor_cell_number = control.censor_cell_number, 
                                     control.censor_time_step   = control.censor_time_step, 
                                     control.censor_real_time   = control.censor_real_time ) )
               )
    flog.info( paste( "End of simulation time step:", pck.env$env$T ) )
    return( list( clones = clones, oncos = onco_clones, ttl_divs = ttl_divs ) )
}


#####   invasion/metastasis   ##################################################
#
.make_IM_clones_oncos <- function( clones, oncos, model_param ) {

   new_clones <- c()
   new_oncos  <- c()
   for ( idx in seq_along( clones ) ) {
      clone1 <- clones[[ idx ]]
      onco1  <- oncos[[  idx ]]
      if ( ! clone1$invasion || is.nan( clone1$im ) ) next

      suc <- clone1$invasion.suc
      while ( suc > 0 ) {
         clone2   <- clone_copy( clone1 )
         onco2    <- onco_copy(  onco1 )
         onco2$id <- clone2$id

         clone2$N_cells <- 1 + rtrunc( 1, spec = model_param$meta.addNcells$spec,
                                       a = -Inf,
                                       b = suc - 1,
                                       lambda = model_param$meta.addNcells$lambda )
         clone2$invasion     <- T
         clone2$invasion.suc <- NaN
         clone2$im           <- NaN   # transformed into metastatic one

         clone2$location.meta[1, ] <- c( round( mvrnorm( 1, mu = rep(0, 3),
                                         Sig = model_param$meta.loc$sig ) ),
                                         clone2$N_cells )

         new_clones <- c( new_clones, clone2 )
         new_oncos  <- c( new_oncos,  onco2 )
         suc <- suc - clone2$N_cells
      }
   }

   return( list( 'clones' = new_clones, 'oncos' = new_oncos ) )
}


#
.initialize_IM_states <- function( clones, oncos ) {

   bag  <- c()
   sack <- c()
   for ( idx in seq_along(clones) ) {
      clone <- clones[[ idx ]]
      onco  <- oncos[[ idx ]]

      if ( ! clone$invasion || is.nan( clone$im ) ) {
         bag  <- c( bag,  clone )
         sack <- c( sack, onco )
         next
      }

      # Initialization. N_cells of <= 0 is updated later altogether
      clone$N_cells <- clone$N_cells - clone$invasion.suc
      clone$invasion.suc <- 0
      clone$invasion     <- F

      bag  <- c( bag,  clone )
      sack <- c( sack, onco )
   }

   return( list( 'clones' = bag, 'oncos' = sack ) )
}


#####   base   ##########################################################
#
.make_BS_clones_oncos <- function ( N_divisions, spec, clones, oncos ) {
   if ( sum( N_divisions ) == 0 ) return( NULL )

   new_clones <- vector()
   new_oncos  <- vector()
   for ( idx in seq_along( clones ) ) {
      if ( is.na( N_divisions[ clones[[ idx ]]$id ] ) ||
                  N_divisions[ clones[[ idx ]]$id ] == 0 ) next

      clone1 <- clones[[ idx ]]
      onco1  <- oncos[[ idx ]]

      N_new_clones <- .get_rtrunc4trial( spec = spec,
                                         max  = 2 * N_divisions[ clones[[ idx ]]$id ],
                                         prob = 1 - (onco1$p0_1 + onco1$p0_2 ) / 2 ) # pchr 1 and 2
      if ( N_new_clones == 0 ) next

      for ( ii in 1:N_new_clones )  {
         new_clone1   <- clone_copy( clone1 )
         new_onco1    <- onco_copy( onco1 )
         new_onco1$id <- new_clone1$id

         # insert mutations
         mut_info <- .make_BS_mut_info( new_onco1 )
         insert_mutations( new_clone1, new_onco1, mut_info )

         new_clones <- c( new_clones, new_clone1 )
         new_oncos  <- c( new_oncos,  new_onco1 )
      }
      clone1$N_cells <- clone1$N_cells - N_new_clones
   }

   return( list( 'clones' = new_clones, 'oncos' = new_oncos ) )
}


#
.make_BS_mut_info <- function( onco1 ) {

   sm <- unlist( ( onco1$sum_prob_1 + onco1$sum_prob_2 ) / 2 )
   num_mut <- rtrunc( 1, spec = 'pois', a = 0, lambda = sm ) # >a (>0)

   types = sample( c('pom', 'del', 'dup' ), size = num_mut,
                   replace = TRUE, prob = (onco1$prob_1 + onco1$prob_2) / 2 )
   mut_info <- data.frame( 'Type' = types )

   return( mut_info )
}


#####   EF   #######################################################
#
.get_EF_data <- function( EF.table, tumbler ) {
   if ( !tumbler ) return( NULL )

   ncol.Rint   <- 12
   ncol.Rother <-  4

   tbl = read.table( file = EF.table, header = TRUE, sep = '\t',
                     stringsAsFactors = FALSE )
   tbl <- .check.Rint.data( tbl, ncol.Rint, ncol.Rother )
   
   if        ( ncol(tbl) == ncol.Rint ) {
      EF.data <- as.data.frame( tbl,
                                colClasses=c( 'character', 'numeric',
                                          rep('character', 5),
                                           rep('numeric', 3), 'character', 'numeric' ) )
   } else if ( ncol(tbl) == ncol.Rother ) {
      EF.data <- as.data.frame( tbl,
                                colClasses=c( 'numeric',
                                          rep('character', 3) ) )
   } else {
      stop( 'Strange case!' )
   }

   return( EF.data )
}


# 
.check.Rint.data  <- function( tbl, ncol.Rint, ncol.Rother ){

   if ( ncol(tbl) == ncol.Rint |  ncol(tbl) == ncol.Rother ) {
       return( tbl )
       } else {
          tbl.new = tbl
          names( tbl.new )[c(1, 10)] <- c( 'When', 'DrvPssRst' )
          tbl.new$TimeDiv <- rep( 'time', nrow(tbl.new) )
          tbl.new$Ratio   <- rep(  NA,    nrow(tbl.new) )
          tbl.new <- tbl.new[ , c(11, 1:10, 12 ) ]
       }

   return( tbl.new )

}

#
.make_EF_clones_oncos <- function( Ndivs, base_counter, EF.data, clones, oncos ) {
   if ( is.null( EF.data ) ) return( NULL )
   switch <- .get_switch( EF.data )

   if        ( switch == 'EF' ) {
      row4slctd_clone_idx <- .get_row4slctd_clone_idx.EF(         base_counter, EF.data, clones )
   } else if ( switch == 'EEL' ) {
      row4slctd_clone_idx <- .get_row4slctd_clone_idx.EEL( Ndivs, base_counter, EF.data, clones )
   }
   if ( is.null( row4slctd_clone_idx ) || length( row4slctd_clone_idx ) == 0 ) return( NULL )

   n_mut4new_clone_idx <- lapply( X = names( row4slctd_clone_idx ), function(X) {
      store <- vector()
      if        ( switch == 'EF' ) {
         cnt <- nrow( row4slctd_clone_idx[[ X ]] )
         while ( cnt > 0 ) {
             # multiple muts may occur in a single cell, rcnt bound >= 1
             rcnt <- rtrunc( 1, spec="pois", a=0, b=cnt,
                             lambda = cnt / clones[[ as.numeric(X) ]]$N_cells )
             store <- c( store, rcnt )
             cnt <- cnt - rcnt
         }
      } else if ( switch == 'EEL' ) {
         count <- table(  row4slctd_clone_idx[[ X ]]$Waiting_division )
         uniq  <- unique( row4slctd_clone_idx[[ X ]]$Waiting_division )
         store <- count[ as.character(uniq) ]
      }
      store
   } )
   names( n_mut4new_clone_idx ) <- names( row4slctd_clone_idx ) # hash


   # insert mutations
   new_clones <- c()
   new_oncos  <- c()
   for ( X in names( n_mut4new_clone_idx ) ) {
      mut_info <- .make_EF_mut_info( row4slctd_clone_idx[[ as.character(X) ]] )

      for ( n_mut in n_mut4new_clone_idx[[ as.character(X) ]] ) {
         clone1   <- clone_copy( clones[[ as.numeric(X) ]] )
         onco1    <- onco_copy( oncos[[ as.numeric(X) ]] )
         onco1$id <- clone1$id

         if        ( switch == 'EF' && clone1$EF.mut_IDs[1] == '' ) {
           clone1$EF.mut_IDs <-                       mut_info[ 1:n_mut, 'Mutation' ]
         } else if ( switch == 'EF' && clone1$EF.mut_IDs[1] != '' ) {
           clone1$EF.mut_IDs <- c( clone1$EF.mut_IDs, mut_info[ 1:n_mut, 'Mutation' ] )
         }

         insert_mutations( clone1, onco1, mut_info = mut_info[ 1:n_mut, ] )
         mut_info <- mut_info[ -(1:n_mut), ]

         new_clones <- c( new_clones, clone1 )
         new_oncos  <- c( new_oncos,  onco1 )
         clones[[ as.numeric(X) ]]$N_cells <-
         clones[[ as.numeric(X) ]]$N_cells - 1 # one new clone
      }
   }

   return( list( 'clones' = new_clones, 'oncos' = new_oncos ) )
}


#
.get_switch <- function( EF.data ) {

   if        ( ! is.null( EF.data$TimeDiv ) ) {
      switch <- 'EF'
   } else if ( ! is.null( EF.data$Waiting_division ) ) {
      switch <- 'EEL'
   } else {
      stop( "Strange case!" )
   }

   return( switch )
}


#
.get_row4slctd_clone_idx.EF <- function( base_counter, EF.data, clones ) {

   slctd_row <- EF.data[ EF.data$When == base_counter, ]
   if ( nrow( slctd_row ) == 0 ) return( NULL )

   row4slctd_clone_idx <- list()
   for ( ii in 1:nrow(slctd_row) ) {
      row1 <- slctd_row[ ii, ]

      if ( row1$Condition == 'NA' || is.na( row1$Condition ) ) {
         slct <- seq_along( clones )
      } else {
         con.muts <- unlist( strsplit( row1[ , 'Condition' ], '[, ]+', perl=T ) )
         slct <- sapply( X=seq_along(clones), function(X) { all( con.muts %in% clones[[ X ]]$EF.mut_IDs ) } )
         slct <- which( slct == T )

         if ( length( slct ) == 0 ) {
             stop( paste('ERROR: for', row1[,'Mutation'],
                   ' no clones have mutations', con.muts, 'at time', base_counter, 'here.', sep=' ') )
         }
      }

      #If the value of all prob is 0, an error occurs in the sample function, so the next line is executed
      prob = sapply( X=slct, function( X ) { clones[[ X ]]$N_cells * clones[[ X ]]$d } )
      if ( all(prob == 0) ) next

      # Forced division with mut, according to clone sizes * div rates
      slct1 <- sample( as.character(slct), 1,
                       prob = sapply( X=slct, function( X )
                          { clones[[ X ]]$N_cells * clones[[ X ]]$d } ) )

      # hash
             row4slctd_clone_idx[[ as.character(slct1) ]] <-
      rbind( row4slctd_clone_idx[[ as.character(slct1) ]], row1 )
   }

   # If the length of row4slctd_clone_idx is 0, return NULL 
   if ( length(row4slctd_clone_idx) == 0 ) {
      return( NULL )
   } 
   return( row4slctd_clone_idx )
}


#
.get_row4slctd_clone_idx.EEL <- function( Ndivs, base_counter, EF.data, clones ) {

   slctd_row <- EF.data[ base_counter < EF.data$Waiting_division &
                         EF.data$Waiting_division <= base_counter + sum( Ndivs ), ]
   if ( nrow( slctd_row ) == 0 ) return( NULL )

   u_wd  <- unique( slctd_row$Waiting_division )
   idx <- which( sapply( X=seq_along( clones ), function(X) {
                    ifelse( is.na( Ndivs[ as.character(clones[[ X ]]$id) ] ) ||
                                   Ndivs[ as.character(clones[[ X ]]$id) ] == 0, F, T ) } ) )

   if ( sum( Ndivs ) < 1e3 ) { # rep = F for small number
      idx_rep <- sapply( X=idx, function(X) { rep( X, times=Ndivs[ as.character(clones[[ X ]]$id) ] ) } )
      slctd_clone_idx <- sample( unlist( idx_rep ), length( u_wd ), rep=F )
   } else {
      slctd_clone_idx <- sample( idx,               length( u_wd ), rep=T,
                                 prob = sapply( X=idx, function(X)
                                    Ndivs[ as.character(clones[[ X ]]$id) ] ) )
   }

   row4slctd_clone_idx <- list()
   for ( ii in seq_along( u_wd ) ) {
      df <- slctd_row[ slctd_row$Waiting_division == u_wd[ ii ], ]
      slct1 <- slctd_clone_idx[ ii ]

      # hash
             row4slctd_clone_idx[[ as.character(slct1) ]] <-
      rbind( row4slctd_clone_idx[[ as.character(slct1) ]], df )
   }

   return( row4slctd_clone_idx )
}


#
.make_EF_mut_info <- function( tbl ) {

   ncol.Rint   <- 12
   ncol.Rother <-  4

   if        ( ncol(tbl) == ncol.Rint ) {
      mut_info <- tbl[ , c( 'Mutation', 'Type', 'Gene',
                            'Chr', 'Chr_stt', 'Chr_end',
                            'Pchr', 'DrvPssRst', 'Ratio' ) ]

   } else if ( ncol(tbl) == ncol.Rother ) {
      mut_info <- tbl[ , c( 'Mutation', 'Type', 'Gene' ) ]

   } else {
      stop( 'Strange case!' )
   }

   return( mut_info )
}


#####   drug intervention   ####################################################
#
.druggalize_clones_oncos <- function( clones, onco_clones, drug_int_param ) {
    if ( is.null(drug_int_param) ) return( list( 'clones' = clones, 'oncos' = onco_clones ) )

    slct <- sapply( X=clones, function(X)
                    { names(X$gene) <- pck.env$onco$name
                      any( X$gene[ drug_int_param[[ 'gene' ]] ] == 1 )
                     } )
    if ( any( is.na(slct) ) ) stop( paste0(drug_int_param[[ 'gene' ]], sep=","), " are out of genes." )

    drugged <- clones[ slct ]
    if ( length( drugged ) == 0 ) return( list( 'clones' = clones, 'oncos' = onco_clones ) )

    N_block <- sapply( X=drugged, function(X) {
        rtrunc( 1, spec="pois", a=-Inf, b=X$N_cells,
                lambda = X$N_cells * drug_int_param[[ 'block_prob' ]] * pck.env$kappa ) } )

    # Drugged
    for ( X in seq_along( drugged ) ) {
        drugged[[ X ]]$N_cells <- drugged[[ X ]]$N_cells - N_block[ X ]
        drugged[[ X ]]$drug.on <- TRUE
    }

    # Relation of drugged with blocked
    blocked <- lapply( X=drugged, function(X) {
        copied <- clone_copy( X );
        X$drug.blocked.child <- copied$id;
        copied; } )

    blocked <- .complete_blocked_clones( N_block, blocked, drug_int_param )

    # Blocked oncos
    blocked_oncos = lapply( X=onco_clones[ slct ], function( X ) onco_copy( X ) )
    sapply( X=seq_along( blocked_oncos ), function( X ) blocked_oncos[[ X ]]$id = blocked[[ X ]]$id )

    # Object update
    clones[ slct ] <- drugged
    clones      = c( clones,      blocked )
    onco_clones = c( onco_clones, blocked_oncos )

    # Update
    if ( length( clones ) > 0 ){
        lapply( 1:length( clones ), FUN = function( id ) {
            update_Hallmarks( clone1 = clones[[ id ]],
                              onco1  = onco_clones[[ id ]] )  
        } )
    }
    pck.env$hall$updateEnviron( pck.env$env, clones )

    return( list( 'clones' = clones, 'oncos' = onco_clones ) )
}


#
.complete_blocked_clones <- function( N_block, blocked, drug_int_param ) {

    # N_block
    sapply( X=seq_along( N_block ), function(X) { blocked[[ X ]]$N_cells = N_block[ X ] } )

    # Drivers into passengers, and flag
    sapply( X=blocked, function(X) {
        slct <- names(X$gene) %in% drug_int_param[[ 'gene' ]];
        X$gene[ slct ]    <- 0;
        X$pasgene[ slct ] <- 1;

        X$drug.blocked.on = TRUE; } )

    # hash
    names(blocked) <- sapply( X=blocked, function( X ) X$id )

    return( blocked )
}


#
.restore_druggalized_clones_oncos <- function( clones, onco_clones, drug_int_param ) {
    if ( is.null(drug_int_param) ) return( list( 'clones' = clones, 'oncos' = onco_clones ) )

    # merge blocked into unblocked
    for ( X in seq_along(clones) ) {
        # unblocked
        if ( clones[[ X ]]$drug.on && ( clones[[ X ]]$drug.blocked.on == F ) ) {
            id <- clones[[ X ]]$drug.blocked.child;
            clones[[ X ]]$N_cells = clones[[ X ]]$N_cells + clones[[ as.character(id) ]]$N_cells;

            clones[[ X ]]$drug.on            <- FALSE;
            clones[[ X ]]$drug.blocked.on    <- FALSE;
            clones[[ X ]]$drug.blocked.child <- -Inf;
        }
    }

    # delete blocked
    TorF        <- !sapply( X=clones, function(X) X$drug.on && ( X$drug.blocked.on == T ) )
    clones      <- clones[ TorF ]
    onco_clones <- onco_clones[ TorF ]

    # Update
    if ( length( clones ) > 0 ){
        lapply( 1:length( clones ), FUN = function( id ) {
            update_Hallmarks( clone1 = clones[[ id ]],
                              onco1  = onco_clones[[ id ]] )  
        } )
    }
    pck.env$hall$updateEnviron( pck.env$env, clones )

    return( list( 'clones' = clones, 'oncos' = onco_clones ) )
}


