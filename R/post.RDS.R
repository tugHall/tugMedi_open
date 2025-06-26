

### Functions to get input data from output files

.update_clone  <-  function( clone1, row_data, rds ){
    
    if ( FALSE ){
            nm_cl = names( clone1 )[ !startsWith( names( clone1 ), '.' )  ]
            names( row_data )
            
            for( nm in names( row_data )[3 : ncol(row_data) ] ){
                print( paste0( nm, '  ', ( nm %in% nm_cl )  )  )
            }
            as.integer( unlist( strsplit(x = '1, 2, 45, 67', split = ', ') ) )
    }
    
    name_row_data  = c( 'ID', 'N_cells', 'ParentID', 'Birth_time', 'ct', 'd', 'i', 'im', 'a',
                        'k', 'K', 'Nmax', 'Ha', 'Him', 'Hi', 'Hd', 'Hb', 'mutden',
                        'genesDysfunc', 'genesWmutsNoDys', 'PointMut_ID', 'CNA_ID' ) 
    
    name_clone     = c( 'id', 'N_cells', 'parent', 'birthday', 'c', 'd', 'i', 'im', 'a', 
                        'k', 'E', 'Nmax', 'Ha', 'Him', 'Hi', 'Hd', 'Hb', 'mutden',
                        'gene', 'pasgene', 'PointMut_ID', 'CNA_ID' )
    
    # names to use non-standard/unique treatment
    nonstd  =  c( 'gene', 'pasgene', 'PointMut_ID', 'CNA_ID' )
    
    for ( i in 1:length( name_row_data ) ){
        nm_row = name_row_data[ i ]
        nm_cl  = name_clone[ i ]
        if ( nm_cl == 'id' ){
            if ( get( nm_cl, clone1 ) !=  row_data[1,nm_row] ){
                print( paste0( 'The clone N ', row_data[1,nm_row], ' in row data and clone1 have different IDs.') )
            }
        } else {
            
            if ( ! nm_cl %in% nonstd ){
                vl = row_data[ 1, nm_row]
                clone1$field( nm_cl, vl )
            } else{
                
                if ( nm_cl == 'gene' | nm_cl == 'pasgene' ){
                    vl = row_data[1, nm_row]
                    vl = unlist( strsplit(x = vl, split = ', ') )
                    vl = as.integer( rds$onco$name %in% vl )
                    clone1$field( nm_cl, vl )
                }
                
                if ( nm_cl == 'PointMut_ID' | nm_cl == 'CNA_ID' ){
                    vl = row_data[1, nm_row ]
                    vl = as.numeric( unlist( strsplit(x = vl, split = ', ') ) )
                    clone1$field( nm_cl, vl )
                }
                
            }
        }
        
    }
    
    return( clone1 )
}

#' Update RDS file from Output/cloneout.txt file
#'
#' @param in.file Name of file with list of input and output files of a simulation
#' @param rds.file name of RDS file with simulation data
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
update_RDS_from_output  <-  function( in.file, rds.file ){
    
    # Get files' names:
    fls = read_files4sim( in.file, check_output = FALSE )
    
    # get data from files
    dataALL  =  get_flow_data( fls$out.clone , fls$out.poms, fls$out.pomsA, fls$out.CNAs, time = NA )
    time_max  = dataALL$time_max
    data_last = dataALL$data_last
    pnt_mut_A = dataALL$pnt_mut_A
    pnt_mut_B = dataALL$pnt_mut
    cna_mut   = dataALL$cna_mut
    
    rds       = readRDS(file = rds.file)
    # Rename file to save updated with the same name
    new_rds.file = paste0( substr(x = rds.file, start = 1, stop = nchar( rds.file) - 4), '_initial.RDS'  )
    
    # Update the clones in rds 
    if ( length( rds$clones) >0 ){
        for( i_cl in 1:length( rds$clones) ){
            w_id = which( data_last[, 'ID'] == rds$clones[[i_cl]]$id )
            row_data = data_last[ w_id , ]
            if ( nrow( row_data ) > 1 ){
                print( 'The rows have the same ID : ' )
                print( row_data )
                print('Only the first row of data will be used.')
                row_data = row_data[1,]
            }
            if ( nrow( row_data ) < 1 ){
                print( paste0( 'There is no data in the file on clone ID = ',  rds$clones[[i_cl]]$id ) )
            }
            rds$clones[[ i_cl ]] = .update_clone( clone1 = rds$clones[[ i_cl ]], 
                                                row_data = row_data, rds = rds )
        }
    }
    
    # Update Hallmarks
    
    file.rename( rds.file, new_rds.file )
    saveRDS(object = rds, file = rds.file )
    
    return( print('File updated') )
}


# TEST
if ( FALSE){
    
    ### TRIALS
    cloneoutfile = './Test1/Output/cloneout.txt'
    pom.file     = './Test1/Output/Mutations/pointMutations_B.txt'
    pomA.file    = './Test1/Output/Mutations/pointMutations_A.txt'
    cna.file     = './Test1/Output/Mutations/CNAs.txt'
    rds.file     = './Test1/Output/Results.sim.01.RDS'
    
    update_RDS_from_output(cloneoutfile = cloneoutfile, pom.file = pom.file, 
                           pomA.file = pomA.file, cna.file = cna.file, rds.file = rds.file )
    
    
}
