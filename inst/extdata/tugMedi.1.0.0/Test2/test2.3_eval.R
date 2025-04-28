#!/usr/bin/Rscript --slave --vanilla


#####################################################
##### For development with source() #################
if (1) {

    # Debugging mode
    files <- Sys.glob( 'R/*.R' )
    files <- files[ ! grepl('R/test|R/_', files, perl=T) ]
    for (file in files) {
       print( paste0( "For development, loading: ", file ) )
       source(file)
    } 
    
} else {
    
    # User mode
    library(tugHall.3)
}
#####################################################


write_evals( input = './Input/ForPosts/evals.parameters.pre.txt' )


