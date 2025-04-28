#!/usr/bin/Rscript --slave --vanilla


#####################################################
##### For development with source() #################
if (0) {

    # Debugging mode
    files <- Sys.glob( 'R/*.R' )
    files <- files[ ! grepl('R/test|R/_', files, perl=T) ]
    for (file in files) {
       print( paste0( "For development, loading: ", file ) )
       source(file)
    } 
    
} else {
    
    # User mode
    library(tugMedi)
}
#####################################################


write_evals( input = './Input/ForPosts/evals.parameters.post.txt' )


