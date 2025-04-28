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
#flog.threshold( 'DEBUG' )
#flog.debug( 'flog level in main: DEBUG' )
#####################################################
print("Simulation starts")

system.time(
start_simulation( input_files = 'Input/DATA/FILES.txt', 
                  saveTo      = 'Output/Results.sim.01.RDS', 
                  seed        = NA
                 )
)

write_VAF( input = './Input/ForPosts/VAF.parameters.txt' )

write_TMB( input = './Input/ForPosts/TMB.parameters.txt' )


