

#####################################################
##### For development with source() #################
if (1) {
    # Debugging mode
    base_dir <- '.'
    files <- Sys.glob( sprintf( '%s/R/*.R', base_dir ) )
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
flog.debug( 'flog level in main: DEBUG' )
#####################################################
flog.info("Generators start")


#write_prm(       input = './Input/ForGenerators/parameters.parameters.txt')

#write_CF(        input = './Input/ForGenerators/CF.parameters.txt')

write_weights(   input = './Input/ForGenerators/hallmark_weights.parameters.txt' )

#write_cloneinit( input  = './Input/ForGenerators/cloneinit.parameters.txt' )

#write_Rint(      input = './Input/ForGenerators/Rint.parameters.txt' )

write_EF.Rint(   input = './Input/ForGenerators/EF.Rint.parameters.txt' )

# Only this depends on write_prm() result
write_EF.Rother( input = './Input/ForGenerators/EF.Rother.parameters.txt' )


##############################################################
flog.info("Simulation starts")


system.time(
start_simulation( input_files = 'Input/DATA/FILES.txt', 
                  saveTo      = 'Output/Results.sim.01.RDS', 
                  seed        = NA
                 )
)


write_VAF(   input = './Input/ForPosts/VAF.parameters.txt' )

write_TMB(   input = './Input/ForPosts/TMB.parameters.txt' )


