

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
    library(tugMedi)
}
#####################################################
#flog.threshold( 'DEBUG' )
flog.debug( 'flog level in main: DEBUG' )
#####################################################
flog.info("Drug Intervention starts")


time_max_add   <-  5
drug_int_param <- list( 'kill_prob' = 0.5,
                        'block_prob'= 1.0,
                        'gene' = c('KRAS', 'APC') ) # Any of them is a target

system.time(
start_simulation( first       = F,
                  input_files = 'Input/DATA/FILES.txt', 
                  loadFrom    = 'Output/Results.sim.01.RDS', 
                  saveTo      = 'Output/Results.sim.02.RDS', 
                  seed        = NA,
                  change_parameters =
                     list( control.censor_cell_number = 1e15,
                           control.censor_time_step   = pck.env$env$T + time_max_add ),
                  drug_int_param = drug_int_param
#                  drug_int_param = NULL
                 )
)

write_realTime_clone( input = './Input/ForPosts/realTime.parameters.txt' )


