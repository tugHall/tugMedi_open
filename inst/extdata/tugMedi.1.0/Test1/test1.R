#!/usr/bin/Rscript --slave --vanilla


library(tugMedi)


#####################################################
print("Generators start")


#write_prm(       input = './Input/ForGenerators/parameters.parameters.txt')

#write_CF(        input = './Input/ForGenerators/CF.parameters.txt')

write_weights(   input = './Input/ForGenerators/hallmark_weights.parameters.txt' )

#write_cloneinit( input  = './Input/ForGenerators/cloneinit.parameters.txt' )

#write_Rint(      input = './Input/ForGenerators/Rint.parameters.txt' )

write_EF.Rint(   input = './Input/ForGenerators/EF.Rint.parameters.txt' )

# Only this depends on write_prm() result
write_EF.Rother( input = './Input/ForGenerators/EF.Rother.parameters.txt' )


##############################################################
print("Simulation starts")


system.time(
start_simulation( input_files = 'Input/DATA/FILES.txt', 
                  saveTo      = 'Output/Results.sim.01.RDS', 
                  seed        = NA
                 )
)


write_VAF(   input = './Input/ForPosts/VAF.parameters.txt' )

write_TMB(   input = './Input/ForPosts/TMB.parameters.txt' )

write_evals( input = './Input/ForPosts/evals.parameters.txt' )

##############################################################
print("Drug Intervention start")


update_RDS_from_output( in.file  = 'Input/DATA/FILES.txt', 
                        rds.file = 'Output/Results.sim.01.RDS' )
                        
drug_int_param <- list( 'kill_prob' = 0.5,
                        'block_prob'= 0.0,
                        'gene' = c('KRAS', 'APC') ) # Any of them is a target


system.time(
    start_simulation( first       = F,
                      input_files = 'Input/DATA/FILES.drgInt.txt', 
                      loadFrom    = 'Output/Results.sim.01.RDS', 
                      saveTo      = 'Output/Results.sim.02.RDS', 
                      seed        = NA,
                      drug_int_param = drug_int_param
                      #                  drug_int_param = NULL
    )
)


write_realTime_clone( input = './Input/ForPosts/realTime.parameters.txt' )




