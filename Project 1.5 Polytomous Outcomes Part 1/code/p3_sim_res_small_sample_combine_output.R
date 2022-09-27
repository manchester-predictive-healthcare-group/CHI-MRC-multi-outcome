### This code is for the small sample validation, so the models are developed on development datasets, and validated on validation datasets

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.5/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")

### Run through scenarios, seed.coef's, and K.sim, and combine all the results from the different paralleised runs of the simulation (defined by seed variable)

for (K.sim in c(3)){
  for (scenario.sim in c(1, 2, 3, 4)){
    for (seed.coef.sim in c(1, 2, 3)){
      for (N.devel.sim in c(2500, 1000, 500, 250)){
        
        print(paste("K.sim = ", K.sim, ", scenario.sim = ", scenario.sim, ", N.devel.sim = ", N.devel.sim, ", seed.coef.sim = ", seed.coef.sim))
        
        ### Load the data for seed = 1
        load(paste("data/sim_run_small_sample_s", scenario.sim, ".", seed.coef.sim, "_K", K.sim, "_N", N.devel.sim, "_seed", 1, ".RData", sep = ""))
        
        ### Assign sim.out.list to a combined .RData file
        sim.out.list.comb <- sim.out.list
        
        ### Remove simulation output to avoid crossover
        rm(sim.out.list)
        
        ### Add the simulation output from the other simulation runs
        for (seed.sim in 2:20){
          
          ### Load the data 
          load(paste("data/sim_run_small_sample_s", scenario.sim, ".", seed.coef.sim, "_K", K.sim, "_N", N.devel.sim, "_seed", seed.sim, ".RData", sep = ""))
          print(N.devel.sim)
          ### Assign sim.out.list to a combined .RData file
          sim.out.list.comb <- c(sim.out.list.comb, sim.out.list)
          
          ### Remove simulation output to avoid crossover
          rm(sim.out.list)
        }
        
        ### Save simulation output
        save.image(paste("data/sim_res_small_sample_s", scenario.sim, ".", seed.coef.sim, "_K", K.sim, "_N", N.devel.sim, ".RData", sep = ""))
      }
    }
  }
}