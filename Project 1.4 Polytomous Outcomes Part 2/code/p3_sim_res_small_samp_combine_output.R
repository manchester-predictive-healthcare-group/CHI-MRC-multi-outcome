### Code to combine and analyse data after simulation has been run

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.4/")

### Load packages
source("code/sim_load_packages.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")

# DGM.vec <- c("mult", "seqlog")
# scenario.vec <- c(1, 2, 3, 4)
# seed.coef.vec <- c(1, 2, 3)
# N.vec <- c(500, 1000, 2500)
# K.vec <- c(3, 5)

DGM.vec <- c("mult", "seqlog")
scenario.vec <- c(1, 2, 3, 4)
seed.coef.vec <- c(1, 2, 3)
N.vec <- c(100, 250, 500)
K.vec <- c(3, 5)

for (a in 1:length(DGM.vec)){
  for (b in 1:length(scenario.vec)){
    for (c in 1:length(seed.coef.vec)){
      for (d in 1:length(N.vec)){
        for (e in 1: length(K.vec)){
          
          DGM.in <- DGM.vec[a]
          scenario.sim <- scenario.vec[b]
          seed.coef.sim <- seed.coef.vec[c]
          N.devel.sim <- N.vec[d]
          K.sim <- K.vec[e]
          
          #           DGM.in <- "seqlog"
          #           scenario.sim <- 1
          #           seed.coef.sim <- 1
          #           N.devel.sim <- 1000
          
          print(Sys.time())
          ### Want to load data and concatenate into common objects
          ## Load for seed.sim = 1
          load(paste("data/s", scenario.sim, "/sim_run_small_samp_DGM", DGM.in, "_scen" , scenario.sim, ".", seed.coef.sim, "_K", 
                     K.sim, "_P", P.sim, "_N", N.devel.sim, "_s", 1, ".RData", sep =""))
          
          ### Put simulation output into combined objects
          ### Create output lists for multinomial and binary logistic modelling approaches for validation in both validation dataset, and development dataset
          ## Validation in validation dataset
          comb.sim.out.calib.data.multinomial <- sim.out.calib.data.multinomial
          #         comb.sim.out.calib.data.multinomial.rcs3 <- sim.out.calib.data.multinomial
          #         comb.sim.out.calib.data.multinomial.rcs4 <- sim.out.calib.data.multinomial
          #         comb.sim.out.calib.data.multinomial.rcs5 <- sim.out.calib.data.multinomial
          comb.sim.out.calib.data.binary <- sim.out.calib.data.binary
          #         comb.sim.out.calib.data.binary.rcs3 <- sim.out.calib.data.binary
          #         comb.sim.out.calib.data.binary.rcs4 <- sim.out.calib.data.binary
          #         comb.sim.out.calib.data.binary.rcs5 <- sim.out.calib.data.binary
          
          comb.sim.out.shrink.multinomial <- sim.out.shrink.multinomial
          #         comb.sim.out.shrink.multinomial.rcs3 <- sim.out.shrink.multinomial.rcs3
          #         comb.sim.out.shrink.multinomial.rcs4 <- sim.out.shrink.multinomial.rcs4
          #         comb.sim.out.shrink.multinomial.rcs5 <- sim.out.shrink.multinomial.rcs5
          comb.sim.out.shrink.binary <- sim.out.shrink.binary
          #         comb.sim.out.shrink.binary.rcs3 <- sim.out.shrink.binary.rcs3
          #         comb.sim.out.shrink.binary.rcs4 <- sim.out.shrink.binary.rcs4
          #         comb.sim.out.shrink.binary.rcs5 <- sim.out.shrink.binary.rcs5
          
          comb.sim.out.Cstat.multinomial <- sim.out.Cstat.multinomial
          #         comb.sim.out.Cstat.multinomial.rcs3 <- sim.out.Cstat.multinomial.rcs3
          #         comb.sim.out.Cstat.multinomial.rcs4 <- sim.out.Cstat.multinomial.rcs4
          #         comb.sim.out.Cstat.multinomial.rcs5 <- sim.out.Cstat.multinomial.rcs5
          comb.sim.out.Cstat.binary <- sim.out.Cstat.binary
          #         comb.sim.out.Cstat.binary.rcs3 <- sim.out.Cstat.binary.rcs3
          #         comb.sim.out.Cstat.binary.rcs4 <- sim.out.Cstat.binary.rcs4
          #         comb.sim.out.Cstat.binary.rcs5 <- sim.out.Cstat.binary.rcs5
          
          
          ### For the 25 different runs of the simulation, load the .RData file and append t the existing objects
          for (seed.sim in 2:20){
            print(seed.sim)
            ## Load the .RData file
            if (file.exists(paste("data/s", scenario.sim, "/sim_run_small_samp_DGM", DGM.in, "_scen" , scenario.sim, ".", seed.coef.sim, "_K", 
                                  K.sim, "_P", P.sim, "_N", N.devel.sim, "_s", seed.sim, ".RData", sep ="")) == FALSE){
              print("ERROR")
            }
            
            ### Load data
            load(paste("data/s", scenario.sim, "/sim_run_small_samp_DGM", DGM.in, "_scen" , scenario.sim, ".", seed.coef.sim, "_K", 
                       K.sim, "_P", P.sim, "_N", N.devel.sim, "_s", seed.sim, ".RData", sep =""))
            
            ### Combine with existing output
            ## Validation in validation dataset
            comb.sim.out.calib.data.multinomial <- c(comb.sim.out.calib.data.multinomial, sim.out.calib.data.multinomial)
            #           comb.sim.out.calib.data.multinomial.rcs3 <- c(comb.sim.out.calib.data.multinomial.rcs3, 
            #                                                    sim.out.calib.data.multinomial.rcs3)
            #           comb.sim.out.calib.data.multinomial.rcs4 <- c(comb.sim.out.calib.data.multinomial.rcs4, 
            #                                                    sim.out.calib.data.multinomial.rcs4)
            #           comb.sim.out.calib.data.multinomial.rcs5 <- c(comb.sim.out.calib.data.multinomial.rcs5, 
            #                                                    sim.out.calib.data.multinomial.rcs5)
            comb.sim.out.calib.data.binary <- c(comb.sim.out.calib.data.binary, sim.out.calib.data.binary)
            #           comb.sim.out.calib.data.binary.rcs3 <- c(comb.sim.out.calib.data.binary.rcs3, 
            #                                                    sim.out.calib.data.binary.rcs3)
            #           comb.sim.out.calib.data.binary.rcs4 <- c(comb.sim.out.calib.data.binary.rcs4, 
            #                                                    sim.out.calib.data.binary.rcs4)
            #           comb.sim.out.calib.data.binary.rcs5 <- c(comb.sim.out.calib.data.binary.rcs5, 
            #                                                    sim.out.calib.data.binary.rcs5)
            
            comb.sim.out.shrink.multinomial <- c(comb.sim.out.shrink.multinomial,
                                                 sim.out.shrink.multinomial)
            #           comb.sim.out.shrink.multinomial.rcs3 <- c(comb.sim.out.shrink.multinomial.rcs3,
            #                                                sim.out.shrink.multinomial.rcs3)
            #           comb.sim.out.shrink.multinomial.rcs4 <- c(comb.sim.out.shrink.multinomial.rcs4,
            #                                                sim.out.shrink.multinomial.rcs4)
            #           comb.sim.out.shrink.multinomial.rcs5 <- c(comb.sim.out.shrink.multinomial.rcs5,
            #                                                sim.out.shrink.multinomial.rcs5)
            comb.sim.out.shrink.binary <- c(comb.sim.out.shrink.binary,
                                            sim.out.shrink.binary)
            #           comb.sim.out.shrink.binary.rcs3 <- c(comb.sim.out.shrink.binary.rcs3,
            #                                                sim.out.shrink.binary.rcs3)
            #           comb.sim.out.shrink.binary.rcs4 <- c(comb.sim.out.shrink.binary.rcs4,
            #                                                sim.out.shrink.binary.rcs4)
            #           comb.sim.out.shrink.binary.rcs5 <- c(comb.sim.out.shrink.binary.rcs5,
            #                                                sim.out.shrink.binary.rcs5)
            
            comb.sim.out.Cstat.multinomial <- c(comb.sim.out.Cstat.multinomial,
                                                sim.out.Cstat.multinomial)
            #           comb.sim.out.Cstat.multinomial.rcs3 <- c(comb.sim.out.Cstat.multinomial.rcs3,
            #                                               sim.out.Cstat.multinomial.rcs3)
            #           comb.sim.out.Cstat.multinomial.rcs4 <- c(comb.sim.out.Cstat.multinomial.rcs4,
            #                                               sim.out.Cstat.multinomial.rcs4)
            #           comb.sim.out.Cstat.multinomial.rcs5 <- c(comb.sim.out.Cstat.multinomial.rcs5,
            #                                               sim.out.Cstat.multinomial.rcs5)
            comb.sim.out.Cstat.binary <- c(comb.sim.out.Cstat.binary,
                                           sim.out.Cstat.binary)
            #           comb.sim.out.Cstat.binary.rcs3 <- c(comb.sim.out.Cstat.binary.rcs3,
            #                                               sim.out.Cstat.binary.rcs3)
            #           comb.sim.out.Cstat.binary.rcs4 <- c(comb.sim.out.Cstat.binary.rcs4,
            #                                               sim.out.Cstat.binary.rcs4)
            #           comb.sim.out.Cstat.binary.rcs5 <- c(comb.sim.out.Cstat.binary.rcs5,
            #                                               sim.out.Cstat.binary.rcs5)
            
          }
          
          save.image(paste("data/results_small_samp_DGM", DGM.in, "_scen" , scenario.sim, ".", seed.coef.sim, "_K", 
                           K.sim, "_P", P.sim, "_N", N.devel.sim, ".RData", sep =""))
        }
      }
    }
  }
}

print("DONE")





