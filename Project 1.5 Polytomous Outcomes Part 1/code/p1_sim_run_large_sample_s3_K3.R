### This code is for the large sample validation, so the models are developed and calibrated on the same dataset

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.5/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions.R")
source("code/sim_functions_results.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")

### Define base scenario number
scenario.sim <- "3"

### Define number of outcomes
K.sim <- 3

## Seed to ensure different generated datasets each time this program is run in parallel
seed.sim <- 1

## Size of development datasets
N.devel.sim <- N.devel.large.sample

### Going to run the simulation multiple times, with differently randomly generated coefficients
for (seed.coef.sim in 1:5){
  
  ### Set seed for loading coefficient simulations
  set.seed(seed.coef.sim)
  
  ## And generate the simulation covariate effects
  coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                    c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)))
  
  ## Multinomial covariate effects
  coef.sim.multinomial <- coef.sim
  beta0.sim.multinomial <- c(-1, -1)*0.35
  
  ## Sequential logistic covariate effects
  coef.sim.seqlog <- coef.sim
  beta0.sim.seqlog <- c(0.5, 0)*2.25
  
#   ## One vs All covariate effects
#   coef.sim.OvA <- rbind(runif(P.sim, -1, 1), runif(P.sim,-1, 1), runif(P.sim,-1, 1))
#   beta0.sim.OvA <- c(-1, -1, -1)*0.35
#   
#   ## One vs One covariate effects
#   coef.sim.OvO.PC <- rbind(runif(P.sim, -1, 1), runif(P.sim,-1, 1), runif(P.sim,-1, 1))
#   beta0.sim.OvO.PC <- c(-1, -1, -1)*0.35
  
  
  ### Generate development datasets
  set.seed(seed.sim)
  
  dat.devel.list <- vector("list", 2)
  names(dat.devel.list) <- c("DGM.multinomial", "DGM.seqlog")
  dat.devel.list[[1]] <- generate.data.DGM.multinomial(K = K.sim, P = P.sim, coef = coef.sim.multinomial, beta0 = beta0.sim.multinomial, N = N.devel.sim)
  dat.devel.list[[2]] <- generate.data.DGM.seqlog(K = K.sim, P = P.sim, coef = coef.sim.seqlog, beta0 = beta0.sim.seqlog, N = N.devel.sim)
  #   dat.devel.list[[3]] <- generate.data.DGM.OvA(K = K.sim, P = P.sim, coef = coef.sim.OvA, beta0 = beta0.sim.OvA, N = N.devel.sim)
  #   dat.devel.list[[4]] <- generate.data.DGM.OvO.PC(K = K.sim, P = P.sim, coef = coef.sim.OvO.PC, beta0 = beta0.sim.OvO.PC, N = N.devel.sim)
  
  ### Check prevalence of each outcome
  prop.table(table(dat.devel.list[[1]]$Y))
  prop.table(table(dat.devel.list[[2]]$Y))
  #   prop.table(table(dat.devel.list[[3]]$Y))
  #   prop.table(table(dat.devel.list[[4]]$Y))
    
  
  ### Run large sample simulation for each dataset and store output
  ## Create vector to store output
  sim.out.list <- vector("list", 2)
  names(sim.out.list) <- c("DGM.multinomial", "DGM.seqlog")
  
  ## Run the simulation
  for (i in 1:2){
    sim.out.list[[i]] <- run.sim.large.sample(dat.devel.list[[i]])
    print(i)
    print(Sys.time())
  }
  rm(i)
  
  
  ### Create tables for each DGM
  sim.res.tables <- vector("list", 2)
  names(sim.res.tables) <- c("DGM.multinomial", "DGM.seqlog")
  for (i in 1:2){
    sim.res.tables[[i]] <- create.table.large.sample(sim.out.list[[i]], dp = 2)
  }
  
  ### Create plots for each DGM using mlr recalibration framework
  sim.res.plot.flex.mlr <- vector("list", 2)
  names(sim.res.plot.flex.mlr) <- c("DGM.multinomial", "DGM.seqlog")
  for (i in 1:2){
    sim.res.plot.flex.mlr[[i]] <- create.plot.flex.mlr(sim.out.list[[i]])
  }
  
  
  ### Create plots for each DGM using binary logistic recalibration framework
  sim.res.plot.flex.po <- vector("list", 2)
  names(sim.res.plot.flex.po) <- c("DGM.multinomial", "DGM.seqlog")
  for (i in 1:2){
    sim.res.plot.flex.po[[i]] <- create.plot.flex.po(sim.out.list[[i]])
  }
  
  ### Save image
  save.image(paste("data/sim_run_large_sample_s", scenario.sim, ".", seed.coef.sim, "_K", K.sim, ".RData", sep = ""))
  
  ### Remove stuff before next simulation run, to ensure no crossover
  rm(dat.devel.list, sim.out.list, sim.res.tables, sim.res.plot.flex.mlr, sim.res.plot.flex.po)
  
  ### Print marker
  print(seed.coef.sim)
  
}

