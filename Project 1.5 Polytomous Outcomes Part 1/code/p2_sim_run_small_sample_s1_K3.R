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

### Load input parameters from qsub file (or assign them manually in a test)
args <- commandArgs(trailingOnly = T)
#args <- c(1, 1, 1000)
## seed for simulation
as.numeric(args[1])
## seed for covariate effects
as.numeric(args[2])
## sample size for development dataset
as.numeric(args[3])

## Seed to ensure different generated datasets each time this program is run in parallel
seed.sim <- as.numeric(args[1])

## Seed for defining the covariate effects
seed.coef.sim <- as.numeric(args[2])

## Size of development datasets
N.devel.sim <- as.numeric(args[3])

### Define base scenario number
scenario.sim <- "1"

## Set number of outcomes
K.sim <- 3

### Set seed for loading coefficient simulations
set.seed(seed.coef.sim)

## And generate the simulation covariate effects
coef.sim <- rbind(c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)))

## Multinomial covariate effects
coef.sim.multinomial <- coef.sim
beta0.sim.multinomial <- c(-1, -1)*0.35

## Sequential logistic covariate effects
coef.sim.seqlog <- coef.sim
beta0.sim.seqlog <- c(0.5, 0)*2.25

# ## One vs All covariate effects
# coef.sim.OvA <- rbind(runif(P.sim, -1, 1), runif(P.sim,-1, 1), runif(P.sim,-1, 1))
# beta0.sim.OvA <- c(-1, -1, -1)*0.35
# 
# ## One vs One covariate effects
# coef.sim.OvO.PC <- rbind(runif(P.sim, -1, 1), runif(P.sim,-1, 1), runif(P.sim,-1, 1))
# beta0.sim.OvO.PC <- c(-1, -1, -1)*0.35

###
###
### Run the simulation N.sim times and store output
###
###

### Create list to store output
sim.out.list <- vector("list", N.sim)

### Set seed for this batch of simulation runs
set.seed(seed.sim)

### Run the sim
for (sim in 1:N.sim){
  
  ### Create data using each DGM
  ## Development data
  dat.devel.list <- vector("list", 2)
  names(dat.devel.list) <- c("DGM.multinomial", "DGM.seqlog")
  dat.devel.list[[1]] <- generate.data.DGM.multinomial(K = K.sim, P = P.sim, coef = coef.sim.multinomial, beta0 = beta0.sim.multinomial, N = N.devel.sim)
  print(Sys.time())
  dat.devel.list[[2]] <- generate.data.DGM.seqlog(K = K.sim, P = P.sim, coef = coef.sim.seqlog, beta0 = beta0.sim.seqlog, N = N.devel.sim)
  print(Sys.time())
#   dat.devel.list[[3]] <- generate.data.DGM.OvA(K = K.sim, P = P.sim, coef = coef.sim.OvA, beta0 = beta0.sim.OvA, N = N.devel.sim)
#   print(Sys.time())
#   dat.devel.list[[4]] <- generate.data.DGM.OvO.PC(K = K.sim, P = P.sim, coef = coef.sim.OvO.PC, beta0 = beta0.sim.OvO.PC, N = N.devel.sim)
#   print(Sys.time())
  
  ## Validation data
  dat.valid.list <- vector("list", 2)
  names(dat.valid.list) <- c("DGM.multinomial", "DGM.seqlog")
  dat.valid.list[[1]] <- generate.data.DGM.multinomial(K = K.sim, P = P.sim, coef = coef.sim.multinomial, beta0 = beta0.sim.multinomial, N = N.valid.sim)
  print(Sys.time())
  dat.valid.list[[2]] <- generate.data.DGM.seqlog(K = K.sim, P = P.sim, coef = coef.sim.seqlog, beta0 = beta0.sim.seqlog, N = N.valid.sim)
  print(Sys.time())
#   dat.valid.list[[3]] <- generate.data.DGM.OvA(K = K.sim, P = P.sim, coef = coef.sim.OvA, beta0 = beta0.sim.OvA, N = N.valid.sim)
#   print(Sys.time())
#   dat.valid.list[[4]] <- generate.data.DGM.OvO.PC(K = K.sim, P = P.sim, coef = coef.sim.OvO.PC, beta0 = beta0.sim.OvO.PC, N = N.valid.sim)
#   print(Sys.time())
  
  ### Run simulation for each dataset and store output
  ## Create vector to store output
  sim.out.list[[sim]] <- vector("list", 2)
  names(sim.out.list[[sim]]) <- c("DGM.multinomial", "DGM.seqlog")
  
  ## Run the simulation
  for (i in 1:2){
    sim.out.list[[sim]][[i]] <- run.sim.small.sample(dat.devel.list[[i]], dat.valid.list[[i]])
    print(i)
    print(Sys.time())
  }
  rm(i)
}


### Save image
save.image(paste("data/sim_run_small_sample_s", scenario.sim, ".", seed.coef.sim, "_K", K.sim, "_N", N.devel.sim, "_seed", seed.sim, ".RData", sep = ""))


