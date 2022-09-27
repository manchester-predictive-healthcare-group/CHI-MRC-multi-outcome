### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.4/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")

### Load input parameters from qsub file
args <- commandArgs(trailingOnly = T)
#args <- c(1,1,1000)
as.numeric(args[1])
as.numeric(args[2])
as.numeric(args[3])

## Seed to ensure different generated datasets each time this program is run in parallel
seed.sim <- as.numeric(args[1])

## Seed for defining the coefficients
seed.coef.sim <- as.numeric(args[2])

## Size of development datasets
N.devel.sim <- as.numeric(args[3])

### Define number of outcomes
K.sim <- 3

### Define base scenario number
scenario.sim <- "2"

### Set seed for loading coef.sim, and generate coef.sim
set.seed(seed.coef.sim)
coef.sim <- rbind(c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
beta0.sim <- c(1, 1)*0.1


### Calculate upper and lower bounds for pred.eval
bounds <- calc.pred.quantiles.DGM.mult(K.sim, P.sim, coef.sim, beta0.sim)
pred.eval <- seq(bounds[1], bounds[2], 0.001)

### Create output lists for multinomial and binary logistic modelling approaches (validated in validation dataset, primary analysis)
sim.out.calib.data.multinomial <- vector("list", N.sim)
# sim.out.calib.data.multinomial.rcs3 <- vector("list", N.sim)
# sim.out.calib.data.multinomial.rcs4 <- vector("list", N.sim)
# sim.out.calib.data.multinomial.rcs5 <- vector("list", N.sim)
sim.out.calib.data.binary <- vector("list", N.sim)
# sim.out.calib.data.binary.rcs3 <- vector("list", N.sim)
# sim.out.calib.data.binary.rcs4 <- vector("list", N.sim)
# sim.out.calib.data.binary.rcs5 <- vector("list", N.sim)

sim.out.shrink.multinomial <- rep(NA, N.sim)
# sim.out.shrink.multinomial.rcs3 <- rep(NA, N.sim)
# sim.out.shrink.multinomial.rcs4 <- rep(NA, N.sim)
# sim.out.shrink.multinomial.rcs5 <- rep(NA, N.sim)
sim.out.shrink.binary <- rep(NA, N.sim)
# sim.out.shrink.binary.rcs3 <- rep(NA, N.sim)
# sim.out.shrink.binary.rcs4 <- rep(NA, N.sim)
# sim.out.shrink.binary.rcs5 <- rep(NA, N.sim)

sim.out.Cstat.multinomial <- rep(NA, N.sim)
# sim.out.Cstat.multinomial.rcs3 <- rep(NA, N.sim)
# sim.out.Cstat.multinomial.rcs4 <- rep(NA, N.sim)
# sim.out.Cstat.multinomial.rcs5 <- rep(NA, N.sim)
sim.out.Cstat.binary <- rep(NA, N.sim)
# sim.out.Cstat.binary.rcs3 <- rep(NA, N.sim)
# sim.out.Cstat.binary.rcs4 <- rep(NA, N.sim)
# sim.out.Cstat.binary.rcs5 <- rep(NA, N.sim)



### Set seet for this run of simulations (this will change for each time this code is run in parallel)
set.seed(seed.sim)

### Run the simulation N.sim times
for (i in 1:N.sim){
  
  ### Generate development daaset
  data.devel <- generate.data.DGM.mult(K = K.sim, P = P.sim, coef = coef.sim, beta0 = beta0.sim, N = N.devel.sim)
  ### Create validation dataset
  data.valid <- generate.data.DGM.mult(K = K.sim, P = P.sim, coef = coef.sim, beta0 = beta0.sim, N = N.valid.sim)
  
  ### Fit model and assign calibration output
  print(Sys.time())
  sim.out.multinomial <- calibrate.model.multinomial.predfix(dat.devel = data.devel, dat.valid = data.valid, pred.eval = pred.eval)
#   print(Sys.time())
#   sim.out.multinomial.rcs3 <- calibrate.model.multinomial.predfix.rcs(dat.devel = data.devel, dat.valid = data.valid, pred.eval = pred.eval, 
#                                                             n.knot.in = 3)
#   print(Sys.time())
#   sim.out.multinomial.rcs4 <- calibrate.model.multinomial.predfix.rcs(dat.devel = data.devel, dat.valid = data.valid, pred.eval = pred.eval, 
#                                                             n.knot.in = 4)
#   print(Sys.time())
#   sim.out.multinomial.rcs5 <- calibrate.model.multinomial.predfix.rcs(dat.devel = data.devel, dat.valid = data.valid, pred.eval = pred.eval, 
#                                                             n.knot.in = 5)
  print(Sys.time())
  sim.out.binary <- calibrate.model.binary.predfix(dat.devel = data.devel, dat.valid = data.valid, pred.eval = pred.eval)
#   print(Sys.time())
#   sim.out.binary.rcs3 <- calibrate.model.binary.predfix.rcs(dat.devel = data.devel, dat.valid = data.valid, pred.eval = pred.eval, 
#                                                        n.knot.in = 3)
#   print(Sys.time())
#   sim.out.binary.rcs4 <- calibrate.model.binary.predfix.rcs(dat.devel = data.devel, dat.valid = data.valid, pred.eval = pred.eval, 
#                                                             n.knot.in = 4)
#   print(Sys.time())
#   sim.out.binary.rcs5 <- calibrate.model.binary.predfix.rcs(dat.devel = data.devel, dat.valid = data.valid, pred.eval = pred.eval, 
#                                                             n.knot.in = 5)
  print(Sys.time())

  ### Put respective simulation output into the correct objects (validated in validation dataset)
  ## Calibration plot data
  sim.out.calib.data.multinomial[[i]] <- sim.out.multinomial[[1]]
#   sim.out.calib.data.multinomial.rcs3[[i]] <- sim.out.multinomial.rcs3[[1]]
#   sim.out.calib.data.multinomial.rcs4[[i]] <- sim.out.multinomial.rcs4[[1]]
#   sim.out.calib.data.multinomial.rcs5[[i]] <- sim.out.multinomial.rcs5[[1]]
  sim.out.calib.data.binary[[i]] <- sim.out.binary[[1]]
#   sim.out.calib.data.binary.rcs3[[i]] <- sim.out.binary.rcs3[[1]]
#   sim.out.calib.data.binary.rcs4[[i]] <- sim.out.binary.rcs4[[1]]
#   sim.out.calib.data.binary.rcs5[[i]] <- sim.out.binary.rcs5[[1]]
  
  ## Shrinkage factor/calibration slope
  ## Cstatistics
  sim.out.shrink.multinomial[i] <- sim.out.multinomial[[2]]
#   sim.out.shrink.multinomial.rcs3[i] <- sim.out.multinomial.rcs3[[2]]
#   sim.out.shrink.multinomial.rcs4[i] <- sim.out.multinomial.rcs4[[2]]
#   sim.out.shrink.multinomial.rcs5[i] <- sim.out.multinomial.rcs5[[2]]
  sim.out.shrink.binary[i] <- sim.out.binary[[2]]
#   sim.out.shrink.binary.rcs3[i] <- sim.out.binary.rcs3[[2]]
#   sim.out.shrink.binary.rcs4[i] <- sim.out.binary.rcs4[[2]]
#   sim.out.shrink.binary.rcs5[i] <- sim.out.binary.rcs5[[2]]
  
  ## Cstatistics
  sim.out.Cstat.multinomial[i] <- sim.out.multinomial[[3]]
#   sim.out.Cstat.multinomial.rcs3[i] <- sim.out.multinomial.rcs3[[3]]
#   sim.out.Cstat.multinomial.rcs4[i] <- sim.out.multinomial.rcs4[[3]]
#   sim.out.Cstat.multinomial.rcs5[i] <- sim.out.multinomial.rcs5[[3]]
  sim.out.Cstat.binary[i] <- sim.out.binary[[3]]
#   sim.out.Cstat.binary.rcs3[i] <- sim.out.binary.rcs3[[3]]
#   sim.out.Cstat.binary.rcs4[i] <- sim.out.binary.rcs4[[3]]
#   sim.out.Cstat.binary.rcs5[i] <- sim.out.binary.rcs5[[3]]
  
  ## Add a variable indicating the simulation run to the calibration plot data
  sim.out.calib.data.multinomial[[i]]$sim.run <- rep(i, nrow(sim.out.calib.data.multinomial[[i]]))
#   sim.out.calib.data.multinomial.rcs3[[i]]$sim.run <- rep(i, nrow(sim.out.calib.data.multinomial.rcs3[[i]]))
#   sim.out.calib.data.multinomial.rcs4[[i]]$sim.run <- rep(i, nrow(sim.out.calib.data.multinomial.rcs4[[i]]))
#   sim.out.calib.data.multinomial.rcs5[[i]]$sim.run <- rep(i, nrow(sim.out.calib.data.multinomial.rcs5[[i]]))
  sim.out.calib.data.binary[[i]]$sim.run <- rep(i, nrow(sim.out.calib.data.binary[[i]]))
#   sim.out.calib.data.binary.rcs3[[i]]$sim.run <- rep(i, nrow(sim.out.calib.data.binary.rcs3[[i]]))
#   sim.out.calib.data.binary.rcs4[[i]]$sim.run <- rep(i, nrow(sim.out.calib.data.binary.rcs4[[i]]))
#   sim.out.calib.data.binary.rcs5[[i]]$sim.run <- rep(i, nrow(sim.out.calib.data.binary.rcs5[[i]]))
  
  print(paste("sim",i))
  print(Sys.time())
}

## Save .RData file, to be combined with others run in parallel for analysis
save.image(paste("data/s", scenario.sim, "/sim_run_small_samp_DGMmult_scen" , scenario.sim, ".", seed.coef.sim, "_K", 
                 K.sim, "_P", P.sim, "_N", N.devel.sim, "_s", seed.sim, ".RData", sep =""))
  
  

