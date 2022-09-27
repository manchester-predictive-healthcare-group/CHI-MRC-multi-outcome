### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.4/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")

### Define base scenario number
scenario.sim <- "2"

### Number of outcome categories
K.sim <- 3

### Size of development datasets
N.devel.sim <- 500000

### Set seed for loading coef.sim, and generate coef.sim
for (seed.coef.sim in 1:5){
  
  ### Set seed
  set.seed(seed.coef.sim)
  
  ### Generate coefficients
  coef.sim <- rbind(c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                    c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))
  beta0.sim <- c(-0.5, 0)*2.25
  
  ### Generate development daaset
  data.devel <- generate.data.DGM.seqlog(K = K.sim, P = P.sim, coef = coef.sim, beta0 = beta0.sim, N = N.devel.sim)
  
  ### Fit model and assign calibration output
  sim.out.multinomial <- calibrate.model.multinomial(dat.devel = data.devel, dat.valid = data.devel)
  sim.out.binary <- calibrate.model.binary(dat.devel = data.devel, dat.valid = data.devel)
  sim.out.multinomial.rcs3 <- calibrate.model.multinomial.rcs(dat.devel = data.devel, dat.valid = data.devel, 
                                                    n.knot.in = 3)
  sim.out.multinomial.rcs4 <- calibrate.model.multinomial.rcs(dat.devel = data.devel, dat.valid = data.devel, 
                                                    n.knot.in = 4)
  sim.out.multinomial.rcs5 <- calibrate.model.multinomial.rcs(dat.devel = data.devel, dat.valid = data.devel, 
                                                    n.knot.in = 5)
  sim.out.multinomial.rcs6 <- calibrate.model.multinomial.rcs(dat.devel = data.devel, dat.valid = data.devel, 
                                                    n.knot.in = 6)
  
  print(seed.coef.sim)
  print(Sys.time())
  
  ## Save .RData file
  save.image(paste("data/s", scenario.sim, "/sim_run_large_sample_DGMseqlog_scen" , scenario.sim, ".", seed.coef.sim, "_K", 
                   K.sim, "_P", P.sim, "_N", N.devel.sim, ".RData", sep =""))
  
  ## Remove stuff to ensure no cross over
  rm(sim.out.multinomial, sim.out.binary, sim.out.multinomial.rcs3, sim.out.multinomial.rcs4, sim.out.multinomial.rcs5, sim.out.multinomial.rcs6)
}


  
  

