### 
### This program will validate the calculated true transition probabilities
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages and functions and input parameters
source("code/z_functions.R")
source("code/z_load_packages.R")
source("code/z_functions_true_transition_probs.R")
source("code/p2/0_input_parameters.R")

### Read in seed/iter
args <- commandArgs(trailingOnly = T)
set <- as.numeric(args[1])
print(paste("set = ", set, sep = ""))

### Generate the x1 and x2 values we want to generate for, according to the set argument
covars <- rbind(c(0, 0, 0, 0, 0, 0, 0, 0),
                c(0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0),
                c(0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5),
                c(1, 1, 1, 1, 1, 1, 1, 1),
                c(-1, -1, -1, -1, -1, -1, -1, -1),
                c(0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5))
x12.eval <- covars[set, 1]
x13.eval <- covars[set, 2]
x15.eval <- covars[set, 3]
x24.eval <- covars[set, 4]
x25.eval <- covars[set, 5]
x34.eval <- covars[set, 6]
x35.eval <- covars[set, 7]
x45.eval <- covars[set, 8]

### Generate some baseline data
n.cohort <- 500000
x.baseline <- data.frame("x12" = rep(x12.eval, n.cohort),
                         "x13" = rep(x13.eval, n.cohort),
                         "x15" = rep(x15.eval, n.cohort),
                         "x24" = rep(x24.eval, n.cohort),
                         "x25" = rep(x25.eval, n.cohort),
                         "x34" = rep(x34.eval, n.cohort),
                         "x35" = rep(x35.eval, n.cohort),
                         "x45" = rep(x45.eval, n.cohort))

### Generate multistate data for these individuals
time.in <- Sys.time()
cohort <- gen.dat.DGM1(n = n.cohort, #number of patients to simulate
                       max.follow = max.follow, #max follow up, informative censoring 1 day after 7 years
                       shape12 = 1, scale12 = scales.sim["12"], #shape and scale for weibull baseline hazard for transition 1 -> 2
                       shape13 = 1, scale13 = scales.sim["13"], #shape and scale for weibull baseline hazard for transition 1 -> 3
                       shape15 = 1, scale15 = scales.sim["15"], #shape and scale for weibull baseline hazard for transition 1 -> 5
                       shape24 = 1, scale24 = scales.sim["24"], #shape and scale for weibull baseline hazard for transition 2 -> 4
                       shape25 = 1, scale25 = scales.sim["25"], #shape and scale for weibull baseline hazard for transition 2 -> 5
                       shape34 = 1, scale34 = scales.sim["34"], #shape and scale for weibull baseline hazard for transition 3 -> 4
                       shape35 = 1, scale35 = scales.sim["35"], #shape and scale for weibull baseline hazard for transition 3 -> 5
                       shape45 = 1, scale45 = scales.sim["45"], #shape and scale for weibull baseline hazard for transition 3 -> 5
                       beta.x12 = beta.x12, #covariate effects for transiion 12
                       beta.x13 = beta.x13, #covariate effects for transiion 13
                       beta.x15 = beta.x15, #covariate effects for transiion 15
                       beta.x24 = beta.x24, #covariate effects for transiion 24
                       beta.x25 = beta.x25, #covariate effects for transiion 25
                       beta.x34 = beta.x34, #covariate effects for transiion 34
                       beta.x35 = beta.x35, #covariate effects for transiion 35
                       beta.x45 = beta.x45, #covariate effects for transiion 45
                       x.in = x.baseline, #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
                       numsteps = max.follow)
time.out <- Sys.time()
time.diff <- time.out - time.in
print(paste("cohort generation took ", time.diff, sep = ""))
time.diff

### Convert into msate format with no censoring
cohort.mstate <- convert.mstate.DGM1.nocens(cohort[["cohort"]], cohort[["max.follow"]])

### Save image
save.image(paste("data/sim/validate_true_transition_probabilities_DGM1_set", set, ".RData", sep = ""))
print("IMAGE SAVED")

