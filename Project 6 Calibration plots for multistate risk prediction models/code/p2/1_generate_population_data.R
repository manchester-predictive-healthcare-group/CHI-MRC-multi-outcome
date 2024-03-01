###
### This program will generate a cohort of size 1,000,000 data according to DGM1
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
library(gems)
library(cubature)
library(dplyr)

### Source functions and input parameters for simulation
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")
source("code/p2/0_input_parameters.R")

### Extract arguments from command line
args <- commandArgs(trailingOnly = T)
seed <- as.numeric(args[1])
print(paste("seed =", seed))

### Choose cohort size
n.cohort <- 1000

### Generate the data, but paralellise the process to improve speed (note I could move to CSF at a later date if this process is highly time consuming,
### but at the moment I think other things require the CSF more (fitting the actual multistate models))
print("START DATA GEN")
## Set seed
set.seed(seed)

## Generate baseline data
x.baseline <- data.frame("x12" = rnorm(n.cohort, 0, 1), 
                         "x13" = rnorm(n.cohort, 0, 1),
                         "x15" = rnorm(n.cohort, 0, 1),
                         "x24" = rnorm(n.cohort, 0, 1),
                         "x25" = rnorm(n.cohort, 0, 1),
                         "x34" = rnorm(n.cohort, 0, 1),
                         "x35" = rnorm(n.cohort, 0, 1),
                         "x45" = rnorm(n.cohort, 0, 1))

Sys.time()
## Generate data
data.out <- gen.dat.DGM1(n = n.cohort, #number of patients to simulate
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
                       
str(data.out)

###
### Calculate true transition probabilities for each individual
###

### Write a function to it for a single individual
print(paste("CALC TRUE RISKS ", Sys.time(), sep = ""))
calc.true.tp.ind <- function(row){
  return(calc.true.transition.probs.DGM1(0, t.eval, row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8],
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
                                         beta.x45 = beta.x45 #covariate effects for transiion 45
  ))
}
Sys.time()

                                 
### Calc true risk for each individual in this dataset
p.true <- apply(dplyr::select(data.out[["cohort"]], x12, x13, x15, x24, x25, x34, x35, x45), 1, calc.true.tp.ind)
p.true <- data.frame(t(p.true))
colnames(p.true) <- paste("p.true", 1:5, sep = "")
print(paste("FINISH TRUE RISKS ", Sys.time(), sep = ""))

### Combine data into one data frame
data.out <- cbind(data.out[[1]], p.true)

### Save data object
saveRDS(data.out, paste("data/sim/popdata/dgm1_seed", seed, ".rds", sep = ""))
