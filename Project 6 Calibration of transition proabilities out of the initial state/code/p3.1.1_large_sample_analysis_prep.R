###
### This program will prep the data ready to calculate calibration for the cohort sie and scenario of interest
### Aalen-Johansen estimator of the entire cohort will be calculated
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
n.cohort <- as.numeric(args[2])
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))

### Load simulated data
load(paste("data/apply_cens_", scen, ".RData", sep = ""))
source("code/z_functions.R")

### Set seed
set.seed(101)

### Extract 100,000 individuals, which is what we're doing the large sample analysis with
data.mstate.reduc <- data.mstate.obj[["data.mstate"]][data.mstate.obj[["data.mstate"]]$patid %in% 1:n.cohort, ]
data.raw.reduc <- data.mstate.obj[["data.raw"]][data.mstate.obj[["data.raw"]]$patid %in% 1:n.cohort, ]
tmat <- data.mstate.obj[["tmat"]]


### 
### Calculate Aalen-Johansen estimator in entire cohort
print(paste("AJ ", Sys.time()))
obs.aj.object <- calc.calib.aj(data.mstate.reduc, 
                        tmat, 
                        t.eval)
print(paste("AJ ", Sys.time()))

### Extract standard errors and estimates seperately
obs.aj <- obs.aj.object[["obs.aj"]]
obs.aj.se <- obs.aj.object[["obs.aj.se"]]

###
### Save image
save.image(paste("data/large_sample_analysis_prep_N", n.cohort, "_", scen, ".RData", sep = ""))
print("IMAGE SAVED")