###
### This program will prepare the data to enable estimation of the pseudo-values in a computationally efficient manner
### 1) Pre-calculate Aalen-Johansen estimator for entire cohort (subgrouped) so this can be used within calculation of pseudo-values
### 2) Pre-calculate the pseudo-value for an individual who doesnt have an event until after t.eval, which will be the same for all individuals
### for which this is the case
###
### Note that because this is for the senstivity analysis, we do not group individuals by predicted risk before calculating pseudo-values.

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]

print(paste("scen = ", scen, sep = ""))
### n.pctls is the number of groups we will split data into before calculating pseudo-values

### Load prepped data
load(paste("data/sim/large_sample_prep_data_", scen, ".RData", sep = ""))

### Load packages
source("code/z_functions.R")
source("code/z_load_packages.R")

###
### Create n.pctls datasets, grouped by predicted risk
###
obs.aj <- calc.calib.aj(data.mstate = data.mstate, 
                        tmat = tmat, 
                        t.eval = t.eval)$obs.aj


###
### Calculate the pseudo value for someone who is uncensored and has not had an event prior to time t
### This pseudo value will be the same for all patients who meet these criteria, so no need to recalculate
### it everytime. This will help with computational time.
###

### Start by identifying an individual for which this is true
patid.pv.same <- data.mstate %>%
  subset(from == 1 & to == 2 & Tstop > t.eval) %>%
  slice(1) %>%
  select(patid) %>%
  as.numeric()

### Calculate the pseudo value for this individual
pv.same <- func.calc.pv.aj(patid.eval = patid.pv.same, 
                           data.mstate = data.mstate, 
                           obs.aj = obs.aj, 
                           tmat = tmat, 
                           n.cohort = nrow(data.raw), 
                           t.eval = t.eval)


### Save image
rm(list=setdiff(ls(), list("data.mstate", "data.raw",
                           "pv.same", "obs.aj", 
                           "p.true", "tmat", "t.eval", "scen", "n.cohort")))
save.image(paste("data/sim/large_sample_prep_data_pv_", scen, "_sens.RData", sep = ""))
print("IMAGE SAVED")


