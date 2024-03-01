###
### This program will prep the data in order to calculate pseudo-values
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")

### Load work space with data in mstate format
load("data/ce/ce_fit_csh_msm_N100000.RData")
source("code/z_functions_ce.R")

### Create a "raw" (not mstate format) dataset
data.raw <- complete.data

### Create data.mstaet formatted data
data.mstate <- complete.data.prep.valid

### Set t.eval
t.eval <- ceiling(10*365.25)

### 
### Calculate Aalen-Johansen estimator in entire cohort
###
print(paste("AJ ", Sys.time()))
obs.aj <- calc.calib.aj.ce(data.mstate, 
                           tmat, 
                           t.eval)
print(paste("AJ ", Sys.time()))

### Save image
rm(complete.data, complete.data.prep.devel, complete.data.prep.valid)
save.image("data/ce/ce_pv_prep.RData")

###
### Calculate the pseudo value for someone who is uncensored and has not had an event prior to time t
### This pseudo value will be the same for all patients who meet these criteria, so no need to recalculate
### it everytime.
###

### Start by finding an id which this is true for
patid.pv.same <- data.mstate %>%
  subset(from == 1 & to == 2 & Tstop > t.eval) %>%
  slice(1) %>%
  select(person_id) %>%
  as.numeric()

### Calculate the pseudo value for this individual
print(paste("START PV SAME", Sys.time()))
pv.same <- func.calc.pv.aj.ce(person_id.eval = patid.pv.same, 
                              data.mstate = data.mstate, 
                              obs.aj = obs.aj, 
                              tmat = tmat, 
                              n.cohort = length(unique(data.mstate$person_id)), 
                              t.eval = t.eval)
print(paste("FINISH PV SAME", Sys.time()))
print(pv.same)

### Save image
save.image("data/ce/ce_pv_prep.RData")


