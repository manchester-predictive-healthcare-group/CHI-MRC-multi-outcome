###
### This program combines the results from the small sample simulation for moderate calibration, ready to be plotted.
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")
source("code/z_functions.R")

### Extract sample size, number of percentiles/groups, and task id from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
n.cohort  <- as.numeric(args[2])
n.pctls <- as.numeric(args[3])

### Print these
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))

### Create an object to store the pseudo-values
### Each list element will contain the pseudo-values for a given state
calib.pv.mod.list.comb <- vector("list", 3)
calib.blr.mod.list.comb <- vector("list", 3)
calib.mlr.mod.list.comb <- vector("list", 3)

### Combine them
for (set in 1:50){
  
  ### Load data
  load(paste("data/sim/small_sample_moderate_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, "_set", set, ".RData", sep = ""))
  
  ### Concatenate
  calib.pv.mod.list.comb <- lapply(c(1,2,3), function(x) {c(calib.pv.mod.list.comb[[x]], calib.pv.mod.list[[x]])})
  calib.blr.mod.list.comb <- lapply(c(1,2,3), function(x) {c(calib.blr.mod.list.comb[[x]], calib.blr.mod.list[[x]])})
  calib.mlr.mod.list.comb <- lapply(c(1,2,3), function(x) {c(calib.mlr.mod.list.comb[[x]], calib.mlr.mod.list[[x]])})
  
}

### Save image
rm(list=setdiff(ls(), list("calib.pv.mod.list.comb", "calib.blr.mod.list.comb",  "calib.mlr.mod.list.comb", 
                           "n.cohort", "scen", "n.pctls")))
save.image(paste("data/sim/small_sample_moderate_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")