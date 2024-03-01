###
### This program combines the results from the small sample simulation sensitivity analyses, ready to be plotted/tabled
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

### Print these
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))

### Create an object to store the pseudo-values
### Each list element will contain the pseudo-values for a given state
calib.aj.mean.list.comb <- vector("list", 3)
calib.blr.mean.list.comb <- vector("list", 3)
calib.mlr.mean.list.comb <- vector("list", 3)
calib.true.mean.list.comb <- vector("list", 3)

### Combine them
for (set in 1:50){
  
  ### Load data
  load(paste("data/sim/small_sample_mean_analysis_sens_", scen, "_n", n.cohort, "_set", set, ".RData", sep = ""))
  
  ### Concatenate into dataset
  calib.aj.mean.list.comb <- lapply(c(1,2,3), function(x) {rbind(calib.aj.mean.list.comb[[x]], calib.aj.mean.list[[x]])})
  calib.blr.mean.list.comb <- lapply(c(1,2,3), function(x) {rbind(calib.blr.mean.list.comb[[x]], calib.blr.mean.list[[x]])})
  calib.mlr.mean.list.comb <- lapply(c(1,2,3), function(x) {rbind(calib.mlr.mean.list.comb[[x]], calib.mlr.mean.list[[x]])})
  calib.true.mean.list.comb <- lapply(c(1,2,3), function(x) {rbind(calib.true.mean.list.comb[[x]], calib.true.mean.list[[x]])})
  
}

### Save image
rm(list=setdiff(ls(), list("calib.aj.mean.list.comb", "calib.blr.mean.list.comb",  "calib.mlr.mean.list.comb", "calib.true.mean.list.comb", 
                           "n.cohort", "scen", "n.pctls")))
save.image(paste("data/sim/small_sample_mean_analysis_sens_", scen, "_n", n.cohort, ".RData", sep = ""))
print("IMAGE SAVED")