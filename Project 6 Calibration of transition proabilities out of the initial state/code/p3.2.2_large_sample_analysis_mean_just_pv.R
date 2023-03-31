###
### This program will calculate mean calibration (or observed risk) for just the pseudo-values in the 
### large sample analysis.
### Pseudo-values were not calculate for M2 or M3 scenarios
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
n.cohort <- as.numeric(args[2])
n.pctls <- as.numeric(args[3])
n.cohort <- 200000
n.pctls <- 20
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.cohort, sep = ""))

### Load workspace
load(paste("data/large_sample_analysis_mean_N", n.cohort, "_", scen, ".RData", sep = ""))
source("code/z_functions.R")

### Remove pv stuff if in workspace to ensure no crossover
rm(pv.comb, calib.pv.list)

###
### Calibration using pseudo values
###
print("START PV")
if (scen %in% c("M1C1", "M1C2", "M1C3")){
  ### Load combined pv data
  load(paste("data/large_sample_analysis_pv_combine_pctls_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
  source("code/z_functions.R")
  
  ### Calculate the columns means of pv.comb, and subtract from p.est.mean 
  ### (note haven't got pseudo values for all patients yet, so just do it for the ones that we do)
  calib.pv.list <- vector("list", 3)
  
  for (i in 1:3){
    calib.pv.list[[i]] <- colMeans(pv.comb[, paste("pv.state", 1:5, sep = "")]) - colMeans(p.est.mean[[i]])
  }
}
print("FINISH PV")


###
### Calibration using Aalen-Johansen estimator (within groups defined by risk)
###

### Load estimated of observed risk using Aalen-Johansen within percentiles
load(paste("data/large_sample_analysis_pv_prep_pctls_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))

### Get vector of Aalen-Johansen estimators for each state, having calculated Aalen-johansen within individuals grouped by predicted risk
obs.aj.pctls.mean <- c(lapply(obs.aj.pctls, function(x) {colMeans(do.call("rbind", x))})[[1]][1],
                       lapply(obs.aj.pctls, function(x) {colMeans(do.call("rbind", x))})[[2]][2],
                       lapply(obs.aj.pctls, function(x) {colMeans(do.call("rbind", x))})[[3]][3],
                       lapply(obs.aj.pctls, function(x) {colMeans(do.call("rbind", x))})[[4]][4],
                       lapply(obs.aj.pctls, function(x) {colMeans(do.call("rbind", x))})[[5]][5])

###
### Calibration using AJ (note I don't have to do this seperately for each estimated risk), we just calculate the observed risk
### and compare with the estimated risks
###
print(paste("START CALIB AJ ", Sys.time(), sep = ""))

### Calculate calibrations
calib.aj.list <- vector("list", 3)

for (i in 1:3){
  calib.aj.list[[i]] <- obs.aj.pctls.mean - colMeans(p.est.mean[[i]])
}

print(paste("FINISH CALIB AJ ", Sys.time(), sep = ""))


### Save image
save.image(paste("data/large_sample_analysis_mean_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")
