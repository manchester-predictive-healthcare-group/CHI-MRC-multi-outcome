###
### This program will calculate moderate alibration (or observed risk) for just the pseudo-values in the 
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
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.cohort, sep = ""))

### Load workspace
load(paste("data/large_sample_analysis_moderate_N", n.cohort, "_", scen, ".RData", sep = ""))
source("code/z_functions.R")

### Remove pv stuff if in workspace to ensure no crossover
rm(pv.comb, calib.pv.list)

###
### Calibration using pseudo values
###
print(paste("START PV", Sys.time()))
if (scen %in% c("M1C1", "M1C2", "M1C3")){
  ### Load combined pv data
  load(paste("data/large_sample_analysis_pv_combine_pctls_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
  source("code/z_functions.R")
  
  ### Calculate the columns means of pv.comb, and subtract from p.est.mean 
  ### (note haven't got pseudo values for all patients yet, so just do it for the ones that we do)
  calib.pv.list <- vector("list", 3)
 
  for (i in 1:3){
    calib.pv.list[[i]] <- calc.calib.pv.moderate(data.raw.reduc, pv.comb, p.est[[i]])
  }
}
print(paste("FINISH PV", Sys.time()))
rm(data.mstate.reduc, data.raw.reduc, p.true, pv.comb, pv.out)
save.image(paste("data/large_sample_analysis_moderate_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")

