###
### This program will calculate calibration (or observed risk) for BLR-IPCW and MLR-IPCW in the 
### large sample analysis
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
# n.cohort <- as.numeric(args[2])
# n.pctls <- as.numeric(args[3])
print(paste("scen = ", scen, sep = ""))
# print(paste("n.cohort = ", n.cohort, sep = ""))
# print(paste("n.pctls = ", n.cohort, sep = ""))

### Load prepped data
load(paste("data/sim/large_sample_prep_data_", scen, ".RData", sep = ""))

###
### Calibration using BLR-IPCW
###

### Define output object
calib.moderate.blr <- vector("list", 3)

### Calculate the calibration plots
for (i in 1:3){
  print(paste("calib blr", i , Sys.time(), sep = " "))
  if (scen == "C1"){
    calib.moderate.blr[[i]] <- calib_msm(data.mstate = data.mstate,
                                     data.raw = data.raw,
                                     j = 1,
                                     s = 0,
                                     t = t.eval,
                                     tp.pred = tp.pred[[i]],
                                     calib.type = "blr",
                                     curve.type = "rcs",
                                     rcs.nk = 5,
                                     CI = 95,
                                     CI.R.boot = 200)
  } else {
    calib.moderate.blr[[i]] <- calib_msm(data.mstate = data.mstate,
                                     data.raw = data.raw,
                                     j = 1,
                                     s = 0,
                                     t = t.eval,
                                     tp.pred = tp.pred[[i]],
                                     calib.type = "blr",
                                     curve.type = "rcs",
                                     rcs.nk = 5,
                                     w.covs = c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45"),
                                     CI = 95,
                                     CI.R.boot = 200)
  }
  
}


### Save image
save.image(paste("data/sim/large_sample_analysis_moderate_blr_rcs_", scen, ".RData", sep = ""))
print("IMAGE SAVED")




