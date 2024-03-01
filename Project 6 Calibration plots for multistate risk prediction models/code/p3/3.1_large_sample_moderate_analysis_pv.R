###
### Estimate moderate calibration for pseudo-value approach.
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

### Load prepped dat afor pv (this contains the appropriate datasets, and also estimate of Aalen-Johansen)
load(paste("data/sim/large_sample_prep_data_", scen, ".RData", sep = ""))

###
### Calibration using pseudo values
###

### Load combined pv data
load(paste("data/sim/large_sample_combine_pv_", scen, "_npctls", 20, ".RData", sep = ""))
source("code/z_functions.R")

### Define output object
calib.moderate.pv <- vector("list", 3)

### Calculate the calibration plots
for (i in 1:3){
  print(paste("calib pv", i , Sys.time(), sep = " "))
  calib.moderate.pv[[i]] <- calib_msm(data.mstate = data.mstate,
                                      data.raw = data.raw,
                                      j = 1,
                                      s = 0,
                                      t = t.eval,
                                      tp.pred = tp.pred[[i]],
                                      calib.type = "pv",
                                      curve.type = "rcs",
                                      rcs.nk = 5,
                                      pv.precalc = pv.comb[,-1],
                                      CI = 95,
                                      CI.type = "parametric")
}

### Remove everything except calib.moderate.pv
rm(list=setdiff(ls(), list("calib.moderate.pv", "scen")))

### Save image
save.image(paste("data/sim/large_sample_analysis_moderate_pv_", scen, ".RData", sep = ""))
print("IMAGE SAVED")
