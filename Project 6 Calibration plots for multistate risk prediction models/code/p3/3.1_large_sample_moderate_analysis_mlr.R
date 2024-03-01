###
### Estimate moderate calibration for MLR-IPCW.
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
### Calibration using MLR-IPCW
###

### Define output object
calib.moderate.mlr <- vector("list", 3)

### Calculate the calibration plots
for (i in 1:3){
  print(paste("calib mlr", i , Sys.time(), sep = " "))
  if (scen == "C1"){
    calib.moderate.mlr[[i]] <- calib_msm(data.mstate = data.mstate,
                                         data.raw = data.raw,
                                         j = 1,
                                         s = 0,
                                         t = t.eval,
                                         tp.pred = tp.pred[[i]],
                                         calib.type = "mlr")
  } else {
    calib.moderate.mlr[[i]] <- calib_msm(data.mstate = data.mstate,
                                         data.raw = data.raw,
                                         j = 1,
                                         s = 0,
                                         t = t.eval,
                                         tp.pred = tp.pred[[i]],
                                         calib.type = "mlr",
                                         w.covs = c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45"))
  }
  
}


### Save image
save.image(paste("data/sim/large_sample_analysis_moderate_mlr_", scen, ".RData", sep = ""))
print("IMAGE SAVED")

