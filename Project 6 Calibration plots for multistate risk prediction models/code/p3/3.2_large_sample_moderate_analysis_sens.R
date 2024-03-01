###
### This program will calculate calibration (or observed risk) for BLR-IPCW and MLR-IPCW and
### assess sensitivity to the specification of the weights
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
# n.cohort <- as.numeric(args[2])
# n.pctls <- as.numeric(args[3])
print(paste("scen = ", scen, sep = ""))
# print(paste("n.cohort = ", n.cohort, sep = ""))
# print(paste("n.pctls = ", n.cohort, sep = ""))

### Load packages
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Load input parameters for that sfcenario (specifically the censoring ones)
source("code/p2/0_input_parameters.R") 

### Load prepped data
load(paste("data/sim/large_sample_prep_data_", scen, ".RData", sep = ""))

###
### Start by obtaining the actual probability of being censored by time t
###
weights.DGM <- calc_weights_DGM(data.raw, 
                                t.eval, 
                                max.weight = 20, 
                                cens_shape = 1, 
                                cens_scale = cens_scale, 
                                cens_beta_x12 = cens_beta_x12,
                                cens_beta_x13 = cens_beta_x13,
                                cens_beta_x15 = cens_beta_x15,
                                cens_beta_x24 = cens_beta_x24,
                                cens_beta_x25 = cens_beta_x25,
                                cens_beta_x34 = cens_beta_x34,
                                cens_beta_x35 = cens_beta_x35,
                                cens_beta_x45 = cens_beta_x45)


###
### Calibration using BLR-IPCW
###

### Define output objects
calib.moderate.blr.DGMw <- vector("list", 3)
calib.moderate.blr.NOw <- vector("list", 3)
calib.moderate.mlr.DGMw <- vector("list", 3)
calib.moderate.mlr.NOw <- vector("list", 3)

### Define predicted risks labels
miscal.labels <- c("cal1", "cal2", "cal3")

### Calculate the calibration plots
for (i in 1:3){
  print(paste(miscal.labels[i] , Sys.time(), sep = " "))
  
  ### BLR-IPCW with perfectly specified weights
  calib.moderate.blr.DGMw[[i]] <- calib_msm(data.mstate = data.mstate,
                                            data.raw = data.raw,
                                            j = 1,
                                            s = 0,
                                            t = t.eval,
                                            tp.pred = tp.pred[[i]],
                                            calib.type = "blr",
                                            curve.type = "rcs",
                                            rcs.nk = 5,
                                            weights = weights.DGM$ipcw.DGM)
  
  ### BLR-IPCW with no weights
  calib.moderate.blr.NOw[[i]] <- calib_msm(data.mstate = data.mstate,
                                           data.raw = data.raw,
                                           j = 1,
                                           s = 0,
                                           t = t.eval,
                                           tp.pred = tp.pred[[i]],
                                           calib.type = "blr",
                                           curve.type = "rcs",
                                           rcs.nk = 5,
                                           weights = rep(1, nrow(data.raw)))
  
  ### MLR-IPCW with perfectly specified weights
  calib.moderate.mlr.DGMw[[i]] <- calib_msm(data.mstate = data.mstate,
                                            data.raw = data.raw,
                                            j = 1,
                                            s = 0,
                                            t = t.eval,
                                            tp.pred = tp.pred[[i]],
                                            calib.type = "mlr",
                                            weights = weights.DGM$ipcw.DGM)
  
  ### MLR-IPCW with no weights
  calib.moderate.mlr.NOw[[i]] <- calib_msm(data.mstate = data.mstate,
                                           data.raw = data.raw,
                                           j = 1,
                                           s = 0,
                                           t = t.eval,
                                           tp.pred = tp.pred[[i]],
                                           calib.type = "mlr",
                                           weights = rep(1, nrow(data.raw)))

  
}

### Save image
save.image(paste("data/sim/large_sample_analysis_moderate_sens_", scen, ".RData", sep = ""))
print("IMAGE SAVED")

