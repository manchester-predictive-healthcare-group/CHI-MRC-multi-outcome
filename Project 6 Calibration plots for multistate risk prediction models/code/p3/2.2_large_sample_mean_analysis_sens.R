###
### This program will estimate mean calibration without grouping individuals by predicted risk (AJ), and using no weights (BLR-IPCW/MLR-IPCW)
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
print(paste("scen = ", scen, sep = ""))

### Load simulated data
load(paste("data/sim/large_sample_prep_data_", scen, ".RData", sep = ""))

### Load packages
source("code/z_functions.R")
source("code/z_load_packages.R")
source("code/p2/0_input_parameters.R")

### Choose number of bootstrap samples
CI.R.boot <- 200

### Define a function which wil calculate the weights wrong by assinging every weight to be 1
weights_manual <- function(data.mstate, data.raw, covs, t, s, landmark.type, j, max.weight, stabilised, max.follow){
  return(data.frame("id" = data.raw$id, "ipcw" = rep(1, nrow(data.raw))))
}

###
### Estimate mean calibration using each of the approaches
###

###
### blr-ipcw
print(paste("START blr-ipcw", Sys.time()))
cl <- makeCluster(3)
registerDoParallel(3)
calib.mean.blr <-(foreach(input=1:3,
                          .packages=c("calibmsm")) %dopar%{
                            
                            if (scen == "C1"){
                              out <- calib_msm(data.mstate = data.mstate,
                                               data.raw = data.raw,
                                               j = 1,
                                               s = 0,
                                               t = t.eval,
                                               tp.pred = tp.pred[[input]],
                                               calib.type = "blr",
                                               w.function = weights_manual,
                                               CI = 95,
                                               CI.R.boot = CI.R.boot,
                                               assess.moderate = FALSE,
                                               assess.mean = TRUE) 
                            } else {
                              out <- calib_msm(data.mstate = data.mstate,
                                               data.raw = data.raw,
                                               j = 1,
                                               s = 0,
                                               t = t.eval,
                                               tp.pred = tp.pred[[input]],
                                               calib.type = "blr",
                                               w.function = weights_manual,
                                               CI = 95,
                                               CI.R.boot = CI.R.boot,
                                               assess.moderate = FALSE,
                                               assess.mean = TRUE) 
                            }
                            
                            ### Return output
                            out
                          })
stopCluster(cl)
str(calib.mean.blr)


### Save output
save.image(paste("data/sim/large_sample_mean_sens_", scen, ".RData", sep = ""))

###
### mlr-ipcw
print(paste("START mlr-ipcw", Sys.time()))
cl <- makeCluster(3)
registerDoParallel(3)
calib.mean.mlr <-(foreach(input=1:3, 
                         .packages=c("calibmsm")) %dopar%{
                           
                           if (scen == "C1"){
                             out <- calib_msm(data.mstate = data.mstate,
                                              data.raw = data.raw,
                                              j = 1,
                                              s = 0,
                                              t = t.eval,
                                              tp.pred = tp.pred[[input]],
                                              calib.type = "mlr",
                                              w.function = weights_manual,
                                              CI = 95,
                                              CI.R.boot = CI.R.boot,
                                              assess.moderate = FALSE,
                                              assess.mean = TRUE) 
                           } else {
                             out <- calib_msm(data.mstate = data.mstate,
                                              data.raw = data.raw,
                                              j = 1,
                                              s = 0,
                                              t = t.eval,
                                              tp.pred = tp.pred[[input]],
                                              calib.type = "mlr",
                                              w.function = weights_manual,
                                              CI = 95,
                                              CI.R.boot = CI.R.boot,
                                              assess.moderate = FALSE,
                                              assess.mean = TRUE) 
                           }
                           
                           ### Return output
                           out
                         })
stopCluster(cl)
str(calib.mean.mlr)


### Save output
save.image(paste("data/sim/large_sample_mean_sens_", scen, ".RData", sep = ""))

###
### Aalen-Johansen
print(paste("START AJ", Sys.time()))
cl <- makeCluster(3)
registerDoParallel(3)
calib.mean.aj <-(foreach(input=1:3, 
                         .packages=c("calibmsm")) %dopar%{
                           
                           if (scen == "C1"){
                             out <- calib_msm(data.mstate = data.mstate,
                                              data.raw = data.raw,
                                              j = 1,
                                              s = 0,
                                              t = t.eval,
                                              tp.pred = tp.pred[[input]],
                                              calib.type = "aj",
                                              CI = 95,
                                              CI.type = "bootstrap",
                                              CI.R.boot = CI.R.boot,
                                              assess.moderate = FALSE,
                                              assess.mean = TRUE) 
                           } else {
                             out <- calib_msm(data.mstate = data.mstate,
                                              data.raw = data.raw,
                                              j = 1,
                                              s = 0,
                                              t = t.eval,
                                              tp.pred = tp.pred[[input]],
                                              calib.type = "aj",
                                              CI = 95,
                                              CI.type = "bootstrap",
                                              CI.R.boot = CI.R.boot,
                                              assess.moderate = FALSE,
                                              assess.mean = TRUE) 
                           }
                           
                           ### Return output
                           out
                         })
stopCluster(cl)
str(calib.mean.aj)

### Save output
save.image(paste("data/sim/large_sample_mean_sens_", scen, ".RData", sep = ""))
print("IMAGE SAVED")

