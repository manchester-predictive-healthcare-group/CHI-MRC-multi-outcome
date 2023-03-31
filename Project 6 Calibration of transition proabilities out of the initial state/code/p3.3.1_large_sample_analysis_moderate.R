###
### This program will calculate calibration (or observed risk) for all methods except pseudo-values in the 
### large sample analysis
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

### Load prepped dat afor pv (this contains the appropriate datasets, and also estimate of Aalen-Johansen)
load(paste("data/large_sample_analysis_prep_N", n.cohort, "_", scen, ".RData", sep = ""))
source("code/z_functions.R")

### Print sample sizes
nrow(data.raw.reduc)
nrow(data.mstate.reduc)

###
### Extract the true risks into their own object
###
p.true <- data.raw.reduc %>% select(paste("p.true", 1:5, sep = ""))

###
### Create miss-calibrated risks (although they still need to sum to 1...)
### Need to get creative for this
###
p.est.perf <- p.true
p.est.over <- exp((log(p.true/(1-p.true)) + 0.5))/(1 + exp((log(p.true/(1-p.true)) + 0.5)))
p.est.under <- exp((log(p.true/(1-p.true)) - 0.5))/(1 + exp((log(p.true/(1-p.true)) - 0.5)))

### ALSO WANT TO HAVE SEPERATE DEVIATIONS WHEN CONSIDERING CHANGES TO MEAN CALIBRATION AND MODERATE CALIBRATION
### Put these into a list
p.est <- list(p.est.perf, p.est.over, p.est.under)


# ###
# ### Calibration using Aalen-Johansen estimator (need to upgrade this to produce a plot)
# ###
# 
# print(paste("START CALIB AJ ", Sys.time(), sep = ""))
# 
# ### Create list
# calib.aj <- vector("list", 3)
# 
# ### Calc calib
# for (ii in 1:3){
#   print(paste("calib aj = ", ii, Sys.time(), sep = " "))
#   calib.aj[[ii]] <- calc.calib.aj.moderate(data.mstate = data.mstate.reduc, 
#                                            data.raw = data.raw.reduc, 
#                                            tmat = tmat, 
#                                            t.eval = t.eval, 
#                                            p.est = p.est[[ii]], 
#                                            num.groups = 10)
#   print(warnings())
# }
# print(paste("FINISH CALIB AJ ", Sys.time(), sep = ""))

###
### Calibration using binary logistic regression at time point t (with and without IPCW)
###

print(paste("START CALIB BLR ", Sys.time(), sep = " "))

### Create list to store calibration of each estimates
calib.blr.list <- vector("list", 3)
calib.blr.ipcw.list <- vector("list", 3)

### Calculate the calibration intecepts and slopes
for (i in 1:3){
  print(paste("calib blr", i , Sys.time(), sep = " "))
  calib.blr.list[[i]] <- calc.calib.blr.mod(data.mstate.reduc,
                                        data.raw.reduc,
                                        t.eval = t.eval,
                                        p.est = p.est[[i]])
  print(warnings())
  print(paste("calib blr ipcw", i , Sys.time(), sep = " "))
  calib.blr.ipcw.list[[i]] <- calc.calib.blr.ipcw.mod(data.mstate.reduc,
                                                  data.raw.reduc,
                                                  t.eval = t.eval,
                                                  p.est = p.est[[i]])
  print(warnings())
}

print(paste("FINISH CALIB BLR ", Sys.time(), sep = " "))

### Save image
save.image(paste("data/large_sample_analysis_moderate_N", n.cohort, "_", scen, ".RData", sep = ""))
print("IMAGE SAVED")

###
### Calibration using nominal recalibration framework at time point t (with and without IPCW)
###

print(paste("START CALIB MLR ", Sys.time(), sep = " "))

### Create list to store calibration of each estimates
calib.mlr.list <- vector("list", 3)
calib.mlr.ipcw.list <- vector("list", 3)

### Calculate the calibration intecepts and slopes
for (i in 1:3){
  
  print(paste("calib mlr", i , Sys.time(), sep = " "))
  calib.mlr.list[[i]] <- calc.calib.mlr.mod(data.mstate.reduc,
                                            data.raw.reduc,
                                            t.eval = t.eval,
                                            p.est = p.est[[i]])
  print(warnings())
  print(paste("calib mlr ipcw", i , Sys.time(), sep = " "))
  calib.mlr.ipcw.list[[i]] <- calc.calib.mlr.ipcw.mod(data.mstate.reduc,
                                                      data.raw.reduc,
                                                      t.eval = t.eval,
                                                      p.est = p.est[[i]])
  print(warnings())
}

print(paste("FINISH CALIB MLR ", Sys.time(), sep = " "))

### Save image
save.image(paste("data/large_sample_analysis_moderate_N", n.cohort, "_", scen, ".RData", sep = ""))
print("IMAGE SAVED")

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

### Save image
rm(list=setdiff(ls(), list("data.mstate.reduc", "data.raw.reduc",
                           "calib.aj", "calib.blr.list", "calib.blr.ipcw.list", 
                           "calib.mlr.list", "calib.mlr.ipcw.list", "calib.pv.list",
                           "p.est", "p.true", "tmat", "t.eval", "scen", "n.cohort", "n.pctls")))
save.image(paste("data/large_sample_analysis_moderate_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")

