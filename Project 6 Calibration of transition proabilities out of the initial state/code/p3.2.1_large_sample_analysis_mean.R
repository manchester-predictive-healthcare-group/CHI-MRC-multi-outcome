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
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))


### Load prepped dat afor pv (this contains the appropriate datasets, and also estimate of Aalen-Johansen)
load(paste("data/large_sample_analysis_prep_N", n.cohort, "_", scen, ".RData", sep = ""))
source("code/z_functions.R")

### Choose number of bootstrap samples
n.boot <- 200

###
### Extract the true risks into their own object
###
p.true <- data.raw.reduc %>% select(paste("p.true", 1:5, sep = ""))

###
### Create miss-calibrated risks (although they still need to sum to 1...)
### Need to get creative for this
###
p.est.perf <- p.true
p.est.over.mean <- exp((log(p.true/(1-p.true)) + 0.5))/(1 + exp((log(p.true/(1-p.true)) + 0.5)))
p.est.under.mean <- exp((log(p.true/(1-p.true)) - 0.5))/(1 + exp((log(p.true/(1-p.true)) - 0.5)))
p.est.4.mean <- exp((log(p.true/(1-p.true)) + c(-0.5, -0.5, 0, 0.5, 0.5)))/(1 + exp((log(p.true/(1-p.true)) + c(-0.5, -0.5, 0, 0.5, 0.5))))
p.est.4.mean <- p.est.4.mean/rowSums(p.est.4.mean)

print(colMeans(p.est.over.mean) - colMeans(p.true))
print(colMeans(p.est.under.mean) - colMeans(p.true))

### ALSO WANT TO HAVE SEPERATE DEVIATIONS WHEN CONSIDERING CHANGES TO MEAN CALIBRATION AND MODERATE CALIBRATION
### Put these into a list
p.est.mean <- list(p.est.perf, p.est.over.mean, p.est.under.mean, p.est.4.mean)

###
### Calibration using binary logistic regression at time point t (with and without IPCW)
###
print(paste("START CALIB BLR ", Sys.time(), sep = ""))

### Create list to store calibration of each estimates
calib.blr.list <- vector("list", 4)
calib.blr.ipcw.list <- vector("list", 4)
calib.blr.ipcw.boot.se.list <- vector("list", 4)

### Calculate the calibration intecepts and slopes
for (i in 1:3){
  print(paste("start blr", i, Sys.time()))
  calib.blr.list[[i]] <- calc.calib.blr(data.mstate.reduc,
                                        data.raw.reduc,
                                        t.eval = t.eval,
                                        p.est = p.est.mean[[i]])
  
  print(paste("start blr ipcw", i, Sys.time()))
  calib.blr.ipcw.list[[i]] <- calc.calib.blr.ipcw(data.mstate.reduc,
                                                  data.raw.reduc,
                                                  t.eval = t.eval,
                                                  p.est = p.est.mean[[i]])
  
  print(paste("start blr boot", i, Sys.time()))
  calib.blr.ipcw.boot.se.list[[i]] <- calc.calib.blr.ipcw.boot.se(data.mstate.reduc,
                                                                  data.raw.reduc,
                                                                  t.eval = t.eval,
                                                                  p.est = p.est.mean[[i]],
                                                                  n.boot = n.boot)
}


print(paste("FINISH CALIB BLR ", Sys.time(), sep = ""))

### Save image
save.image(paste("data/large_sample_analysis_mean_N", n.cohort, "_", scen, ".RData", sep = ""))
print("IMAGE SAVED")


###
### Calibration using nominal recalibration framework of van Hoorde, at time point t
###
print(paste("START CALIB MLR ", Sys.time(), sep = ""))

### Create list to store calibration of each estimates
calib.mlr.list <- vector("list", 4)
calib.mlr.ipcw.list <- vector("list", 4)
calib.mlr.ipcw.boot.se.list <- vector("list", 4)

### Calculate the calibration intecepts and slopes
for (i in 1:3){
  print(paste("start mlr", i, Sys.time()))
  calib.mlr.list[[i]] <- calc.calib.mlr(data.mstate.reduc,
                                        data.raw.reduc,
                                        t.eval = t.eval,
                                        p.est = p.est.mean[[i]])
  
  print(paste("start mlr ipcw", i, Sys.time()))
  calib.mlr.ipcw.list[[i]] <- calc.calib.mlr.ipcw(data.mstate.reduc,
                                                  data.raw.reduc,
                                                  t.eval = t.eval,
                                                  p.est = p.est.mean[[i]])
  
  print(paste("start mlr boot", i, Sys.time()))
  calib.mlr.ipcw.boot.se.list[[i]] <- calc.calib.mlr.ipcw.boot.se(data.mstate.reduc,
                                                                  data.raw.reduc,
                                                                  t.eval = t.eval,
                                                                  p.est = p.est.mean[[i]],
                                                                  n.boot = n.boot)
}

print(paste("FINISH CALIB MLR ", Sys.time(), sep = ""))

### Save image
save.image(paste("data/large_sample_analysis_mean_N", n.cohort, "_", scen, ".RData", sep = ""))
print("IMAGE SAVED")


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
rm(list = setdiff(ls(), list("calib.aj.list", 
                             "calib.mlr.ipcw.boot.se.list", "calib.mlr.ipcw.list", "calib.mlr.list",
                             "calib.blr.ipcw.boot.se.list", "calib.blr.ipcw.list", "calib.blr.list",
                             "obs.aj", "obs.aj.se", "obs.aj.pctls.mean",
                             "p.true", "p.est.mean",
                             "scen", "n.boot", "n.cohort")))

save.image(paste("data/large_sample_analysis_mean_N", n.cohort, "_", scen, ".RData", sep = ""))
print("IMAGE SAVED")
  
