###
### This program will estimate the true calibration curve based off the true risks and the predicted risks using restricted cubic splines.
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### True calibration is the same irrespective of the censoring that has been applied, 
### so just need to load data for one scenario
load(paste("data/sim/large_sample_prep_data_C1.RData", sep = ""))
source("code/z_functions.R")

###
### Estimating the true calibration curves
###

### This will be done by regressing the predicted risks on the true risks, and generating predicted-observed risks using this model
### We can treat the true risks as pseudo-values and insert into calib_msm using pv.precalc to achieve this

### Repeat for rcs.nk = 3, 4 and 5
for (rcs.nk in c(3,4,5)){
  print(paste("rcs.nk", rcs.nk, Sys.time()))
  ### Define output object
  true.calib <- vector("list", 3)
  for (i in 1:3){
    ### Estimate true calibration plot
    true.calib[[i]] <- calib_msm(data.mstate = data.mstate, 
                                 data.raw = data.raw,
                                 j = 1,
                                 s = 0,
                                 t = t.eval,
                                 calib.type = "pv",
                                 tp.pred = tp.pred[[i]],
                                 pv.precalc = tp.pred[[1]],
                                 curve.type = "rcs",
                                 rcs.nk = rcs.nk)[["plotdata"]]
  }
  
  saveRDS(true.calib, paste("data/sim/true.calib.lineplot.rcs.nk", rcs.nk, ".rds", sep = ""))
  
}



rcs.nk = 3
### Define output object
true.calib <- vector("list", 3)
for (i in 1:3){
  ### Estimate true calibration plot
  true.calib[[i]] <- calib_msm(data.mstate = data.mstate, 
                               data.raw = data.raw,
                               j = 1,
                               s = 0,
                               t = t.eval,
                               calib.type = "pv",
                               tp.pred = tp.pred[[i]],
                               pv.precalc = tp.pred[[1]],
                               curve.type = "rcs",
                               rcs.nk = rcs.nk)[["plotdata"]]
}



