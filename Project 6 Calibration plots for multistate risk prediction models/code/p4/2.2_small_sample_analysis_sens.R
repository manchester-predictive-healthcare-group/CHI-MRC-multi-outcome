###
### This program runs the assessment of mean calibration in the small sample simulation, using no weights, and no percentile groups
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
n.cohort <- as.numeric(args[2])
set <- as.numeric(args[3])
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("set = ", set, sep = ""))

### Load packages
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Load prepped data
load(paste("data/sim/small_sample_prep_data_", scen, ".RData", sep = ""))

### Choose number of simulations
n.sim <- 1000/50

### Set seed
set.seed(set)

###
### Calibration using binary logistic regression at time point t (with and without IPCW)
###
print(paste("START CALIB BLR ", Sys.time(), sep = ""))

### Create list to store calibration of each estimates
### (include a list to store true calibration in each simulation iteration)
calib.blr.mean.list <- vector("list", 3)
calib.mlr.mean.list <- vector("list", 3)
calib.aj.mean.list <- vector("list", 3)
calib.true.mean.list <- vector("list", 3)

### Calculate the calibration intecepts and slopes
for (i in 1:3){
  
  ### BLR-IPCW
  calib.blr.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
  colnames(calib.blr.mean.list[[i]]) <- paste("state", 1:5, sep = "")
  
  ### MLR-IPCW
  calib.mlr.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
  colnames(calib.mlr.mean.list[[i]]) <- paste("state", 1:5, sep = "")
  
  ### Pseudo-value
  calib.aj.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
  colnames(calib.aj.mean.list[[i]]) <- paste("state", 1:5, sep = "")
  
  ### True
  calib.true.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
  colnames(calib.true.mean.list[[i]]) <- paste("state", 1:5, sep = "")
  
}


###
### Run the simulation loop
###
for (sim in 1:n.sim){
  
  print(paste("cohort size = ", n.cohort, ", sim =", sim, Sys.time(), sep = " "))
  
  ### Choose patids at random
  patids.vec <- sample((1:1000000), n.cohort, replace = FALSE)
  # We exclude individual 390375, as this person was removed from M1C2 dataset as they caused an error
  
  ### Extract individuals
  data.mstate.reduc <- data.mstate[data.mstate$id %in% patids.vec, ]
  data.raw.reduc <- data.raw[data.raw$id %in% patids.vec, ]

  ###
  ### Calculate calibration for each type of predicted risk
  ###
  for (i in 1:3){
 
    ### Extract predicted risks
    tp.pred[[i]]$id <- 1:nrow(tp.pred[[i]])
    tp.pred.reduc <- tp.pred[[i]][tp.pred[[i]]$id %in% patids.vec, ] |>
      dplyr::select(!id)
    
    if (scen == "C1"){
      
      ###
      ### Estimate mean calibration with BLR-IPCW
      calib.mean <- calib_msm(data.mstate = data.mstate.reduc,
                              data.raw = data.raw.reduc,
                              j = 1,
                              s = 0,
                              t = t.eval,
                              tp.pred = tp.pred.reduc,
                              calib.type = "blr",
                              weights = rep(1, nrow(data.raw.reduc)),
                              assess.moderate = FALSE,
                              assess.mean = TRUE)[["mean"]] 
      
      ### Assign output
      calib.blr.mean.list[[i]][sim, ] <- calib.mean
      
      ###
      ### Estimate mean calibration with MLR-IPCW
      calib.mean <- calib_msm(data.mstate = data.mstate.reduc,
                              data.raw = data.raw.reduc,
                              j = 1,
                              s = 0,
                              t = t.eval,
                              tp.pred = tp.pred.reduc,
                              calib.type = "mlr",
                              weights = rep(1, nrow(data.raw.reduc)),
                              assess.moderate = FALSE,
                              assess.mean = TRUE)[["mean"]] 
      
      ### Assign output
      calib.mlr.mean.list[[i]][sim, ] <- calib.mean
      
      ###
      ### Estimate mean calibration with Aalen-Johansen
      calib.mean <- calib_msm(data.mstate = data.mstate.reduc,
                              data.raw = data.raw.reduc,
                              j = 1,
                              s = 0,
                              t = t.eval,
                              tp.pred = tp.pred.reduc,
                              calib.type = "aj",
                              assess.moderate = FALSE,
                              assess.mean = TRUE)[["mean"]]
      
      ### Assign output
      calib.aj.mean.list[[i]][sim, ] <- calib.mean
      
    } else if (scen %in% c("C2", "C3")){
      
      ###
      ### Estimate mean calibration with BLR-IPCW
      calib.mean <- calib_msm(data.mstate = data.mstate.reduc,
                              data.raw = data.raw.reduc,
                              j = 1,
                              s = 0,
                              t = t.eval,
                              tp.pred = tp.pred.reduc,
                              calib.type = "blr",
                              weights = rep(1, nrow(data.raw.reduc)),
                              assess.moderate = FALSE,
                              assess.mean = TRUE)[["mean"]] 
      
      ### Assign output
      calib.blr.mean.list[[i]][sim, ] <- calib.mean
      
      ###
      ### Estimate mean calibration with MLR-IPCW
      calib.mean <- calib_msm(data.mstate = data.mstate.reduc,
                              data.raw = data.raw.reduc,
                              j = 1,
                              s = 0,
                              t = t.eval,
                              tp.pred = tp.pred.reduc,
                              calib.type = "mlr",
                              weights = rep(1, nrow(data.raw.reduc)),
                              assess.moderate = FALSE,
                              assess.mean = TRUE)[["mean"]] 
      
      ### Assign output
      calib.mlr.mean.list[[i]][sim, ] <- calib.mean
      
      ###
      ### Estimate mean calibration with Aalen-Johansen
      calib.mean <- calib_msm(data.mstate = data.mstate.reduc,
                              data.raw = data.raw.reduc,
                              j = 1,
                              s = 0,
                              t = t.eval,
                              tp.pred = tp.pred.reduc,
                              calib.type = "aj",
                              assess.moderate = FALSE,
                              assess.mean = TRUE)[["mean"]]
      
      ### Assign output
      calib.aj.mean.list[[i]][sim, ] <- calib.mean
      
    }
    
    ###
    ### Save true calibration
    calib.true.mean.list[[i]][sim, ] <- colMeans(data.raw.reduc[, paste("p.true", 1:5, sep = "")] - tp.pred.reduc)
    
  }
  
}


###
### Save image
###
rm(list = setdiff(ls(), list("scen", "n.cohort", "n.sim", "set",
                             "calib.blr.mean.list",
                             "calib.mlr.mean.list",
                             "calib.aj.mean.list",
                             "calib.true.mean.list")))

save.image(paste("data/sim/small_sample_mean_analysis_sens_", scen, "_n", n.cohort, "_set", set, ".RData", sep = ""))
print("IMAGE SAVED")



