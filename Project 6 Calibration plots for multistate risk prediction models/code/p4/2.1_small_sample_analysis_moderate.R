###
### This program runs the assessment of moderate calibration in the small sample simulation
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
n.cohort <- as.numeric(args[2])
n.pctls <- as.numeric(args[3])
set <- as.numeric(args[4])
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))
print(paste("set = ", set, sep = ""))

### Load packages
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Load prepped data
load(paste("data/sim/small_sample_prep_data_", scen, ".RData", sep = ""))

### Choose number of simulations
n.sim <- 200/50

### Set seed
set.seed(set)

### Create list to store calibration of each estimates
### (include a list to store true calibration in each simulation iteration)
calib.blr.mod.list <- vector("list", 3)
calib.mlr.mod.list <- vector("list", 3)
calib.pv.mod.list <- vector("list", 3)

### Calculate the calibration intecepts and slopes
for (i in 1:3){
  
  ### BLR-IPCW
  calib.blr.mod.list[[i]] <- vector("list", n.sim)
  
  ### MLR-IPCW
  calib.mlr.mod.list[[i]] <- vector("list", n.sim)
  
  ### Pseudo-value
  calib.pv.mod.list[[i]] <- vector("list", n.sim)
  
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
      ### Estimate mod calibration with BLR-IPCW
      calib.mod <- calib_msm(data.mstate = data.mstate.reduc,
                             data.raw = data.raw.reduc,
                             j = 1,
                             s = 0,
                             t = t.eval,
                             tp.pred = tp.pred.reduc,
                             calib.type = "blr",
                             curve.type = "rcs",
                             rcs.nk = 4,
                             assess.moderate = TRUE,
                             assess.mean = FALSE)[["plotdata"]] 
      
      ### Assign output
      calib.blr.mod.list[[i]][[sim]] <- calib.mod
      
      ###
      ### Estimate mod calibration with MLR-IPCW
      calib.mod <- calib_msm(data.mstate = data.mstate.reduc,
                             data.raw = data.raw.reduc,
                             j = 1,
                             s = 0,
                             t = t.eval,
                             tp.pred = tp.pred.reduc,
                             calib.type = "mlr",
                             assess.moderate = TRUE,
                             assess.mean = FALSE)[["plotdata"]] 
      
      ### Assign output
      calib.mlr.mod.list[[i]][[sim]] <- calib.mod
      
      ###
      ### Estimate mod calibration with Aalen-Johansen
      calib.mod <- calib_msm(data.mstate = data.mstate.reduc,
                             data.raw = data.raw.reduc,
                             j = 1,
                             s = 0,
                             t = t.eval,
                             tp.pred = tp.pred.reduc,
                             calib.type = "pv",
                             curve.type = "rcs",
                             rcs.nk = 4,
                             assess.moderate = TRUE,
                             assess.mean = FALSE)[["plotdata"]]
      
      ### Assign output
      calib.pv.mod.list[[i]][[sim]] <- calib.mod
      
    } else if (scen %in% c("C2", "C3")){
      
      ###
      ### Estimate mod calibration with BLR-IPCW
      calib.mod <- calib_msm(data.mstate = data.mstate.reduc,
                             data.raw = data.raw.reduc,
                             j = 1,
                             s = 0,
                             t = t.eval,
                             tp.pred = tp.pred.reduc,
                             calib.type = "blr",
                             curve.type = "rcs",
                             rcs.nk = 4,
                             w.covs = c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45"),
                             assess.moderate = TRUE,
                             assess.mean = FALSE)[["plotdata"]] 
      
      ### Assign output
      calib.blr.mod.list[[i]][[sim]] <- calib.mod
      
      ###
      ### Estimate mod calibration with MLR-IPCW
      calib.mod <- calib_msm(data.mstate = data.mstate.reduc,
                             data.raw = data.raw.reduc,
                             j = 1,
                             s = 0,
                             t = t.eval,
                             tp.pred = tp.pred.reduc,
                             calib.type = "mlr",
                             w.covs = c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45"),
                             assess.moderate = TRUE,
                             assess.mean = FALSE)[["plotdata"]] 
      
      ### Assign output
      calib.mlr.mod.list[[i]][[sim]] <- calib.mod
      
      ###
      ### Estimate mod calibration with Aalen-Johansen
      calib.mod <- calib_msm(data.mstate = data.mstate.reduc,
                             data.raw = data.raw.reduc,
                             j = 1,
                             s = 0,
                             t = t.eval,
                             tp.pred = tp.pred.reduc,
                             calib.type = "pv",
                             curve.type = "rcs",
                             rcs.nk = 4,
                             pv.n.pctls = n.pctls,
                             assess.moderate = TRUE,
                             assess.mean = FALSE)[["plotdata"]]
      
      ### Assign output
      calib.pv.mod.list[[i]][[sim]] <- calib.mod
      
    }
  }
}

str(calib.blr.mod.list)
str(calib.mlr.mod.list)
str(calib.pv.mod.list)

###
### Save image
###
rm(list = setdiff(ls(), list("scen", "n.cohort", "n.sim", "n.pctls", "set",
                             "calib.blr.mod.list",
                             "calib.mlr.mod.list",
                             "calib.pv.mod.list")))

save.image(paste("data/sim/small_sample_moderate_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, "_set", set, ".RData", sep = ""))
print("IMAGE SAVED")



