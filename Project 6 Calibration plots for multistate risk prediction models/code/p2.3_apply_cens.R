###
### This program will generate a cohort of size 1,000,000 data according to DGM1
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen.DGM <- args[1]
scen.cens <- args[2]
print(paste("DGM =", scen.DGM, ", cens =", scen.cens, sep = " "))

### Load simulated data
if (scen.DGM == "M1"){
  load("data/generate_population_data_DGM1.RData")
} else if (scen.DGM == "M2"){
  load("data/generate_population_data_DGM2_2.RData")
} else if (scen.DGM == "M3"){
  load("data/generate_population_data_DGM2_4.RData")
}

### Source latest functions
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")


### Set the coefficients for the censoring mechanism
if (scen.cens == "C1"){
  cens_beta_x1 <- 0
  cens_beta_x2 <- 0
} else if (scen.cens == "C2"){
  cens_beta_x1 <- 0.25
  cens_beta_x2 <- -0.25
} else if (scen.cens == "C3"){
  cens_beta_x1 <- 1
  cens_beta_x2 <- -1
} else if (scen.cens == "C4"){
  cens_beta_x1 <- 3
  cens_beta_x2 <- 3
} else if (scen.cens == "C5"){
  cens_beta_x1 <- 1
  cens_beta_x2 <- 1
}


### Pick a t to evaluate at
t.eval <- ceiling(365.25*7)

### Extract the raw data for analysis
data.pop.raw <- data.pop[["cohort"]]

### Set eed for censoring which happens at random
set.seed(101)

###
### Apply censoring and convert to format that can be analysed using mstate
###

### Censoring is applied at random, at a scale resulting in 20% censored by 7 years follow up
print(paste("START CONVERT MSTATE ", Sys.time(), sep = ""))
cens_scale <- calc.scale(0.6)

### Load simulated data
if (scen.DGM == "M1"){
  data.mstate.obj <- convert.mstate.DGM1.cens(cohort.in = data.pop.raw, 
                                              max.follow = data.pop[["max.follow"]], 
                                              cens_shape = 1,
                                              cens_scale = cens_scale,
                                              cens_beta_x1 = cens_beta_x1,
                                              cens_beta_x2 = cens_beta_x2)
} else if (scen.DGM %in% c("M2", "M3")){
  data.mstate.obj <- convert.mstate.DGM2.combine45.cens(cohort.in = data.pop.raw, 
                                              max.follow = data.pop[["max.follow"]], 
                                              cens_shape = 1,
                                              cens_scale = cens_scale,
                                              cens_beta_x1 = cens_beta_x1,
                                              cens_beta_x2 = cens_beta_x2)
} 


print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))


###
### Save image
save.image(paste("data/apply_cens_", scen.DGM, scen.cens, ".RData", sep = ""))
print("IMAGE SAVED")

### There is an individual who is censored at the same time as an event (by complete chance), causing errors, need to remove
### Happens for "M1C2"
if (scen.DGM == "M1" & scen.cens == "C2"){
  data.mstate.obj[["data.mstate"]] <- data.mstate.obj[["data.mstate"]][data.mstate.obj[["data.mstate"]]$patid != 390375, ]
  data.mstate.obj[["data.raw"]] <- data.mstate.obj[["data.raw"]][data.mstate.obj[["data.raw"]]$patid != 390375, ]
  
  ###
  ### Save image
  save.image(paste("data/apply_cens_", scen.DGM, scen.cens, ".RData", sep = ""))
  print("IMAGE SAVED")
  
}
