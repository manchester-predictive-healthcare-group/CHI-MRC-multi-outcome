###
### This program apply censoring.
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
library(dplyr)
library(simsurv)
library(mstate)

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
print(paste("cens = ", scen, sep = " "))

### Source functions and input parameters for simulation
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")
source("code/p2/0_input_parameters.R") 

### Read in population data

### Write function to read in data and assign patids
extract_function <- function(seed){
  data <- readRDS(paste("data/sim/popdata/dgm1_seed", seed, ".rds", sep = ""))
  data$patid <- ((seed-1)*1000+1):(seed*1000)
  return(data)
}

### Extract data
data.pop.raw <- do.call("rbind", lapply(1:1000, extract_function))
str(data.pop.raw)

### Set seed for censoring which happens at random
set.seed(101)

###
### Apply censoring and convert to format that can be analysed using mstate
###
data.mstate.obj <- convert.mstate.integer.DGM1.cens(cohort.in = data.pop.raw, 
                                                    max.follow = max.follow, 
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

print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))


###
### Save image
saveRDS(data.mstate.obj[["data.mstate"]], paste("data/sim/popdata/pop.mstate.dgm1.", scen, ".RData", sep = ""))
saveRDS(data.mstate.obj[["data.raw"]], paste("data/sim/popdata/pop.raw.dgm1.", scen, ".RData", sep = ""))
saveRDS(data.mstate.obj[["tmat"]], paste("data/sim/popdata/tmat.dgm1.RData", sep = ""))
print("FINISHED")


