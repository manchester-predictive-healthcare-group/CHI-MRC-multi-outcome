###
### This program will load all the pseudo values and combine into a datset to be loaded
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")
source("code/z_functions.R")

### Define n.cohort
n.cohort <- 200000

### Extract sample size, number of percentiles/groups, and task id from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]

### Print these
print(paste("scen = ", scen, sep = ""))

### Assign set.size
set.size <- 500

### Assign n.cohort
n.cohort <- 200000

### We calculated the pseudo-values for sets of 500 individuals at a time, so need to calculate the number of sets of individuals
n.set <- (n.cohort/set.size)

### Create an object to store the pseudo-values
### Each list element will contain the pseudo-values for a given state
pv.comb <- data.frame(matrix(NA, nrow = 0, ncol = 6))

### Combine them
for (set in 1:n.set){
  ### Load data
  load(paste("data/sim/pv/large_sample_calc_pv_", scen, "_sens_set", set, ".RData", sep = ""))
  ### Concatenate into dataset
  pv.comb <- rbind(pv.comb, pv.out.set)
}

### Assign column names
colnames(pv.comb) <- c("patid", paste("pv.state", 1:5, sep = ""))

### Get str
str(pv.comb)

### Order each of the data frames and concatenate
pv.comb <- arrange(pv.comb, patid)

### Save image
rm(list=setdiff(ls(), list("pv.comb", "n.cohort", "scen")))
save.image(paste("data/sim/large_sample_combine_pv_sens_", scen, ".RData", sep = ""))
