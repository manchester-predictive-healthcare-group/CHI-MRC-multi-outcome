###
### This program will combine the pseudo-values into one workspace
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")
source("code/z_functions_ce.R")

### Extract sample size, number of percentiles/groups, and task id from command line
args <- commandArgs(trailingOnly = T)
n.devel <- as.numeric(args[1])
n.pctls <- as.numeric(args[2])

### Print these
print(paste("n.devel = ", n.devel, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))

### We will calculate the pseudo-values for 500 individuals at a time, and need to calculate the number of sets of individuals there will be in each quantile
n.set <- (100000/n.pctls)/500

### Define a dataset that includes all combinations of simulation parameters (i.e. simulation cases)
sims_parameters <- crossing(
  scen = factor(c("M1C1", "M1C2", "M1C3"), levels = c("M1C1", "M1C2", "M1C3")),
  state = 1:6,
  pctl = 1:n.pctls,
  set = 1:n.set)

### Assign max.state
max.state <- max(sims_parameters$state)

### Create an object to store the pseudo-values
### Each list element will contain the pseudo-values for a given state
pv.out <- vector("list", max.state)

for (state in 1:max.state){
  
  ### Print progress
  print(paste("state =", state, Sys.time()))
  
  ### Create data frame to store pseudo-values
  pv.out[[state]] <- data.frame(matrix(NA, nrow = 0, ncol = 2))
  
  ### Combine them
  for (pctl in 1:n.pctls){
    for (set in 1:n.set){
      ### Load data
      load(paste("data/ce/ce_pv_calc_set_N", n.devel, "_S", state, "_npctls", n.pctls, "_pctl", pctl, "_set", set, ".RData", sep = ""))
      ### Concatenate into dataset
      pv.out[[state]] <- rbind(pv.out[[state]], pv.out.set)
    }
  }
  
  ### Assign column names
  colnames(pv.out[[state]]) <- c("person_id", paste("pv.state", state, sep = ""))
}

### Order each of the data frames and concatenate
for (state in 1:max.state){
  pv.out[[state]] <- arrange(pv.out[[state]], person_id)
}

### Remove person_id from all but first data frame
for (state in 2:max.state){
  pv.out[[state]] <- pv.out[[state]][,2]
}

### Combine into a single dataframe
pv.comb <- do.call("cbind", pv.out)
colnames(pv.comb) <- c("person_id", paste("pv.state", 1:max.state, sep = ""))
colnames(pv.comb)
head(pv.comb)

### Save image
rm(list=setdiff(ls(), list("pv.comb", "n.devel", "n.pctls")))
save.image(paste("data/ce/ce_pv_combine_N", n.devel, "_npctls", n.pctls, ".RData", sep = ""))



