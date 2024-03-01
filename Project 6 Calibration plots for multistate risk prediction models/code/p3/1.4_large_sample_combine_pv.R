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
n.pctls <- as.numeric(args[2])

### Print these
print(paste("scen = ", scen, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))

### Define n.cohort
n.cohort <- 200000

### We calculated the pseudo-values for sets of 500 individuals at a time, so need to calculate the number of sets of individuals
n.set <- (n.cohort/n.pctls)/500

### Assign max.state
max.state <- max(5)

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
      load(paste("data/sim/pv/large_sample_calc_pv_", scen, 
                 "_S", state,"_npctls", n.pctls, "_pctl", pctl, "_set", set, ".RData", sep = ""))
      ### Concatenate into dataset
      pv.out[[state]] <- rbind(pv.out[[state]], pv.out.set)
    }
  }
  
  ### Assign column names
  colnames(pv.out[[state]]) <- c("patid", paste("pv.state", state, sep = ""))
}

### Get str
str(pv.out)

### Order each of the data frames and concatenate
for (state in 1:max.state){
  pv.out[[state]] <- arrange(pv.out[[state]], patid)
}

### Remove patid from all but first data frame
for (state in 2:max.state){
  pv.out[[state]] <- pv.out[[state]][,2]
}

### Combine into a single dataframe
pv.comb <- do.call("cbind", pv.out)
colnames(pv.comb) <- c("patid", paste("pv.state", 1:max.state, sep = ""))
colnames(pv.comb)
head(pv.comb)

### Save image
rm(list=setdiff(ls(), list("pv.comb", "n.cohort", "scen", "n.pctls")))
save.image(paste("data/sim/large_sample_combine_pv_", scen, "_npctls", n.pctls, ".RData", sep = ""))
