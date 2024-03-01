###
### This program will calculate the pseudo-values for a set of individuals. This is done seperately based on a number of input parameters, 
### and should be parallelised over the variable set (from 1 to n.cohort/set.size)
### For the sensitivity analysis, we DO NOT group individuals by their predicted risk before calculating pseudo-values.
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions.R")
source("code/z_load_packages.R")

### Extract sample size, number of percentiles/groups, and task id from command line
args <- commandArgs(trailingOnly = T)
task.id <- as.numeric(args[1])

### Assign set.size
set.size <- 500

### Assign n.cohort
n.cohort <- 200000

### We will calculate the pseudo-values for 500 individuals at a time, and need to calculate the number of sets of individuals
n.set <- n.cohort/set.size

### Define a dataset that includes all combinations of simulation parameters (i.e. simulation cases)
sims_parameters <- crossing(
  scen = factor(c("C1", "C2", "C3"), levels = c("C1", "C2", "C3")),
  set = 1:n.set)

### Assign variables based on table
scen <- as.character(sims_parameters$scen[task.id])
set <- sims_parameters$set[task.id]

### Print these
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("set = ", set, sep = ""))

### Load the prepped Aalen-Johasen estimates, and pv.same estimates
load(paste("data/sim/large_sample_prep_data_pv_", scen, "_sens.RData", sep = ""))

### Extract the patids for individuals in the desired pctl (sorted by state of interest)
patids <- data.raw$patid[((set-1)*set.size+1):(set*set.size)]

### Calculate a pseudo-value for each individual
print(paste("START PV ", Sys.time()))
pv.temp <- lapply(patids, func.calc.pv.aj.efficient,
                  data.mstate = data.mstate,
                  obs.aj = obs.aj,
                  tmat = tmat,
                  n.cohort = nrow(data.raw),
                  t.eval = t.eval,
                  pv.same = pv.same)
print(paste("END PV ", Sys.time()))

### Only interested in pseudo-values for the state of interest (which we sorted by)
pv.temp <- do.call("rbind", pv.temp)

### Combine with patids, so we have a record of who these pseudo-values are for
pv.out.set <- data.frame(cbind(patids, pv.temp))
colnames(pv.out.set) <- c("patid", paste("pv.", 1:5, sep = ""))

### Remove everything except the pseudo-values
rm(list=setdiff(ls(), list("pv.out.set", "scen", "set", "n.cohort")))

### Save image
save.image(paste("data/sim/pv/large_sample_calc_pv_", scen, "_sens_set", set, ".RData", sep = ""))

