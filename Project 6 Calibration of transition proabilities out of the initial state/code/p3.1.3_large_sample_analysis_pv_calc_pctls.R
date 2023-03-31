###
### This program will calculate the pseudo-values for a set of individuals. This is done seperately based on a number of input parameters, 
### and should be parallelised over the variable set (from 1 to n.cohort/set.size)
### 1) Scenario
### 2) Size of cohort
### 3) Number of percentiles we grouped individuals by
### 4) State of interest we can to calculate the pseudo-values for
### 5) Which subgroup (pctl) we are calculating the pseudo-values for
### 6) Which set of individuals we are calculating for
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")
source("code/z_functions.R")

### Extract sample size, number of percentiles/groups, and task id from command line
args <- commandArgs(trailingOnly = T)
n.cohort <- as.numeric(args[1])
n.pctls <- as.numeric(args[2])
task.id <- as.numeric(args[3])

### Assign set.size
set.size <- 500

### We will calculate the pseudo-values for 500 individuals at a time, and need to calculate the number of sets of individuals there will be in each quantile
n.set <- (n.cohort/n.pctls)/set.size

### Define a dataset that includes all combinations of simulation parameters (i.e. simulation cases)
sims_parameters <- crossing(
  scen = factor(c("M1C1", "M1C2", "M1C3"), levels = c("M1C1", "M1C2", "M1C3")),
  state = 1:5,
  pctl = 1:n.pctls,
  set = 1:n.set)

### Assign variables based on table
scen <- as.character(sims_parameters$scen[task.id])
state <- sims_parameters$state[task.id]
pctl <- sims_parameters$pctl[task.id]
set <- sims_parameters$set[task.id]

### Print these
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))
print(paste("state = ", state, sep = ""))
print(paste("pctl = ", pctl, sep = ""))
print(paste("set = ", set, sep = ""))

### Load the prepped Aalen-Johasen estimates, and pv.same estimates
load(paste("data/large_sample_analysis_pv_prep_pctls_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))

### Extract the patids for individuals in the desired pctl (sorted by state of interest)
patids <- data.sort.pctls[[state]][[pctl]]$patid[((set-1)*set.size+1):(set*set.size)]

### Calculate a pseudo-value for each individual
print(paste("START PV ", Sys.time()))
pv.temp <- lapply(patids, func.calc.pv.aj.efficient,
                  data.mstate = subset(data.mstate.reduc, patid %in% data.sort.pctls[[state]][[pctl]]$patid),
                  obs.aj = obs.aj.pctls[[state]][[pctl]],
                  tmat = tmat,
                  n.cohort = nrow(data.sort.pctls[[state]][[pctl]]),
                  t.eval = t.eval,
                  pv.same = pv.same.pctls[[state]][[pctl]])
print(paste("END PV ", Sys.time()))

### Only interested in pseudo-values for the state of interest (which we sorted by)
pv.temp <- do.call("rbind", pv.temp)[, state]

### Combine with patids, so we have a record of who these pseudo-values are for
pv.out.set <- data.frame(cbind(patids, pv.temp))
colnames(pv.out.set) <- c("patid", paste("pv.", state, sep = ""))

### Remove everything except the pseudo-values
rm(list=setdiff(ls(), list("pv.out.set", "n.cohort", "scen", "state", "n.pctls", "pctl", "set")))

### Save image
save.image(paste("data/large_sample_analysis_pv_calc_pctls_N", n.cohort, "_", scen, 
                 "_S", state,"_npctls", n.pctls, "_pctl", pctl, "_set", set, ".RData", sep = ""))

