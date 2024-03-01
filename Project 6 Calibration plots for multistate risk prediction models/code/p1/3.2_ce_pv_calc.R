###
### This program will calculate pseudo-values for each patient in the validation cohort
###

### Note that this has to be run for a range of values of 'set', which should be parallelised
### i.e. the program calculates transition probabilities for groups of set.size individual.
### 'set' should take the values 1:100000/set.size
### In this simulation we used set.size = 100

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
n.devel <- as.numeric(args[1])
n.pctls <- as.numeric(args[2])
task.id <- as.numeric(args[3])

### Load prepped data to calculate pseudo values with
load(paste("data/ce/ce_pv_prep_N", n.devel, "_npctls", n.pctls, ".RData", sep = ""))
source("code/z_functions_ce.R")

### Define number of states
max.state <- max(data.mstate$to)

### Assign set.size
set.size <- 500

### We will calculate the pseudo-values for 500 individuals at a time, and need to calculate the number of sets of individuals there will be in each quantile
n.set <- (100000/n.pctls)/set.size

### Define a dataset that includes all combinations of simulation parameters (i.e. simulation cases)
sims_parameters <- crossing(
  state = 1:max.state,
  pctl = 1:n.pctls,
  set = 1:n.set)

### Assign variables based on table
state <- sims_parameters$state[task.id]
pctl <- sims_parameters$pctl[task.id]
set <- sims_parameters$set[task.id]

### Print these
print(paste("n.devel = ", n.devel, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))
print(paste("state = ", state, sep = ""))
print(paste("pctl = ", pctl, sep = ""))
print(paste("set = ", set, sep = ""))

### Extract the person_ids for individuals in the desired pctl (sorted by state of interest)
person_ids <- data.sort.pctls[[state]][[pctl]]$person_id[((set-1)*set.size+1):(set*set.size)]

### Calculate a pseudo-value for each individual
print(paste("START PV ", Sys.time()))
pv.temp <- lapply(person_ids, func.calc.pv.aj.efficient.ce,
                  data.mstate = subset(data.mstate, person_id %in% data.sort.pctls[[state]][[pctl]]$person_id),
                  obs.aj = obs.aj.pctls[[state]][[pctl]],
                  tmat = tmat,
                  n.cohort = nrow(data.sort.pctls[[state]][[pctl]]),
                  t.eval = t.eval,
                  pv.same = pv.same.pctls[[state]][[pctl]])
print(paste("END PV ", Sys.time()))

### Only interested in pseudo-values for the state of interest (which we sorted by)
pv.temp <- do.call("rbind", pv.temp)[, state]

### Combine with person_ids, so we have a record of who these pseudo-values are for
pv.out.set <- data.frame(cbind(person_ids, pv.temp))
colnames(pv.out.set) <- c("person_id", paste("pv.", state, sep = ""))

### Remove everything except the pseudo-values
rm(list=setdiff(ls(), list("pv.out.set", "n.devel", "state", "n.pctls", "pctl", "set")))

### Save image
save.image(paste("data/ce/ce_pv_calc_set_N", n.devel, "_S", state,"_npctls", n.pctls, "_pctl", pctl, "_set", set, ".RData", sep = ""))