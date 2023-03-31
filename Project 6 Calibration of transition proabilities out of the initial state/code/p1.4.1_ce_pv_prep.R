###
### This program will prep the data in order to calculate pseudo-values
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")

### Extract number of pctls from command line
args <- commandArgs(trailingOnly = T)
n.devel <- as.numeric(args[1])
n.pctls <- as.numeric(args[2])
print(paste("n.devel = ", n.devel, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))

### Load work space with data in mstate format
load(paste("data/ce_combine_probtrans_N", n.devel, ".RData", sep = ""))
source("code/z_functions_ce.R")

### For ndevel <- 100000, the dataframe  was saved at pt. instead of pstate, so change this
if (n.devel == 100000){
  colnames(data.raw)[colnames(data.raw) %in% paste("pt.", 1:6, sep = "")] <- paste("pstate", 1:6, sep = "")
}

### Define number of states
max.state <- max(data.mstate$to)

### Create a list to store the sorted data, the Aalen-Johansen estimator for each group, 
### and the pseudo-value for an individual who doesnt have an event prior to time t.eval
data.sort.pctls <- vector("list", max.state)
obs.aj.pctls <- vector("list", max.state)
pv.same.pctls <- vector("list", max.state)

### Loop through the states of interest
for (state in 1:max.state){
  
  print(paste("state = ", state, Sys.time()))
  ###
  ### Create n.pctls datasets, grouped by predicted risk
  ###
  
  ### Arrange data by risk of category of interest
  data.raw.sort.temp <- arrange(data.raw, !!sym(paste("pstate", state, sep = "")))
  
  ### Create a list of length n.pctls
  data.sort.pctls[[state]] <- vector("list", n.pctls)
  group.size <- nrow(data.raw.sort.temp)/n.pctls
  for (i in 1:n.pctls){
    data.sort.pctls[[state]][[i]] <- slice(data.raw.sort.temp, (((i-1)*group.size)+1):(i*group.size))
  }
  
  ###
  ### Calculate the Aalen-Johansen estimator within each group
  ###
  
  ### Create a list to store them
  obs.aj.pctls[[state]] <- vector("list", n.pctls)
  
  ### Calculate obs.aj for each group
  for (i in 1:n.pctls){
    print(paste("Calc obs.aj = ", state, "pctl = ", i, Sys.time()))
    obs.aj.pctls[[state]][[i]] <- calc.calib.aj.ce(data.mstate = subset(data.mstate, person_id %in% data.sort.pctls[[state]][[i]]$person_id), 
                                                tmat = tmat, t.eval = t.eval)$obs.aj
  }
  
  
  ###
  ### Calculate the pseudo value for someone who is uncensored and has not had an event prior to time t
  ### This pseudo value will be the same for all patients who meet these criteria, so no need to recalculate
  ### it everytime. This will help with computational time.
  ###
  
  ### Create a vector to store the person_ids
  pv.same.pctls[[state]] <- vector("list", n.pctls)
  
  ### Calculate the person_ids
  for (i in 1:n.pctls){
    
    ### Start by identifying an individual for which this is true
    person_id.pv.same <- subset(data.mstate, person_id %in% data.sort.pctls[[state]][[i]]$person_id) %>%
      subset(from == 1 & to == 2 & Tstop > t.eval) %>%
      slice(1) %>%
      select(person_id) %>%
      as.numeric()
    
    ### Calculate the pseudo value for this individual
    print(paste("PV SAME STATE = ", state, "pctl = ", i, Sys.time()))
    pv.same.pctls[[state]][[i]] <- func.calc.pv.aj.ce(person_id.eval = person_id.pv.same, 
                                                   data.mstate = subset(data.mstate, person_id %in% data.sort.pctls[[state]][[i]]$person_id), 
                                                   obs.aj = obs.aj.pctls[[state]][[i]], 
                                                   tmat = tmat, 
                                                   n.cohort = nrow(data.sort.pctls[[state]][[i]]), 
                                                   t.eval = t.eval)
  }
  
}

### Save image
rm(list=setdiff(ls(), list("data.mstate", "data.raw", 
                           "data.sort.pctls", "pv.same.pctls", "obs.aj.pctls", "n.pctls", "max.state",
                           "tmat", "t.eval", "scen", "n.devel")))
save.image(paste("data/ce_pv_prep_N", n.devel, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")



