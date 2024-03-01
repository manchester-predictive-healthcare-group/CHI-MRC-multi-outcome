###
### This program will prepare the data to enable estimation of the pseudo-values in a computationally efficient manner
### 1) Pre-calculate Aalen-Johansen estimator for entire cohort (subgrouped) so this can be used within calculation of pseudo-values
### 2) Pre-calculate the pseudo-value for an individual who doesnt have an event until after t.eval, which will be the same for all individuals
### for which this is the case
###
### Note that we do this for each state seperately, as we want to order and subgroup the individuals differently depending on which
### state we can generating transition probabilities/pseudo-values for
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
n.pctls <- as.numeric(args[2])

print(paste("scen = ", scen, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))
### n.pctls is the number of groups we will split data into before calculating pseudo-values

### Load prepped data
load(paste("data/sim/large_sample_prep_data_", scen, ".RData", sep = ""))

### Load packages
source("code/z_functions.R")
source("code/z_load_packages.R")

### Define number of states
max.state <- max(data.mstate$to)

### Create a list to store the sorted data, the Aalen-Johansen estimator for each group, 
### and the pseudo-value for an individual who doesnt have an event prior to time t.eval
data.sort.pctls <- vector("list", max.state)
obs.aj.pctls <- vector("list", max.state)
pv.same.pctls <- vector("list", max.state)

### Loop through the states of interest
for (state in 1:max.state){
  
  ###
  ### Create n.pctls datasets, grouped by predicted risk
  ###
  
  ### Arrange data by risk of category of interest
  p.est.data.frame <- arrange(data.raw, !!sym(paste("p.true", state, sep = "")))
  
  ### Create a list of length n.pctls
  data.sort.pctls[[state]] <- vector("list", n.pctls)
  group.size <- nrow(p.est.data.frame)/n.pctls
  for (i in 1:n.pctls){
    data.sort.pctls[[state]][[i]] <- slice(p.est.data.frame, (((i-1)*group.size)+1):(i*group.size))
  }
  
  ###
  ### Calculate the Aalen-Johansen estimator within each group
  ###

  ### Create a list to store them
  obs.aj.pctls[[state]] <- vector("list", n.pctls)
  
  ### Calculate obs.aj for each group
  for (i in 1:n.pctls){
    print(paste("Calc obs.aj = ", state, "pctl = ", i, Sys.time()))
    obs.aj.pctls[[state]][[i]] <- calc.calib.aj(data.mstate = subset(data.mstate, patid %in% data.sort.pctls[[state]][[i]]$patid), 
                                       tmat = tmat, t.eval = t.eval)$obs.aj
  }
  
  
  ###
  ### Calculate the pseudo value for someone who is uncensored and has not had an event prior to time t
  ### This pseudo value will be the same for all patients who meet these criteria, so no need to recalculate
  ### it everytime. This will help with computational time.
  ###
  
  ### Create a vector to store the patids
  pv.same.pctls[[state]] <- vector("list", n.pctls)
  
  ### Calculate the patids
  for (i in 1:n.pctls){
    
    ### Start by identifying an individual for which this is true
    patid.pv.same <- subset(data.mstate, patid %in% data.sort.pctls[[state]][[i]]$patid) %>%
      subset(from == 1 & to == 2 & Tstop > t.eval) %>%
      slice(1) %>%
      select(patid) %>%
      as.numeric()
    
    ### Calculate the pseudo value for this individual
    print(paste("PV SAME STATE = ", state, "pctl = ", i, Sys.time()))
    pv.same.pctls[[state]][[i]] <- func.calc.pv.aj(patid.eval = patid.pv.same, 
                                          data.mstate = subset(data.mstate, patid %in% data.sort.pctls[[state]][[i]]$patid), 
                                          obs.aj = obs.aj.pctls[[state]][[i]], 
                                          tmat = tmat, 
                                          n.cohort = nrow(data.sort.pctls[[state]][[i]]), 
                                          t.eval = t.eval)
  }
  
}


### Save image
rm(list=setdiff(ls(), list("data.mstate", "data.raw", "p.est.data.frame",
                           "data.sort.pctls", "pv.same.pctls", "obs.aj.pctls", "n.pctls", "max.state",
                           "p.true", "tmat", "t.eval", "scen", "n.cohort")))
save.image(paste("data/sim/large_sample_prep_data_pv_", scen, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")


