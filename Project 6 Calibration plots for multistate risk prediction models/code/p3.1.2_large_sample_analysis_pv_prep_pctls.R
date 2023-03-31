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

### Load packages
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
n.cohort <- as.numeric(args[2])
n.pctls <- as.numeric(args[3])

print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))
### n.pctls is the number of groups we will split data into before calculating pseudo-values

### Load prepped data
load(paste("data/large_sample_analysis_prep_N", n.cohort, "_", scen, ".RData", sep = ""))
source("code/z_functions.R")

###
### Convert mstate data to integer format (note this function will only work for M1 scenarios) for computational reasons
### I have checked and this has minor impact on results (can compare obs.aj with obs.aj.integer)
###
data.mstate.reduc <- convert.mstate.integer.DGM1.cens.preexist(cohort.in = data.raw.reduc,
                                                               max.follow = ceiling(10*365.25) + 1)

### Print sample sizes
nrow(data.raw.reduc)
nrow(data.mstate.reduc)

###
### Extract the true risks into their own object
###
p.true <- data.raw.reduc %>% select(paste("p.true", 1:5, sep = ""))

###
### Create miss-calibrated risks (although they still need to sum to 1...)
### Need to get creative for this
###
p.est.perf <- p.true
p.est.over <- exp((log(p.true/(1-p.true)) + 0.5))/(1 + exp((log(p.true/(1-p.true)) + 0.5)))
p.est.under <- exp((log(p.true/(1-p.true)) - 0.5))/(1 + exp((log(p.true/(1-p.true)) - 0.5)))

### Create p.est object and assign colnames
p.est <- p.est.perf
colnames(p.est) <- paste("p.est", 1:5, sep = "")

### Note we want to use the order of individuals risk to group them up before calculating pseudo-values
### p.est.over and p.est.under do not change the order of patients, therefore we don't have to recalculate the pseudo-values each time
### Add the predicted risks to a dataset
p.est.data.frame <- cbind(data.raw.reduc, p.est)

### Define number of states
max.state <- max(data.mstate.reduc$to)

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
  p.est.data.frame <- arrange(p.est.data.frame, !!sym(paste("p.est", state, sep = "")))
  
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
    obs.aj.pctls[[state]][[i]] <- calc.calib.aj(data.mstate = subset(data.mstate.reduc, patid %in% data.sort.pctls[[state]][[i]]$patid), 
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
    patid.pv.same <- subset(data.mstate.reduc, patid %in% data.sort.pctls[[state]][[i]]$patid) %>%
      subset(from == 1 & to == 2 & Tstop > t.eval) %>%
      slice(1) %>%
      select(patid) %>%
      as.numeric()
    
    ### Calculate the pseudo value for this individual
    print(paste("PV SAME STATE = ", state, "pctl = ", i, Sys.time()))
    pv.same.pctls[[state]][[i]] <- func.calc.pv.aj(patid.eval = patid.pv.same, 
                                          data.mstate = subset(data.mstate.reduc, patid %in% data.sort.pctls[[state]][[i]]$patid), 
                                          obs.aj = obs.aj.pctls[[state]][[i]], 
                                          tmat = tmat, 
                                          n.cohort = nrow(data.sort.pctls[[state]][[i]]), 
                                          t.eval = t.eval)
  }
  
}


### Save image
rm(list=setdiff(ls(), list("data.mstate.reduc", "data.raw.reduc", "p.est.data.frame",
                           "data.sort.pctls", "pv.same.pctls", "obs.aj.pctls", "n.pctls", "max.state",
                           "p.true", "tmat", "t.eval", "scen", "n.cohort")))
save.image(paste("data/large_sample_analysis_pv_prep_pctls_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")


