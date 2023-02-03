###
### Report event rates for scenario s1_nocorr under each DGM and sample size ###
###

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Source packages
source("code/sim_function_load_packages.R")

### Load functions for the simulation
source("code/sim_function_generate_data_all_2cont.R")
source("code/sim_function_fit_models_all_2cont.R")
source("code/sim_function_calc_true_risk_all_2cont.R")
source("code/sim_function_run_simulation_2cont.R")

### Set global parameters, which will be used whichever DGM is chosen

## Set number of times we run the simultion
k.sim <- 1000

## Size of development dataset
n.devel.sim <- as.numeric(1000)

## Maximum follow up before censoring
max.follow.sim <- 10

## Followup time at which to evaluate
t.eval.sim <- 10

## Shape/scale/covariate effects for censoring mechanism
baseline_cens.sim <- c(1, 95)
COV_beta_cens.sim <- c(0.1, 0.1)


#####################################################################
### Create list to store event rates and other output of interest ###
#####################################################################
output.object <- data.frame(matrix(NA, nrow = k.sim, ncol = 12))
colnames(output.object) <- c("follow.up.A", "follow.up.B", "follow.up.AB", "n.events.A", "n.events.B", "n.events.AB", 
                             "rate.A", "rate.B", "rate.AB", "risk.A", "risk.B", "risk.AB")


##########################
##########################
### Run the simulation ###
##########################
##########################
time.start <- Sys.time()

### Assign DGM.in
DGM.in <- "c.clay"

### The following paramters are specific to the data generating mechanism
## Shape and scale paramters for baseline hazard of each transition
if (DGM.in == "msm"){
  
} else if(DGM.in %in% c("c.clay", "c.gumb", "c.joe", "c.fgm", "c.frank", "c.clay.rot", "c.gumb.rot", "c.joe.rot")){
  ## Baseline hazard parameters for marginal distributions
  baselineA.sim <- c(1,115)
  baselineB.sim <- c(1,120)
  ## Covariate effects
  COV_betaA.sim <- c(0.25, 0.6)
  COV_betaB.sim <- c(0.73, 0.15)
  ## Name of copula and rotation
  if (DGM.in == "c.clay"){
    copula.sim <- "clayton"
    rotate.cop.sim <- 0
    ## Copula association parameter
    copula.param.sim <- 0.000001
  } else if (DGM.in == "c.clay.rot"){
    
  } else if (DGM.in == "c.gumb"){
    
  } else if (DGM.in == "c.clay.rot"){
    
  } else if (DGM.in == "c.joe"){
    
  } else if (DGM.in == "c.joe.rot"){
    
  } else if (DGM.in == "c.fgm"){
    
  } else if (DGM.in == "c.frank"){
    
  }
} else if(DGM.in %in% c("f.normal", "f.gamma")){
  if (DGM.in == "f.normal"){
    
  } else if (DGM.in == "f.gamma"){
    
  }
}

for (k in 1:k.sim){
  
  ###
  ### Set seed and print progress
  ###
  set.seed(k)
  print(paste("DGM.in = ", DGM.in, "n.devel.sim = ", n.devel.sim, "iter = ", k, Sys.time()))
  
  ##############################################################
  ##############################################################
  ### Step 1: Generate development and validation datasets #####
  ##############################################################
  ##############################################################
  
  ### Note that we still generate both, so the exact same development data is generated (based off the seed) as was generated 
  ### in the simulation used in the paper.
  
  ### First generate the predictor variables datasets
  x.baseline.devel <- data.frame("x1" = rnorm(n.devel.sim, 0, 1), "x2" = rnorm(n.devel.sim, 0, 1))
  
  ### Generate the full datasets
  ## Development
  if (DGM.in == "msm"){
    ## Development
    dat.devel <- gen.dat.msm(n = n.devel.sim,
                             max.follow = max.follow.sim,
                             shape12 = shape12.sim, scale12 = scale12.sim, #shape and scale for weibull baseline hazard for transition 1 -> 2
                             shape13 = shape13.sim, scale13 = scale13.sim, #shape and scale for weibull baseline hazard for transition 1 -> 3
                             shape24 = shape24.sim, scale24 = scale24.sim, #shape and scale for weibull baseline hazard for transition 2 -> 4
                             shape34 = shape34.sim, scale34 = scale34.sim, #shape and scale for weibull baseline hazard for transition 3 -> 4
                             beta12.cont = beta12.cont.sim, beta12.cont2 = beta12.cont2.sim, #covariate effects for transiion 12
                             beta13.cont = beta13.cont.sim, beta13.cont2 = beta13.cont2.sim,  #covariate effects for transiion 13
                             beta24.cont = beta24.cont.sim, beta24.cont2 = beta24.cont2.sim,  #covariate effects for transiion 24
                             beta34.cont = beta34.cont.sim, beta34.cont2 = beta34.cont2.sim,  #covariate effects for transiion 34
                             baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim, #shape/scale and covariate effects for censoring dist
                             x.in = x.baseline.devel,
                             numsteps = numsteps.sim)
    
  } else if (DGM.in %in% c("c.clay", "c.gumb", "c.joe", "c.fgm", "c.frank", "c.clay.rot", "c.gumb.rot", "c.joe.rot")){
    ## Development
    dat.devel <- gen.dat.copula(n = n.devel.sim, 
                                max.follow = max.follow.sim,
                                copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
                                baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
                                baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim,
                                x.in = x.baseline.devel)
    
  } else if(DGM.in %in% c("f.normal", "f.gamma")){
    ## Development
    dat.devel <- gen.dat.frailty(n = n.devel.sim, 
                                 max.follow = max.follow.sim,
                                 frail.dist = frail.dist.sim, frail.eff = 1, frail.var = frail.var.sim,
                                 baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                 baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
                                 baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim,
                                 x.in = x.baseline.devel)
    
  }
  
  ### Also want to create formatted versions of the development and validation datasets, which will be used at model validation stage
  ### This will have one row per individual, with a joint outcome time until A+B
  dat.devel.format <- format.dataset.for.validation(dat.devel)
  
  ### Calculate #events and event rates for A, B and AB seperately
  dat.devel.A <-subset(dat.devel, outcome_char == "A")
  survfit.A <- survfit(Surv(time, status) ~ 1, data= dat.devel.A)
  
  dat.devel.B <-subset(dat.devel, outcome_char == "B")
  survfit.B <- survfit(Surv(time, status) ~ 1, data= dat.devel.B)
  
  survfit.AB <- survfit(Surv(time, status) ~ 1, data= dat.devel.format)
  
  ### Extract into output dataset
  output.object[k, "follow.up.A"] <- sum(dat.devel.A$time)
  output.object[k, "follow.up.B"] <- sum(dat.devel.B$time)
  output.object[k, "follow.up.AB"] <- sum(dat.devel.format$time)
  
  output.object[k, "n.events.A"] <- sum(dat.devel.A$status)
  output.object[k, "n.events.B"] <- sum(dat.devel.B$status)
  output.object[k, "n.events.AB"] <- sum(dat.devel.format$status)
  
  output.object[k, "rate.A"] <- sum(dat.devel.A$status)/sum(dat.devel.A$time)
  output.object[k, "rate.B"] <- sum(dat.devel.B$status)/sum(dat.devel.B$time)
  output.object[k, "rate.AB"] <- sum(dat.devel.format$status)/sum(dat.devel.format$time)
  
  output.object[k, "risk.A"] <- min(summary(survfit.A)$surv)
  output.object[k, "risk.B"] <- min(summary(survfit.B)$surv)
  output.object[k, "risk.AB"] <- min(summary(survfit.AB)$surv)
  
}
time.end <- Sys.time()

### Turn output objects into means and standard errors
rm(list = setdiff(ls(), list("output.object", "n.devel.sim")))

save.image(paste("data/event_rates_s1_nocorr.RData", sep = ""))
