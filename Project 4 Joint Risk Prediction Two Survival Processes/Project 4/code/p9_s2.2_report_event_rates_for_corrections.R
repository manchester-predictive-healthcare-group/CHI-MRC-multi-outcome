###
### Report event rates for scenario s2.2 under each DGM and sample size ###
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
output.object <- vector("list", 6)
names(output.object) <- c("msm", "c.clay", "c.gumb", "c.frank", "f.normal", "f.gamma")
for (DGM.num in 1:6){
  output.object[[DGM.num]]  <- data.frame(matrix(NA, nrow = k.sim, ncol = 12))
  colnames(output.object[[DGM.num]]) <- c("follow.up.A", "follow.up.B", "follow.up.AB", "n.events.A", "n.events.B", "n.events.AB", 
                                          "rate.A", "rate.B", "rate.AB", "risk.A", "risk.B", "risk.AB")
}

##########################
##########################
### Run the simulation ###
##########################
##########################
time.start <- Sys.time()

### Assign counters for assinging output
DGM.num <- 1

for (DGM.in in c("msm", "c.clay", "c.gumb", "c.frank", "f.normal", "f.gamma")){
  
  ### The following paramters are specific to the data generating mechanism
  ## Shape and scale paramters for baseline hazard of each transition
  if (DGM.in == "msm"){
    ## Rates
    shape12.sim <- 1
    scale12.sim <- 42
    shape13.sim <- 1
    scale13.sim <- 43
    shape24.sim <- 1
    scale24.sim <- scale13.sim*0.33
    shape34.sim <- 1 
    scale34.sim <- scale12.sim*0.33
    
    ## Covariate effects for each transition
    beta12.cont.sim <- 0.22
    beta12.cont2.sim <- 0.7
    beta13.cont.sim <- 0.9
    beta13.cont2.sim <- 0.15
    beta24.cont.sim <- 0.9
    beta24.cont2.sim <- 0.15
    beta34.cont.sim <- 0.22
    beta34.cont2.sim <- 0.7
    
    ## Number of sampling steps when simulating data using multistate model
    numsteps.sim <- 50000
    
  } else if(DGM.in %in% c("c.clay", "c.gumb", "c.joe", "c.fgm", "c.frank", "c.clay.rot", "c.gumb.rot", "c.joe.rot")){
    ## Baseline hazard parameters for marginal distributions
    baselineA.sim <- c(1,35)
    baselineB.sim <- c(1,35)
    ## Covariate effects
    COV_betaA.sim <- c(0.3, 0.85)
    COV_betaB.sim <- c(0.9, 0.22)
    ## Name of copula and rotation
    if (DGM.in == "c.clay"){
      copula.sim <- "clayton"
      rotate.cop.sim <- 0
      ## Copula association parameter
      copula.param.sim <- 2
    } else if (DGM.in == "c.clay.rot"){
      copula.sim <- "clayton"
      rotate.cop.sim <- 90
    } else if (DGM.in == "c.gumb"){
      copula.sim <- "gumbel"
      rotate.cop.sim <- 0
      ## Copula association parameter
      copula.param.sim <- 1.54
    } else if (DGM.in == "c.clay.rot"){
      copula.sim <- "gumbel"
      rotate.cop.sim <- 90
    } else if (DGM.in == "c.joe"){
      copula.sim <- "joe"
      rotate.cop.sim <- 0
    } else if (DGM.in == "c.joe.rot"){
      copula.sim <- "joe"
      rotate.cop.sim <- 90
    } else if (DGM.in == "c.fgm"){
      copula.sim <- "fgm"
      rotate.cop.sim <- 0
    } else if (DGM.in == "c.frank"){
      copula.sim <- "frank"
      rotate.cop.sim <- 0
      ## Copula association parameter
      copula.param.sim <- 4.2
    }
  } else if(DGM.in %in% c("f.normal", "f.gamma")){
    if (DGM.in == "f.normal"){
      ## Baseline parameters
      baselineA.sim <- c(1,50.5)
      baselineB.sim <- c(1,51.5)
      ## Covariate effects
      COV_betaA.sim <- c(0.4, 1)
      COV_betaB.sim <- c(1.1, 0.35)
      ## Frailty association parameter
      frail.var.sim <- 1.38
      ## Assign distribution
      frail.dist.sim <- "normal"
      int.method.sim <- "pcubature"
    } else if (DGM.in == "f.gamma"){
      ## Baseline hazard parameters for marginal distributions
      baselineA.sim <- c(1,25)
      baselineB.sim <- c(1,25)
      ## Covariate effects
      COV_betaA.sim <- c(0.4, 1)
      COV_betaB.sim <- c(1.1, 0.35)
      ## Frailty association parameter
      frail.var.sim <- 1.45
      ## Assigm distribution
      frail.dist.sim <- "gamma"
      int.method.sim <- "hcubature"
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
      output.object[[DGM.num]][k, "follow.up.A"] <- sum(dat.devel.A$time)
      output.object[[DGM.num]][k, "follow.up.B"] <- sum(dat.devel.B$time)
      output.object[[DGM.num]][k, "follow.up.AB"] <- sum(dat.devel.format$time)
      
      output.object[[DGM.num]][k, "n.events.A"] <- sum(dat.devel.A$status)
      output.object[[DGM.num]][k, "n.events.B"] <- sum(dat.devel.B$status)
      output.object[[DGM.num]][k, "n.events.AB"] <- sum(dat.devel.format$status)
      
      output.object[[DGM.num]][k, "rate.A"] <- sum(dat.devel.A$status)/sum(dat.devel.A$time)
      output.object[[DGM.num]][k, "rate.B"] <- sum(dat.devel.B$status)/sum(dat.devel.B$time)
      output.object[[DGM.num]][k, "rate.AB"] <- sum(dat.devel.format$status)/sum(dat.devel.format$time)
      
      output.object[[DGM.num]][k, "risk.A"] <- min(summary(survfit.A)$surv)
      output.object[[DGM.num]][k, "risk.B"] <- min(summary(survfit.B)$surv)
      output.object[[DGM.num]][k, "risk.AB"] <- min(summary(survfit.AB)$surv)
      
    }
  ## Add to DGM counter
  DGM.num <- DGM.num + 1
}
time.end <- Sys.time()

### Turn output objects into means and standard errors
rm(list = setdiff(ls(), list("output.object", "n.devel.sim")))

save.image(paste("data/event_rates_s2.2.RData", sep = ""))
