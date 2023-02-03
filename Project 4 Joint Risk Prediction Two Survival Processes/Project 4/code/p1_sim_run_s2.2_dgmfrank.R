### Simulation file for scenario s2.1, frank copula DGM ###

### Clear workspace
rm(list=ls())

### Pull in arguments from qsub file
### arg1 = seed, arg2 = number of simulation runs in the loop, arg3 = development cohort size, arg4 = validation cohort size
args <- commandArgs(trailingOnly = T)
str(args)
as.numeric(args[1])
as.numeric(args[2])
as.numeric(args[3])
as.numeric(args[4])

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

## Set seed for simulation
set.seed.sim <- as.numeric(args[1])

## Set number of times we run the simultion
k.sim <- as.numeric(args[2])

## Define a vector which indicates which analysis methods should be used
## OPTIONS ARE: c("product", "joint", "msm", "c.clay", "c.clay.rot", "c.gumb", "c.gumb.rot", "c.joe", "c.joe.rot", "c.fgm", "c.frank", 
##                "f.normal_weib", "f.gamma_weib", "f.normal", "f.gamma")
anal.methods.sim <- c("msm", "product", "joint", "c.clay", "c.gumb", "c.frank", "f.normal_weib", "f.gamma_weib")

## Define what values of x1 we will evaluate strong calibration at
str.calib.eval.sim <- c(-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2)

## Size of development/validation dataset
n.devel.sim <- as.numeric(args[3])
n.valid.sim <- as.numeric(args[4])

## Maximum follow up before censoring
max.follow.sim <- 10

## Followup time at which to evaluate
t.eval.sim <- 10

## Number of iterations we run the stan model for, for fraily models
n.iter.sim <- 4000

## Shape/scale/covariate effects for censoring mechanism
baseline_cens.sim <- c(1, 95)
COV_beta_cens.sim <- c(0.1, 0.1)


##############################################################
### Define input specific to the data generating mechanism ###
##############################################################

### Choice of data generating mechanism
# DGM.sim <- "msm"
# DGM.sim <- "f.normal"
# DGM.sim <- "f.gamma"
# DGM.sim <- "c.clay"
# DGM.sim <- "c.gumb"
# DGM.sim <- "c.joe"
# DGM.sim <- "c.fgm"
DGM.sim <- "c.frank"
# DGM.sim <- "c.clay.rot"
# DGM.sim <- "c.gumb.rot"
# DGM.sim <- "c.joe.rot"

### The following paramters are specific to the data generating mechanism
## Shape and scale paramters for baseline hazard of each transition
if (DGM.sim == "msm"){
  
} else if(DGM.sim %in% c("c.clay", "c.gumb", "c.joe", "c.fgm", "c.frank", "c.clay.rot", "c.gumb.rot", "c.joe.rot")){
  ## Baseline hazard parameters for marginal distributions
  baselineA.sim <- c(1,35)
  baselineB.sim <- c(1,35)
  ## Covariate effects
  COV_betaA.sim <- c(0.3, 0.85)
  COV_betaB.sim <- c(0.9, 0.22)
  ## Copula association parameter
  copula.param.sim <- 4.2
  ## Name of copula and rotation
  if (DGM.sim == "c.clay"){
    copula.sim <- "clayton"
    rotate.cop.sim <- 0
  } else if (DGM.sim == "c.clay.rot"){
    copula.sim <- "clayton"
    rotate.cop.sim <- 90
  } else if (DGM.sim == "c.gumb"){
    copula.sim <- "gumbel"
    rotate.cop.sim <- 0
  } else if (DGM.sim == "c.clay.rot"){
    copula.sim <- "gumbel"
    rotate.cop.sim <- 90
  } else if (DGM.sim == "c.joe"){
    copula.sim <- "joe"
    rotate.cop.sim <- 0
  } else if (DGM.sim == "c.joe.rot"){
    copula.sim <- "joe"
    rotate.cop.sim <- 90
  } else if (DGM.sim == "c.fgm"){
    copula.sim <- "fgm"
    rotate.cop.sim <- 0
  } else if (DGM.sim == "c.frank"){
    copula.sim <- "frank"
    rotate.cop.sim <- 0
  }
} else if(DGM.sim %in% c("f.normal", "f.gamma")){

}


################################################
### Create output object for the simulations ###
################################################

## Harrels C
sim.out.C.Har <- data.frame(matrix(0, nrow = k.sim, ncol = length(anal.methods.sim)))
colnames(sim.out.C.Har) <- anal.methods.sim


## Uno's C
sim.out.C.Uno <- data.frame(matrix(0, nrow = k.sim, ncol = length(anal.methods.sim)))
colnames(sim.out.C.Uno) <- anal.methods.sim


## Strong calibration at specific values of X
# For each combination of x1 and x2 there will be a seperate data.frame
# Going to just evaluate with x2 = 0 and x2 = 1, which doesn't really make sense, but I don't think we are going to use this output anyway now
str.calib.dat <- data.frame(matrix(0, nrow = k.sim, ncol = length(anal.methods.sim)))
colnames(str.calib.dat) <- anal.methods.sim

# Create output object for when x2 = 0
sim.out.calib.str <- vector("list", length(str.calib.eval.sim))
for (i in 1:length(sim.out.calib.str)){
  sim.out.calib.str[[i]] <- str.calib.dat
}
names(sim.out.calib.str) <- paste("x1=",str.calib.eval.sim, sep = "")

# The output object for x2 = 1 has an identical structure

# Create a list with two identical elements of this format
sim.out.calib.str <- list(sim.out.calib.str, sim.out.calib.str)
names(sim.out.calib.str) <- c("x2=0", "x2=1")


## Calibration plots
# The data for each calibration plot will be stored in a data.frame, as list elements
# This is the validation dataset, with predicted values of risk for each individual/method
sim.out.calib.plot <- vector("list", k.sim)


####################
### Set the seed ###
####################
set.seed(set.seed.sim)
time.start <- Sys.time()


##########################
##########################
### Run the simulation ###
##########################
##########################

for (k in 1:k.sim){
  
  ##############################################################
  ##############################################################
  ### Step 1: Generate development and validation datasets #####
  ##############################################################
  ##############################################################
  
  ### First generate the predictor variables datasets
  x.baseline.devel <- data.frame("x1" = rnorm(n.devel.sim, 0, 1), "x2" = rnorm(n.devel.sim, 0, 1))
  x.baseline.valid <- data.frame("x1" = rnorm(n.valid.sim, 0, 1), "x2" = rnorm(n.valid.sim, 0, 1))
  
  ### Generate the full datasets
  ## Development
  if (DGM.sim == "msm"){
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
    
    ## Validation
    dat.valid <- gen.dat.msm(n = n.valid.sim,
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
                             x.in = x.baseline.valid,
                             numsteps = numsteps.sim)
  } else if (DGM.sim %in% c("c.clay", "c.gumb", "c.joe", "c.fgm", "c.frank", "c.clay.rot", "c.gumb.rot", "c.joe.rot")){
    ## Development
    dat.devel <- gen.dat.copula(n = n.devel.sim, 
                                max.follow = max.follow.sim,
                                copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
                                baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
                                baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim,
                                x.in = x.baseline.devel)
    
    ## Validation
    dat.valid <- gen.dat.copula(n = n.valid.sim, 
                                max.follow = max.follow.sim,
                                copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
                                baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
                                baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim,
                                x.in = x.baseline.valid)
    
  } else if(DGM.sim %in% c("f.normal", "f.gamma")){
    ## Development
    dat.devel <- gen.dat.frailty(n = n.devel.sim, 
                                 max.follow = max.follow.sim,
                                 frail.dist = frail.dist.sim, frail.eff = 1, frail.var = frail.var.sim,
                                 baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                 baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
                                 baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim,
                                 x.in = x.baseline.devel)
    
    ## Validation
    dat.valid <- gen.dat.frailty(n = n.valid.sim, 
                                 max.follow = max.follow.sim,
                                 frail.dist = frail.dist.sim, frail.eff = 1, frail.var = frail.var.sim,
                                 baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                 baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
                                 baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim,
                                 x.in = x.baseline.valid)
  }
  
  
  
  ### Also want to create formatted versions of the development and validation datasets, which will be used at model validation stage
  ### This will have one row per individual, with a joint outcome time until A+B
  
  ## NB: We only need to format the development dataset for calulation of Uno's C statistic
  dat.devel.format <- format.dataset.for.validation(dat.devel)
  dat.valid.format <- format.dataset.for.validation(dat.valid)
  
  
  ###########################################################
  ###########################################################
  ### Step 2: Fit each model to the development dataset #####
  ###########################################################
  ###########################################################
  
  ### Fit each model for methods that have been specified in the anal.methods.sim vector
  ### Also want add into a list for easier use later
  
  ## Create list
  fit.list <- vector("list", length(anal.methods.sim))
  names(fit.list) <- anal.methods.sim
  
  ## Now fit the models and add to list
  if ("msm" %in% anal.methods.sim){
    fit.msm <- fit.model.msm(dat.devel)
    fit.list[["msm"]] <- fit.msm
  }
  if ("product" %in% anal.methods.sim){
    fit.product <- fit.model.product(dat.devel)
    fit.list[["product"]] <- fit.product
  }
  if ("joint" %in% anal.methods.sim){
    fit.joint <- fit.model.joint(dat.devel)
    fit.list[["joint"]] <- fit.joint
  }
  if ("c.clay" %in% anal.methods.sim){
    fit.c.clay <- fit.model.copula(dat.devel, copula = "clayton", rotate.cop = 0)
    fit.list[["c.clay"]] <- fit.c.clay
  }
  if ("c.clay.rot" %in% anal.methods.sim){
    fit.c.clay.rot <- fit.model.copula(dat.devel, copula = "clayton", rotate.cop = 90)
    fit.list[["c.clay.rot"]] <- fit.c.clay.rot
  }
  if ("c.gumb" %in% anal.methods.sim){
    fit.c.gumb <- fit.model.copula(dat.devel, copula = "gumbel", rotate.cop = 0)
    fit.list[["c.gumb"]] <- fit.c.gumb
  }
  if ("c.gumb.rot" %in% anal.methods.sim){
    fit.c.gumb.rot <- fit.model.copula(dat.devel, copula = "gumbel", rotate.cop = 90)
    fit.list[["c.gumb.rot"]] <- fit.c.gumb.rot
  }
  if ("c.joe" %in% anal.methods.sim){
    fit.c.joe <- fit.model.copula(dat.devel, copula = "joe", rotate.cop = 0)
    fit.list[["c.joe"]] <- fit.c.joe
  }
  if ("c.joe.rot" %in% anal.methods.sim){
    fit.c.joe.rot <- fit.model.copula(dat.devel, copula = "joe", rotate.cop = 90)
    fit.list[["c.joe.rot"]] <- fit.c.joe.rot
  }
  if ("c.fgm" %in% anal.methods.sim){
    fit.c.fgm <- fit.model.copula(dat.devel, copula = "fgm", rotate.cop = 0)
    fit.list[["c.fgm"]] <- fit.c.fgm
  }
  if ("c.frank" %in% anal.methods.sim){
    fit.c.frank <- fit.model.copula(dat.devel, copula = "frank", rotate.cop = 0)
    fit.list[["c.frank"]] <- fit.c.frank
  }
  if ("f.normal" %in% anal.methods.sim){
    fit.f.normal <- fit.model.frailty(dat.devel, frail.dist = "normal", baseline.dist = "exp", n.iter = n.iter.sim)
    fit.list[["f.normal"]] <- fit.f.normal
  }
  if ("f.normal_weib" %in% anal.methods.sim){
    fit.f.normal_weib <- fit.model.frailty(dat.devel, frail.dist = "normal", baseline.dist = "weibull", n.iter = n.iter.sim)
    fit.list[["f.normal_weib"]] <- fit.f.normal_weib
  }
  if ("f.gamma" %in% anal.methods.sim){
    fit.f.gamma <- fit.model.frailty(dat.devel, frail.dist = "gamma", baseline.dist = "exp", n.iter = n.iter.sim)
    fit.list[["f.gamma"]] <- fit.f.gamma
  }
  if ("f.gamma_weib" %in% anal.methods.sim){
    fit.f.gamma_weib <- fit.model.frailty(dat.devel, frail.dist = "gamma", baseline.dist = "weibull", n.iter = n.iter.sim)
    fit.list[["f.gamma_weib"]] <- fit.f.gamma_weib
  }
  
  
  
  ################################################################################
  ################################################################################
  ### Step 3: Generate predicted risks and true risks for validation dataset #####
  ################################################################################
  ################################################################################
  
  ### Generate predicted risks and true risks, this will be used to calculate discrimination and calibration
  dat.valid.pred <- generate.predicted.and.true.risks(dat.valid.format)
  
  ### Also save this to output object, as it will be used for calibration plots
  sim.out.calib.plot[[k]] <- dat.valid.pred
  
  
  ########################################################
  ########################################################
  ### Step 4: Calculate discrimination of each model #####
  ########################################################
  ########################################################
  
  
  ##################
  ### Harrel's C ###
  ##################
  
  ### Specify an output vector
  C.Har <- vector("numeric", length = length(anal.methods.sim))
  
  ### Name the elements of the vector
  names(C.Har) <- anal.methods.sim
  
  ### Fit each model for methods that have been specified in the anal.methods.sim vector
  for (i in 1:length(anal.methods.sim)){
    C.Har[anal.methods.sim[i]] <- calc.C.Har(dat.valid.pred, paste("predrisk.", anal.methods.sim[i], sep = ""))
  }
  
  ### Store in output object
  sim.out.C.Har[k, ] <- C.Har
  
  
  ###############
  ### Uno's C ###
  ###############
  
  ### Specify an output vector
  C.Uno <- vector("numeric", length = length(anal.methods.sim))
  
  ### Name the elements of the vector
  names(C.Uno) <- anal.methods.sim
  
  ### Fit each model for methods that have been specified in the anal.methods.sim vector
  for (i in 1:length(anal.methods.sim)){
    C.Uno[anal.methods.sim[i]] <- calc.C.Uno(dat.devel.format, dat.valid.pred, paste("predrisk.", anal.methods.sim[i], sep = ""))
  }
  
  ### Store in output object
  sim.out.C.Uno[k, ] <- C.Uno
  
  
  ####################################################################################
  ####################################################################################
  ### Step 5: Calculate calibration at specific values of X (strong calibration) #####
  ####################################################################################
  ####################################################################################
  
  ### Want to compre predicted risk with true risk, for a range of values across predictor space x1 and x2
  ### Note that we do not require the validation dataset for this
  
  ### Define the list we will store the output in
  str.calib <- vector("list", 2)
  names(str.calib) <- c("x2=0", "x2=1")
  
  ### We want to consider a range of x1 values, and both x2 values
  x1.eval <- str.calib.eval.sim
  x2.eval <- c(0, 1)
  
  ### Create a list within each entry of str.calib
  str.calib[[1]] <- vector("list", length(x1.eval))
  str.calib[[2]] <- vector("list", length(x1.eval))
  
  ### Assign names
  names(str.calib[[1]]) <- paste("x1=",x1.eval, sep = "")
  names(str.calib[[2]]) <- paste("x1=",x1.eval, sep = "")
  
  ### Now do the calculation for each entity
  for (j in 1:length(x2.eval)){
    for (i in 1:length(x1.eval)){
      ### Calculate truerisk
      if (DGM.sim == "msm"){
        truerisk <- as.numeric(calc.true.risk.msm(t.eval = t.eval.sim,  
                                                  x1.eval = x1.eval[i], 
                                                  x2.eval = x2.eval[j],
                                                  shape12 = shape12.sim, scale12 = scale12.sim,
                                                  shape13 = shape13.sim, scale13 = scale13.sim,
                                                  shape24 = shape24.sim, scale24 = scale24.sim,
                                                  shape34 = shape34.sim, scale34 = scale34.sim,
                                                  beta12.cont = beta12.cont.sim, beta12.cont2 = beta12.cont2.sim,
                                                  beta13.cont = beta13.cont.sim, beta13.cont2 = beta13.cont2.sim,
                                                  beta24.cont = beta24.cont.sim, beta24.cont2 = beta24.cont2.sim,
                                                  beta34.cont = beta34.cont.sim, beta34.cont2 = beta34.cont2.sim)["risk.joint.true"])
        
      } else if(DGM.sim %in% c("c.clay", "c.gumb", "c.joe", "c.fgm", "c.frank", "c.clay.rot", "c.gumb.rot", "c.joe.rot")){
        truerisk <- as.numeric(calc.true.risk.copula(t.eval = t.eval.sim,  
                                                     x1.eval = x1.eval[i], 
                                                     x2.eval = x2.eval[j],
                                                     copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
                                                     baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                                     baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.joint.true"])
        
      } else if(DGM.sim %in% c("f.gamma","f.normal")){
        truerisk <- as.numeric(calc.true.risk.frailty(t.eval = t.eval.sim, 
                                                      x1.eval = x1.eval[i], 
                                                      x2.eval = x2.eval[j],
                                                      frail.dist = frail.dist.sim, frail.var = frail.var.sim,
                                                      baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                                      baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
                                                      int.method = int.method.sim)["risk.joint.true"])
      }
      
      ### Calculate predicted risk for each analysis method, calculate calibration, and put into vector
      
      ## Create object where results will be stored
      calib.output <- vector("numeric", length = length(anal.methods.sim))
      
      ## Assign names
      names(calib.output) <- anal.methods.sim
      
      ## Calculate calibration for each method 
      if ("msm" %in% anal.methods.sim){
        calib.output["msm"] <- as.numeric(calc.risk.msm(fit.in = fit.msm, 
                                                        t.eval = t.eval.sim, 
                                                        x1.eval = x1.eval[i], 
                                                        x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("product" %in% anal.methods.sim){
        calib.output["product"] <- as.numeric(calc.risk.product(fit.in = fit.product, 
                                                                t.eval = t.eval.sim, 
                                                                x1.eval = x1.eval[i], 
                                                                x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("joint" %in% anal.methods.sim){
        calib.output["joint"] <- as.numeric(calc.risk.joint(fit.in = fit.joint, 
                                                            t.eval = t.eval.sim, 
                                                            x1.eval = x1.eval[i], 
                                                            x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("c.clay" %in% anal.methods.sim){
        calib.output["c.clay"] <- as.numeric(calc.risk.copula(fit.in = fit.c.clay, 
                                                              t.eval = t.eval.sim, 
                                                              x1.eval = x1.eval[i], 
                                                              x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("c.clay.rot" %in% anal.methods.sim){
        calib.output["c.clay.rot"] <- as.numeric(calc.risk.copula(fit.in = fit.c.clay.rot, 
                                                                  t.eval = t.eval.sim,  
                                                                  x1.eval = x1.eval[i], 
                                                                  x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("c.gumb" %in% anal.methods.sim){
        calib.output["c.gumb"] <- as.numeric(calc.risk.copula(fit.in = fit.c.gumb, 
                                                              t.eval = t.eval.sim,  
                                                              x1.eval = x1.eval[i], 
                                                              x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("c.gumb.rot" %in% anal.methods.sim){
        calib.output["c.gumb.rot"] <- as.numeric(calc.risk.copula(fit.in = fit.c.gumb.rot, 
                                                                  t.eval = t.eval.sim,  
                                                                  x1.eval = x1.eval[i], 
                                                                  x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("c.joe" %in% anal.methods.sim){
        calib.output["c.joe"] <- as.numeric(calc.risk.copula(fit.in = fit.c.joe, 
                                                             t.eval = t.eval.sim,  
                                                             x1.eval = x1.eval[i], 
                                                             x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("c.joe.rot" %in% anal.methods.sim){
        calib.output["c.joe.rot"] <- as.numeric(calc.risk.copula(fit.in = fit.c.joe.rot, 
                                                                 t.eval = t.eval.sim,  
                                                                 x1.eval = x1.eval[i], 
                                                                 x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("c.fgm" %in% anal.methods.sim){
        calib.output["c.fgm"] <- as.numeric(calc.risk.copula(fit.in = fit.c.fgm, 
                                                             t.eval = t.eval.sim, 
                                                             x1.eval = x1.eval[i], 
                                                             x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("c.frank" %in% anal.methods.sim){
        calib.output["c.frank"] <- as.numeric(calc.risk.copula(fit.in = fit.c.frank, 
                                                             t.eval = t.eval.sim, 
                                                             x1.eval = x1.eval[i], 
                                                             x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("f.normal" %in% anal.methods.sim){
        calib.output["f.normal"] <- as.numeric(calc.risk.frailty(fit.in = fit.f.normal, 
                                                                 t.eval = t.eval.sim, 
                                                                 x1.eval = x1.eval[i], 
                                                                 x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("f.normal_weib" %in% anal.methods.sim){
        calib.output["f.normal_weib"] <- as.numeric(calc.risk.frailty(fit.in = fit.f.normal_weib, 
                                                                      t.eval = t.eval.sim, 
                                                                      x1.eval = x1.eval[i], 
                                                                      x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("f.gamma" %in% anal.methods.sim){
        calib.output["f.gamma"] <- as.numeric(calc.risk.frailty(fit.in = fit.f.gamma, 
                                                                t.eval = t.eval.sim, 
                                                                x1.eval = x1.eval[i], 
                                                                x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      if ("f.gamma_weib" %in% anal.methods.sim){
        calib.output["f.gamma_weib"] <- as.numeric(calc.risk.frailty(fit.in = fit.f.gamma_weib, 
                                                                     t.eval = t.eval.sim, 
                                                                     x1.eval = x1.eval[i], 
                                                                     x2.eval = x2.eval[j])["risk.joint.est"]) - truerisk
      }
      
      ### Assign calibration vector to main object
      sim.out.calib.str[[j]][[i]][k, ] <- calib.output
      
    }
  }
  
  print(k)
  Sys.time()
}

time.end <- Sys.time()

rm(x.baseline.devel, x.baseline.valid, dat.devel, dat.valid, dat.devel.format, dat.valid.format, dat.valid.pred, calib.output,
   fit.list, fit.product, fit.joint, fit.msm, fit.c.clay, fit.c.gumb, fit.c.frank, fit.f.normal_weib, fit.f.gamma_weib)
   
save.image(paste("data/s2.2/sim_run_s2.2_dgmfrank_n", n.devel.sim, "v", n.valid.sim, "s", set.seed.sim, ".RData", sep = ""))

print("ALL DONE IMAGE SAVED")