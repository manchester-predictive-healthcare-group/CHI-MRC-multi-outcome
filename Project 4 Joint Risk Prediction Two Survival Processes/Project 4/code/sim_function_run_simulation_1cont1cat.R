#####################################################################################################################
### This program will create functions to calculate joint risk by fitting each analysis method to a simulated dataset
### Seperate functions will be made to fit the model, and then calculate a predicted risk for a given set of covariates
### This means we only have to fit the initial model once, to reduce computation time (although for msm its the 
### generation of risk scores which is the lengthy part of the process)
#####################################################################################################################

### Clear workspace
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Load packages
# library(survival)
# library(gems)
# library(mstate)
# library(dplyr)
# library(tidyr)
# library(simsurv)
# library(frailtypack)
# library(CopulaCenR)
# library(copula)
# library(cubature)
# library(GJRM)
# library(rstan)


#############################################################################
### A few misc functions that will be used in the run_simulation function ###
#############################################################################


##################
### Function 1 ###
##################

### This function will create formatted versions of the development and validation datasets, which will be used at model validation stage
### This will have one row per individual, with a joint outcome time until A+B
format.dataset.for.validation <- function(data.in){
  
  ## Create wide version of dataset
  data.anal <- tidyr::pivot_wider(data.in, id_cols = id, names_from = outcome_char, values_from = c(time, status, x1, x2))
  
  ## Remove duplicate predictors variables and rename
  data.anal <- data.anal %>% 
    dplyr::select(id, time_A, time_B, status_A, status_B, x1_A, x2_A) %>%
    rename(x1 = x1_A, x2 = x2_A)
  
  ## Create a variable for until until both A and B, and the appropriate status indicator (both must be uncensored to observe event)
  data.anal <- data.anal %>%
    mutate(time = pmax(time_A, time_B),
           status = pmin(status_A, status_B))
  
  ## Remove variables we dont need
  data.anal <- data.anal %>% 
    dplyr::select(-time_A, -time_B, -status_A, -status_B)
  #head(data.anal)
  ## Turn into dataframe
  data.anal <- data.frame(data.anal)
  
  ## Return output
  return(data.anal)}



##################
### Function 2 ###
##################

### This function will calculate a predicted risk and true risk for each individual in the input dataset (this will be validation dataset)
### The inputted dataset must be formatted using the format.dataset.for.validation function (function 1)

generate.predicted.and.true.risks <- function(data.in){
  
  ### Create new output object
  data.out <- data.in
  
  ### Create empty vectors in which to store the predicted risks, 
  ### and run a loop to calculate the predicted risks for each individual for each analysis method
  if ("msm" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.msm <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.msm"] <- as.numeric(calc.risk.msm(fit.in = fit.msm, 
                                                              t.eval = t.eval.sim, 
                                                              x1.eval = data.out[i, "x1"], 
                                                              x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("product" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.product <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.product"] <- as.numeric(calc.risk.product(fit.in = fit.product, 
                                                                      t.eval = t.eval.sim, 
                                                                      x1.eval = data.out[i, "x1"], 
                                                                      x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("joint" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.joint <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.joint"] <- as.numeric(calc.risk.joint(fit.in = fit.joint, 
                                                                  t.eval = t.eval.sim, 
                                                                  x1.eval = data.out[i, "x1"], 
                                                                  x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("c.clay" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.c.clay <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.c.clay"] <- as.numeric(calc.risk.copula(fit.in = fit.c.clay, 
                                                                    t.eval = t.eval.sim, 
                                                                    x1.eval = data.out[i, "x1"], 
                                                                    x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("c.clay.rot" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.c.clay.rot <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.c.clay.rot"] <- as.numeric(calc.risk.copula(fit.in = fit.c.clay.rot, 
                                                                        t.eval = t.eval.sim, 
                                                                        x1.eval = data.out[i, "x1"], 
                                                                        x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("c.gumb" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.c.gumb <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.c.gumb"] <- as.numeric(calc.risk.copula(fit.in = fit.c.gumb, 
                                                                    t.eval = t.eval.sim, 
                                                                    x1.eval = data.out[i, "x1"], 
                                                                    x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("c.gumb.rot" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.c.gumb.rot <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.c.gumb.rot"] <- as.numeric(calc.risk.copula(fit.in = fit.c.gumb.rot, 
                                                                        t.eval = t.eval.sim, 
                                                                        x1.eval = data.out[i, "x1"], 
                                                                        x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("c.joe" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.c.joe <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.c.joe"] <- as.numeric(calc.risk.copula(fit.in = fit.c.joe, 
                                                                   t.eval = t.eval.sim, 
                                                                   x1.eval = data.out[i, "x1"], 
                                                                   x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("c.joe.rot" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.c.joe.rot <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.c.joe.rot"] <- as.numeric(calc.risk.copula(fit.in = fit.c.joe.rot, 
                                                                       t.eval = t.eval.sim, 
                                                                       x1.eval = data.out[i, "x1"], 
                                                                       x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("c.fgm" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.c.fgm <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.c.fgm"] <- as.numeric(calc.risk.copula(fit.in = fit.c.fgm, 
                                                                   t.eval = t.eval.sim, 
                                                                   x1.eval = data.out[i, "x1"], 
                                                                   x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("f.normal" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.f.normal <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.f.normal"] <- as.numeric(calc.risk.frailty(fit.in = fit.f.normal, 
                                                                       t.eval = t.eval.sim, 
                                                                       x1.eval = data.out[i, "x1"], 
                                                                       x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("f.normal_weib" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.f.normal_weib <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.f.normal_weib"] <- as.numeric(calc.risk.frailty(fit.in = fit.f.normal_weib, 
                                                                            t.eval = t.eval.sim, 
                                                                            x1.eval = data.out[i, "x1"], 
                                                                            x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("f.gamma" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.f.gamma <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.f.gamma"] <- as.numeric(calc.risk.frailty(fit.in = fit.f.gamma, 
                                                                      t.eval = t.eval.sim, 
                                                                      x1.eval = data.out[i, "x1"], 
                                                                      x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  if ("f.gamma_weib" %in% anal.methods.sim){
    # Create empty vector
    data.out$predrisk.f.gamma_weib <- rep(0, nrow(data.out))
    # Run the loop
    for (i in 1:nrow(data.out)){
      data.out[i, "predrisk.f.gamma_weib"] <- as.numeric(calc.risk.frailty(fit.in = fit.f.gamma_weib, 
                                                                           t.eval = t.eval.sim, 
                                                                           x1.eval = data.out[i, "x1"], 
                                                                           x2.eval = as.numeric(data.out[i, "x2"] == 1))["risk.joint.est"])
    }
  }
  
  print("predicted risks generated")
  print(Sys.time())
  
  ### Also want to do a loop to calculate all the true risks, which is different depending on the data generating mechanism
  
  ### First create variableto store these in
  data.out$truerisk <- rep(0, nrow(data.out))
  
  ### Calculate the true risk for each individual
  if (DGM.sim == "msm"){
    for (i in 1:nrow(data.out)){
      data.out[i, "truerisk"] <- as.numeric(calc.true.risk.msm(t.eval = t.eval.sim, 
                                                               x1.eval = data.out[i, "x1"], 
                                                               x2.eval = as.numeric(data.out[i, "x2"] == 1),
                                                               shape12 = shape12.sim, scale12 = scale12.sim,
                                                               shape13 = shape13.sim, scale13 = scale13.sim,
                                                               shape24 = shape24.sim, scale24 = scale24.sim,
                                                               shape34 = shape34.sim, scale34 = scale34.sim,
                                                               beta12.cont = beta12.cont.sim, beta12.cat = beta12.cat.sim,
                                                               beta13.cont = beta13.cont.sim, beta13.cat = beta13.cat.sim,
                                                               beta24.cont = beta24.cont.sim, beta24.cat = beta24.cat.sim,
                                                               beta34.cont = beta34.cont.sim, beta34.cat = beta34.cat.sim)["risk.joint.true"])
    }
  } else if(DGM.sim %in% c("c.clay", "c.gumb", "c.joe", "c.fgm", "c.clay.rot", "c.gumb.rot", "c.joe.rot")){
    for (i in 1:nrow(data.out)){
      data.out[i, "truerisk"] <- as.numeric(calc.true.risk.copula(t.eval = t.eval.sim, 
                                                                  x1.eval = data.out[i, "x1"], 
                                                                  x2.eval = as.numeric(data.out[i, "x2"] == 1),
                                                                  copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
                                                                  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                                                  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.joint.true"])
    }
  } else if(DGM.sim %in% c("f.gamma","f.normal")){
    for (i in 1:nrow(data.out)){
      data.out[i, "truerisk"] <- as.numeric(calc.true.risk.frailty(t.eval = t.eval.sim, 
                                                                   x1.eval = data.out[i, "x1"], 
                                                                   x2.eval = as.numeric(data.out[i, "x2"] == 1),
                                                                   frail.dist = frail.dist.sim, frail.var = frail.var.sim,
                                                                   baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
                                                                   baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
                                                                   int.method = int.method.sim)["risk.joint.true"])
    }
  }
  
  print("true risks generated")
  print(Sys.time())
  
  return(data.out)
  
}


##################
### Function 3 ###
##################

### Function to calculate Harrel's C for input dataset
calc.C.Har <- function(data.in, predrisk){
  return(as.numeric(rcorr.cens(1-data.in[ , predrisk], Surv(data.in$time, data.in$status))["C Index"]))}


##################
### Function 4 ###
##################

### Specify a function to calculate Uno's C for given input development and validation datasets, and predicted risks in validation dataset
calc.C.Uno <- function(data.devel.in, data.valid.in, predrisk){
  return(UnoC(Surv.rsp = Surv(data.devel.in$time, data.devel.in$status), 
              Surv.rsp.new = Surv(data.valid.in$time, data.valid.in$status), 
              lpnew = data.valid.in[ , predrisk]))
}



######################
######################
### Run simulation ###
######################
######################

run_simulation <- function(){
  
  for (k in 1:k.sim){
    
    ##############################################################
    ##############################################################
    ### Step 1: Generate development and validation datasets #####
    ##############################################################
    ##############################################################
    
    ### First generate the predictor variables datasets
    x.baseline.devel <- data.frame("x1" = rnorm(n.devel.sim, 0, 1), "x2" = rbinom(n.devel.sim, 1, 0.5))
    x.baseline.valid <- data.frame("x1" = rnorm(n.valid.sim, 0, 1), "x2" = rbinom(n.valid.sim, 1, 0.5))
    
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
                               beta12.cont = beta12.cont.sim, beta12.cat = beta12.cat.sim, #covariate effects for transiion 12
                               beta13.cont = beta13.cont.sim, beta13.cat = beta13.cat.sim, #covariate effects for transiion 13
                               beta24.cont = beta24.cont.sim, beta24.cat = beta24.cat.sim, #covariate effects for transiion 24
                               beta34.cont = beta34.cont.sim, beta34.cat = beta34.cat.sim, #covariate effects for transiion 34
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
                               beta12.cont = beta12.cont.sim, beta12.cat = beta12.cat.sim, #covariate effects for transiion 12
                               beta13.cont = beta13.cont.sim, beta13.cat = beta13.cat.sim, #covariate effects for transiion 13
                               beta24.cont = beta24.cont.sim, beta24.cat = beta24.cat.sim, #covariate effects for transiion 24
                               beta34.cont = beta34.cont.sim, beta34.cat = beta34.cat.sim, #covariate effects for transiion 34
                               baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim, #shape/scale and covariate effects for censoring dist
                               x.in = x.baseline.valid,
                               numsteps = numsteps.sim)
    } else if (DGM.sim %in% c("c.clay", "c.gumb", "c.joe", "c.fgm", "c.clay.rot", "c.gumb.rot", "c.joe.rot")){
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
    if ("f.normal" %in% anal.methods.sim){
      fit.f.normal <- fit.model.frailty(dat.devel, frail.dist = "normal", baseline.dist = "exp", n.iter = 5000)
      fit.list[["f.normal"]] <- fit.f.normal
    }
    if ("f.normal_weib" %in% anal.methods.sim){
      fit.f.normal_weib <- fit.model.frailty(dat.devel, frail.dist = "normal", baseline.dist = "weibull", n.iter = 5000)
      fit.list[["f.normal_weib"]] <- fit.f.normal_weib
    }
    if ("f.gamma" %in% anal.methods.sim){
      fit.f.gamma <- fit.model.frailty(dat.devel, frail.dist = "gamma", baseline.dist = "exp", n.iter = 5000)
      fit.list[["f.gamma"]] <- fit.f.gamma
    }
    if ("f.gamma_weib" %in% anal.methods.sim){
      fit.f.gamma_weib <- fit.model.frailty(dat.devel, frail.dist = "gamma", baseline.dist = "weibull", n.iter = 5000)
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
    
    # if ("msm" %in% anal.methods.sim){
    #   C.Har["msm"] <- calc.C.Har(dat.valid.pred, "predrisk.msm")
    # }
    # if ("product" %in% anal.methods.sim){
    #   C.Har["product"] <- calc.C.Har(dat.valid.pred, "predrisk.product")
    # }
    # if ("joint" %in% anal.methods.sim){
    #   C.Har["joint"] <- calc.C.Har(dat.valid.pred, "predrisk.joint")
    # }
    # if ("c.clay" %in% anal.methods.sim){
    #   C.Har["c.clay"] <- calc.C.Har(dat.valid.pred, "predrisk.c.clay")
    # }
    # if ("c.clay.rot" %in% anal.methods.sim){
    #   C.Har["c.clay.rot"] <- calc.C.Har(dat.valid.pred, "predrisk.c.clay.rot")
    # }
    # if ("c.gumb" %in% anal.methods.sim){
    #   C.Har["c.gumb"] <- calc.C.Har(dat.valid.pred, "predrisk.c.gumb")
    # }
    # if ("c.gumb.rot" %in% anal.methods.sim){
    #   C.Har["c.gumb.rot"] <- calc.C.Har(dat.valid.pred, "predrisk.c.gumb.rot")
    # }
    # if ("c.joe" %in% anal.methods.sim){
    #   C.Har["c.joe"] <- calc.C.Har(dat.valid.pred, "predrisk.c.joe")
    # }
    # if ("c.joe.rot" %in% anal.methods.sim){
    #   C.Har["c.joe.rot"] <- calc.C.Har(dat.valid.pred, "predrisk.c.joe.rot")
    # }
    # if ("c.fgm" %in% anal.methods.sim){
    #   C.Har["c.fgm"] <- calc.C.Har(dat.valid.pred, "predrisk.c.fgm")
    # }
    
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
    
    # if ("msm" %in% anal.methods.sim){
    #   C.Uno["msm"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.msm")
    # }
    # if ("product" %in% anal.methods.sim){
    #   C.Uno["product"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.product")
    # }
    # if ("joint" %in% anal.methods.sim){
    #   C.Uno["joint"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.joint")
    # }
    # if ("c.clay" %in% anal.methods.sim){
    #   C.Uno["c.clay"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.c.clay")
    # }
    # if ("c.clay.rot" %in% anal.methods.sim){
    #   C.Uno["c.clay.rot"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.c.clay.rot")
    # }
    # if ("c.gumb" %in% anal.methods.sim){
    #   C.Uno["c.gumb"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.c.gumb")
    # }
    # if ("c.gumb.rot" %in% anal.methods.sim){
    #   C.Uno["c.gumb.rot"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.c.gumb.rot")
    # }
    # if ("c.joe" %in% anal.methods.sim){
    #   C.Uno["c.joe"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.c.joe")
    # }
    # if ("c.joe.rot" %in% anal.methods.sim){
    #   C.Uno["c.joe.rot"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.c.joe.rot")
    # }
    # if ("c.fgm" %in% anal.methods.sim){
    #   C.Uno["c.fgm"] <- calc.C.Uno(dat.devel.format, dat.valid.pred, "predrisk.c.fgm")
    # }
    
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
                                                    beta12.cont = beta12.cont.sim, beta12.cat = beta12.cat.sim,
                                                    beta13.cont = beta13.cont.sim, beta13.cat = beta13.cat.sim,
                                                    beta24.cont = beta24.cont.sim, beta24.cat = beta24.cat.sim,
                                                    beta34.cont = beta34.cont.sim, beta34.cat = beta34.cat.sim)["risk.joint.true"])
          
        } else if(DGM.sim %in% c("c.clay", "c.gumb", "c.joe", "c.fgm", "c.clay.rot", "c.gumb.rot", "c.joe.rot")){
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
  
}

save.image("data/sim_function_run_simulation_1cont1cat.RData")
