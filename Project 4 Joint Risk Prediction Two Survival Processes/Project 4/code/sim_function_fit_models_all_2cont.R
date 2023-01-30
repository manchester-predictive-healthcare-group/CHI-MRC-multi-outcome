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

##################
##################
### FIT MSM ###
##################
##################

### Function to fit the multistate model
fit.model.msm <- function(data.in){
  
  ## Get data in format for analysis with MSM model
  data.anal.pre <- tidyr::pivot_wider(data.in, id_cols = id, names_from = outcome_char, values_from = c(time, status, x1, x2))
  
  ## Remove duplicate predictors variables and rename
  data.anal.pre <- data.anal.pre %>% 
    dplyr::select(id, time_A, time_B, status_A, status_B, x1_A, x2_A) %>%
    rename(x1 = x1_A, x2 = x2_A)
  
  ## Assign event times
  data.anal.pre$state1 <- rep(0, nrow(data.anal.pre))
  data.anal.pre$state2 <- data.anal.pre$time_A
  data.anal.pre$state3 <- data.anal.pre$time_B
  data.anal.pre$state4 <- pmax(data.anal.pre$time_A, data.anal.pre$time_B)
  
  
  ## Assign censoring indicator  
  ## First assign all 0 (censored)
  data.anal.pre$state2.s <- rep(0, nrow(data.anal.pre))
  data.anal.pre$state3.s <- rep(0, nrow(data.anal.pre))
  data.anal.pre$state4.s <- rep(0, nrow(data.anal.pre))
  
  ## Entry to state2 will be observed if time_A < time_B and status_A = 1 
  ## Entry to state3 will be observed if time_B < time_A and status_B = 1 
  ## Entry to state4 will be observed if status_A and status_B = 1
  data.anal.pre <- data.anal.pre %>% 
    mutate(state2.s = case_when(time_A < time_B & status_A == 1 ~ 1, TRUE ~ as.numeric(as.character(state2.s))),
           state3.s = case_when(time_B < time_A & status_B == 1 ~ 1, TRUE ~ as.numeric(as.character(state3.s))),
           state4.s = case_when(status_A == 1 & status_B == 1 ~ 1, TRUE ~ as.numeric(as.character(state4.s)))
    )
  
  ## Select only variables of interest
  data.anal.pre <- dplyr::select(data.anal.pre, id, state1, state2, state3, state4, state2.s, state3.s, state4.s, x1, x2)
  
  ## Note that we only observe the initial transition into A/B if it happens first 
  ## (i.e. the event time for entry into state 2 and state 3 should be the same, 
  ## with censoring on whichever event doesnt happen first)
  ## This will be dealt with using the msprep functionality to deal with this
  
  ## First create a transition matrix for the allowed transitions
  tmat <- transMat(x = list(c(2,3), c(4), c(4), c()), 
                   names = c("state1", "state2", "state3", "state4"))
  
  ## Change from tibble into data.frame
  data.anal.pre <- data.frame(data.anal.pre)
  #str(data.anal.pre)
  
  ## Now can prepr the data into format which can be analysed uing mstate
  data.anal <- msprep(data.anal.pre, trans = tmat, time = c(NA, "state2", "state3", "state4"),
                      status = c(NA, "state2.s", "state3.s", "state4.s"), keep = c("x1","x2"))
  
  ## Note how the censoring transitions into state 2/3 are now done so at the correct time
  #head(data.anal, n = 20)
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("x1", "x2")
  data.anal <- expand.covs(data.anal, covs, longnames = FALSE)
  head(data.anal)
  
  
  ### Now to fit the cause-specific hazard models MSM and store for output
  msm.model <- coxph(Surv(Tstart, Tstop, status) ~ x1.1 + x1.2 + x1.3 + x1.4 + 
                        x2.1 + x2.2 + x2.3 + x2.4 + strata(trans), 
                      data = data.anal, method = "breslow")
  
  ### Also save the development data, the format of which will be used to create "newdata" arguments for when we generate risk scores
  msm.data <- data.anal.pre
  
  ### Also save transition matrix
  msm.tmat <- tmat
  
  ### Create output object
  output.obj <- vector("list",3)
  names(output.obj) <- c("msm.model", "msm.data", "msm.tmat")
  
  output.obj[[1]] <- msm.model
  output.obj[[2]] <- msm.data
  output.obj[[3]] <- msm.tmat
  
  return(output.obj)
}


####################
####################
### PREDRISK MSM ###
####################
####################

### Function to generate predicted risks, using a fitted multistate model

### fit.in = the output object from fit.model.msm function
### t.eval is time to predict risks for
### x1.eval = x1 value to assess risks at
### x2.eval = x2 value to assess risks at
calc.risk.msm <- function(fit.in, t.eval, x1.eval, x2.eval){

### Create a patient for which we run probtrans
## Start by creating a base using the formatted data used to fit the msm
pat <- fit.in[["msm.data"]][1, ]

# Assign x1 and x2 to what we want it to be
pat$x1 <- x1.eval
pat$x2 <- x2.eval

# Turn into an msprep object
pat <- msprep(pat, trans = fit.in[["msm.tmat"]], time = c(NA, "state2", "state3", "state4"),
               status = c(NA, "state2.s", "state3.s", "state4.s"), keep = c("x1","x2"))
# Make it four rows (one for each model)
pat <- pat[rep(1, 4), c("x1", "x2") ]
# Assign trans variable and associate transition matrix
pat$trans <- 1:4
attr(pat, "trans") <- fit.in[["msm.tmat"]]
# Expand covariates to allow different effects per transition, after creating covs object
covs <- c("x1", "x2")
pat <- expand.covs(pat, covs, longnames = FALSE)
pat$strata <- pat$trans


## Apply the multistate model
msf <- msfit(fit.in[["msm.model"]], pat, trans = fit.in[["msm.tmat"]])

## Generate risks for patients from time t = 0
pt <- probtrans(msf, predt = 0)

## Calculate joint risks, joint survival, and marginal survival probabilities, at time t.eval

# Joint risk is probability of being in state 4
risk.joint.est <- pt[[1]]$pstate4[max(which(pt[[1]]$time < t.eval))]

# Risk of A is probbility of being in state 2 or 4
risk.marg.est.A <- (pt[[1]]$pstate2[max(which(pt[[1]]$time < t.eval))] + pt[[1]]$pstate4[max(which(pt[[1]]$time < t.eval))])

# Marginal survival probability of B is probability of not being in state 2 or 4
risk.marg.est.B <- (pt[[1]]$pstate3[max(which(pt[[1]]$time < t.eval))] + pt[[1]]$pstate4[max(which(pt[[1]]$time < t.eval))])

# Joint survival calculated through math (similar eqn used for copulas)
surv.joint.est <- 1 - risk.marg.est.A - risk.marg.est.B + risk.joint.est

## Create output object
output.obj <- c("surv.marg.est.A" = 1 - risk.marg.est.A,
                "surv.marg.est.B" = 1 - risk.marg.est.B,
                "surv.joint.est" = surv.joint.est,
                "risk.marg.est.A" = risk.marg.est.A,
                "risk.marg.est.B" = risk.marg.est.B,
                "risk.joint.est" = risk.joint.est)

return(output.obj)
}



##################
##################
### FIT COPULA ###
##################
##################

### Function for fitting the copula model

fit.model.copula <- function(data.in, copula, rotate.cop){
  
  stopifnot(copula %in% c("clayton","gumbel","joe","fgm","galambos","amh","frank"), local = TRUE)
  
  ## Get data in format for analysis with copula model using GJRM
  data.anal <- tidyr::pivot_wider(data.in, id_cols = id, names_from = outcome_char, values_from = c(time, status, x1, x2))
  
  ## Remove duplicate predictors variables and rename
  data.anal <- data.anal %>% 
    dplyr::select(id, time_A, time_B, status_A, status_B, x1_A, x2_A) %>%
    rename(x1 = x1_A, x2 = x2_A)
  
  ## Fit the copula model
  eq1 <- time_A ~ s(log(time_A), bs = "mpi") + x1 + x2
  eq2 <- time_B ~ s(log(time_B), bs = "mpi") + x1 + x2
  
  ## Assign the copula Bivariate distribution
  if (copula == "normal" & rotate.cop == 0){BivD.in <- "N"}
  if (copula == "clayton" & rotate.cop == 0){BivD.in <- "C0"}
  if (copula == "gumbel" & rotate.cop == 0){BivD.in <- "G0"}
  if (copula == "joe" & rotate.cop == 0){BivD.in <- "J0"}
  if (copula == "galambos" & rotate.cop == 0){BivD.in <- "GAL0"}
  if (copula == "fgm" & rotate.cop == 0){BivD.in <- "FGM"}
  if (copula == "frank" & rotate.cop == 0){BivD.in <- "F"}
  if (copula == "amh" & rotate.cop == 0){BivD.in <- "AMH"}
  if (copula == "clayton" & rotate.cop == 90){BivD.in <- "C90"}
  if (copula == "gumbel" & rotate.cop == 90){BivD.in <- "G90"}
  if (copula == "joe" & rotate.cop == 90){BivD.in <- "J90"}
  if (copula == "galambos" & rotate.cop == 90){BivD.in <- "GAL90"}
  
  gjrm.fit <- gjrm(list(eq1, eq2), data = data.anal, surv = TRUE, 
                   margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = BivD.in)
  
  # summary(gjrm.fit)$tableP1
  # summary(gjrm.fit)$tableP2
  # summary(gjrm.fit)$theta.a
  
  ## Extract coefficients and estimate of association parameter in copula
  betaA.1.est <- data.frame(summary(gjrm.fit)$tableP1)["x1", "Estimate"]
  betaA.2.est <- data.frame(summary(gjrm.fit)$tableP1)["x2", "Estimate"]
  
  betaB.1.est <- data.frame(summary(gjrm.fit)$tableP2)["x1", "Estimate"]
  betaB.2.est <- data.frame(summary(gjrm.fit)$tableP2)["x2", "Estimate"]
  
  eta.est <- summary(gjrm.fit)$theta.a
  
  
  ### Note hat the GJRM is a proportional hazards model, and does not give the baseline hazard, need to calculate
  ### this for each marginal distribution separately.
  
  ## Fit a coxph model
  coxph.cop.A <- coxph(Surv(time_A, status_A) ~ x1 + x2, data = data.anal)
  coxph.cop.B <- coxph(Surv(time_B, status_B) ~ x1 + x2, data = data.anal)
  
  ## Extract baseline hazards
  basehaz.A <- basehaz(coxph.cop.A, centered = FALSE)
  basehaz.B <- basehaz(coxph.cop.B, centered = FALSE)
  
  ### Create output object
  output.obj <- vector("list", 10)
  names(output.obj) <- c("gjrm.fit", "betaA.1.est", "betaA.2.est", "betaB.1.est", "betaB.2.est", 
                         "eta.est", "basehaz.A", "basehaz.B", "copula", "rotate.cop")
  
  output.obj[[1]] <- gjrm.fit
  output.obj[[2]] <- betaA.1.est
  output.obj[[3]] <- betaA.2.est
  output.obj[[4]] <- betaB.1.est
  output.obj[[5]] <- betaB.2.est
  output.obj[[6]] <- eta.est
  output.obj[[7]] <- basehaz.A
  output.obj[[8]] <- basehaz.B
  output.obj[[9]] <- copula
  output.obj[[10]] <- rotate.cop
     
  return(output.obj)
}


#######################
#######################
### PREDRISK COPULA ###
#######################
#######################

### fit.in = the output object from fit.model.copula function
### t.eval is time to predict risks for
### x1.eval = x1 value to assess risks at
### x2.eval = x2 value to assess risks at

calc.risk.copula <- function(fit.in, t.eval, x1.eval, x2.eval){
  
  ## Extract baseline hazards for marginal distributions from fitted copula model object
  basehaz.A <- fit.in[["basehaz.A"]]
  basehaz.B <- fit.in[["basehaz.B"]]
  
  ## Extract baseline hazards at the time point of interest
  basehaz.A.t <- basehaz.A$hazard[max(which(basehaz.A$time < t.eval))]
  basehaz.B.t <- basehaz.B$hazard[max(which(basehaz.B$time < t.eval))]
  
  ## Calculate marginal risks
  surv.marg.est.A <- exp(-basehaz.A.t*exp(fit.in[["betaA.1.est"]]*x1.eval + fit.in[["betaA.2.est"]]*x2.eval))
  surv.marg.est.B <- exp(-basehaz.B.t*exp(fit.in[["betaB.1.est"]]*x1.eval + fit.in[["betaB.2.est"]]*x2.eval))
  
  risk.marg.est.A <- 1 - exp(-basehaz.A.t*exp(fit.in[["betaA.1.est"]]*x1.eval + fit.in[["betaA.2.est"]]*x2.eval))
  risk.marg.est.B <- 1 - exp(-basehaz.B.t*exp(fit.in[["betaB.1.est"]]*x1.eval + fit.in[["betaB.2.est"]]*x2.eval))
  
  ## Calculate copula using association parameter estimated from the data
  if (fit.in[["copula"]] %in% c("fgm")){
    cl.est <- fgmCopula(param = fit.in[["eta.est"]], dim = 2)
  } else {
    cl.est <- archmCopula(family = tolower(fit.in[["copula"]]), param = fit.in[["eta.est"]], dim = 2)
  }
  
  ## Rotate if required
  if (fit.in[["rotate.cop"]] == 90){
    cl.est <- rotCopula(cl.est, flip = c(TRUE, FALSE))
  } 
  
  ## Calculate joint survival probability
  surv.joint.est <- pCopula(c(surv.marg.est.A,surv.marg.est.B), cl.est)
  
  ## Calculate joint risk
  risk.joint.est <- 1 - surv.marg.est.A - surv.marg.est.B + surv.joint.est 
  
  ## Create output object
  output.obj <- c("surv.marg.est.A" = surv.marg.est.A,
                  "surv.marg.est.B" = surv.marg.est.B,
                  "surv.joint.est" = surv.joint.est,
                  "risk.marg.est.A" = risk.marg.est.A,
                  "risk.marg.est.B" = risk.marg.est.B,
                  "risk.joint.est" = risk.joint.est)
  
  return(output.obj)
}


################################################
################################################
### Fit models for product of marginal risks ###
################################################
################################################

### Write a function to fit the marginal coxph models
fit.model.product <- function(data.in){
  
  ## Get data in format for analysis with copula model using GJRM
  data.anal <- tidyr::pivot_wider(data.in, id_cols = id, names_from = outcome_char, values_from = c(time, status, x1, x2))
  
  ## Remove duplicate predictors variables and rename
  data.anal <- data.anal %>% 
    dplyr::select(id, time_A, time_B, status_A, status_B, x1_A, x2_A) %>%
    rename(x1 = x1_A, x2 = x2_A)
  
  ## Make into dataframe
  data.anal <- data.frame(data.anal)
  
  ## Fit a coxph model for each outcome
  coxph.A <- coxph(Surv(time_A, status_A) ~ x1 + x2, data = data.anal)
  coxph.B <- coxph(Surv(time_B, status_B) ~ x1 + x2, data = data.anal)
  
  ## Extract baseline hazards
  basehaz.A <- basehaz(coxph.A, centered = FALSE)
  basehaz.B <- basehaz(coxph.B, centered = FALSE)
  
  ## Create output object
  output.object <- vector("list", 4)
  names(output.object) <- c("coef.A", "coef.B", "basehaz.A", "basehaz.B")
  
  ## Assign required output to output object
  output.object[[1]] <- coxph.A$coefficients
  output.object[[2]] <- coxph.B$coefficients
  output.object[[3]] <- basehaz.A
  output.object[[4]] <- basehaz.B
  
  return(output.object)
}


###################################################
###################################################
### Predict risks for product of marginal risks ###
###################################################
###################################################

calc.risk.product <- function(fit.in, t.eval, x1.eval, x2.eval){
  
  ## Extract baseline hazards for marginal distributions from fitted copula model object
  basehaz.A <- fit.in[["basehaz.A"]]
  basehaz.B <- fit.in[["basehaz.B"]]
  
  ## Extract baseline hazards at the time point of interest
  basehaz.A.t <- basehaz.A$hazard[max(which(basehaz.A$time < t.eval))]
  basehaz.B.t <- basehaz.B$hazard[max(which(basehaz.B$time < t.eval))]
  
  ## Calculate marginal survival probabilities and risks
  surv.marg.est.A <- as.numeric(exp(-basehaz.A.t*exp(fit.in[["coef.A"]][1]*x1.eval + fit.in[["coef.A"]][2]*x2.eval)))
  surv.marg.est.B <- as.numeric(exp(-basehaz.B.t*exp(fit.in[["coef.B"]][1]*x1.eval + fit.in[["coef.B"]][2]*x2.eval)))
  
  risk.marg.est.A <- 1 - as.numeric(exp(-basehaz.A.t*exp(fit.in[["coef.A"]][1]*x1.eval + fit.in[["coef.A"]][2]*x2.eval)))
  risk.marg.est.B <- 1 - as.numeric(exp(-basehaz.B.t*exp(fit.in[["coef.B"]][1]*x1.eval + fit.in[["coef.B"]][2]*x2.eval)))
  
  ## Calculate joint survival probability
  surv.joint.est <- surv.marg.est.A*surv.marg.est.B
  
  ## Calculate joint risk
  risk.joint.est <- risk.marg.est.A*risk.marg.est.B
  
  ## Create output object
  output.obj <- c("surv.marg.est.A" = surv.marg.est.A,
                  "surv.marg.est.B" = surv.marg.est.B,
                  "surv.joint.est" = surv.joint.est,
                  "risk.marg.est.A" = risk.marg.est.A,
                  "risk.marg.est.B" = risk.marg.est.B,
                  "risk.joint.est" = risk.joint.est)
  
  ## Return output object
  return(output.obj)
}
  

#######################################
#######################################
### Fit model for joint risk method ###
#######################################
#######################################

### Write a function to fit the marginal coxph models
fit.model.joint <- function(data.in){
  
  ## Get data in format for analysis with copula model using GJRM
  data.anal <- tidyr::pivot_wider(data.in, id_cols = id, names_from = outcome_char, values_from = c(time, status, x1, x2))
  
  ## Remove duplicate predictors variables and rename
  data.anal <- data.anal %>% 
    dplyr::select(id, time_A, time_B, status_A, status_B, x1_A, x2_A) %>%
    rename(x1 = x1_A, x2 = x2_A)
  
  ## Create a variable for until until both A and B, and the appropriate status indicator (both must be uncensored to observe event)
  data.anal <- data.anal %>%
    mutate(time_joint = pmax(time_A, time_B),
           status_joint = pmin(status_A, status_B))
  
  ## Make into dataframe
  data.anal <- data.frame(data.anal)
  
  ## Fit a coxph model
  coxph.joint <- coxph(Surv(time_joint, status_joint) ~ x1 + x2, data = data.anal)
  
  ## Extract baseline hazards
  basehaz.joint <- basehaz(coxph.joint, centered = FALSE)
  
  ## Create output object
  output.object <- vector("list", 2)
  names(output.object) <- c("coef.joint", "basehaz.joint")
  
  ## Assign required output to output object
  output.object[[1]] <- coxph.joint$coefficients
  output.object[[2]] <- basehaz.joint
  
  return(output.object)
}


###########################################
###########################################
### Predict risks for joint risk method ###
###########################################
###########################################

calc.risk.joint <- function(fit.in, t.eval, x1.eval, x2.eval){

  ## Extract baseline hazards for marginal distributions from fitted copula model object
  basehaz.joint <- fit.in[["basehaz.joint"]]
  
  ## Extract baseline hazards at the time point of interest
  basehaz.joint.t <- basehaz.joint$hazard[max(which(basehaz.joint$time < t.eval))]
  
  ## Calculate joint risks
  risk.joint.est <- 1 - as.numeric(exp(-basehaz.joint.t*exp(fit.in[["coef.joint"]][1]*x1.eval + fit.in[["coef.joint"]][2]*x2.eval)))
  
  ## Create output object
  output.obj <- c("risk.joint.est" = risk.joint.est)
  
  ## Return output object
  return(output.obj)

}



####################################
####################################
### Fit model for frailty method ###
####################################
####################################

### Function for fitting the copula model

fit.model.frailty <- function(data.in, frail.dist, baseline.dist, n.iter){
  
  stopifnot(frail.dist %in% c("gamma","normal"), local = TRUE)
  
#   data.in <- dat.frailty.gamma
#   frail.dist <- "gamma"
#   baseline.dist <- "exp"
  
    ## Get data in format for analysis with rstan
    data.in.uncens.A <- filter(data.in, status == 1 & outcome == 1)
    data.in.uncens.B <- filter(data.in, status == 1 & outcome == 2)
    data.in.cens.A <- filter(data.in, status == 0 & outcome == 1)
    data.in.cens.B <- filter(data.in, status == 0 & outcome == 2)

    stan_data <- list(N = length(unique(data.in$id)),
                         N_cens_A = length(unique(data.in.cens.A$id)),
                         N_uncens_A = length(unique(data.in.uncens.A$id)),
                         N_cens_B = length(unique(data.in.cens.B$id)),
                         N_uncens_B = length(unique(data.in.uncens.B$id)),
                         #stacked_N = nrow(data.in.cens),
                         
                         P = 2,
                         
                         IDs_cens_A = data.in.cens.A$id,
                         IDs_uncens_A = data.in.uncens.A$id,
                         IDs_cens_B = data.in.cens.B$id,
                         IDs_uncens_B = data.in.uncens.B$id,
                         
                         X_cens_A = data.in.cens.A %>% dplyr::select(x1, x2) %>% data.matrix(),
                         X_uncens_A = data.in.uncens.A %>% dplyr::select(x1, x2) %>% data.matrix(),
                         X_cens_B = data.in.cens.B %>% dplyr::select(x1, x2) %>% data.matrix(),
                         X_uncens_B = data.in.uncens.B %>% dplyr::select(x1, x2) %>% data.matrix(),
                         
                         times_cens_A = data.in.cens.A$time,
                         times_uncens_A = data.in.uncens.A$time,
                         times_cens_B = data.in.cens.B$time,
                         times_uncens_B = data.in.uncens.B$time)
  
  
  ## Assign the appropriate stan model depending on what frailty distribution
  if (frail.dist == "normal" & baseline.dist == "exp"){
    sm <- stan_model("code/ExpSurvModel_multivariate_vectorized_normal_exp.stan")
    fit <- sampling(sm, data=stan_data, seed=42, chains=1, cores=1, iter=n.iter)
  } else if (frail.dist == "gamma" & baseline.dist == "exp"){
    sm <- stan_model("code/ExpSurvModel_multivariate_vectorized_gamma_exp.stan")
    fit <- sampling(sm, data=stan_data, seed=42, chains=1, cores=1, iter=n.iter)
  } else if (frail.dist == "normal" & baseline.dist == "weibull"){
    sm <- stan_model("code/ExpSurvModel_multivariate_vectorized_normal_weibull_inverted.stan")
    fit <- sampling(sm, data=stan_data, seed=42, chains=1, cores=1, iter=n.iter)
  } else if (frail.dist == "gamma" & baseline.dist == "weibull"){
    sm <- stan_model("code/ExpSurvModel_multivariate_vectorized_gamma_weibull_inverted.stan")
    fit <- sampling(sm, data=stan_data, seed=42, chains=1, cores=1, iter=n.iter)
  }
  
  ## Assign estimates to an object and output
  if (baseline.dist == "exp"){
    model.estimates <- summary(fit,  par = c("betas_A[1]", "betas_A[2]",
                                             "betas_B[1]", "betas_B[2]",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  } else if (baseline.dist == "weibull"){
    model.estimates <- summary(fit,  par = c("betas_A[1]", "betas_A[2]",
                                             "betas_B[1]", "betas_B[2]",
                                             "shape_A", "shape_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  }
  print("output assigned")
  ### Create output object
  output.obj <- vector("list", 11)
  names(output.obj) <- c("betaA.1", "betaA.2", "betaB.1", "betaB.2", 
                         "bh.shape.A", "bh.int.A", "bh.shape.B", "bh.int.B",
                         "frail.var.est", 
                         "frail.dist", "bh.dist")
  
  
  if (baseline.dist == "weibull"){
    output.obj[[1]] <- model.estimates$summary["betas_A[1]","mean"]
    output.obj[[2]] <- model.estimates$summary["betas_A[2]","mean"]
    output.obj[[3]] <- model.estimates$summary["betas_B[1]","mean"]
    output.obj[[4]] <- model.estimates$summary["betas_B[2]","mean"]
    output.obj[[5]] <- model.estimates$summary["shape_A","mean"]
    output.obj[[6]] <- model.estimates$summary["intercept_A","mean"]
    output.obj[[7]] <- model.estimates$summary["shape_B","mean"]
    output.obj[[8]] <- model.estimates$summary["intercept_B","mean"]
  } else {
    output.obj[[1]] <- model.estimates$summary["betas_A[1]","mean"]
    output.obj[[2]] <- model.estimates$summary["betas_A[2]","mean"]
    output.obj[[3]] <- model.estimates$summary["betas_B[1]","mean"]
    output.obj[[4]] <- model.estimates$summary["betas_B[2]","mean"]
    output.obj[[5]] <- 1
    output.obj[[6]] <- model.estimates$summary["intercept_A","mean"]
    output.obj[[7]] <- 1
    output.obj[[8]] <- model.estimates$summary["intercept_B","mean"]
  }
  output.obj[[9]] <- model.estimates$summary["frail_param","mean"]
  output.obj[[10]] <- frail.dist
  output.obj[[11]] <- baseline.dist
  
  ## Return output object
  return(output.obj)
}



########################
########################
### PREDRISK FRAILTY ###
########################
########################

### fit.in = the output object from fit.model.copula function
### t.eval is time to predict risks for
### x1.eval = x1 value to assess risks at
### x2.eval = x2 value to assess risks at

### First I need to calculate the baseline hazard at time t.eval.
### I want to do it seperately
### This needs to be done first an
calc.risk.frailty <- function(fit.in, t.eval, x1.eval, x2.eval){
  
#   fit.in <- model.frailty.gamma.exponential
#   t.eval <- 10
#   x1.eval <- 0.1
#   x2.eval <- 0
  
  ## Create functions for the conditional survival probabilities on frailty term. Note that it is also technically conditional
  ## on x1 and x2, but these have already been inputted at the start of the calc.risk.frailty function, and we want to keep these
  ## constant, as we aren't integrating over them
  if (fit.in[["frail.dist"]] == "gamma"){
    surv.marg.A.func <- function(frail.term){
      return(pweibull(t.eval, 
                      shape = fit.in[["bh.shape.A"]], 
                      scale = 1/(frail.term*exp(fit.in[["bh.int.A"]] + fit.in[["betaA.1"]]*x1.eval + fit.in[["betaA.2"]]*x2.eval)),
                      lower.tail = FALSE))
    }
    surv.marg.B.func <- function(frail.term){
      return(pweibull(t.eval, 
                      shape = fit.in[["bh.shape.B"]], 
                      scale = 1/(frail.term*exp(fit.in[["bh.int.B"]] + fit.in[["betaB.1"]]*x1.eval + fit.in[["betaB.2"]]*x2.eval)),
                      lower.tail = FALSE))
    }
    risk.marg.A.func <- function(frail.term){
      return(1-pweibull(t.eval, 
                      shape = fit.in[["bh.shape.A"]], 
                      scale = 1/(frail.term*exp(fit.in[["bh.int.A"]] + fit.in[["betaA.1"]]*x1.eval + fit.in[["betaA.2"]]*x2.eval)),
                      lower.tail = FALSE))
    }
    risk.marg.B.func <- function(frail.term){
      return(1-pweibull(t.eval, 
                      shape = fit.in[["bh.shape.B"]], 
                      scale = 1/(frail.term*exp(fit.in[["bh.int.B"]] + fit.in[["betaB.1"]]*x1.eval + fit.in[["betaB.2"]]*x2.eval)),
                      lower.tail = FALSE))
    }
  } else if (fit.in[["frail.dist"]] == "normal"){
    surv.marg.A.func <- function(frail.term){
      return(pweibull(t.eval, 
                      shape = fit.in[["bh.shape.A"]], 
                      scale = 1/(exp(frail.term + fit.in[["bh.int.A"]] + fit.in[["betaA.1"]]*x1.eval + fit.in[["betaA.2"]]*x2.eval)),
                      lower.tail = FALSE))
    }
    surv.marg.B.func <- function(frail.term){
      return(pweibull(t.eval, 
                      shape = fit.in[["bh.shape.B"]], 
                      scale = 1/(exp(frail.term + fit.in[["bh.int.B"]] + fit.in[["betaB.1"]]*x1.eval + fit.in[["betaB.2"]]*x2.eval)),
                      lower.tail = FALSE))
    }
    risk.marg.A.func <- function(frail.term){
      return(1-pweibull(t.eval, 
                        shape = fit.in[["bh.shape.A"]], 
                        scale = 1/(exp(frail.term + fit.in[["bh.int.A"]] + fit.in[["betaA.1"]]*x1.eval + fit.in[["betaA.2"]]*x2.eval)),
                        lower.tail = FALSE))
    }
    risk.marg.B.func <- function(frail.term){
      return(1-pweibull(t.eval, 
                        shape = fit.in[["bh.shape.B"]], 
                        scale = 1/(exp(frail.term + fit.in[["bh.int.B"]] + fit.in[["betaB.1"]]*x1.eval + fit.in[["betaB.2"]]*x2.eval)),
                        lower.tail = FALSE))
    }
  }
  
  ### Integrate these functions over the frailty distributions to get the marginal risks/surival probabilities,
  ### and Integrate the product of these functions over the frailty distributions to get the joint risks/survival probabilities
  
  ## First need to create functions which are the product of the entity we want to integrate over, and the density function of the 
  ## random effect (with variance estimated from the data)
  if (fit.in[["frail.dist"]] == "gamma"){
    surv.marg.est.A.for.int <- function(z.in){
      surv.marg.A.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    surv.marg.est.B.for.int <- function(z.in){
      surv.marg.B.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    risk.marg.est.A.for.int <- function(z.in){
      risk.marg.A.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    risk.marg.est.B.for.int <- function(z.in){
      risk.marg.B.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    surv.joint.est.for.int <- function(z.in){
      surv.marg.A.func(z.in)*surv.marg.B.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    risk.joint.est.for.int <- function(z.in){
      risk.marg.A.func(z.in)*risk.marg.B.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
  } else if (fit.in[["frail.dist"]] == "normal"){
    surv.marg.est.A.for.int <- function(z.in){
      surv.marg.A.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    surv.marg.est.B.for.int <- function(z.in){
      surv.marg.B.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    risk.marg.est.A.for.int <- function(z.in){
      risk.marg.A.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    risk.marg.est.B.for.int <- function(z.in){
      risk.marg.B.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    surv.joint.est.for.int <- function(z.in){
      surv.marg.A.func(z.in)*surv.marg.B.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    risk.joint.est.for.int <- function(z.in){
      risk.marg.A.func(z.in)*risk.marg.B.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
  }
  
  ## Now do the integration
  if (fit.in[["frail.dist"]] == "gamma"){
    ## Define upper limit for integration at the 99.999th percentile of the ditribution
    upper.lim <- qgamma(0.99999, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    
    ## Now run the imputation
    surv.marg.est.A <- cubintegrate(f = surv.marg.est.A.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    surv.marg.est.B <- cubintegrate(f = surv.marg.est.B.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    risk.marg.est.A <- cubintegrate(f = risk.marg.est.A.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    risk.marg.est.B <- cubintegrate(f = risk.marg.est.B.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    surv.joint.est <- cubintegrate(f = surv.joint.est.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    risk.joint.est <- cubintegrate(f = risk.joint.est.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
  } else if (fit.in[["frail.dist"]] == "normal"){
    ## Define lower and upper limit for integration at the 0.01th 99.99th percentile of the ditribution
    lower.lim <- qnorm(0.00001, mean = 0, sd = fit.in[["frail.var.est"]])
    upper.lim <- qnorm(0.99999, mean = 0, sd = fit.in[["frail.var.est"]])
    
    ## Now run the imputation
    surv.marg.est.A <- cubintegrate(f = surv.marg.est.A.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    surv.marg.est.B <- cubintegrate(f = surv.marg.est.B.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    risk.marg.est.A <- cubintegrate(f = risk.marg.est.A.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    risk.marg.est.B <- cubintegrate(f = risk.marg.est.B.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    surv.joint.est <- cubintegrate(f = surv.joint.est.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    risk.joint.est <- cubintegrate(f = risk.joint.est.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
  }
  
  ## Create output object
  output.obj <- c("surv.marg.est.A" = surv.marg.est.A,
                  "surv.marg.est.B" = surv.marg.est.B,
                  "surv.joint.est" = surv.joint.est,
                  "risk.marg.est.A" = risk.marg.est.A,
                  "risk.marg.est.B" = risk.marg.est.B,
                  "risk.joint.est" = risk.joint.est)
  
  return(output.obj)
}

save.image("data/sim_function_fit_models_all_2cont.RData")
# 
# # 
# load("data/sim_function_fit_models_all.RData")
# load("data/sim_function_generate_data_all.RData")
# load("data/sim_function_calc_true_risk_all.RData")
# # 
# # 
# ### Define the number of individuals we want to generate
# n.sim <- 1000
# 
# 
# 
# ## Set number of times we run the simultion
# k.sim <- 50
# 
# ## Define a vector which indicates which analysis methods should be used
# ## OPTIONS ARE: c("product", "joint", "msm", "c.clay", "c.clay.rot", "c.gumb", "c.gumb.rot", "c.joe", "c.joe.rot", "c.fgm")
# anal.methods.sim <- c("msm", "product", "joint", "c.clay", "c.gumb", "c.fgm", "f.normal", "f.gamma", "f.normal_weib", "f.gamma_weib")
# 
# ## Define what values of x1 we will evaluate strong calibration at
# str.calib.eval.sim <- c(-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2)
# 
# 
# 
# 
# 
# 
# 
# 
# ## Size of development/validation dataset
# n.devel.sim <- 1000
# n.valid.sim <- 200
# 
# ## Maximum follow up before censoring
# max.follow.sim <- 100000
# 
# ## Followup time at which to evaluate
# t.eval.sim <- 10
# 
# ### The following paramters are specific to the data generating mechanism
# ## Shape and scale paramters for baseline hazard of each transition
#   shape12.sim <- 1
#   scale12.sim <- 15
#   shape13.sim <- 1
#   scale13.sim <- 10
#   shape24.sim <- 1
#   scale24.sim <- 10
#   shape34.sim <- 1 
#   scale34.sim <- 15
#   
#   ## Covariate effects for each transition
#   beta12.cont.sim <- 0.25
#   beta12.cat.sim <- 0.5
#   beta13.cont.sim <- 0.3
#   beta13.cat.sim <- 0.25
#   beta24.cont.sim <- 0.1
#   beta24.cat.sim <- 0.1
#   beta34.cont.sim <- 0.25
#   beta34.cat.sim <- 0.5
# 
# beta12.cont.sim <- 0
# beta12.cat.sim <- 0
# beta13.cont.sim <- 0
# beta13.cat.sim <- 0
# beta24.cont.sim <- 0
# beta24.cat.sim <- 0
# beta34.cont.sim <- 0
# beta34.cat.sim <- 0
#   
#   ## Number of sampling steps when simulating data using multistate model
#   numsteps.sim <- 50000
#   
# 
#   ## Baseline hazard parameters for marginal distributions
#   baselineA.sim <- c(1,1.5)
#   baselineB.sim <- c(1,1.0)
#   ## Covariate effects
#   COV_betaA.sim <- c(0.5, 0.5)
#   COV_betaB.sim <- c(0.25, 0.75)
#   ## Copula association parameter
#   copula.param.sim <- 3
#   ## Name of copula and rotation
#   if (DGM.sim == "c.clay"){
#     copula.sim <- "clayton"
#     rotate.cop.sim <- 0
#   } else if (DGM.sim == "c.clay.rot"){
#     copula.sim <- "clayton"
#     rotate.cop.sim <- 90
#   } else if (DGM.sim == "c.gumb"){
#     copula.sim <- "gumbel"
#     rotate.cop.sim <- 0
#   } else if (DGM.sim == "c.clay.rot"){
#     copula.sim <- "gumbel"
#     rotate.cop.sim <- 90
#   } else if (DGM.sim == "c.joe"){
#     copula.sim <- "joe"
#     rotate.cop.sim <- 0
#   } else if (DGM.sim == "c.joe.rot"){
#     copula.sim <- "joe"
#     rotate.cop.sim <- 90
#   } else if (DGM.sim == "c.fgm"){
#     copula.sim <- "fgm"
#     rotate.cop.sim <- 0
#   }
# 
#   ## Baseline hazard parameters for marginal distributions
#   baselineA.sim <- c(1,0.15)
#   baselineB.sim <- c(1,0.1)
#   ## Covariate effects
#   COV_betaA.sim <- c(0.5, 0.5)
#   COV_betaB.sim <- c(0.25, 0.75)
#   ## Frailty association parameter
#   frail.var.sim <- 0.3
#   ## Assigm distribution
#   if (DGM.sim == "f.normal"){
#     frail.dist.sim <- "normal"
#     int.method.sim <- "pcubature"
#   } else if (DGM.sim == "f.gamma"){
#     frail.dist.sim <- "gamma"
#     int.method.sim <- "hcubature"
#   }
# 
# ## Shape/scale/covariate effects for censoring mechanism
# baseline_cens.sim <- c(1, 15000)
# COV_beta_cens.sim <- c(0, 0)
# 
# ### Set.seed
# set.seed(505)
# 
# ### Create the dataset of baseline characteristics
# x.baseline.devel <- data.frame("x1" = rnorm(n.devel.sim, 0, 1), "x2" = rbinom(n.devel.sim, 1, 0.5))
# 
# dat.devel <- gen.dat.copula(n = n.devel.sim, 
#                             max.follow = max.follow.sim,
#                             copula = "clayton", copula.param = copula.param.sim, rotate.cop = 0,
#                             baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
#                             baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
#                             baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim,
#                             x.in = x.baseline.devel)
# 
# str(dat.devel)
# 1/365
# min(dat.devel$time)
# min(dat.devel$cens.time)
# 1/365
# dat.devel$time[dat.devel$time < 0.00005] <- 0.00005
# 0.00001
# 
# fit.c.clay <- fit.model.copula(dat.devel, copula = "clayton", rotate.cop = 0)
# fit.c.fgm <- fit.model.copula(dat.devel, copula = "fgm", rotate.cop = 0)
# fit.c.gumbel<- fit.model.copula(dat.devel, copula = "gumbel", rotate.cop = 0)
# str(fit.c.clay)
# 
# dat.devel.error <- dat.devel
# min(dat.devel.error$time)
# 
# dat.devel.error2 <- dat.devel
# min(dat.devel.error2$time)
# 
# dat.devel.error3 <- dat.devel
# min(dat.devel.error2$time)
# 
# dat.devel <- gen.dat.frailty(n = n.devel.sim, 
#                              max.follow = max.follow.sim,
#                              frail.dist = "gamma", frail.eff = 1, frail.var = frail.var.sim,
#                              baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
#                              baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
#                              baseline_cens = baseline_cens.sim, COV_beta_cens = COV_beta_cens.sim,
#                              x.in = x.baseline.devel)
# 
# fit.c.clay <- fit.model.copula(dat.devel, copula = "clayton", rotate.cop = 0)
# fit.c.fgm <- fit.model.copula(dat.devel, copula = "fgm", rotate.cop = 0)
# fit.c.gumbel<- fit.model.copula(dat.devel, copula = "gumbel", rotate.cop = 0)
# 
# fit.c.clay[[1]]
# 
# ################################################################
# ### FRAILTY ###
# ################################################################
# # dat.frailty.gamma <- gen.dat.frailty(n = n.sim, 
# #                                      max.follow = 10000,
# #                                      frail.dist = "gamma", frail.eff = 1, frail.var = 0.7,
# #                                      baselineA = c(1,15), COV_betaA = c(0.2, 0.5), 
# #                                      baselineB = c(1,10), COV_betaB = c(0.1, 0.4), 
# #                                      baseline_cens = c(1, 20), COV_beta_cens = c(0.25, 0.25),
# #                                      x.in = x.baseline)
# 
# dat.frailty.normal <- gen.dat.frailty(n = n.sim, 
#                                       max.follow = 10000,
#                                       frail.dist = "normal", frail.eff = 1, frail.var = 0.3,
#                                       baselineA = c(1,15), COV_betaA = c(0.2, 0.5), 
#                                       baselineB = c(1,10), COV_betaB = c(0.1, 0.4), 
#                                       baseline_cens = c(1, 20), COV_beta_cens = c(0.25, 0.25),
#                                       x.in = x.baseline)
# 
# ## Fit models
# model.frailty.gamma.weibull <- fit.model.frailty(data.in = dat.frailty.gamma, 
#                                            frail.dist = "gamma", 
#                                            baseline.dist = "weibull")
# 
# model.frailty.gamma.exponential <- fit.model.frailty(data.in = dat.frailty.gamma, 
#                                                frail.dist = "gamma", 
#                                                baseline.dist = "exp")
# 
# 
# 
# 
# model.frailty.normal.weibull <- fit.model.frailty(data.in = dat.frailty.normal, 
#                                                 frail.dist = "normal", 
#                                                 baseline.dist = "weibull")
# 
# model.frailty.normal.exponential <- fit.model.frailty(data.in = dat.frailty.normal, 
#                                                     frail.dist = "normal", 
#                                                     baseline.dist = "exp")
# 
# 
# risks.gamma.weibull <- calc.risk.frailty(model.frailty.gamma.weibull, t.eval = 10, x1.eval = 0.1, x2.eval = 0)
# risks.gamma.exponential <- calc.risk.frailty(model.frailty.gamma.exponential, t.eval = 10, x1.eval = 0.1, x2.eval = 0)
# 
# risks.normal.weibull <- calc.risk.frailty(model.frailty.normal.weibull, t.eval = 10, x1.eval = 0.1, x2.eval = 0)
# risks.normal.exponential <- calc.risk.frailty(model.frailty.normal.exponential, t.eval = 10, x1.eval = 0.1, x2.eval = 0)
# 
# 
# 
# 
# model.frailty.gamma.weibull
# model.frailty.gamma.exponential
# model.frailty.normal.weibull
# model.frailty.norml.exponential
# 
# 
# 1/exp(-2.91)
# true.risk.frailty.gamma <- calc.true.risk.frailty(
#   t.eval = 10, x1.eval = 0.1, x2.eval = 0,
#   frail.dist = "gamma", frail.var = 0.7,
#   baselineA <- c(1,15), COV_betaA = c(0.2, 0.5), 
#   baselineB <- c(1,10), COV_betaB = c(0.1, 0.4),
#   int.lower = 0.0001, int.upper = 20, int.method = "hcubature")
# 
# 
# 
# true.risk.frailty.normal <- calc.true.risk.frailty(
#   t.eval = 10, x1.eval = 0.1, x2.eval = 0,
#   frail.dist = "normal", frail.var = 0.3,
#   baselineA <- c(1,15), COV_betaA = c(0.2, 0.5), 
#   baselineB <- c(1,10), COV_betaB = c(0.1, 0.4),
#   int.lower = -10, int.upper = 10, int.method = "pcubature")
# 
# risks.gamma.weibull
# risks.gamma.exponential
# risks.normal.weibull
# risks.normal.exponential
# 
# 
# true.risk.frailty.gamma
# true.risk.frailty.normal
# 
# 
# load("data/test_frailty_model.RData")
# 
# 
# 
# 
# # 
# #load("data/summary_simulation_so_far.RData")
# 
# dat.devel <- dat.msm[1:2000, ]
# dat.valid <- dat.msm[5001: 5500, ]
#   
# fit.product <- fit.model.product(dat.devel)
# pred.risk.product <- calc.risk.product(fit.product, t.eval = 10, x1.eval = 0.1, x2.eval = 0)  
# 
# fit.joint <- fit.model.joint(dat.devel)
# pred.risk.joint <- calc.risk.joint(fit.joint, t.eval = 10, x1.eval = 0.1, x2.eval = 0)  
# pred.risk.product
# pred.risk.joint
# 
# fit.copula.clayton <- fit.model.copula(dat.devel, copula = "clayton", rotate.cop = 0)
# 
# summary(fit.copula.clayton[[1]])
# fit.copula.clayton[[9]]
# fit.copula.clayton[[10]]
# 
# pred.risk.copula.clayton <- calc.risk.copula(fit.in = fit.copula.clayton, t.eval = 10, x1.eval = 0.1, x2.eval = 0)
# pred.risk.copula.clayton
# 
# fit.msm <- fit.model.msm(dat.devel)
# pred.risk.msm <- calc.risk.msm(fit.msm, t.eval = 10, x1.eval = 0.1, x2.eval = 0)
# 
# pred.risk.copula.clayton
# pred.risk.msm
# 
# 
# 
# 
# 
# ###########################################
# ### MSM
# ################################################
# ### MSM
# n.devel.sim <- 30000
# 
# ### Set.seed
# set.seed(505)
# 
# ### Create the dataset of baseline characteristics
# x.baseline.devel <- data.frame("x1" = rnorm(n.devel.sim, 0, 1), "x2" = rbinom(n.devel.sim, 1, 0.5))
# 
# dat.msm <- gen.dat.msm(n = n.devel.sim,
#                        max.follow = 100000,
#                        shape12 = 1, scale12 = 15, #shape and scale for weibull baseline hazard for transition 1 -> 2
#                        shape13 = 1, scale13 = 10, #shape and scale for weibull baseline hazard for transition 1 -> 3
#                        shape24 = 1, scale24 = 8, #shape and scale for weibull baseline hazard for transition 2 -> 4
#                        shape34 = 1, scale34 = 13,
#                        beta12.cont = 0.25, beta12.cat = 0.5,
#                        beta13.cont = 0.3, beta13.cat = 0.25,
#                        beta24.cont = 0.3, beta24.cat = 0.25,
#                        beta34.cont = 0.25, beta34.cat = 0.5, #covariate effects for transiion 34
#                        baseline_cens = c(1,15000), COV_beta_cens = c(0, 0),
#                        x.in = x.baseline.devel,
#                        numsteps = 50000)
# 
# 
# true.risk.msm <- calc.true.risk.msm(t.eval = 10, x1.eval = 1.5, x2.eval = 1,
#                                     shape12 = 1, scale12 = 15,
#                                     shape13 = 1, scale13 = 10,
#                                     shape24 = 1, scale24 = 8,
#                                     shape34 = 1, scale34 = 13,
#                                     beta12.cont = 0.25, beta12.cat = 0.5,
#                                     beta13.cont = 0.3, beta13.cat = 0.25,
#                                     beta24.cont = 0.3, beta24.cat = 0.25,
#                                     beta34.cont = 0.25, beta34.cat = 0.5)
# 
# 
# calc.true.risk.msm(t.eval = 10, x1.eval = 1.5, x2.eval = 1,
#                    shape12 = 1, scale12 = 15,
#                    shape13 = 1, scale13 = 10,
#                    shape24 = 1, scale24 = 8,
#                    shape34 = 1, scale34 = 13,
#                    beta12.cont = 0.25, beta12.cat = 0.5,
#                    beta13.cont = 0.3, beta13.cat = 0.25,
#                    beta24.cont = 0.3, beta24.cat = 0.25,
#                    beta34.cont = 0.25, beta34.cat = 0.5)
# 
# 
# calc.true.risk.msm(t.eval = 10, x1.eval = 0, x2.eval = 0,
#                    shape12 = 1, scale12 = 15,
#                    shape13 = 1, scale13 = 10,
#                    shape24 = 1, scale24 = 8,
#                    shape34 = 1, scale34 = 13,
#                    beta12.cont = 0.25, beta12.cat = 0.5,
#                    beta13.cont = 0.3, beta13.cat = 0.25,
#                    beta24.cont = 0.3, beta24.cat = 0.25,
#                    beta34.cont = 0.25, beta34.cat = 0.5)
# 
# fit.msm <- fit.model.msm(dat.msm)
# pred.risk.msm <- calc.risk.msm(fit.msm, t.eval = 10, x1.eval = 1.5, x2.eval = 1)
# 
# true.risk.msm
# pred.risk.msm
# 
# 
# 
# 
# n.devel.sim <- 30000
# 
# ### Set.seed
# set.seed(505)
# 
# ### Create the dataset of baseline characteristics
# x.baseline.devel.fix <- data.frame("x1" = rep(1.5, n.devel.sim), "x2" = rep(1, n.devel.sim))
# 
# dat.msm.fix <- gen.dat.msm(n = n.devel.sim,
#                        max.follow = 100000,
#                        shape12 = 1, scale12 = 15, #shape and scale for weibull baseline hazard for transition 1 -> 2
#                        shape13 = 1, scale13 = 10, #shape and scale for weibull baseline hazard for transition 1 -> 3
#                        shape24 = 1, scale24 = 8, #shape and scale for weibull baseline hazard for transition 2 -> 4
#                        shape34 = 1, scale34 = 13,
#                        beta12.cont = 0.25, beta12.cat = 0.5,
#                        beta13.cont = 0.3, beta13.cat = 0.25,
#                        beta24.cont = 0.3, beta24.cat = 0.25,
#                        beta34.cont = 0.25, beta34.cat = 0.5, #covariate effects for transiion 34
#                        baseline_cens = c(1,30000), COV_beta_cens = c(0, 0),
#                        x.in = x.baseline.devel.fix,
#                        numsteps = 50000)
# 
# true.risk.msm <- calc.true.risk.msm(t.eval = 10, x1.eval = 1.5, x2.eval = 1,
#                                     shape12 = 1, scale12 = 15,
#                                     shape13 = 1, scale13 = 10,
#                                     shape24 = 1, scale24 = 8,
#                                     shape34 = 1, scale34 = 13,
#                                     beta12.cont = 0.25, beta12.cat = 0.5,
#                                     beta13.cont = 0.3, beta13.cat = 0.25,
#                                     beta24.cont = 0.3, beta24.cat = 0.25,
#                                     beta34.cont = 0.25, beta34.cat = 0.5)
# 
# 
# dat.msm.fix.wide <- tidyr::pivot_wider(dat.msm.fix, id_cols = id, names_from = outcome_char, values_from = c(time, status, x1, x2))
# 
# ## Remove duplicate predictors variables and rename
# dat.msm.fix.wide <- dat.msm.fix.wide %>% 
#   dplyr::select(id, time_A, time_B, status_A, status_B, x1_A, x2_A) %>%
#   rename(x1 = x1_A, x2 = x2_A)
# 
# str(dat.msm.fix.wide)
# sum(dat.msm.fix.wide$time_A < 10)/nrow(dat.msm.fix.wide)
# sum(dat.msm.fix.wide$time_B < 10)/nrow(dat.msm.fix.wide)
# sum(dat.msm.fix.wide$time_A < 10 & dat.msm.fix.wide$time_B < 10)/nrow(dat.msm.fix.wide)
# sum(dat.msm.fix$time < 10 & dat.msm.fix$outcome == 2)/sum(dat.msm.fix$outcome == 2)
# 
# true.risk.msm
# pred.risk.msm
# 
# 
# save.image("data/TEST_MSM.RData")
