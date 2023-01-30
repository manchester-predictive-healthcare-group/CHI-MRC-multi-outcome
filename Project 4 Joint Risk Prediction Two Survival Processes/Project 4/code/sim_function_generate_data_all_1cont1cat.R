#######################################################################################
### This workspace will contain all the functions for generating data in the simulation
### DGM's are multistate model, copula, and frailty approches
### For each DGM, data is reformatted into a common data model
#######################################################################################

### Clear workspace
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Load packages
library(survival)
library(gems)
library(mstate)
library(dplyr)
library(tidyr)
library(simsurv)
library(frailtypack)
library(CopulaCenR)
library(copula)


##################
##################
### MSM ###
##################
##################

gen.dat.msm <- function(n, #number of patients to simulate
                        max.follow, #maximum follow up
                        shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                        shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                        shape24, scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                        shape34, scale34, #shape and scale for weibull baseline hazard for transition 3 -> 4
                        beta12.cont, beta12.cat, #covariate effects for transiion 12
                        beta13.cont, beta13.cat, #covariate effects for transiion 13
                        beta24.cont, beta24.cat, #covariate effects for transiion 24
                        beta34.cont, beta34.cat, #covariate effects for transiion 34
                        baseline_cens, COV_beta_cens, #weibull baseline parameters, and covariate effects for censoring mechanism
                        x.in, #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
                        numsteps) #number of sampler steps in gems data generation process
{
  
  #     n.sim <- 50
  #     x.baseline <- data.frame("x1" = rnorm(n.sim, 0, 1), "x2" = rbinom(n.sim, 1, 0.5))
  #     n <- n.sim
  #     shape12 <- 1.3
  #     scale12 <- 1
  #     shape13 <- 0.8
  #     scale13 <- 1
  #     shape24 <- 0.8
  #     scale24 <- 1.5
  #     shape34 <- 1.3
  #     scale34 <- 1.5
  #     beta12.cont <- 0
  #     beta13.cont <- 0
  #     beta24.cont <- 0
  #     beta34.cont <- 0
  #     beta12.cat <- 0.5
  #     beta13.cat <- 0.25
  #     beta24.cat <- 0.25
  #     beta34.cat <- 0.5
  #     x.in <- x.baseline
  #     numsteps <- 1000
  
  
  ## Generate a baseline covariate data frame
  bl <- x.in
  
  ## Generate an empty hazard matrix
  hf <- generateHazardMatrix(4)
  hf
  
  ## Change the entries of the transitions we want to allow
  ## Define the transitions as weibull
  hf[[1, 2]] <- function(t, shape, scale, beta.cont, beta.cat) {
    exp(bl["x1"]*beta.cont + bl["x2"]*beta.cat)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[1, 3]] <- function(t, shape, scale, beta.cont, beta.cat) {
    exp(bl["x1"]*beta.cont + bl["x2"]*beta.cat)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[2, 4]] <- function(t, shape, scale, beta.cont, beta.cat) {
    exp(bl["x1"]*beta.cont + bl["x2"]*beta.cat)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[3, 4]] <- function(t, shape, scale, beta.cont, beta.cat) {
    exp(bl["x1"]*beta.cont + bl["x2"]*beta.cat)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf
  
  
  ## Note that by using (t + sum(history)) we are implementing a clock forward approach
  
  
  ## Generate an empty parameter matrix
  par <- generateParameterMatrix(hf)
  
  ## Use the vector of scales in each transition hazard
  par[[1, 2]] <- list(shape = shape12, scale = scale12, 
                      beta.cont = beta12.cont, beta.cat = beta12.cat)
  par[[1, 3]] <- list(shape = shape13, scale = scale13, 
                      beta.cont = beta13.cont, beta.cat = beta13.cat)
  par[[2, 4]] <- list(shape = shape24, scale = scale24, 
                      beta.cont = beta24.cont, beta.cat = beta24.cat)
  par[[3, 4]] <- list(shape = shape34, scale = scale34, 
                      beta.cont = beta34.cont, beta.cat = beta34.cat)
  
  
  ## Generate the cohort
  time.in <- Sys.time()
  cohort <- simulateCohort(transitionFunctions = hf, parameters = par,
                           cohortSize = n, baseline = bl, to = 10000, sampler.steps = numsteps)
  time.out <- Sys.time()
  time.diff <- time.out - time.in
  
  
  ## Get data into the common data model
  
  ## Extract event times
  dat.temp <- cohort@time.to.state
  
  ## Condition A is the minimum of time until state 2 or state 4 is entered
  ## Condition B is the minimum of time until state 3 or state 4 is entered
  ## To take the minimum, we need to assign the NA's to something, if we set the NA's to be equal to state 4, then it will not affected the min function
  dat.temp <- data.frame(t(apply(dat.temp, 1, function(x) {ifelse(is.na(x), max(x, na.rm = TRUE), x)})))
  
  
  ## Take the min
  dat.temp$t.A <- pmin(dat.temp$"State.2", dat.temp$"State.4")
  dat.temp$t.B <- pmin(dat.temp$"State.3", dat.temp$"State.4")
  
  ## Create the censoring indicator
  dat.temp$cens.A <- rep(1, nrow(dat.temp))
  dat.temp$cens.B <- rep(1, nrow(dat.temp))
  
  ## Add an id variable
  dat.temp$id <- 1:nrow(dat.temp)
  
  ## Add the predictor variables to the output matrix
  dat.temp <- cbind(dat.temp, x.in)
  
  ## Remove unnecary variables
  dat.temp <- dplyr::select(dat.temp, t.A, t.B, cens.A, cens.B, id, x1, x2)
  
  ## Merge the columns into a single outcome variable in long format
  ## Do two seperate pivot's, to get the outcome variable, and the censoring variable both in long format
  ## (should be able to do it in one go? but can't figure out code)
  dat.temp.time <- dat.temp %>% 
    pivot_longer(cols = c("t.A", "t.B"), names_to = c("time.ind"), values_to = "time") %>%
    dplyr::select(time.ind, time, id, x1, x2)
  head(dat.temp.time)
  
  dat.temp.cens <- dat.temp %>% 
    pivot_longer(cols = c("cens.A", "cens.B"), names_to = c("cens.ind"), values_to = "status") %>%
    dplyr::select(status)
  head(dat.temp.cens)
  
  ## Combine the two into one dataset
  dat.temp <- cbind(dat.temp.time, dat.temp.cens)
  rm(dat.temp.time, dat.temp.cens)
  head(dat.temp)
  
  ## Create a variable for the margin indicator
  dat.temp$outcome <- rep(0, nrow(dat.temp))
  dat.temp$outcome[dat.temp$time.ind == "t.A"] <- 1
  dat.temp$outcome[dat.temp$time.ind == "t.B"] <- 2
  
  ## Create a factor version, and a character version (character version will contain variable names for when changing to wide format)
  dat.temp$outcome_fac <- as.factor(dat.temp$outcome)
  dat.temp$outcome_char <- dat.temp$outcome
  dat.temp$outcome_char[dat.temp$outcome == 1] <- "A"
  dat.temp$outcome_char[dat.temp$outcome == 2] <- "B"
  
  ## Drope the time.ind variable
  dat.temp <- dplyr::select(dat.temp, -time.ind)
  
  ## Order the variables 
  dat.temp <- dat.temp[,c("id","x1","x2","time","status","outcome","outcome_fac","outcome_char")]
  
  ## Turn x2 into a factor variables
  dat.temp$x2 <- as.factor(dat.temp$x2)
  
  ###
  ### Create and add the censoring indicator
  ###
  
  ## Generate the censoring times for each individual
  cens.times <- simsurv("weibull", lambdas = 1/baseline_cens[2], gammas = baseline_cens[1], x = x.in, 
                        betas = c("x1" = COV_beta_cens[1], "x2" = COV_beta_cens[2]))
  
  ## Create two rows for each censoring time to add to the dataset
  cens.times <- rbind(cens.times, cens.times)
  cens.times <- arrange(cens.times, id)
  
  ## Add the censoring time for each individual to the dataset
  dat.temp$cens.time <- cens.times$eventtime
  
  ## Change the censoring indicator to be 0 if the event happens after the censoring time
  dat.temp$status[dat.temp$cens.time < dat.temp$time] <- 0
  
  ## If the event is censored, change the event time to be equal to the censoring time
  dat.temp$time[dat.temp$status == 0] <- dat.temp$cens.time[dat.temp$status == 0]
  
  ## Finally, censor all observations that happen after some maximum follow up time
  dat.temp <- dat.temp %>% 
    mutate(time = case_when(time > max.follow ~ max.follow, TRUE ~ as.numeric(as.character(time))))
  
  dat.temp <- dat.temp %>% 
    mutate(status = case_when(time == max.follow ~ 0, TRUE ~ as.numeric(as.character(status))))
  
  
  ###
  ### Finally make the minimum possible event time (censored or not) to be 1 day = 1/365.
  ### This should not alter results, is realistic, and will stop us getting errors when event time too close to zero
  ###
  dat.temp$time[dat.temp$time < 1/365] <- 1/365
  
  ## Return dataset
  return(dat.temp)
}



##################
##################
### COPULA ###
##################
##################

gen.dat.copula <- function (n, # number of observations
                            max.follow, #maximum follow up
                            copula, copula.param, rotate.cop, #type of copula, association parameter, and whether copula should be rotated
                            baselineA, COV_betaA, #outcome1: baselineA = shape and scale of weibull distribution, COV_betaA = covariate effects
                            baselineB, COV_betaB, #outcome2: baselineB = shape and scale of weibull distribution, COV_betaB = covariate effects
                            baseline_cens, COV_beta_cens, #weibull baseline parameters, and covariate effects for censoring mechanism
                            x.in) #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
{

#   n <- n.sim
#    copula <- "clayton"
#   copula.param <- 3
#    baselineA <- c(1,15)
#   COV_betaA <- c(0.5, 0.5)
#    baselineB <- c(1,10)
#   COV_betaB <- c(0.25, 0.75)
#    COV_beta_cens <- c(0,0)
# baseline_cens = c(1,20)
#    x.in <- x.baseline
  
  ## First turn input dataframe into a matrix, which is what is required for this function to work
  ## (save a copy of x as a dataframe, which is required for generating the censoring indicator)
  x.data.frame <- x.in
  x.in <- as.matrix(x.in)
  
  ## Define var_list (variable names)
  var_list <- c("x1", "x2")
  
  ## Assign distribution of baseline hazards
  dist1 <- "Weibull"
  dist2 <- "Weibull"
  
  ## The rest of this code is edited from the CopulaCenR::data_sim_copula function
  cl <- archmCopula(family = tolower(copula), param = copula.param, 
                    dim = 2)
  if (rotate.cop == 90){cl <- rotCopula(cl, flip = c(TRUE, FALSE))}
  if (rotate.cop == 180){cl <- rotCopula(cl, flip = c(TRUE, TRUE))}
  
  Cop <- rCopula(n, cl)
  u <- Cop[, 1]
  v <- Cop[, 2]
  if (dist1 != "Gompertz") {
    k1 <- baselineA[1]
    lambda1 <- baselineA[2]
  }
  if (dist1 == "Gompertz") {
    a1 <- baselineA[1]
    b1 <- baselineA[2]
  }
  if (dist2 != "Gompertz") {
    k2 <- baselineB[1]
    lambda2 <- baselineB[2]
  }
  if (dist2 == "Gompertz") {
    a2 <- baselineB[1]
    b2 <- baselineB[2]
  }
  COV_betaA <- matrix(COV_betaA, ncol = 1)
  COV_betaB <- matrix(COV_betaB, ncol = 1)
  colnames(x.in) = var_list
  dat1 <- data.frame(id = seq(1:n), outcome = rep(1, n), x.in)
  dat2 <- data.frame(id = seq(1:n), outcome = rep(2, n), x.in)
  if (dist1 == "Weibull") {
    t1 <- as.vector(((-log(u))/(exp(matrix(unlist(dat1[,var_list]),
                                           ncol = length(var_list)) %*% COV_betaA)))^(1/k1)*lambda1)
    t2 <- as.vector(((-log(v))/(exp(matrix(unlist(dat2[, var_list]), 
                                           ncol = length(var_list)) %*% COV_betaB)))^(1/k2)*lambda2)
  }
  if (dist1 == "Loglogistic") {
    t1 <- as.vector(((-1 + 1/u)/(exp(matrix(unlist(dat1[, var_list]), 
                                            ncol = length(var_list)) %*% COV_betaA)))^(1/k1)*lambda1)
    t2 <- as.vector(((-1 + 1/v)/(exp(matrix(unlist(dat2[, var_list]), 
                                            ncol = length(var_list)) %*% COV_betaB)))^(1/k2)*lambda2)
  }
  if (dist1 == "Gompertz") {
    t1 <- as.vector(1/a1 * log(1 - a1 * log(u) * exp(-1 * matrix(unlist(dat1[, var_list]), 
                                                                 ncol = length(var_list)) %*% COV_betaA)/b1))
    t2 <- as.vector(1/a2 * log(1 - a2 * log(v) * exp(-1 * matrix(unlist(dat2[, var_list]), 
                                                                 ncol = length(var_list)) %*% COV_betaB)/b2))
  }
  dat1 <- cbind(dat1, data.frame(time = t1))
  dat2 <- cbind(dat2, data.frame(time = t2))
  
  dat <- rbind(dat1, dat2)
  dat <- dat[order(dat$id), ]
  
  ## Create a factor version, and a character version (character version will contain variable names for when changing to wide format)
  dat$outcome_fac <- as.factor(dat$outcome)
  dat$outcome_char <- dat$outcome
  dat$outcome_char[dat$outcome == 1] <- "A"
  dat$outcome_char[dat$outcome == 2] <- "B"
  
  ## Add a status indicator for all individuals
  dat$status <- rep(1, nrow(dat))
  
  ## Remove rownames
  rownames(dat) <- NULL
  
  ## Order columns
  dat <- dat[,c("id","x1","x2","time","status","outcome","outcome_fac","outcome_char")]
  
  ## Turn x2 into a factor variables
  dat$x2 <- as.factor(dat$x2)
  
  ###
  ### Create and add the censoring indicator
  ###
  
  ## Generate the censoring times for each individual
  cens.times <- simsurv("weibull", lambdas = 1/baseline_cens[2], gammas = baseline_cens[1], x = x.data.frame, 
                        betas = c("x1" = COV_beta_cens[1], "x2" = COV_beta_cens[2]))
  
  ## Create two rows for each censoring time to add to the dataset
  cens.times <- rbind(cens.times, cens.times)
  cens.times <- arrange(cens.times, id)
  
  ## Add the censoring time for each individual to the dataset
  dat$cens.time <- cens.times$eventtime
  
  ## Change the censoring indicator to be 0 if the event happens after the censoring time
  dat$status[dat$cens.time < dat$time] <- 0
  
  ## If the event is censored, change the event time to be equal to the censoring time
  dat$time[dat$status == 0] <- dat$cens.time[dat$status == 0]
  
  ## Finally, censor all observations that happen after some maximum follow up time
  dat <- dat %>% 
    mutate(time = case_when(time > max.follow ~ max.follow, TRUE ~ as.numeric(as.character(time))))
  
  dat <- dat %>% 
    mutate(status = case_when(time == max.follow ~ 0, TRUE ~ as.numeric(as.character(status))))
  
  
  ###
  ### Finally make the minimum possible event time (censored or not) to be 1 day = 1/365.
  ### This should not alter results, is realistic, and will stop us getting errors when event time too close to zero
  ###
  dat$time[dat$time < 1/365] <- 1/365
  
  ## Return dataset
  return(dat)
}


##################
##################
### FRAILTY ###
##################
##################

### NB: Variance for the gamma is applied on the scale of the hazard, 
### whereas we take the exponent of the normal random effect and apply it on the exp() scale

gen.dat.frailty <- function(n, #number of observations to simulate
                            max.follow, #maximum follow up
                            frail.dist, frail.eff, frail.var, #frailty distribution, frail.eff = 0/1 to remove/include frailty term, variance of random effect, 
                            baselineA, COV_betaA, #1) shape and scale of weibull baseline hazard, 2) covariate effects
                            baselineB, COV_betaB, #1) shape and scale of weibull baseline hazard, 2) covariate effects
                            baseline_cens, COV_beta_cens, #weibull baseline parameters, and covariate effects for censoring mechanism
                            x.in) #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
{
  
  #   n = n.sim
  #   var = 1/3
  #   baselineA = c(1,2)
  #   COV_betaA = c(0.5, 1.5)
  #   baselineB = c(1,0.5)
  #   COV_betaB = c(0.5, 1.5)
  #   x.in = x.baseline
  
  ## Simulate random effects
  if (frail.dist == "normal"){rfrail <- rnorm(n, 0, frail.var)}
  
  if (frail.dist == "gamma"){rfrail <- log(rgamma(n, shape = 1/frail.var, scale = frail.var))}
  
  if (!(frail.dist %in% c("normal", "gamma"))){print("frail.type mut be either normal or gamma")}
  
  stopifnot(frail.dist %in% c("normal", "gamma"))
  
  ## Add these to baseline characteritics/input data
  x.in$rfrail <-rfrail
  
  ## Use simsurv to generate survival times
  o.dat.frail.outcomes.1 <- simsurv("weibull", lambdas = 1/baselineA[2], gammas = baselineA[1], x = x.in, 
                                    betas = c("x1" = COV_betaA[1], "x2" = COV_betaA[2], "rfrail" = frail.eff))
  o.dat.frail.outcomes.2 <- simsurv("weibull", lambdas = 1/baselineB[2], gammas = baselineB[1], x = x.in, 
                                    betas = c("x1" = COV_betaB[1], "x2" = COV_betaB[2], "rfrail" = frail.eff))
  
  ## Combine survival times with baseline data
  o.dat.frail.1 <- cbind(o.dat.frail.outcomes.1, x.in)
  o.dat.frail.2 <- cbind(o.dat.frail.outcomes.2, x.in)
  
  ## Add a variable to indicate which outcome the row corresponds to
  o.dat.frail.1$outcome <- rep(1, nrow(o.dat.frail.1))
  o.dat.frail.2$outcome <- rep(2, nrow(o.dat.frail.2))
  
  ## Combine the two outcomes
  o.dat.frail <- rbind(o.dat.frail.1, o.dat.frail.2)
  colnames(o.dat.frail)[colnames(o.dat.frail) == "eventtime"] <- "time"
  
  ## Arrange by patient
  o.dat.frail <- arrange(o.dat.frail, id)
  
  ## Create a factor version, and a character version (character version will contain variable names for when changing to wide format)
  o.dat.frail$outcome_fac <- as.factor(o.dat.frail$outcome)
  o.dat.frail$outcome_char <- o.dat.frail$outcome
  o.dat.frail$outcome_char[o.dat.frail$outcome == 1] <- "A"
  o.dat.frail$outcome_char[o.dat.frail$outcome == 2] <- "B"
  
  ## Order columns
  o.dat.frail <- o.dat.frail[,c("id","x1","x2","rfrail","time","status","outcome","outcome_fac","outcome_char")]
  
  ## Turn x2 into a factor variables
  o.dat.frail$x2 <- as.factor(o.dat.frail$x2)
  
  ###
  ### Create and add the censoring indicator
  ###
  
  ## Generate the censoring times for each individual
  cens.times <- simsurv("weibull", lambdas = 1/baseline_cens[2], gammas = baseline_cens[1], x = x.in, 
                        betas = c("x1" = COV_beta_cens[1], "x2" = COV_beta_cens[2]))
  
  ## Create two rows for each censoring time to add to the dataset
  cens.times <- rbind(cens.times, cens.times)
  cens.times <- arrange(cens.times, id)
  
  ## Add the censoring time for each individual to the dataset
  o.dat.frail$cens.time <- cens.times$eventtime
  
  ## Change the censoring indicator to be 0 if the event happens after the censoring time
  o.dat.frail$status[o.dat.frail$cens.time < o.dat.frail$time] <- 0
  
  ## If the event is censored, change the event time to be equal to the censoring time
  o.dat.frail$time[o.dat.frail$status == 0] <- o.dat.frail$cens.time[o.dat.frail$status == 0]
  
  ## Finally, censor all observations that happen after some maximum follow up time
  o.dat.frail <- o.dat.frail %>% 
    mutate(time = case_when(time > max.follow ~ max.follow, TRUE ~ as.numeric(as.character(time))))
  
  o.dat.frail <- o.dat.frail %>% 
    mutate(status = case_when(time == max.follow ~ 0, TRUE ~ as.numeric(as.character(status))))
  
  
  ###
  ### Finally make the minimum possible event time (censored or not) to be 1 day = 1/365.
  ### This should not alter results, is realistic, and will stop us getting errors when event time too close to zero
  ###
  o.dat.frail$time[o.dat.frail$time < 1/365] <- 1/365
  
  ## Return dataset
  return(o.dat.frail)
}


##############
### Save image
##############
save.image("data/sim_function_generate_data_all_1cont1cat.RData")


# ###
# ### Generate some small datasets to look at and check all under common data model
# ###
# 
# ### Define the number of individuals we want to generate
# n.sim <- 50
# 
# ### Create the dataset of baseline characteristics
# x.baseline <- data.frame("x1" = rnorm(n.sim, 0, 1), "x2" = rbinom(n.sim, 1, 0.5))
# 
# 
# ### Frailty
# dat.frailty.50 <- gen.dat.frailty(n = n.sim, 
#                                frail.dist = "gamma", frail.eff = 1, frail.var = 0.001,
#                                baselineA = c(1,1.5), COV_betaA = c(0, 0), 
#                                baselineB = c(1,1), COV_betaB = c(0, 0), 
#                                x.in = x.baseline)
# 
# 
# ### Copula
# dat.copula.50 <- gen.dat.copula(n = n.sim, 
#                              copula = "clayton", copula.param = 3,
#                              baselineA = c(1,1.5), COV_betaA = c(0, 0), 
#                              baselineB = c(1,1), COV_betaB = c(0, 0), 
#                              x.in = x.baseline)
# 
# 
# ### MSM
# dat.msm.50 <- gen.dat.msm(n = n.sim,
#                           shape12 = 1, scale12 = 1.5, #shape and scale for weibull baseline hazard for transition 1 -> 2
#                           shape13 = 1, scale13 = 1, #shape and scale for weibull baseline hazard for transition 1 -> 3
#                           shape24 = 1, scale24 = 1, #shape and scale for weibull baseline hazard for transition 2 -> 4
#                           shape34 = 1, scale34 = 1.5, #shape and scale for weibull baseline hazard for transition 3 -> 4
#                           beta12.cont = 0, beta12.cat = 0, #covariate effects for transiion 12
#                           beta13.cont = 0, beta13.cat = 0, #covariate effects for transiion 13
#                           beta24.cont = 0, beta24.cat = 0, #covariate effects for transiion 24
#                           beta34.cont = 0, beta34.cat = 0, #covariate effects for transiion 34
#                           x.in = x.baseline,
#                           numsteps = 50000)
# 
# 
# ###
# ### Generate some test data with no association, and see if the outcome prevalence matches output we expect from weibull
# ###
# 
# ### Define the number of individuals we want to generate
# n.sim <- 10000
# 
# ### Create the dataset of baseline characteristics
# x.baseline <- data.frame("x1" = rnorm(n.sim, 0, 1), "x2" = rbinom(n.sim, 1, 0.5))
# 
# 
# ### Frailty
# dat.frailty <- gen.dat.frailty(n = n.sim, 
#                                frail.dist = "gamma", frail.eff = 1, frail.var = 0.001,
#                                baselineA = c(1,1.5), COV_betaA = c(0, 0), 
#                                baselineB = c(1,1), COV_betaB = c(0, 0), 
#                                x.in = x.baseline)
# 
# mean(dat.frailty$time[dat.frailty$outcome == 1])
# mean(dat.frailty$time[dat.frailty$outcome == 2])
# 
# 
# ### Copula
# dat.copula <- gen.dat.copula(n = n.sim, 
#                              copula = "clayton", copula.param = 3,
#                              baselineA = c(1,1.5), COV_betaA = c(0, 0), 
#                              baselineB = c(1,1), COV_betaB = c(0, 0), 
#                              x.in = x.baseline)
# 
# mean(dat.copula$time[dat.copula$outcome == 1])
# mean(dat.copula$time[dat.copula$outcome == 2])
# 
# 
# ### MSM
# dat.msm <- gen.dat.msm(n = n.sim,
#                        shape12 = 1, scale12 = 1.5, #shape and scale for weibull baseline hazard for transition 1 -> 2
#                        shape13 = 1, scale13 = 1, #shape and scale for weibull baseline hazard for transition 1 -> 3
#                        shape24 = 1, scale24 = 1, #shape and scale for weibull baseline hazard for transition 2 -> 4
#                        shape34 = 1, scale34 = 1.5, #shape and scale for weibull baseline hazard for transition 3 -> 4
#                        beta12.cont = 0, beta12.cat = 0, #covariate effects for transiion 12
#                        beta13.cont = 0, beta13.cat = 0, #covariate effects for transiion 13
#                        beta24.cont = 0, beta24.cat = 0, #covariate effects for transiion 24
#                        beta34.cont = 0, beta34.cat = 0, #covariate effects for transiion 34
#                        x.in = x.baseline,
#                        numsteps = 50000)
# 
# mean(dat.msm$time[dat.msm$outcome == 1])
# mean(dat.msm$time[dat.msm$outcome == 2])

