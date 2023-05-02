###################################################################################################
###################################################################################################
###
### Functions in section 3 are  functions to calculate observed risks/apply validation
###
###################################################################################################
###################################################################################################

###
### 3.1) Calculate Aalen-Johansen estimator for a cohort at time t.eval
###
calc.calib.aj.ce <- function(data.mstate, tmat, t.eval){
#   str(data.mstate)
#   events(data.mstate)
#   data.mstate <- data.mstate.deciles[[i]][[j]]
#   tmat <- tmat 
#   t.eval <- t.eval
#   events(data.mstate)
  ### Assign max state number
  max.state <- max(data.mstate$to)
  
  ### Fit csh's with no predictors
  csh.aj <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data.mstate)
  
  ### Calculate cumulative incidence functions
  msfit.aj <- msfit(csh.aj, trans = tmat)
  
  ### Calculate Aalen-Johansen estimator
  pt.aj <- probtrans(msfit.aj, predt = 0)
  
  ### Extract the closest time in the data to the time we want to evaluate at
  t.eval.dat <- pt.aj[[1]]$time[max(which(pt.aj[[1]]$time <= t.eval))]
  
  ### Extract AJ estimator at this time point
  obs.aj <- pt.aj[[1]][pt.aj[[1]]$time == t.eval.dat, paste("pstate", 1:max.state, sep = "")]
  
  ### Extract AJ standard error  at this time point
  obs.aj.se <- pt.aj[[1]][pt.aj[[1]]$time == t.eval.dat, paste("se", 1:max.state, sep = "")]
  
  ### Create output object
  output.object <- list("obs.aj" = obs.aj, "obs.aj.se" = obs.aj.se)
  
  return(output.object)
}


###
### 3.2) Calculate landmark Aalen-Johansen estimator for a cohort at time t.eval
###
calc.calib.lmaj.ce <- function(data.mstate, t.eval){
  
  #   t.eval <- ceiling(7*365.25)
  #     data.mstate <- data.mstate.obj[["data.mstate"]]
  #     tmat <- data.mstate.obj[["tmat"]]
  ### Assign max state number
  max.state <- max(data.mstate$to)
  
  ### Calculate observed risks using landmark Aalen-Johansen
  pt.lmaj <- LMAJ(msdata = data.mstate, s = 0, from = 1)
  
  ### Extract the closest time in the data to the time we want to evaluate at
  t.eval.dat <- pt.lmaj$time[max(which(pt.lmaj$time <= t.eval))]
  
  ### Get the LMAJ estimator at this time point
  obs.lmaj <- pt.lmaj[pt.lmaj$time == t.eval.dat, paste("pstate", 1:max.state, sep = "")]
  
  return(obs.lmaj)
}


###
### 3.3) Function to extract all individuals in state j, at time t.eval, from a dataset in mstate format
### Used in other functions which assess calibration using blr/mlr at specific time points (3.4 and 3.5)
### 
extract.ids.states.ce <- function(data.mstate, j, t.eval){
  
  ### For non-absorbing states, to be in state j at time t, you must have an observations from state j, where Tstart <= t.eval < Tstop
  if (j < max(data.mstate$to)){
    ## Extract ids
    ids.state.j <- subset(data.mstate, from == j & Tstart <= t.eval & t.eval < Tstop) %>%
      select(person_id) %>%
      distinct(person_id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$person_id)
  } else if (j == max(data.mstate$to)){
    ### For absorbind state, just have to have moved into it
    ids.state.j <- subset(data.mstate, to == max(data.mstate$to) & t.eval >= Tstop & status == 1) %>%
      select(person_id) %>%
      distinct(person_id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$person_id)
  }
  
  return(ids.state.j)
}

###
### EXPLORE
###
# data.raw.noNA %>%
#   group_by(state2.bin) %>%
#   dplyr::summarize(mean = mean(pstate2))
# colnames(data.raw)
# temp.data <- data.raw[, c(1:20, 39)]
# table(data.raw.noNA$state.poly)
# msm.cox.fit


###
### 3.4A) Assess calibration using binary logistic regression
### Note this does not apply ipcw, and therefore will be incorrect in the presence of censoring
###


###
### 3.4B) Assess calibration using binary logistic regression and IPCW weights
###
calc.calib.blr.ipcw.ce <- function(data.mstate, data.raw, t.eval, p.est, max.weight = 10){
  
  ## Let's do state 1 first
  ## Want to know which individual are in state j at time t
#   data.mstate <- data.mstate, 
#   data.raw <- data.raw, 
#   t.eval <- t.eval, 
#   p.est <- select(data.raw, paste("pstate", 1:6, sep = ""))
#   max.weight <- 10
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states.ce(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #         intersect(ids.state.list[[1]], ids.state.list[[3]])
  #       intersect(ids.state.list[[1]], ids.state.list[[4]])
  #       intersect(ids.state.list[[1]], ids.state.list[[5]])
  #       intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(person_id %in% ids.state.list[[1]] ~ 1,
                                  person_id %in% ids.state.list[[2]] ~ 2,
                                  person_id %in% ids.state.list[[3]] ~ 3,
                                  person_id %in% ids.state.list[[4]] ~ 4,
                                  person_id %in% ids.state.list[[5]] ~ 5,
                                  person_id %in% ids.state.list[[6]] ~ 6),
           state1.bin = case_when(state.poly == 1 ~ 1,
                                  state.poly != 1 ~ 0),
           state2.bin = case_when(state.poly == 2 ~ 1,
                                  state.poly != 2 ~ 0),
           state3.bin = case_when(state.poly == 3 ~ 1,
                                  state.poly != 3 ~ 0),
           state4.bin = case_when(state.poly == 4 ~ 1,
                                  state.poly != 4 ~ 0),
           state5.bin = case_when(state.poly == 5 ~ 1,
                                  state.poly != 5 ~ 0),
           state6.bin = case_when(state.poly == 6 ~ 1,
                                  state.poly != 6 ~ 0))
  
  ### Add the predicted risks, and the logit transormation of the predicted risks to the dataset
  p.est.logit <- log(p.est/(1-p.est))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est, p.est.logit)
  
  ### This data could be used to fit invalid logistic regression calibration models
  
  ###
  ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
  ###
  weights <- calc.weights.ce(data.raw, t.eval, max.weight = max.weight)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights)
  
  ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
  ### (to ensure weights/offsets and individuals being modelled are consistent)
  data.raw.noNA <- data.raw %>% subset(!is.na(state.poly))
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1 <- glm(state1.bin ~ p.est.logit1, weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm2 <- glm(state2.bin ~ p.est.logit2, weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm3 <- glm(state3.bin ~ p.est.logit3, weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm4 <- glm(state4.bin ~ p.est.logit4, weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm5 <- glm(state5.bin ~ p.est.logit5, weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm6 <- glm(state6.bin ~ p.est.logit6, weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  
  slopes <- c("state1" = as.numeric(coefficients(lrm1)[2]),
              "state2" = as.numeric(coefficients(lrm2)[2]),
              "state3" = as.numeric(coefficients(lrm3)[2]),
              "state4" = as.numeric(coefficients(lrm4)[2]),
              "state5" = as.numeric(coefficients(lrm5)[2]),
              "state6" = as.numeric(coefficients(lrm6)[2]))
  
  slopes.se <- c("state1" = as.numeric(summary(lrm1)$coefficients[2,2]),
                 "state2" = as.numeric(summary(lrm2)$coefficients[2,2]),
                 "state3" = as.numeric(summary(lrm3)$coefficients[2,2]),
                 "state4" = as.numeric(summary(lrm4)$coefficients[2,2]),
                 "state5" = as.numeric(summary(lrm5)$coefficients[2,2]),
                 "state6" = as.numeric(summary(lrm6)$coefficients[2,2]))
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1.offset <- glm(state1.bin ~ offset(p.est.logit1), weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm2.offset <- glm(state2.bin ~ offset(p.est.logit2), weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm3.offset <- glm(state3.bin ~ offset(p.est.logit3), weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm4.offset <- glm(state4.bin ~ offset(p.est.logit4), weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm5.offset <- glm(state5.bin ~ offset(p.est.logit5), weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm6.offset <- glm(state6.bin ~ offset(p.est.logit6), weights = ipcw, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  
  intercepts <- c("state1" = as.numeric(coefficients(lrm1.offset)),
                  "state2" = as.numeric(coefficients(lrm2.offset)),
                  "state3" = as.numeric(coefficients(lrm3.offset)),
                  "state4" = as.numeric(coefficients(lrm4.offset)),
                  "state5" = as.numeric(coefficients(lrm5.offset)),
                  "state6" = as.numeric(coefficients(lrm6.offset)))
  
  intercepts.se <- c("state1" = as.numeric(summary(lrm1.offset)$coefficients[2]),
                     "state2" = as.numeric(summary(lrm2.offset)$coefficients[2]),
                     "state3" = as.numeric(summary(lrm3.offset)$coefficients[2]),
                     "state4" = as.numeric(summary(lrm4.offset)$coefficients[2]),
                     "state5" = as.numeric(summary(lrm5.offset)$coefficients[2]),
                     "state6" = as.numeric(summary(lrm6.offset)$coefficients[2]))
  
  ###
  ### Create predicted observed risk for each individual
  ###
  pred.obs1 <- predict(lrm1.offset, newdata = data.raw.noNA, type = "response")
  pred.obs2 <- predict(lrm2.offset, newdata = data.raw.noNA, type = "response")
  pred.obs3 <- predict(lrm3.offset, newdata = data.raw.noNA, type = "response")
  pred.obs4 <- predict(lrm4.offset, newdata = data.raw.noNA, type = "response")
  pred.obs5 <- predict(lrm5.offset, newdata = data.raw.noNA, type = "response")
  pred.obs6 <- predict(lrm6.offset, newdata = data.raw.noNA, type = "response")
  
  ### Get difference in risk
  diff.pred.obs1 <- mean(pred.obs1 - data.raw.noNA$p.est1)
  diff.pred.obs2 <- mean(pred.obs2 - data.raw.noNA$p.est2)
  diff.pred.obs3 <- mean(pred.obs3 - data.raw.noNA$p.est3)
  diff.pred.obs4 <- mean(pred.obs4 - data.raw.noNA$p.est4)
  diff.pred.obs5 <- mean(pred.obs5 - data.raw.noNA$p.est5)
  diff.pred.obs6 <- mean(pred.obs6 - data.raw.noNA$p.est6)
  diff.pred.obs <- c(diff.pred.obs1, diff.pred.obs2, diff.pred.obs3, diff.pred.obs4, diff.pred.obs5, diff.pred.obs6)
  
  ### Create output object
  output.obj <- list("int" = intercepts, "int.se" = intercepts.se,
                     "slopes" = slopes, "slopes.se" = slopes.se,
                     "diff.pred.obs" = diff.pred.obs)
  
  return(output.obj)
  
}



###
### 3.4C) Assess calibration using binary logistic regression and IPCW weights
### Standard errors calculated through bootstrapping
###
calc.calib.blr.ipcw.boot.ce <- function(data.mstate, data.raw, t.eval, p.est, max.weight = 10, n.boot = 200){
  
  ## Let's do state 1 first
  ## Want to know which individual are in state j at time t
#   data.mstate = data.mstate
#   data.raw = data.raw
#   t.eval <- t.eval
#   p.est <- select(data.raw, paste("pstate", 1:6, sep = ""))
#   max.weight <- 10
#   n.boot <- 5
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states.ce(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #         intersect(ids.state.list[[1]], ids.state.list[[3]])
  #       intersect(ids.state.list[[1]], ids.state.list[[4]])
  #       intersect(ids.state.list[[1]], ids.state.list[[5]])
  #       intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  ###
  
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(person_id %in% ids.state.list[[1]] ~ 1,
                                  person_id %in% ids.state.list[[2]] ~ 2,
                                  person_id %in% ids.state.list[[3]] ~ 3,
                                  person_id %in% ids.state.list[[4]] ~ 4,
                                  person_id %in% ids.state.list[[5]] ~ 5,
                                  person_id %in% ids.state.list[[6]] ~ 6),
           state1.bin = case_when(state.poly == 1 ~ 1,
                                  state.poly != 1 ~ 0),
           state2.bin = case_when(state.poly == 2 ~ 1,
                                  state.poly != 2 ~ 0),
           state3.bin = case_when(state.poly == 3 ~ 1,
                                  state.poly != 3 ~ 0),
           state4.bin = case_when(state.poly == 4 ~ 1,
                                  state.poly != 4 ~ 0),
           state5.bin = case_when(state.poly == 5 ~ 1,
                                  state.poly != 5 ~ 0),
           state6.bin = case_when(state.poly == 6 ~ 1,
                                  state.poly != 6 ~ 0))
  
  ### Add the predicted risks, and the logit transormation of the predicted risks to the dataset
  p.est.logit <- log(p.est/(1-p.est))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est, p.est.logit)
  
  ### This data could be used to fit invalid logistic regression calibration models
  
  ### Now write a function that will be put into a boot function to calculate bootstrapped standard errors for intercept and slopes
  get_int_slope_boot <- function(data, indices){
    
    ### Create bootstrap sample based on indices
    data.boot <- data[indices, ]
    
    ###
    ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
    ###
    weights <- calc.weights.ce(data.boot, t.eval, max.weight = max.weight)
    
    ### Add these to dataset
    data.boot <- cbind(data.boot, weights)
    
    ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
    ### (to ensure weights/offsets and individuals being modelled are consistent)
    data.boot.noNA <- data.boot %>% subset(!is.na(state.poly))
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1 <- glm(state1.bin ~ p.est.logit1, weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm2 <- glm(state2.bin ~ p.est.logit2, weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm3 <- glm(state3.bin ~ p.est.logit3, weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm4 <- glm(state4.bin ~ p.est.logit4, weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm5 <- glm(state5.bin ~ p.est.logit5, weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm6 <- glm(state6.bin ~ p.est.logit6, weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    
    slopes <- c("state1" = as.numeric(coefficients(lrm1)[2]),
                "state2" = as.numeric(coefficients(lrm2)[2]),
                "state3" = as.numeric(coefficients(lrm3)[2]),
                "state4" = as.numeric(coefficients(lrm4)[2]),
                "state5" = as.numeric(coefficients(lrm5)[2]),
                "state6" = as.numeric(coefficients(lrm6)[2]))
    
    slopes.se <- c("state1" = as.numeric(summary(lrm1)$coefficients[2,2]),
                   "state2" = as.numeric(summary(lrm2)$coefficients[2,2]),
                   "state3" = as.numeric(summary(lrm3)$coefficients[2,2]),
                   "state4" = as.numeric(summary(lrm4)$coefficients[2,2]),
                   "state5" = as.numeric(summary(lrm5)$coefficients[2,2]),
                   "state6" = as.numeric(summary(lrm6)$coefficients[2,2]))
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1.offset <- glm(state1.bin ~ offset(p.est.logit1), weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm2.offset <- glm(state2.bin ~ offset(p.est.logit2), weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm3.offset <- glm(state3.bin ~ offset(p.est.logit3), weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm4.offset <- glm(state4.bin ~ offset(p.est.logit4), weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm5.offset <- glm(state5.bin ~ offset(p.est.logit5), weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm6.offset <- glm(state6.bin ~ offset(p.est.logit6), weights = ipcw, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    
    intercepts <- c("state1" = as.numeric(coefficients(lrm1.offset)),
                    "state2" = as.numeric(coefficients(lrm2.offset)),
                    "state3" = as.numeric(coefficients(lrm3.offset)),
                    "state4" = as.numeric(coefficients(lrm4.offset)),
                    "state5" = as.numeric(coefficients(lrm5.offset)),
                    "state6" = as.numeric(coefficients(lrm6.offset)))
    
    intercepts.se <- c("state1" = as.numeric(summary(lrm1.offset)$coefficients[2]),
                       "state2" = as.numeric(summary(lrm2.offset)$coefficients[2]),
                       "state3" = as.numeric(summary(lrm3.offset)$coefficients[2]),
                       "state4" = as.numeric(summary(lrm4.offset)$coefficients[2]),
                       "state5" = as.numeric(summary(lrm5.offset)$coefficients[2]),
                       "state6" = as.numeric(summary(lrm6.offset)$coefficients[2]))
    
    ###
    ### Create predicted observed risk for each individual
    ###
    pred.obs1 <- predict(lrm1.offset, newdata = data.boot.noNA, type = "response")
    pred.obs2 <- predict(lrm2.offset, newdata = data.boot.noNA, type = "response")
    pred.obs3 <- predict(lrm3.offset, newdata = data.boot.noNA, type = "response")
    pred.obs4 <- predict(lrm4.offset, newdata = data.boot.noNA, type = "response")
    pred.obs5 <- predict(lrm5.offset, newdata = data.boot.noNA, type = "response")
    pred.obs6 <- predict(lrm6.offset, newdata = data.boot.noNA, type = "response")
    
    ### Get difference in risk
    diff.pred.obs1 <- mean(pred.obs1 - data.boot.noNA$p.est1)
    diff.pred.obs2 <- mean(pred.obs2 - data.boot.noNA$p.est2)
    diff.pred.obs3 <- mean(pred.obs3 - data.boot.noNA$p.est3)
    diff.pred.obs4 <- mean(pred.obs4 - data.boot.noNA$p.est4)
    diff.pred.obs5 <- mean(pred.obs5 - data.boot.noNA$p.est5)
    diff.pred.obs6 <- mean(pred.obs6 - data.boot.noNA$p.est6)
    diff.pred.obs <- c(diff.pred.obs1, diff.pred.obs2, diff.pred.obs3, diff.pred.obs4, diff.pred.obs5, diff.pred.obs6)
    
    ### Create output object
    output.obj <- c(intercepts, slopes, diff.pred.obs)
    
    return(output.obj)
    
  }
  
  ### Run the bootstrapping
  boot.obj <- boot(data.raw, statistic = get_int_slope_boot, R = n.boot)
  colnames(boot.obj$t) <- c(paste("int", 1:6, sep = ""), 
                            paste("slope", 1:6, sep = ""), 
                            paste("diff.pred.obs", 1:6, sep = ""))  
  
  ### Calculate standard errors
  se <- sqrt(apply(boot.obj$t, 2, var))
  names(se) <- c(paste("int", 1:6, sep = ""), 
                 paste("slope", 1:6, sep = ""), 
                 paste("diff.pred.obs", 1:6, sep = ""))  
  
  ### Return output
  output.object <- list("boot.obj" = boot.obj$t, "se" = se)
  return(output.object)
  
}


###
### 3.5A) Function to calculate calibration after applying the nominal calibration framework of van Hoorde et al,
### No weights
###


###
### 3.5B) Function to calculate calibration after applying the nominal calibration framework of van Hoorde et al,
### weighted using IPCW weights
###
calc.calib.mlr.ipcw.ce <- function(data.mstate, data.raw, t.eval, p.est, max.weight = 10){
  
  #   ## Let's do state 1 first
  #   ## Want to know which individual are in state j at time t
  #   data.mstate <- data.mstate.reduc
  #   data.raw <- data.raw.reduc
  #   t.eval <- ceiling(7*365.25)
  #   p.est <- p.true
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states.ce(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #         intersect(ids.state.list[[1]], ids.state.list[[3]])
  #       intersect(ids.state.list[[1]], ids.state.list[[4]])
  #       intersect(ids.state.list[[1]], ids.state.list[[5]])
  #       intersect(ids.state.list[[4]], ids.state.list[[5]])
  

  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(person_id %in% ids.state.list[[1]] ~ 1,
                                  person_id %in% ids.state.list[[2]] ~ 2,
                                  person_id %in% ids.state.list[[3]] ~ 3,
                                  person_id %in% ids.state.list[[4]] ~ 4,
                                  person_id %in% ids.state.list[[5]] ~ 5,
                                  person_id %in% ids.state.list[[6]] ~ 6))
  
  
  ### Add p.est to dataset
  data.raw <- data.frame(data.raw, p.est)
  
  ### Add linear predictors from a multinomial framework
  data.raw <- data.raw %>%
    mutate(mlr.lp1 = log(p.est2/p.est1),
           mlr.lp2 = log(p.est3/p.est1),
           mlr.lp3 = log(p.est4/p.est1),
           mlr.lp4 = log(p.est5/p.est1),
           mlr.lp5 = log(p.est6/p.est1),
           state.poly.fac = as.factor(state.poly))
  
  ###
  ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
  ###
  weights <- calc.weights.ce(data.raw, t.eval, max.weight = max.weight)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights)
  
  ###
  ### Now fit the nominal calibration framework using the weights
  ###
  
  ## Add constraints
  i <- diag(5)
  i1 <- rbind(1, 0, 0, 0, 0)
  i2 <- rbind(0, 1, 0, 0, 0)
  i3 <- rbind(0, 0, 1, 0, 0)
  i4 <- rbind(0, 0, 0, 1, 0)
  i5 <- rbind(0, 0, 0, 0, 1)
  clist <- list("(Intercept)" = i, "mlr.lp1" = i1, "mlr.lp2" = i2, "mlr.lp3" = i3, "mlr.lp4" = i4, "mlr.lp5" = i5)
  clist
  
  ### Apply nominal recalibration framework
  ###
  
  ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
  ### (to ensure weights and individuals being modelled are consistent)
  data.raw.noNA <- data.raw %>% subset(!is.na(state.poly))
  
  ### Fit models
  calib.model <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4 + mlr.lp5, constraints = clist, weights = ipcw,
                            data = data.raw.noNA, family = multinomial(refLevel = "1"))
  
  calib.model.offset <- vgam(data.raw.noNA$state.poly.fac ~ 1, offset = as.matrix(data.raw.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4", "mlr.lp5")]),
                                   weights = data.raw.noNA$ipcw,
                                   family = multinomial(refLevel = "1"))

  
  ###
  ### Create predicted observed risk for each individual
  ###
  pred.obs <- data.frame("mlr.lp1" = coefficients(calib.model.offset)[1] + data.raw.noNA$mlr.lp1,
                               "mlr.lp2" = coefficients(calib.model.offset)[2] + data.raw.noNA$mlr.lp2,
                               "mlr.lp3" = coefficients(calib.model.offset)[3] + data.raw.noNA$mlr.lp3,
                               "mlr.lp4" = coefficients(calib.model.offset)[4] + data.raw.noNA$mlr.lp4,
                               "mlr.lp5" = coefficients(calib.model.offset)[5] + data.raw.noNA$mlr.lp5)
  pred.obs <- mutate(pred.obs, 
                           p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                           p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                           p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                           p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                           p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                           p6 = exp(mlr.lp5)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)))
  
  ### Get difference in risk
  diff.pred.obs1 <- mean(pred.obs[,"p1"] - data.raw.noNA$p.est1)
  diff.pred.obs2 <- mean(pred.obs[,"p2"] - data.raw.noNA$p.est2)
  diff.pred.obs3 <- mean(pred.obs[,"p3"] - data.raw.noNA$p.est3)
  diff.pred.obs4 <- mean(pred.obs[,"p4"] - data.raw.noNA$p.est4)
  diff.pred.obs5 <- mean(pred.obs[,"p5"] - data.raw.noNA$p.est5)
  diff.pred.obs6 <- mean(pred.obs[,"p6"] - data.raw.noNA$p.est6)
  diff.pred.obs <- c(diff.pred.obs1, diff.pred.obs2, diff.pred.obs3, 
                     diff.pred.obs4, diff.pred.obs5, diff.pred.obs6)
  
  
  ### Create output object
  output.object <- list("int" = calib.model.offset@coefficients, 
                        "int.se" = sqrt(diag(vcov(calib.model.offset))), 
                        "slopes" = calib.model@coefficients[paste("mlr.lp", 1:5, sep = "")], 
                        "slopes.se" = sqrt(diag(vcov(calib.model))[paste("mlr.lp", 1:5, sep = "")]),
                        "diff.pred.obs" = diff.pred.obs)
  
  return(output.object)
  
}



###
### 3.5C) Function to calculate calibration after applying the nominal calibration framework of van Hoorde et al,
### weighted using IPCW weights.
### Calculate standard errors using bootstrapping
###
calc.calib.mlr.ipcw.boot.ce <- function(data.mstate, data.raw, t.eval, p.est, max.weight = 10, n.boot = 200){
  
  #   ## Let's do state 1 first
  #   ## Want to know which individual are in state j at time t
  #   data.mstate <- data.mstate.reduc
  #   data.raw <- data.raw.reduc
  #   t.eval <- ceiling(7*365.25)
  #   p.est <- p.true
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states.ce(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #         intersect(ids.state.list[[1]], ids.state.list[[3]])
  #       intersect(ids.state.list[[1]], ids.state.list[[4]])
  #       intersect(ids.state.list[[1]], ids.state.list[[5]])
  #       intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(person_id %in% ids.state.list[[1]] ~ 1,
                                  person_id %in% ids.state.list[[2]] ~ 2,
                                  person_id %in% ids.state.list[[3]] ~ 3,
                                  person_id %in% ids.state.list[[4]] ~ 4,
                                  person_id %in% ids.state.list[[5]] ~ 5,
                                  person_id %in% ids.state.list[[6]] ~ 6))
  
  
  ### Add p.est to dataset
  data.raw <- data.frame(data.raw, p.est)
  
  ### Add linear predictors from a multinomial framework
  data.raw <- data.raw %>%
    mutate(mlr.lp1 = log(p.est2/p.est1),
           mlr.lp2 = log(p.est3/p.est1),
           mlr.lp3 = log(p.est4/p.est1),
           mlr.lp4 = log(p.est5/p.est1),
           mlr.lp5 = log(p.est6/p.est1),
           state.poly.fac = as.factor(state.poly))
  
  ### Now write a function that will be put into a boot function to calculate bootstrapped standard errors for intercept and slopes
  get_int_slope_boot <- function(data, indices){
    
    ### Create bootstrap sample based on indices
    data.boot <- data[indices, ]
    
    ###
    ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
    ###
    weights <- calc.weights.ce(data.boot, t.eval, max.weight = max.weight)
    
    ### Add these to dataset
    data.boot <- cbind(data.boot, weights)
    
    ###
    ### Now fit the nominal calibration framework using the weights
    ###
    
    ## Add constraints
    i <- diag(5)
    i1 <- rbind(1, 0, 0, 0, 0)
    i2 <- rbind(0, 1, 0, 0, 0)
    i3 <- rbind(0, 0, 1, 0, 0)
    i4 <- rbind(0, 0, 0, 1, 0)
    i5 <- rbind(0, 0, 0, 0, 1)
    clist <- list("(Intercept)" = i, "mlr.lp1" = i1, "mlr.lp2" = i2, "mlr.lp3" = i3, "mlr.lp4" = i4, "mlr.lp5" = i5)
    clist
    
    ### Apply nominal recalibration framework
    ###
    
    ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
    ### (to ensure weights and individuals being modelled are consistent)
    data.boot.noNA <- data.boot %>% subset(!is.na(state.poly))
    
    ### Fit models
    calib.model <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4 + mlr.lp5, constraints = clist, weights = ipcw,
                        data = data.boot.noNA, family = multinomial(refLevel = "1"))
    
    calib.model.offset <- vgam(data.boot.noNA$state.poly.fac ~ 1, offset = as.matrix(data.boot.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4", "mlr.lp5")]),
                               weights = data.boot.noNA$ipcw,
                               family = multinomial(refLevel = "1"))
    
    
    ###
    ### Create predicted observed risk for each individual
    ###
    pred.obs <- data.frame("mlr.lp1" = coefficients(calib.model.offset)[1] + data.boot.noNA$mlr.lp1,
                           "mlr.lp2" = coefficients(calib.model.offset)[2] + data.boot.noNA$mlr.lp2,
                           "mlr.lp3" = coefficients(calib.model.offset)[3] + data.boot.noNA$mlr.lp3,
                           "mlr.lp4" = coefficients(calib.model.offset)[4] + data.boot.noNA$mlr.lp4,
                           "mlr.lp5" = coefficients(calib.model.offset)[5] + data.boot.noNA$mlr.lp5)
    
    pred.obs <- mutate(pred.obs, 
                       p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                       p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                       p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                       p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                       p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)),
                       p6 = exp(mlr.lp5)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4) + exp(mlr.lp5)))
    
    ### Get difference in risk
    diff.pred.obs1 <- mean(pred.obs[,"p1"] - data.boot.noNA$p.est1)
    diff.pred.obs2 <- mean(pred.obs[,"p2"] - data.boot.noNA$p.est2)
    diff.pred.obs3 <- mean(pred.obs[,"p3"] - data.boot.noNA$p.est3)
    diff.pred.obs4 <- mean(pred.obs[,"p4"] - data.boot.noNA$p.est4)
    diff.pred.obs5 <- mean(pred.obs[,"p5"] - data.boot.noNA$p.est5)
    diff.pred.obs6 <- mean(pred.obs[,"p6"] - data.boot.noNA$p.est6)
    diff.pred.obs <- c(diff.pred.obs1, diff.pred.obs2, diff.pred.obs3, 
                       diff.pred.obs4, diff.pred.obs5, diff.pred.obs6)
    
    
    ### Create output object
    output.object <- c(calib.model.offset@coefficients, calib.model@coefficients[paste("mlr.lp", 1:5, sep = "")], diff.pred.obs)
    
    return(output.object)
  }
  
  ### Run the bootstrapping
  boot.obj <- boot(data.raw, statistic = get_int_slope_boot, R = n.boot)
  colnames(boot.obj$t) <- c(paste("int", 2:6, sep = ""), 
                            paste("slope", 2:6, sep = ""),
                            paste("diff.pred.obs", 1:6, sep = ""))
  
  ### Calculate standard errors
  se <- sqrt(apply(boot.obj$t, 2, var))
  names(se) <- c(paste("int", 2:6, sep = ""), 
                 paste("slope", 2:6, sep = ""), 
                 paste("diff.pred.obs", 1:6, sep = ""))
  
  ### Return output
  output.object <- list("boot.obj" = boot.obj$t, "se" = se)
  return(output.object)
  
}



###
### 3.6A) Define a function to calculate pseudo-value for an individual, using the Aalen-Johansen estimator
### 

### obs.aj is pre-calculated Aalen-Johansen estimator in entire cohort
func.calc.pv.aj.ce <- function(person_id.eval, data.mstate, obs.aj, tmat, n.cohort, t.eval){
  
  ### Calculate AJ estimate without patient in dataset
  est.drop.pat <- calc.calib.aj.ce(subset(data.mstate, person_id != person_id.eval), 
                                tmat = tmat, t.eval = t.eval)
  
  ### Retain just the estimate (not the standard error)
  est.drop.pat <- est.drop.pat[["obs.aj"]]
  
  ### Calculate the pseudo-value
  pv.pat <- n.cohort*obs.aj - (n.cohort-1)*est.drop.pat
  pv.pat <- as.numeric(pv.pat[1, ])
  
  return(pv.pat)
}



###
### 3.6B) Define function to calculate pseudo value, but there is an extra argument containing 
### the pseudo-value for someone who has an event after time t.eval.
### This pseudo value will be the same for all patients who meet these criteria, so no need to recalculate
### it everytime.
###
### obs.aj is pre-calculated Aalen-Johansen estimator in entire cohort data.mstate
### pv.same is the pre-calculated pseudo-valus
func.calc.pv.aj.efficient.ce <- function(person_id.eval, data.mstate, obs.aj, tmat, n.cohort, t.eval, pv.same){
  
  ### Create temporary dataset of the first row of the data of the individuals of interest
  data.mstate.temp <- subset(data.mstate, person_id == person_id.eval & from == 1 & to == 2)
  
  if (data.mstate.temp$Tstop > t.eval){
    ### Assign pv.pat to be the pre-calculated value of pseudo value
    pv.pat <- pv.same
    
  } else if (data.mstate.temp$Tstop <= t.eval){
    ### Calculate pseudo-value for individual
    pv.pat <- func.calc.pv.aj.ce(person_id.eval = person_id.eval, 
                                 data.mstate = data.mstate, 
                                 obs.aj = obs.aj, 
                                 tmat = tmat, 
                                 n.cohort = n.cohort, 
                                 t.eval = t.eval)
  }
  
  return(pv.pat)
}



###
### 4.1B) Assess moderate calibration using binary logistic regression and IPCW
###
calc.calib.blr.ipcw.mod.ce <- function(data.mstate, data.raw, t.eval, p.est, max.weight = 10){
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states.ce(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #         intersect(ids.state.list[[1]], ids.state.list[[3]])
  #       intersect(ids.state.list[[1]], ids.state.list[[4]])
  #       intersect(ids.state.list[[1]], ids.state.list[[5]])
  #       intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  ### 
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(person_id %in% ids.state.list[[1]] ~ 1,
                                  person_id %in% ids.state.list[[2]] ~ 2,
                                  person_id %in% ids.state.list[[3]] ~ 3,
                                  person_id %in% ids.state.list[[4]] ~ 4,
                                  person_id %in% ids.state.list[[5]] ~ 5,
                                  person_id %in% ids.state.list[[6]] ~ 6),
           state1.bin = case_when(state.poly == 1 ~ 1,
                                  state.poly != 1 ~ 0),
           state2.bin = case_when(state.poly == 2 ~ 1,
                                  state.poly != 2 ~ 0),
           state3.bin = case_when(state.poly == 3 ~ 1,
                                  state.poly != 3 ~ 0),
           state4.bin = case_when(state.poly == 4 ~ 1,
                                  state.poly != 4 ~ 0),
           state5.bin = case_when(state.poly == 5 ~ 1,
                                  state.poly != 5 ~ 0),
           state6.bin = case_when(state.poly == 6 ~ 1,
                                  state.poly != 6 ~ 0),
           state.poly.fac = as.factor(state.poly))
  
  ### Add the predicted risks, and the logit transormation of the predicted risks to the dataset
  p.est.logit <- log(p.est/(1-p.est))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est, p.est.logit)
  
  ###
  ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
  ###
  weights <- calc.weights.ce(data.raw, t.eval, max.weight = max.weight)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights)
  
  ###
  ### Fit the models
  ###
  
  ### Fit loess recalibration models to the uncensored observations at time t to calculate intercepts
  loess1 <- loess(state1.bin ~ p.est1, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw"])
  loess2 <- loess(state2.bin ~ p.est2, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw"])
  loess3 <- loess(state3.bin ~ p.est3, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw"])
  loess4 <- loess(state4.bin ~ p.est4, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw"])
  loess5 <- loess(state5.bin ~ p.est5, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw"])
  loess6 <- loess(state6.bin ~ p.est6, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw"])
  
  ### Created 'predicted observed' probabilities for each individual
  data.raw$loess.pred.obs1 <- predict(loess1, newdata = data.raw)
  data.raw$loess.pred.obs2 <- predict(loess2, newdata = data.raw)
  data.raw$loess.pred.obs3 <- predict(loess3, newdata = data.raw)
  data.raw$loess.pred.obs4 <- predict(loess4, newdata = data.raw)
  data.raw$loess.pred.obs5 <- predict(loess5, newdata = data.raw)
  data.raw$loess.pred.obs6 <- predict(loess6, newdata = data.raw)
  
  ### Produce plots for each and store in a list
  
  ### Create plot titles that will be used
  plot.titles <- c("Healthy", "CVD", "T2D", "CKD", "Multimorbidity", "Death")
  plots.list <- vector("list", 6)
  for (i in 1:6){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("loess.pred.obs", i, sep = "")]
    
    ### Create the plots
    plots.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>% 
                                      arrange(pred) %>% select(person_id, pred, obs)) +
      geom_line(aes(x = pred, y = obs), color = "red") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, " (", plot.titles[i], ")", sep = ""))
  }
  
  ### Calculate the ECI
  ECI <- calc.ECI.ce(data.raw, p.est, data.raw[, paste("loess.pred.obs", 1:6, sep = "")])

  
  ### Return output object
  output.object <- list("plots.list" = plots.list,
                        "ECI" = ECI)
  return(output.object)
  
}


###
### 4.2B) Assess moderate calibration using multinomial nominal recalibration framework and IPCW
###
calc.calib.mlr.ipcw.mod.ce <- function(data.mstate, data.raw, t.eval, p.est, ps.int = 4, degree = 3, max.weight = 10){
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states.ce(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #         intersect(ids.state.list[[1]], ids.state.list[[3]])
  #       intersect(ids.state.list[[1]], ids.state.list[[4]])
  #       intersect(ids.state.list[[1]], ids.state.list[[5]])
  #       intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  ### 
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(person_id %in% ids.state.list[[1]] ~ 1,
                                  person_id %in% ids.state.list[[2]] ~ 2,
                                  person_id %in% ids.state.list[[3]] ~ 3,
                                  person_id %in% ids.state.list[[4]] ~ 4,
                                  person_id %in% ids.state.list[[5]] ~ 5,
                                  person_id %in% ids.state.list[[6]] ~ 6),
           state.poly.fac = as.factor(state.poly))
  
  ### Add p.est to dataset
  data.raw <- data.frame(data.raw, p.est)
  
  ### Add linear predictors from a multinomial framework
  data.raw <- data.raw %>%
    mutate(mlr.lp1 = log(p.est2/p.est1),
           mlr.lp2 = log(p.est3/p.est1),
           mlr.lp3 = log(p.est4/p.est1),
           mlr.lp4 = log(p.est5/p.est1),
           mlr.lp5 = log(p.est6/p.est1))
  
  ###
  ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
  ###
  weights <- calc.weights.ce(data.raw, t.eval, max.weight = max.weight)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights)

  ###
  ### Miss-specified weighted model
  ###
  
  ### Apply nominal recalibration framework with vector spline smoothers
  calib.model <- vgam(state.poly.fac ~ sm.ps(mlr.lp1, ps.int = ps.int, degree = degree) + 
                        sm.ps(mlr.lp2, ps.int = ps.int, degree = degree) + 
                        sm.ps(mlr.lp3, ps.int = ps.int, degree = degree) + 
                        sm.ps(mlr.lp4, ps.int = ps.int, degree = degree) + 
                        sm.ps(mlr.lp5, ps.int = ps.int, degree = degree), weights = data.raw[!is.na(data.raw$state.poly), "ipcw"],
                            data = data.raw[!is.na(data.raw$state.poly), ], family = multinomial(refLevel = "1"))
  
  
  ###
  ### Generate predicted-observed risks and add to data.raw
  ###
  
  ### For all other functions, I just generate predicted risks for all individuals, and then just plot for those who were uncensored
  ### However, some of the censored individuals are causing an error, therefore I must generate NA vectors, then assign the
  ### predicted observed probabilities to the crrect individuals
  
  ## Create dataframe to store
  dat.mlr.pred.obs <- data.frame(matrix(NA, ncol = 6, nrow = nrow(data.raw)))
  ## Assign colnames
  colnames(dat.mlr.pred.obs) <- paste("mlr.pred.obs", 1:ncol(dat.mlr.pred.obs), sep = "")
  ## Calc pred.obs for those who are uncesored
  mlr.pred.obs <- predict(calib.model, newdata = data.raw[!is.na(data.raw$state.poly), ], type = "response")
  ## Assign to appropriate individuals
  dat.mlr.pred.obs[!is.na(data.raw$state.poly), ] <- mlr.pred.obs
  
  ### Then add it to data.raw
  data.raw <- cbind(data.raw, dat.mlr.pred.obs)
  
  ### Create plots
  
  ### Create plot titles that will be used
  plot.titles <- c("Healthy", "CVD", "T2D", "CKD", "Multimorbidity", "Death")
  plots.list <- vector("list", 6)
  for (i in 1:6){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("mlr.pred.obs", i, sep = "")]
    
    ### Create the plots
    plots.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>%
                                      arrange(pred) %>%  select(person_id, pred, obs)) +
      geom_point(aes(x = pred, y = obs), color = "red", size = 0.5, alpha = 0.05) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, " (", plot.titles[i], ")", sep = ""))
  }
  
  ### Calculate the ECI
  ECI <- calc.ECI.ce(data.raw, p.est, dat.mlr.pred.obs)
  
  ### Return output object
  output.object <- list("plots.list" = plots.list, 
                        "ECI" = ECI)
  return(output.object)
  
}


###
### 4.3) Calculate calibration using pseudo values
###
calc.calib.pv.moderate.ce <- function(data.raw, pv.comb, p.est){
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Combine data.raw and predicted probabilities
  data.pv <- cbind(pv.comb, p.est)
  
  ### Fit loess recalibration models to the uncensored observations at time t to calculate intercepts
  loess1 <- loess(pv.state1 ~ p.est1, data = data.pv)
  loess2 <- loess(pv.state2 ~ p.est2, data = data.pv)
  loess3 <- loess(pv.state3 ~ p.est3, data = data.pv)
  loess4 <- loess(pv.state4 ~ p.est4, data = data.pv)
  loess5 <- loess(pv.state5 ~ p.est5, data = data.pv)
  loess6 <- loess(pv.state6 ~ p.est6, data = data.pv)
  
  ### Created 'predicted observed' probabilities for each individual
  data.pv$loess.pred.obs1 <- predict(loess1, newdata = data.pv)
  data.pv$loess.pred.obs2 <- predict(loess2, newdata = data.pv)
  data.pv$loess.pred.obs3 <- predict(loess3, newdata = data.pv)
  data.pv$loess.pred.obs4 <- predict(loess4, newdata = data.pv)
  data.pv$loess.pred.obs5 <- predict(loess5, newdata = data.pv)
  data.pv$loess.pred.obs6 <- predict(loess6, newdata = data.pv)
  
  ### Calculate the ECI
  ECI <- calc.ECI.ce(data.raw, p.est, data.pv[, paste("loess.pred.obs", 1:6, sep = "")])
  
  ### Create plot titles that will be used
  plot.titles <- c("Healthy", "CVD", "T2D", "CKD", "Multimorbidity", "Death")
  plots.list <- vector("list", 6)
  for (i in 1:6){
    
    ### Creaet variables to plot 
    data.pv$pred <- data.pv[, paste("p.est", i, sep = "")]
    data.pv$obs <- data.pv[, paste("loess.pred.obs", i, sep = "")]
    
    ### Create the plots
    plots.list[[i]] <- ggplot(data = data.pv %>% 
                                arrange(pred) %>% select(pred, obs), 
                              aes(x = pred, y = obs, color = "red")) +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.pv$pred))) + 
      ylim(c(0, max(data.pv$pred))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, " (", plot.titles[i], ")", sep = ""))
  }
  
  
  ### Create output object and return it
  output.object <- list("plots.list" = plots.list, "ECI" = ECI)
  return(output.object)
  
}



###
### Misc 1) Function to calculate IPC weights estimating the model from the data
###
calc.weights.ce <- function(data.raw, t.eval, max.weight = 10){
  
  ### First need to estimate the censoring mechanism from the data
  
  ### We add a censoring time that would be observed from the data (this would be the min of censoring time, and death)
  ### We take the min of these as the event time
  ### It's an "event" if censoring happens, and "censored" if death happens (as we want to estimate the rate of censoring)
    data.raw <- data.raw %>%
      mutate(dtcens_var.s = case_when(c6 == 1 ~ 0,
                                      c6 == 0 ~ 1))
  
  ###
  ### Create models for censoring in order to calculate the IPCW weights
  ###
  
  ### A model where we do not adjust for predictors covariates (i.e. miss-specified ipcw weights)
  cens.model <- coxph(as.formula(paste("Surv(dtcens_var, dtcens_var.s) ~ ", 
                                       paste(covs, collapse = "+"), 
                                       sep = "")), 
                      data = data.raw)
  
  ### Calculate a data frame containing probability of censored and uncenosred at each time point
  ### The weights will be the probability of being uncensored, at the time of the event for each individual

  ## Extract baseline hazard
  data.weights <- basehaz(cens.model, centered = TRUE)
  ## Add lp to data.raw
  data.raw$lp <- predict(cens.model, type = "lp", reference = "sample")
  
  
  ### We do not convert hazards to probabilities yet (like for mspec), as weights are different for each individual, so we do this
  ### within prob.uncens.func.pspec below
  
  ### Create weights for the cohort at time t.eval
  ### Note for individuals who died, we take the probability of them being uncensored at the time of death
  ### For individuals still alive, we take the probability of being uncensored at time t.eval
  
  ### If individual has death prior to the time we are evaluating, assign weight at time of death
  ### Get location of individuals who had death prior to evaluation time (or censoring)
  obs.death.prior <- !is.na(data.raw$state.poly) & data.raw$o6 < t.eval
  
  ###
  ### Now create probability of censoring weights
  ###
  
  ### First assign all individuals a weight of the probability of being uncesored at time t.eval
  ### This is the linear predictor times the cumulative hazard at time t.eval, and appropriate transformation to get a risk
  data.raw$pcw <- as.numeric(exp(-exp(data.raw$lp)*data.weights$hazard[max(which(data.weights$time <= t.eval))]))
  
  ## Write a function which will extract the uncensored probability for an individual at a given time
  prob.uncens.func <- function(input){
    ## Assign t and person_id
    t <- input[1]
    lp <- input[2]
    
    ## Get hazard at appropriate time
    bhaz.t <- data.weights$hazard[max(which(data.weights$time <= t))]
    
    ## Return risk
    return(exp(-exp(lp)*bhaz.t))
  }
  
  ## Apply this function to all the times at which individuals have died prior to censoring, and assign to the appropriate individuals
  data.raw$pcw[obs.death.prior] <- apply(data.raw[obs.death.prior, c("o6", "lp")], 1, FUN = prob.uncens.func)
 
  ### Invert these
  data.raw$ipcw <- 1/data.raw$pcw
  
  ### Finally cap these at 10
  data.raw$ipcw <- pmin(data.raw$ipcw, max.weight)
  
  ### Create output object
  output.object <- data.frame("ipcw" = data.raw$ipcw, "pcw" = data.raw$pcw)
  
  return(output.object)
  
}


###
### Misc 2) Function to calculate IPC weights estimating the model from the data
###
calc.weights.stab.ce <- function(data.raw, t.eval, max.weight = 10){
  
  ### First need to estimate the censoring mechanism from the data
  
  ### We add a censoring time that would be observed from the data (this would be the min of censoring time, and death)
  ### We take the min of these as the event time
  ### It's an "event" if censoring happens, and "censored" if death happens (as we want to estimate the rate of censoring)
  data.raw <- data.raw %>%
    mutate(dtcens_var.s = case_when(c6 == 1 ~ 0,
                                    c6 == 0 ~ 1))
  
  ###
  ### Create models for censoring in order to calculate the IPCW weights
  ###
  
  ### Fully adjusted model to estimate probability of being censored
  cens.model <- coxph(as.formula(paste("Surv(dtcens_var, dtcens_var.s) ~ ", 
                                       paste(covs, collapse = "+"), 
                                       sep = "")), 
                      data = data.raw)
  
  ### Intercept only model to estimate marginal risk
  cens.model.int <- coxph(as.formula(paste("Surv(dtcens_var, dtcens_var.s) ~ 1", 
                                           sep = "")), 
                          data = data.raw)
  
  ### Calculate a data frame containing probability of censored and uncenosred at each time point
  ### The weights will be the probability of being uncensored, at the time of the event for each individual
  
  ## Extract baseline hazard
  data.weights <- basehaz(cens.model, centered = TRUE)
  data.weights.int <- basehaz(cens.model.int, centered = TRUE)
  
  ## Add lp to data.raw
  data.raw$lp <- predict(cens.model, type = "lp", reference = "sample")
  
  
  ### We do not convert hazards to probabilities yet (like for mspec), as weights are different for each individual, so we do this
  ### within prob.uncens.func.pspec below
  
  ### Create weights for the cohort at time t.eval
  ### Note for individuals who died, we take the probability of them being uncensored at the time of death
  ### For individuals still alive, we take the probability of being uncensored at time t.eval
  
  ### If individual has death prior to the time we are evaluating, assign weight at time of death
  ### Get location of individuals who had death prior to evaluation time (or censoring)
  obs.death.prior <- !is.na(data.raw$state.poly) & data.raw$o6 < t.eval
  
  ###
  ### Now create probability of censoring weights
  ###
  
  ### First assign all individuals a weight of the probability of being uncesored at time t.eval
  ### This is the linear predictor times the cumulative hazard at time t.eval, and appropriate transformation to get a risk
  data.raw$pcw <- as.numeric(exp(-exp(data.raw$lp)*data.weights$hazard[max(which(data.weights$time <= t.eval))]))
  
  ## Write a function which will extract the uncensored probability for an individual at a given time
  prob.uncens.func <- function(input){
    ## Assign t and person_id
    t <- input[1]
    lp <- input[2]
    
    ## Get hazard at appropriate time
    bhaz.t <- data.weights$hazard[max(which(data.weights$time <= t))]
    
    ## Return risk
    return(exp(-exp(lp)*bhaz.t))
  }
  
  ## Apply this function to all the times at which individuals have died prior to censoring, and assign to the appropriate individuals
  data.raw$pcw[obs.death.prior] <- apply(data.raw[obs.death.prior, c("o6", "lp")], 1, FUN = prob.uncens.func)
  
  ### Invert these
  data.raw$ipcw <- 1/data.raw$pcw
  
  ### Finally calculate numerator, marginal probability of being uncensored at time t.eval
  data.raw$ipcw.numer <- as.numeric(exp(-data.weights.int$hazard[max(which(data.weights.int$time <= t.eval))]))
  
  ### Calculate stabilised weights
  data.raw$ipcw.stab <- data.raw$ipcw.numer*data.raw$ipcw
  
  ### Finally cap these at 10
  data.raw$ipcw <- pmin(data.raw$ipcw, max.weight)
  data.raw$ipcw.stab <- pmin(data.raw$ipcw.stab, max.weight)
  
  ### Create output object
  output.object <- data.frame("ipcw" = data.raw$ipcw, "ipcw.stab" = data.raw$ipcw.stab, "pcw" = data.raw$pcw)
  
  return(output.object)
  
}

###
### Misc 2) A program to calculate ECI
###
calc.ECI.ce <- function(data.raw.in, p.est.in, p.obs.in){
  
  ### Reduvce to individuals who are uncensored
  p.est.in <- p.est.in[!is.na(data.raw.in$state.poly.fac), ]
  p.obs.in <-p.obs.in[!is.na(data.raw.in$state.poly.fac), ]
  data.raw.in <- data.raw.in[!is.na(data.raw.in$state.poly.fac), ]
  
  ### Calculate estimated calibration index
  # Calc event rates
  K1 <- sum(data.raw.in$state.poly.fac == 1)/nrow(data.raw.in)
  K2 <- sum(data.raw.in$state.poly.fac == 2)/nrow(data.raw.in)
  K3 <- sum(data.raw.in$state.poly.fac == 3)/nrow(data.raw.in)
  K4 <- sum(data.raw.in$state.poly.fac == 4)/nrow(data.raw.in)
  K5 <- sum(data.raw.in$state.poly.fac == 5)/nrow(data.raw.in)
  K6 <- sum(data.raw.in$state.poly.fac == 6)/nrow(data.raw.in)
  
  # Calc numerator
  ECI.numer <- (sum((p.est.in[,1] - p.obs.in[,1])^2) + sum((p.est.in[,2] - p.obs.in[,2])^2) +
                  sum((p.est.in[,3] - p.obs.in[,3])^2) + sum((p.est.in[,4] - p.obs.in[,4])^2) +
                  sum((p.est.in[,5] - p.obs.in[,5])^2) + sum((p.est.in[,6] - p.obs.in[,6])^2))
  # Calc denominator
  ECI.denom <- (sum((p.est.in[,1] - K1)^2) + sum((p.est.in[,2] - K2)^2) +
                  sum((p.est.in[,3] - K3)^2) + sum((p.est.in[,4] - K4)^2) +
                  sum((p.est.in[,5] - K5)^2) + sum((p.est.in[,6] - K6)^2))
  # Calc ECI
  ECI <- ECI.numer/ECI.denom
  
  return(ECI)
}


###
### Misc 4) Same as 4.1B), but using stabilised weights
###
calc.calib.blr.ipcw.stab.mod.ce <- function(data.mstate, data.raw, t.eval, p.est, max.weight = 10){
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states.ce(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #         intersect(ids.state.list[[1]], ids.state.list[[3]])
  #       intersect(ids.state.list[[1]], ids.state.list[[4]])
  #       intersect(ids.state.list[[1]], ids.state.list[[5]])
  #       intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  ### 
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(person_id %in% ids.state.list[[1]] ~ 1,
                                  person_id %in% ids.state.list[[2]] ~ 2,
                                  person_id %in% ids.state.list[[3]] ~ 3,
                                  person_id %in% ids.state.list[[4]] ~ 4,
                                  person_id %in% ids.state.list[[5]] ~ 5,
                                  person_id %in% ids.state.list[[6]] ~ 6),
           state1.bin = case_when(state.poly == 1 ~ 1,
                                  state.poly != 1 ~ 0),
           state2.bin = case_when(state.poly == 2 ~ 1,
                                  state.poly != 2 ~ 0),
           state3.bin = case_when(state.poly == 3 ~ 1,
                                  state.poly != 3 ~ 0),
           state4.bin = case_when(state.poly == 4 ~ 1,
                                  state.poly != 4 ~ 0),
           state5.bin = case_when(state.poly == 5 ~ 1,
                                  state.poly != 5 ~ 0),
           state6.bin = case_when(state.poly == 6 ~ 1,
                                  state.poly != 6 ~ 0),
           state.poly.fac = as.factor(state.poly))
  
  ### Add the predicted risks, and the logit transormation of the predicted risks to the dataset
  p.est.logit <- log(p.est/(1-p.est))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est, p.est.logit)
  
  ###
  ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
  ###
  weights <- calc.weights.stab.ce(data.raw, t.eval, max.weight = max.weight)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights)
  
  ###
  ### Fit the models
  ###
  
  ### Fit loess recalibration models to the uncensored observations at time t to calculate intercepts
  loess1 <- loess(state1.bin ~ p.est1, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.stab"])
  loess2 <- loess(state2.bin ~ p.est2, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.stab"])
  loess3 <- loess(state3.bin ~ p.est3, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.stab"])
  loess4 <- loess(state4.bin ~ p.est4, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.stab"])
  loess5 <- loess(state5.bin ~ p.est5, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.stab"])
  loess6 <- loess(state6.bin ~ p.est6, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.stab"])
  
  ### Created 'predicted observed' probabilities for each individual
  data.raw$loess.pred.obs1 <- predict(loess1, newdata = data.raw)
  data.raw$loess.pred.obs2 <- predict(loess2, newdata = data.raw)
  data.raw$loess.pred.obs3 <- predict(loess3, newdata = data.raw)
  data.raw$loess.pred.obs4 <- predict(loess4, newdata = data.raw)
  data.raw$loess.pred.obs5 <- predict(loess5, newdata = data.raw)
  data.raw$loess.pred.obs6 <- predict(loess6, newdata = data.raw)
  
  ### Create plot titles that will be used
  plot.titles <- c("Healthy", "CVD", "T2D", "CKD", "Multimorbidity", "Death")
  plots.list <- vector("list", 6)
  for (i in 1:6){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("loess.pred.obs", i, sep = "")]
    
    ### Create the plots
    plots.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>% 
                                arrange(pred) %>% select(person_id, pred, obs)) +
      geom_line(aes(x = pred, y = obs), color = "red") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, " (", plot.titles[i], ")", sep = ""))
  }
  
  ### Calculate the ECI
  ECI <- calc.ECI.ce(data.raw, p.est, data.raw[, paste("loess.pred.obs", 1:6, sep = "")])
  
  
  ### Return output object
  output.object <- list("plots.list" = plots.list,
                        "ECI" = ECI)
  return(output.object)
  
}
