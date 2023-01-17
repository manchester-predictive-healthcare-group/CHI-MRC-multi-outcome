### Program with all functions to be used in this project
###

###################################################################################################
###################################################################################################
###
### Functions in section 1 arefunctions to generate data from DGM1, applying censoring, 
### and transform into a format that can be analysed by mstate
###
###################################################################################################
###################################################################################################

###
### 1.1) Generate data according to DGM 1
###
gen.dat.DGM1 <- function(n, #number of patients to simulate
                         max.follow, #maximum follow up
                         shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                         shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                         shape15, scale15, #shape and scale for weibull baseline hazard for transition 1 -> 5
                         shape24, scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                         shape25, scale25, #shape and scale for weibull baseline hazard for transition 2 -> 5
                         shape34, scale34, #shape and scale for weibull baseline hazard for transition 3 -> 4
                         shape35, scale35, #shape and scale for weibull baseline hazard for transition 3 -> 5
                         shape45, scale45, #shape and scale for weibull baseline hazard for transition 3 -> 5
                         beta12.x1, beta12.x2, #covariate effects for transiion 12
                         beta13.x1, beta13.x2, #covariate effects for transiion 13
                         beta15.x1, beta15.x2, #covariate effects for transiion 15
                         beta24.x1, beta24.x2, #covariate effects for transiion 24
                         beta25.x1, beta25.x2, #covariate effects for transiion 25
                         beta34.x1, beta34.x2, #covariate effects for transiion 34
                         beta35.x1, beta35.x2, #covariate effects for transiion 35
                         beta45.x1, beta45.x2, #covariate effects for transiion 45
                         x.in, #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
                         numsteps) #number of sampler steps in gems data generation process
{
  
#   n.cohort <- 100
#   x.baseline <- data.frame("x1" = rnorm(n.cohort, 0, 1), "x2" = rnorm(n.cohort, 0, 1))
#   n <- n.cohort
#   max.follow <- ceiling(365.25*7)
#   
#   ### Baseline hazards
#   shape12 <- 1
#   scale12 <- 1588.598
#   
#   shape13 <- 1
#   scale13 <- 0.5*1588.598
#   
#   shape15 <- 1
#   scale15 <- 5*1588.598
#   
#   shape24 <- 1
#   scale24 <- 1588.598
#   
#   shape25 <- 1
#   scale25 <- 5*1588.598
#   
#   shape34 <- 1
#   scale34 <- 0.5*1588.598
#   
#   shape35 <- 1
#   scale35 <- 5*1588.598
#   
#   shape45 <- 1
#   scale45 <- 5*1588.598
#   
#   #qweibull(0.8, 1, 1588.598)
#   
#   ## Covariate effects
#   beta12.x1 <- 1
#   beta12.x2 <- 1
#   beta13.x1 <- 0.5
#   beta13.x2 <- 0.5
#   beta15.x1 <- 1
#   beta15.x2 <- 0.5
#   beta24.x1 <- 0.5
#   beta24.x2 <- 1
#   beta25.x1 <- 1
#   beta25.x2 <- 1
#   beta34.x1 <- 0.5
#   beta34.x2 <- 0.5
#   beta35.x1 <- 1
#   beta35.x2 <- 0.5
#   beta45.x1 <- 0.5
#   beta45.x2 <- 1
#   
#   x.in <- x.baseline
#   numsteps <- max.follow
  
  ## Generate a baseline covariate data frame
  bl <- x.in
  
  ## Generate an empty hazard matrix
  hf <- generateHazardMatrix(5)
  #hf
  
  ## Change the entries of the transitions we want to allow
  ## Define the transitions as weibull
  hf[[1, 2]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[1, 3]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[1, 5]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[2, 4]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[2, 5]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[3, 4]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[3, 5]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[4, 5]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  #print(hf)
  
  
  ## Note that by using (t + sum(history)) we are implementing a clock forward approach
  
  
  ## Generate an empty parameter matrix
  par <- generateParameterMatrix(hf)
  
  ## Use the vector of scales in each transition hazard
  par[[1, 2]] <- list(shape = shape12, scale = scale12, 
                      beta.x1 = beta12.x1, beta.x2 = beta12.x2)
  par[[1, 3]] <- list(shape = shape13, scale = scale13, 
                      beta.x1 = beta13.x1, beta.x2 = beta13.x2)
  par[[1, 5]] <- list(shape = shape15, scale = scale15, 
                      beta.x1 = beta15.x1, beta.x2 = beta15.x2)
  par[[2, 4]] <- list(shape = shape24, scale = scale24, 
                      beta.x1 = beta24.x1, beta.x2 = beta24.x2)
  par[[2, 5]] <- list(shape = shape25, scale = scale25, 
                      beta.x1 = beta25.x1, beta.x2 = beta25.x2)
  par[[3, 4]] <- list(shape = shape34, scale = scale34, 
                      beta.x1 = beta34.x1, beta.x2 = beta34.x2)
  par[[3, 5]] <- list(shape = shape35, scale = scale35, 
                      beta.x1 = beta35.x1, beta.x2 = beta35.x2)
  par[[4, 5]] <- list(shape = shape45, scale = scale45, 
                      beta.x1 = beta45.x1, beta.x2 = beta45.x2)
  
  ## Generate the cohort
  time.in <- Sys.time()
  cohort <- simulateCohort(transitionFunctions = hf, parameters = par,
                           cohortSize = n, baseline = bl, to = max.follow, sampler.steps = numsteps)
  time.out <- Sys.time()
  time.diff <- time.out - time.in

  ## Get data into the common data model
  cohort.out <- data.frame(cohort@time.to.state, cohort@baseline, patid = 1:nrow(cohort@time.to.state))
  
  return(list("cohort" = cohort.out, "max.follow" = max.follow))
}


###
### 1.2A) Convert data into a format that can be analysed using the mstate package
### This will normally be the step where we apply censoring indicator (see 1.1C), but here we do not.
### This is so we can get raw event rates, and compare with the calculated true transition probabilities, to verify both those formula,
### and the data generating mechanism
###
convert.mstate.DGM1.nocens <- function(cohort.in, max.follow){
  
  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat.mstate.temp <- select(cohort.in, paste("State.", 1:5, sep = ""))
  colnames(dat.mstate.temp) <- paste0("state", 1:5)
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat.mstate.temp.noNA <- dat.mstate.temp
  dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), max.follow, x)})))
  #head(dat.mstate.temp.noNA)
  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
  colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
    paste0("state", 1:ncol(dat.mstate.temp), ".s")
  #head(dat.mstate.temp.noNA)
  
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat.mstate.temp))){
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
  }
  #head(dat.mstate.temp.noNA)
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat.mstate.temp <- dat.mstate.temp.noNA
  
  ## Now need to add baseline data
  dat.mstate.temp$x1 <- cohort.in$x1
  dat.mstate.temp$x2 <- cohort.in$x2
  dat.mstate.temp$patid <- cohort.in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3,5), c(4,5), c(4,5), c(5), c()),
                   names = paste0("state", 1:5))
  
  
  ##############################################################################
  # Found a bug possibly??
  # When there are not categorical covariates, the msprep "keep =" functionality
  # does not appear to work.
  
  # For example if I run the following code, which does not keep x2 (the categorical 
  # variable), we get all the x1 values jumbled up and in the wrong order
  # cohort.dat.wide.clockreset.test <- msprep(cohort.dat2, trans = tmat, time = c(NA, "state2", "state3"),
  #                                      status = c(NA, "state2.s", "state3.s"), keep = c("x1"))
  # head(cohort.dat.wide.clockreset.test)
  
  # However when we include x2, as in the example below, it is fine
  ##############################################################################
  #
  # THIS BUG NO LONGER APPEARS TO BE HAPPENING - FOUND IN PREVIOUS CODE WHICH IS WHERE ABOVE TEXT IS FROM
  # LEAVING FOR FUTURE REFERENCE TO CHECK FOR THIS ISSUE
  #
  ##############################################################################
  
  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:5)),
                                 status = c(NA, paste0("state", 2:5, ".s")), 
                                 keep = c("x1","x2","patid"))
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("x1", "x2")
  dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
  head(dat.mstate.temp.wide)
  
  return(list("data.mstate" = dat.mstate.temp.wide,
              "tmat" = tmat,
              "data.raw" = cohort.in))
         
}


###
### 1.2B) Convert data into a format that can be analysed using the mstate package
### This is what will be used in the actual simulation, is it includes censoring
###
###
convert.mstate.DGM1.cens <- function(cohort.in,
                                     max.follow,
                                     cens_shape,
                                     cens_scale,
                                     cens_beta_x1,
                                     cens_beta_x2){

#   cohort.in <- cohort[["cohort"]]
#   max.follow <- cohort[["max.follow"]]
#   cens_shape <- 1
#   cens_scale <- 4000
#   cens_beta_x1 <- 0
#   cens_beta_x2 <- 0
  
  ## Generate the censoring times for each individual
  cens.times <- simsurv("weibull", lambdas = 1/cens_scale, gammas = cens_shape, x = select(cohort.in, c("x1", "x2")), 
                        betas = c("x1" = cens_beta_x1, "x2" = cens_beta_x2))
  ## Make maximum censoring time the time of max follow
  cens.times$eventtime <- pmin(cens.times$eventtime, rep(max.follow, nrow(cens.times)))
  
  ## Create raw dataset for output (the data that was inputted, plus cens.time and patid) and ignore this for rest of function
  data.raw <- data.frame(cohort.in, "cens.times" = cens.times$eventtime)
  
  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat.mstate.temp <- select(cohort.in, paste("State.", 1:5, sep = ""))
  colnames(dat.mstate.temp) <- paste0("state", 1:5)
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat.mstate.temp.noNA <- dat.mstate.temp
  dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), max.follow, x)})))
  #head(dat.mstate.temp.noNA)
  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
  colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
    paste0("state", 1:ncol(dat.mstate.temp), ".s")
  #head(dat.mstate.temp.noNA)
  
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat.mstate.temp))){
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
  }
  #head(dat.mstate.temp.noNA)
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat.mstate.temp <- dat.mstate.temp.noNA
  
  ## Add the censoring times
  dat.mstate.temp$cens.time <- cens.times$eventtime
  
  ## If any events happen after censoring, reduce event time to the censoring time, and set the event indicator to 0
  dat.mstate.temp <- dat.mstate.temp %>%
    mutate(state2 = case_when(state2 < cens.time ~ state2,
                              state2 >= cens.time ~ cens.time),
           state3 = case_when(state3 < cens.time ~ state3,
                              state3 >= cens.time ~ cens.time),
           state4 = case_when(state4 < cens.time ~ state4,
                              state4 >= cens.time ~ cens.time),
           state5 = case_when(state5 < cens.time ~ state5,
                              state5 >= cens.time ~ cens.time)) %>%
    mutate(state2.s = case_when(state2 < cens.time ~ state2.s,
                                state2 == cens.time ~ 0),
           state3.s = case_when(state3 < cens.time ~ state3.s,
                                state3 == cens.time ~ 0),
           state4.s = case_when(state4 < cens.time ~ state4.s,
                                state4 == cens.time ~ 0),
           state5.s = case_when(state5 < cens.time ~ state5.s,
                                state5 == cens.time ~ 0))
  
  ## Now need to add baseline data
  dat.mstate.temp$x1 <- cohort.in$x1
  dat.mstate.temp$x2 <- cohort.in$x2
  dat.mstate.temp$patid <- cohort.in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3,5), c(4,5), c(4,5), c(5), c()),
                   names = paste0("state", 1:5))
  
  
  ##############################################################################
  # Found a bug possibly??
  # When there are not categorical covariates, the msprep "keep =" functionality
  # does not appear to work.
  
  # For example if I run the following code, which does not keep x2 (the categorical 
  # variable), we get all the x1 values jumbled up and in the wrong order
  # cohort.dat.wide.clockreset.test <- msprep(cohort.dat2, trans = tmat, time = c(NA, "state2", "state3"),
  #                                      status = c(NA, "state2.s", "state3.s"), keep = c("x1"))
  # head(cohort.dat.wide.clockreset.test)
  
  # However when we include x2, as in the example below, it is fine
  ##############################################################################
  #
  # THIS BUG NO LONGER APPEARS TO BE HAPPENING - FOUND IN PREVIOUS CODE WHICH IS WHERE ABOVE TEXT IS FROM
  # LEAVING FOR FUTURE REFERENCE TO CHECK FOR THIS ISSUE
  #
  ##############################################################################
  
  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:5)),
                                 status = c(NA, paste0("state", 2:5, ".s")), 
                                 keep = c("x1","x2","cens.time","patid"))
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("x1", "x2")
  dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
  head(dat.mstate.temp.wide)
  
  return(list("data.mstate" = dat.mstate.temp.wide,
              "tmat" = tmat,
              "data.raw" = data.raw))
  
}


###
### 1.2C) and 1.2D) are same as 1.2A) and 1.2B), but all event times are rounded to the nearest 1, to reduce computational time
### Some event times have to have 1 added onto them, to avoid simultaneous transitions
###

###
### 1.2C): Equivalent to 1.2A), but with integer event times
###
convert.mstate.integer.DGM1.nocens <- function(cohort.in, max.follow){
  
  ## Turn all event times to integers and create a gap of at least one between all events
  cohort.in <- cohort.in %>%
    mutate(State.1 = round(State.1),
           State.2 = round(State.2),
           State.3 = round(State.3),
           State.4 = round(State.4),
           State.5 = round(State.5)) %>%
    mutate(State.2 = case_when(!is.na(State.2) & !is.na(State.1) & State.1 == State.2 ~ State.2 + 1,
                               TRUE ~ State.2),
           State.3 = case_when(!is.na(State.3) & !is.na(State.1) & State.3 == State.1 ~ State.3 + 1,
                               TRUE ~ State.3),
           State.4 = case_when(!is.na(State.4) & !is.na(State.2) & State.4 == State.2 ~ State.4 + 1,
                               !is.na(State.4) & !is.na(State.3) & State.4 == State.3 ~ State.4 + 1,
                               TRUE ~ State.4),
           State.5 = case_when(!is.na(State.5) & !is.na(State.1) & State.5 == State.1 ~ State.5 + 1,
                               !is.na(State.5) & !is.na(State.2) & State.5 == State.2 ~ State.5 + 1,
                               !is.na(State.5) & !is.na(State.3) & State.5 == State.3 ~ State.5 + 1,
                               !is.na(State.5) & !is.na(State.4) & State.5 == State.4 ~ State.5 + 1,
                               TRUE ~ State.5))
  
  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat.mstate.temp <- select(cohort.in, paste("State.", 1:5, sep = ""))
  colnames(dat.mstate.temp) <- paste0("state", 1:5)
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat.mstate.temp.noNA <- dat.mstate.temp
  dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), max.follow, x)})))
  #head(dat.mstate.temp.noNA)
  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
  colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
    paste0("state", 1:ncol(dat.mstate.temp), ".s")
  #head(dat.mstate.temp.noNA)
  
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat.mstate.temp))){
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
  }
  #head(dat.mstate.temp.noNA)
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat.mstate.temp <- dat.mstate.temp.noNA
  
  ## Now need to add baseline data
  dat.mstate.temp$x1 <- cohort.in$x1
  dat.mstate.temp$x2 <- cohort.in$x2
  dat.mstate.temp$patid <- cohort.in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3,5), c(4,5), c(4,5), c(5), c()),
                   names = paste0("state", 1:5))
  
  
  ##############################################################################
  # Found a bug possibly??
  # When there are not categorical covariates, the msprep "keep =" functionality
  # does not appear to work.
  
  # For example if I run the following code, which does not keep x2 (the categorical 
  # variable), we get all the x1 values jumbled up and in the wrong order
  # cohort.dat.wide.clockreset.test <- msprep(cohort.dat2, trans = tmat, time = c(NA, "state2", "state3"),
  #                                      status = c(NA, "state2.s", "state3.s"), keep = c("x1"))
  # head(cohort.dat.wide.clockreset.test)
  
  # However when we include x2, as in the example below, it is fine
  ##############################################################################
  #
  # THIS BUG NO LONGER APPEARS TO BE HAPPENING - FOUND IN PREVIOUS CODE WHICH IS WHERE ABOVE TEXT IS FROM
  # LEAVING FOR FUTURE REFERENCE TO CHECK FOR THIS ISSUE
  #
  ##############################################################################
  
  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:5)),
                                 status = c(NA, paste0("state", 2:5, ".s")), 
                                 keep = c("x1","x2","patid"))
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("x1", "x2")
  dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
  head(dat.mstate.temp.wide)
  
  return(list("data.mstate" = dat.mstate.temp.wide,
              "tmat" = tmat,
              "data.raw" = cohort.in))
  
}

###
### 1.2D): Equivalent to 1.2B), but with integer event times
###
convert.mstate.integer.DGM1.cens <- function(cohort.in,
                                             max.follow,
                                             cens_shape,
                                             cens_scale,
                                             cens_beta_x1,
                                             cens_beta_x2){
  
  #   cohort.in <- cohort[["cohort"]]
  #   max.follow <- cohort[["max.follow"]]
  #   cens_shape <- 1
  #   cens_scale <- 4000
  #   cens_beta_x1 <- 0
  #   cens_beta_x2 <- 0
  
  ## Turn all event times to integers and create a gap of at least one between all events
  cohort.in <- cohort.in %>%
    mutate(State.1 = round(State.1),
           State.2 = round(State.2),
           State.3 = round(State.3),
           State.4 = round(State.4),
           State.5 = round(State.5)) %>%
    mutate(State.2 = case_when(!is.na(State.2) & !is.na(State.1) & State.1 == State.2 ~ State.2 + 1,
                               TRUE ~ State.2),
           State.3 = case_when(!is.na(State.3) & !is.na(State.1) & State.3 == State.1 ~ State.3 + 1,
                               TRUE ~ State.3),
           State.4 = case_when(!is.na(State.4) & !is.na(State.2) & State.4 == State.2 ~ State.4 + 1,
                               !is.na(State.4) & !is.na(State.3) & State.4 == State.3 ~ State.4 + 1,
                               TRUE ~ State.4),
           State.5 = case_when(!is.na(State.5) & !is.na(State.1) & State.5 == State.1 ~ State.5 + 1,
                               !is.na(State.5) & !is.na(State.2) & State.5 == State.2 ~ State.5 + 1,
                               !is.na(State.5) & !is.na(State.3) & State.5 == State.3 ~ State.5 + 1,
                               !is.na(State.5) & !is.na(State.4) & State.5 == State.4 ~ State.5 + 1,
                               TRUE ~ State.5))
  
  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat.mstate.temp <- select(cohort.in, paste("State.", 1:5, sep = ""))
  colnames(dat.mstate.temp) <- paste0("state", 1:5)
  
  ## Generate the censoring times for each individual
  cens.times <- simsurv("weibull", lambdas = 1/cens_scale, gammas = cens_shape, x = select(cohort.in, c("x1", "x2")), 
                        betas = c("x1" = cens_beta_x1, "x2" = cens_beta_x2))
  
  ## If censoring time < 1, set to 1
  cens.times$eventtime[cens.times$eventtime < 1] <- 1
  
  ## Make maximum censoring time the time of max follow
  cens.times$eventtime <- pmin(round(cens.times$eventtime), rep(max.follow, nrow(cens.times)))
  
  ## Create raw dataset for output (the data that was inputted, plus cens.time and patid) and ignore this for rest of function
  data.raw <- data.frame(cohort.in, "cens.times" = cens.times$eventtime)
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat.mstate.temp.noNA <- dat.mstate.temp
  dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), max.follow, x)})))
  #head(dat.mstate.temp.noNA)
  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
  colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
    paste0("state", 1:ncol(dat.mstate.temp), ".s")
  #head(dat.mstate.temp.noNA)
  
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat.mstate.temp))){
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
  }
  #head(dat.mstate.temp.noNA)
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat.mstate.temp <- dat.mstate.temp.noNA
  
  ## Add the censoring times
  dat.mstate.temp$cens.time <- cens.times$eventtime
  
  ## If any events happen after censoring, reduce event time to the censoring time, and set the event indicator to 0
  dat.mstate.temp <- dat.mstate.temp %>%
    mutate(state2 = case_when(state2 < cens.time ~ state2,
                              state2 >= cens.time ~ cens.time),
           state3 = case_when(state3 < cens.time ~ state3,
                              state3 >= cens.time ~ cens.time),
           state4 = case_when(state4 < cens.time ~ state4,
                              state4 >= cens.time ~ cens.time),
           state5 = case_when(state5 < cens.time ~ state5,
                              state5 >= cens.time ~ cens.time)) %>%
    mutate(state2.s = case_when(state2 < cens.time ~ state2.s,
                                state2 == cens.time ~ 0),
           state3.s = case_when(state3 < cens.time ~ state3.s,
                                state3 == cens.time ~ 0),
           state4.s = case_when(state4 < cens.time ~ state4.s,
                                state4 == cens.time ~ 0),
           state5.s = case_when(state5 < cens.time ~ state5.s,
                                state5 == cens.time ~ 0))
  
  ## Now need to add baseline data
  dat.mstate.temp$x1 <- cohort.in$x1
  dat.mstate.temp$x2 <- cohort.in$x2
  dat.mstate.temp$patid <- cohort.in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3,5), c(4,5), c(4,5), c(5), c()),
                   names = paste0("state", 1:5))
  
  
  ##############################################################################
  # Found a bug possibly??
  # When there are not categorical covariates, the msprep "keep =" functionality
  # does not appear to work.
  
  # For example if I run the following code, which does not keep x2 (the categorical 
  # variable), we get all the x1 values jumbled up and in the wrong order
  # cohort.dat.wide.clockreset.test <- msprep(cohort.dat2, trans = tmat, time = c(NA, "state2", "state3"),
  #                                      status = c(NA, "state2.s", "state3.s"), keep = c("x1"))
  # head(cohort.dat.wide.clockreset.test)
  
  # However when we include x2, as in the example below, it is fine
  ##############################################################################
  #
  # THIS BUG NO LONGER APPEARS TO BE HAPPENING - FOUND IN PREVIOUS CODE WHICH IS WHERE ABOVE TEXT IS FROM
  # LEAVING FOR FUTURE REFERENCE TO CHECK FOR THIS ISSUE
  #
  ##############################################################################
  
  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:5)),
                                 status = c(NA, paste0("state", 2:5, ".s")), 
                                 keep = c("x1","x2","cens.time","patid"))
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("x1", "x2")
  dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
  head(dat.mstate.temp.wide)
  
  return(list("data.mstate" = dat.mstate.temp.wide,
              "tmat" = tmat,
              "data.raw" = data.raw))
  
}


###
### 1.2E) Convert a dataset into mstate format, but when the censoring mechanism has already been generated, and is a variable within cohort.in
###
convert.mstate.integer.DGM1.cens.preexist <- function(cohort.in,
                                                      max.follow){
  
  ## Turn all event times to integers and create a gap of at least one between all events
  cohort.in <- cohort.in %>%
    mutate(State.1 = round(State.1),
           State.2 = round(State.2),
           State.3 = round(State.3),
           State.4 = round(State.4),
           State.5 = round(State.5)) %>%
    mutate(State.2 = case_when(!is.na(State.2) & !is.na(State.1) & State.1 == State.2 ~ State.2 + 1,
                               TRUE ~ State.2),
           State.3 = case_when(!is.na(State.3) & !is.na(State.1) & State.3 == State.1 ~ State.3 + 1,
                               TRUE ~ State.3),
           State.4 = case_when(!is.na(State.4) & !is.na(State.2) & State.4 == State.2 ~ State.4 + 1,
                               !is.na(State.4) & !is.na(State.3) & State.4 == State.3 ~ State.4 + 1,
                               TRUE ~ State.4),
           State.5 = case_when(!is.na(State.5) & !is.na(State.1) & State.5 == State.1 ~ State.5 + 1,
                               !is.na(State.5) & !is.na(State.2) & State.5 == State.2 ~ State.5 + 1,
                               !is.na(State.5) & !is.na(State.3) & State.5 == State.3 ~ State.5 + 1,
                               !is.na(State.5) & !is.na(State.4) & State.5 == State.4 ~ State.5 + 1,
                               TRUE ~ State.5))
  
  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat.mstate.temp <- select(cohort.in, paste("State.", 1:5, sep = ""))
  colnames(dat.mstate.temp) <- paste0("state", 1:5)
  
  ## If censoring time < 1, set to 1
  cohort.in$cens.times[cohort.in$cens.times < 1] <- 1
  
  ## Make maximum censoring time the time of max follow
  cohort.in$cens.times <- pmin(round(cohort.in$cens.times), rep(max.follow, nrow(cohort.in)))
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat.mstate.temp.noNA <- dat.mstate.temp
  dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), max.follow, x)})))
  #head(dat.mstate.temp.noNA)
  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
  colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
    paste0("state", 1:ncol(dat.mstate.temp), ".s")
  #head(dat.mstate.temp.noNA)
  
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat.mstate.temp))){
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
  }
  #head(dat.mstate.temp.noNA)
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat.mstate.temp <- dat.mstate.temp.noNA
  
  ## Add the censoring times
  dat.mstate.temp$cens.time <- cohort.in$cens.times
  
  ## If any events happen after censoring, reduce event time to the censoring time, and set the event indicator to 0
  dat.mstate.temp <- dat.mstate.temp %>%
    mutate(state2 = case_when(state2 < cens.time ~ state2,
                              state2 >= cens.time ~ cens.time),
           state3 = case_when(state3 < cens.time ~ state3,
                              state3 >= cens.time ~ cens.time),
           state4 = case_when(state4 < cens.time ~ state4,
                              state4 >= cens.time ~ cens.time),
           state5 = case_when(state5 < cens.time ~ state5,
                              state5 >= cens.time ~ cens.time)) %>%
    mutate(state2.s = case_when(state2 < cens.time ~ state2.s,
                                state2 == cens.time ~ 0),
           state3.s = case_when(state3 < cens.time ~ state3.s,
                                state3 == cens.time ~ 0),
           state4.s = case_when(state4 < cens.time ~ state4.s,
                                state4 == cens.time ~ 0),
           state5.s = case_when(state5 < cens.time ~ state5.s,
                                state5 == cens.time ~ 0))
  
  ## Now need to add baseline data
  dat.mstate.temp$x1 <- cohort.in$x1
  dat.mstate.temp$x2 <- cohort.in$x2
  dat.mstate.temp$patid <- cohort.in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3,5), c(4,5), c(4,5), c(5), c()),
                   names = paste0("state", 1:5))
  
  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:5)),
                                 status = c(NA, paste0("state", 2:5, ".s")), 
                                 keep = c("x1","x2","cens.time","patid"))
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("x1", "x2")
  dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
  
  return(dat.mstate.temp.wide)
  
}

###################################################################################################
###################################################################################################
###
### Functions in section 2 are functions to generate data from DGM2, applying censoring, 
### and transform into a format that can be analysed by mstate
###
###################################################################################################
###################################################################################################

###
### 2.1) Generate data according to DGM 2
###
gen.dat.DGM2 <- function(n, #number of patients to simulate
                         max.follow, #maximum follow up
                         shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                         shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                         shape16, scale16, #shape and scale for weibull baseline hazard for transition 1 -> 6
                         shape24, scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                         shape26, scale26, #shape and scale for weibull baseline hazard for transition 2 -> 6
                         shape35, scale35, #shape and scale for weibull baseline hazard for transition 3 -> 5
                         shape36, scale36, #shape and scale for weibull baseline hazard for transition 3 -> 6
                         shape46, scale46, #shape and scale for weibull baseline hazard for transition 4 -> 6
                         shape56, scale56, #shape and scale for weibull baseline hazard for transition 5 -> 6
                         beta12.x1, beta12.x2, #covariate effects for transiion 1 -> 2
                         beta13.x1, beta13.x2, #covariate effects for transiion 1 -> 3
                         beta16.x1, beta16.x2, #covariate effects for transiion 1 -> 6
                         beta24.x1, beta24.x2, #covariate effects for transiion 2 -> 4
                         beta26.x1, beta26.x2, #covariate effects for transiion 2 -> 6
                         beta35.x1, beta35.x2, #covariate effects for transiion 3 -> 5
                         beta36.x1, beta36.x2, #covariate effects for transiion 3 -> 6
                         beta46.x1, beta46.x2, #covariate effects for transiion 4 -> 6
                         beta56.x1, beta56.x2, #covariate effects for transiion 5 -> 6
                         x.in, #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
                         numsteps) #number of sampler steps in gems data generation process
{
  
#     n.cohort <- 100
#     x.baseline <- data.frame("x1" = rnorm(n.cohort, 0, 1), "x2" = rnorm(n.cohort, 0, 1))
#     n <- n.cohort
#     max.follow <- ceiling(365.25*7)
#     
#   ### Baseline hazards
#   shape12 <- 1
#   scale12 <- 1588.598
#   
#   shape13 <- 1
#   scale13 <- 0.5*1588.598
#   
#   shape16 <- 1
#   scale16 <- 5*1588.598
#   
#   shape24 <- 1
#   scale24 <- 1588.598
#   
#   shape26 <- 1
#   scale26 <- 5*1588.598
#   
#   shape35 <- 1
#   scale35 <- 0.5*1588.598
#   
#   shape36 <- 1
#   scale36 <- 5*1588.598
#   
#   shape46 <- 1
#   scale46 <- 5*1588.598
#   
#   shape56 <- 1
#   scale56 <- 5*1588.598
#   
#   #qweibull(0.8, 1, 1588.598)
#   
#   ## Covariate effects
#   beta12.x1 <- 1
#   beta12.x2 <- 1
#   beta13.x1 <- 0.5
#   beta13.x2 <- 0.5
#   beta16.x1 <- 1
#   beta16.x2 <- 0.5
#   beta24.x1 <- 0.5
#   beta24.x2 <- 1
#   beta26.x1 <- 1
#   beta26.x2 <- 1
#   beta35.x1 <- 0.5
#   beta35.x2 <- 0.5
#   beta36.x1 <- 1
#   beta36.x2 <- 0.5
#   beta46.x1 <- 0.5
#   beta46.x2 <- 1
#   beta56.x1 <- 0.5
#   beta56.x2 <- 1
#     
#     x.in <- x.baseline
#     numsteps <- max.follow
    
  
  ## Generate a baseline covariate data frame
  bl <- x.in
  
  ## Generate an empty hazard matrix
  hf <- generateHazardMatrix(6)
  #hf
  
  ## Change the entries of the transitions we want to allow
  ## Define the transitions as weibull
  hf[[1, 2]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[1, 3]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[1, 6]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[2, 4]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[2, 6]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[3, 5]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[3, 6]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[4, 6]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[5, 6]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  #print(hf)
  
  
  ## Note that by using (t + sum(history)) we are implementing a clock forward approach
  
  
  ## Generate an empty parameter matrix
  par <- generateParameterMatrix(hf)
  
  ## Use the vector of scales in each transition hazard
  par[[1, 2]] <- list(shape = shape12, scale = scale12, 
                      beta.x1 = beta12.x1, beta.x2 = beta12.x2)
  par[[1, 3]] <- list(shape = shape13, scale = scale13, 
                      beta.x1 = beta13.x1, beta.x2 = beta13.x2)
  par[[1, 6]] <- list(shape = shape16, scale = scale16, 
                      beta.x1 = beta16.x1, beta.x2 = beta16.x2)
  par[[2, 4]] <- list(shape = shape24, scale = scale24, 
                      beta.x1 = beta24.x1, beta.x2 = beta24.x2)
  par[[2, 6]] <- list(shape = shape26, scale = scale26, 
                      beta.x1 = beta26.x1, beta.x2 = beta26.x2)
  par[[3, 5]] <- list(shape = shape35, scale = scale35, 
                      beta.x1 = beta35.x1, beta.x2 = beta35.x2)
  par[[3, 6]] <- list(shape = shape36, scale = scale36, 
                      beta.x1 = beta36.x1, beta.x2 = beta36.x2)
  par[[4, 6]] <- list(shape = shape46, scale = scale46, 
                      beta.x1 = beta46.x1, beta.x2 = beta46.x2)
  par[[5, 6]] <- list(shape = shape56, scale = scale56, 
                      beta.x1 = beta56.x1, beta.x2 = beta56.x2)
  
  ## Generate the cohort
  time.in <- Sys.time()
  cohort <- simulateCohort(transitionFunctions = hf, parameters = par,
                           cohortSize = n, baseline = bl, to = max.follow, sampler.steps = numsteps)
  time.out <- Sys.time()
  time.diff <- time.out - time.in
  
  ## Get data into the common data model
  cohort.out <- data.frame(cohort@time.to.state, cohort@baseline, patid = 1:nrow(cohort@time.to.state))
  
  return(list("cohort" = cohort.out, "max.follow" = max.follow))
}



###
### 2.2A) Convert data into a format that can be analysed using the mstate package
### Leaving states 4 and 5 seperate, because for this function we want to generate data, simply to check the true transition probabilities are correct
### Therefore we want to know how many people are seperately in states 4 and 5
### This will also normally be the step where we apply censoring indicator (see 1.1C), but here we do not.
### This is so we can get raw event rates, and compare with the calculated true transition probabilities, to verify both those formula,
### and the data generating mechanism
###
convert.mstate.DGM2.seperate45.nocens <- function(cohort.in,
                                                  max.follow){
  
  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat.mstate.temp <- select(cohort.in, paste("State.", 1:6, sep = ""))
  colnames(dat.mstate.temp) <- paste0("state", 1:6)
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat.mstate.temp.noNA <- dat.mstate.temp
  dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), max.follow, x)})))
  #head(dat.mstate.temp.noNA)
  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
  colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
    paste0("state", 1:ncol(dat.mstate.temp), ".s")
  #head(dat.mstate.temp.noNA)
  
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat.mstate.temp))){
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
  }
  #head(dat.mstate.temp.noNA)
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat.mstate.temp <- dat.mstate.temp.noNA
  
  ## Now need to add baseline data
  dat.mstate.temp$x1 <- cohort.in$x1
  dat.mstate.temp$x2 <- cohort.in$x2
  dat.mstate.temp$patid <- cohort.in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3,6), c(4,6), c(5,6), c(6), c(6), c()),
                   names = paste0("state", 1:6))
  
  
  ##############################################################################
  # Found a bug possibly??
  # When there are not categorical covariates, the msprep "keep =" functionality
  # does not appear to work.
  
  # For example if I run the following code, which does not keep x2 (the categorical 
  # variable), we get all the x1 values jumbled up and in the wrong order
  # cohort.dat.wide.clockreset.test <- msprep(cohort.dat2, trans = tmat, time = c(NA, "state2", "state3"),
  #                                      status = c(NA, "state2.s", "state3.s"), keep = c("x1"))
  # head(cohort.dat.wide.clockreset.test)
  
  # However when we include x2, as in the example below, it is fine
  ##############################################################################
  #
  # THIS BUG NO LONGER APPEARS TO BE HAPPENING - FOUND IN PREVIOUS CODE WHICH IS WHERE ABOVE TEXT IS FROM
  # LEAVING FOR FUTURE REFERENCE TO CHECK FOR THIS ISSUE
  #
  ##############################################################################
  
  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:6)),
                                 status = c(NA, paste0("state", 2:6, ".s")), 
                                 keep = c("x1","x2","patid"))
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("x1", "x2")
  dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
  head(dat.mstate.temp.wide)
  
  return(list("data.mstate" = dat.mstate.temp.wide,
              "tmat" = tmat,
              "data.raw" = cohort.in))
  
}



###
### 2.2B) Convert data into a format that can be analysed using the mstate package
### Comining states 4 and 5 seperate, because we want to analyse data by falsely assuming a structure where ordering of developing states 2 and 3
### does not impact transition rates

### We also apply censoring
###
convert.mstate.DGM2.combine45.cens <- function(cohort.in,
                                               max.follow,
                                               cens_shape,
                                               cens_scale,
                                               cens_beta_x1,
                                               cens_beta_x2){
  
  #   cens_shape <- 1
  #   cens_scale <- 5*1588.598
  #   cens_beta_x1 <- 1
  #   cens_beta_x2 <- 0.5
  #   data.in <- cohort
  #   str(cohort.in@baseline)
  
  ## Generate the censoring times for each individual
  cens.times <- simsurv("weibull", lambdas = 1/cens_scale, gammas = cens_shape, x = select(cohort.in, c("x1", "x2")), 
                        betas = c("x1" = cens_beta_x1, "x2" = cens_beta_x2))
  ## Make maximum censoring time the time of max follow
  cens.times$eventtime <- pmin(cens.times$eventtime, rep(max.follow, nrow(cens.times)))
  
  ## Create raw dataset for output (the data that was inputted, plus cens.time and patid) and ignore this for rest of function
  data.raw <- data.frame(cohort.in, "cens.times" = cens.times$eventtime)
  
  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat.mstate.temp <- select(cohort.in, paste("State.", 1:6, sep = ""))
  colnames(dat.mstate.temp) <- paste0("state", 1:6)
  
  ## Combine states 4 and 5 (they are mutually exclusive, so just take the number which isn't NA, if either are not NA)
  dat.mstate.temp$state45 <- rep(NA, nrow(dat.mstate.temp))
  dat.mstate.temp$state45[!is.na(dat.mstate.temp$state4)] <- dat.mstate.temp$state4[!is.na(dat.mstate.temp$state4)]
  dat.mstate.temp$state45[!is.na(dat.mstate.temp$state5)] <- dat.mstate.temp$state5[!is.na(dat.mstate.temp$state5)]
  
  ## Turn state 45 into state 4, and turn state 6 (death) into a new state 5
  dat.mstate.temp$state4 <- dat.mstate.temp$state45
  dat.mstate.temp$state5 <- dat.mstate.temp$state6
  
  ## Remove states 45 and 6
  dat.mstate.temp <- dat.mstate.temp[, 1:5]
  
  ## Now it has the structure of DGM1, but was generated through DGM2
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat.mstate.temp.noNA <- dat.mstate.temp
  dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), max.follow, x)})))
  #head(dat.mstate.temp.noNA)
  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
  colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
    paste0("state", 1:ncol(dat.mstate.temp), ".s")
  #head(dat.mstate.temp.noNA)
  
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat.mstate.temp))){
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
  }
  #head(dat.mstate.temp.noNA)
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat.mstate.temp <- dat.mstate.temp.noNA
  
  ## Add the censoring times
  dat.mstate.temp$cens.time <- cens.times$eventtime
  
  ## If any events happen after censoring, reduce event time to the censoring time, and set the event indicator to 0
  dat.mstate.temp <- dat.mstate.temp %>%
    mutate(state2 = case_when(state2 < cens.time ~ state2,
                              state2 >= cens.time ~ cens.time),
           state3 = case_when(state3 < cens.time ~ state3,
                              state3 >= cens.time ~ cens.time),
           state4 = case_when(state4 < cens.time ~ state4,
                              state4 >= cens.time ~ cens.time),
           state5 = case_when(state5 < cens.time ~ state5,
                              state5 >= cens.time ~ cens.time)) %>%
    mutate(state2.s = case_when(state2 < cens.time ~ state2.s,
                                state2 == cens.time ~ 0),
           state3.s = case_when(state3 < cens.time ~ state3.s,
                                state3 == cens.time ~ 0),
           state4.s = case_when(state4 < cens.time ~ state4.s,
                                state4 == cens.time ~ 0),
           state5.s = case_when(state5 < cens.time ~ state5.s,
                                state5 == cens.time ~ 0))
  
  ## Now need to add baseline data
  dat.mstate.temp$x1 <- cohort.in$x1
  dat.mstate.temp$x2 <- cohort.in$x2
  dat.mstate.temp$patid <- cohort.in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3,5), c(4,5), c(4,5), c(5), c()),
                   names = paste0("state", 1:5))
  
  
  ##############################################################################
  # Found a bug possibly??
  # When there are not categorical covariates, the msprep "keep =" functionality
  # does not appear to work.
  
  # For example if I run the following code, which does not keep x2 (the categorical 
  # variable), we get all the x1 values jumbled up and in the wrong order
  # cohort.dat.wide.clockreset.test <- msprep(cohort.dat2, trans = tmat, time = c(NA, "state2", "state3"),
  #                                      status = c(NA, "state2.s", "state3.s"), keep = c("x1"))
  # head(cohort.dat.wide.clockreset.test)
  
  # However when we include x2, as in the example below, it is fine
  ##############################################################################
  #
  # THIS BUG NO LONGER APPEARS TO BE HAPPENING - FOUND IN PREVIOUS CODE WHICH IS WHERE ABOVE TEXT IS FROM
  # LEAVING FOR FUTURE REFERENCE TO CHECK FOR THIS ISSUE
  #
  ##############################################################################
  
  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:5)),
                                 status = c(NA, paste0("state", 2:5, ".s")), 
                                 keep = c("x1","x2","cens.time","patid"))
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("x1", "x2")
  dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
  head(dat.mstate.temp.wide)
  
  return(list("data.mstate" = dat.mstate.temp.wide,
              "tmat" = tmat,
              "data.raw" = data.raw))
  
}


###
### 2.2C) and 2.2D) are same as 2.2A) and 2.2B), but all event times are rounded to the nearest 1, to reduce computational time
### Some event times have to have 1 added onto them, to avoid simultaneous transitions
###

### NB HAVENT ACTUALLY DONE THESE YET, AS IT DIDNT APPEAR TO REDUCE COMPUTATIONAL TIMES MUCH

# ###
# ### 2.2C): Equivalent to 2.2A), but with integer event times
# ###
# convert.mstate.integer.DGM2.seperate45.nocens <- function(cohort.in,
#                                                   max.follow){
#   
#   ## Turn event times into a dataframe and make the colnames not have any spaces in them
#   dat.mstate.temp <- select(cohort.in, paste("State.", 1:6, sep = ""))
#   colnames(dat.mstate.temp) <- paste0("state", 1:6)
#   
#   ## Now set any transitions that didn't happen to the maximum value of follow up
#   ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
#   ## at this point in time (when follow up stops)
#   dat.mstate.temp.noNA <- dat.mstate.temp
#   dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), max.follow, x)})))
#   #head(dat.mstate.temp.noNA)
#   
#   ## Add censoring variables
#   dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
#   colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
#     paste0("state", 1:ncol(dat.mstate.temp), ".s")
#   #head(dat.mstate.temp.noNA)
#   
#   
#   ## If it is not an NA value (from original dataset), set the censoring indicator to 1
#   for (i in (2:ncol(dat.mstate.temp))){
#     dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
#   }
#   #head(dat.mstate.temp.noNA)
#   
#   ## Rename dataset to what it was before, and remove excess dataset
#   dat.mstate.temp <- dat.mstate.temp.noNA
#   
#   ## Now need to add baseline data
#   dat.mstate.temp$x1 <- cohort.in$x1
#   dat.mstate.temp$x2 <- cohort.in$x2
#   dat.mstate.temp$patid <- cohort.in$patid
#   
#   ### Now we can use msprep from the mstate package to turn into wide format
#   ## First create a transition matrix corresponding to the columns
#   tmat <- transMat(x = list(c(2,3,6), c(4,6), c(5,6), c(6), c(6), c()),
#                    names = paste0("state", 1:6))
#   
#   
#   ##############################################################################
#   # Found a bug possibly??
#   # When there are not categorical covariates, the msprep "keep =" functionality
#   # does not appear to work.
#   
#   # For example if I run the following code, which does not keep x2 (the categorical 
#   # variable), we get all the x1 values jumbled up and in the wrong order
#   # cohort.dat.wide.clockreset.test <- msprep(cohort.dat2, trans = tmat, time = c(NA, "state2", "state3"),
#   #                                      status = c(NA, "state2.s", "state3.s"), keep = c("x1"))
#   # head(cohort.dat.wide.clockreset.test)
#   
#   # However when we include x2, as in the example below, it is fine
#   ##############################################################################
#   #
#   # THIS BUG NO LONGER APPEARS TO BE HAPPENING - FOUND IN PREVIOUS CODE WHICH IS WHERE ABOVE TEXT IS FROM
#   # LEAVING FOR FUTURE REFERENCE TO CHECK FOR THIS ISSUE
#   #
#   ##############################################################################
#   
#   ## Now can prepare the data into wide format
#   dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
#                                  time = c(NA, paste0("state", 2:6)),
#                                  status = c(NA, paste0("state", 2:6, ".s")), 
#                                  keep = c("x1","x2","patid"))
#   
#   ## Want to expand the covariates to allow different covariate effects per transition
#   covs <- c("x1", "x2")
#   dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
#   head(dat.mstate.temp.wide)
#   
#   return(list("data.mstate" = dat.mstate.temp.wide,
#               "tmat" = tmat,
#               "data.raw" = cohort.in))
#   
# }
# 
# 
# 
# ###
# ### 2.2D): Equivalent to 2.2B), but with integer event times
# ###
# convert.mstate.integer.DGM2.combine45.cens <- function(cohort.in,
#                                                max.follow,
#                                                cens_shape,
#                                                cens_scale,
#                                                cens_beta_x1,
#                                                cens_beta_x2){
#   
#   #   cens_shape <- 1
#   #   cens_scale <- 5*1588.598
#   #   cens_beta_x1 <- 1
#   #   cens_beta_x2 <- 0.5
#   #   data.in <- cohort
#   #   str(cohort.in@baseline)
#   
#   ## Generate the censoring times for each individual
#   cens.times <- simsurv("weibull", lambdas = 1/cens_scale, gammas = cens_shape, x = select(cohort.in, c("x1", "x2")), 
#                         betas = c("x1" = cens_beta_x1, "x2" = cens_beta_x2))
#   ## Make maximum censoring time the time of max follow
#   cens.times$eventtime <- pmin(cens.times$eventtime, rep(max.follow, nrow(cens.times)))
#   
#   ## Create raw dataset for output (the data that was inputted, plus cens.time and patid) and ignore this for rest of function
#   data.raw <- data.frame(cohort.in, "cens.times" = cens.times$eventtime)
#   
#   ## Turn event times into a dataframe and make the colnames not have any spaces in them
#   dat.mstate.temp <- select(cohort.in, paste("State.", 1:6, sep = ""))
#   colnames(dat.mstate.temp) <- paste0("state", 1:6)
#   
#   ## Combine states 4 and 5 (they are mutually exclusive, so just take the number which isn't NA, if either are not NA)
#   dat.mstate.temp$state45 <- rep(NA, nrow(dat.mstate.temp))
#   dat.mstate.temp$state45[!is.na(dat.mstate.temp$state4)] <- dat.mstate.temp$state4[!is.na(dat.mstate.temp$state4)]
#   dat.mstate.temp$state45[!is.na(dat.mstate.temp$state5)] <- dat.mstate.temp$state5[!is.na(dat.mstate.temp$state5)]
#   
#   ## Turn state 45 into state 4, and turn state 6 (death) into a new state 5
#   dat.mstate.temp$state4 <- dat.mstate.temp$state45
#   dat.mstate.temp$state5 <- dat.mstate.temp$state6
#   
#   ## Remove states 45 and 6
#   dat.mstate.temp <- dat.mstate.temp[, 1:5]
#   
#   ## Now it has the structure of DGM1, but was generated through DGM2
#   
#   ## Now set any transitions that didn't happen to the maximum value of follow up
#   ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
#   ## at this point in time (when follow up stops)
#   dat.mstate.temp.noNA <- dat.mstate.temp
#   dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), max.follow, x)})))
#   #head(dat.mstate.temp.noNA)
#   
#   ## Add censoring variables
#   dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
#   colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
#     paste0("state", 1:ncol(dat.mstate.temp), ".s")
#   #head(dat.mstate.temp.noNA)
#   
#   
#   ## If it is not an NA value (from original dataset), set the censoring indicator to 1
#   for (i in (2:ncol(dat.mstate.temp))){
#     dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
#   }
#   #head(dat.mstate.temp.noNA)
#   
#   ## Rename dataset to what it was before, and remove excess dataset
#   dat.mstate.temp <- dat.mstate.temp.noNA
#   
#   ## Add the censoring times
#   dat.mstate.temp$cens.time <- cens.times$eventtime
#   
#   ## If any events happen after censoring, reduce event time to the censoring time, and set the event indicator to 0
#   dat.mstate.temp <- dat.mstate.temp %>%
#     mutate(state2 = case_when(state2 < cens.time ~ state2,
#                               state2 >= cens.time ~ cens.time),
#            state3 = case_when(state3 < cens.time ~ state3,
#                               state3 >= cens.time ~ cens.time),
#            state4 = case_when(state4 < cens.time ~ state4,
#                               state4 >= cens.time ~ cens.time),
#            state5 = case_when(state5 < cens.time ~ state5,
#                               state5 >= cens.time ~ cens.time)) %>%
#     mutate(state2.s = case_when(state2 < cens.time ~ state2.s,
#                                 state2 == cens.time ~ 0),
#            state3.s = case_when(state3 < cens.time ~ state3.s,
#                                 state3 == cens.time ~ 0),
#            state4.s = case_when(state4 < cens.time ~ state4.s,
#                                 state4 == cens.time ~ 0),
#            state5.s = case_when(state5 < cens.time ~ state5.s,
#                                 state5 == cens.time ~ 0))
#   
#   ## Now need to add baseline data
#   dat.mstate.temp$x1 <- cohort.in$x1
#   dat.mstate.temp$x2 <- cohort.in$x2
#   dat.mstate.temp$patid <- cohort.in$patid
#   
#   ### Now we can use msprep from the mstate package to turn into wide format
#   ## First create a transition matrix corresponding to the columns
#   tmat <- transMat(x = list(c(2,3,5), c(4,5), c(4,5), c(5), c()),
#                    names = paste0("state", 1:5))
#   
#   
#   ##############################################################################
#   # Found a bug possibly??
#   # When there are not categorical covariates, the msprep "keep =" functionality
#   # does not appear to work.
#   
#   # For example if I run the following code, which does not keep x2 (the categorical 
#   # variable), we get all the x1 values jumbled up and in the wrong order
#   # cohort.dat.wide.clockreset.test <- msprep(cohort.dat2, trans = tmat, time = c(NA, "state2", "state3"),
#   #                                      status = c(NA, "state2.s", "state3.s"), keep = c("x1"))
#   # head(cohort.dat.wide.clockreset.test)
#   
#   # However when we include x2, as in the example below, it is fine
#   ##############################################################################
#   #
#   # THIS BUG NO LONGER APPEARS TO BE HAPPENING - FOUND IN PREVIOUS CODE WHICH IS WHERE ABOVE TEXT IS FROM
#   # LEAVING FOR FUTURE REFERENCE TO CHECK FOR THIS ISSUE
#   #
#   ##############################################################################
#   
#   ## Now can prepare the data into wide format
#   dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
#                                  time = c(NA, paste0("state", 2:5)),
#                                  status = c(NA, paste0("state", 2:5, ".s")), 
#                                  keep = c("x1","x2","cens.time","patid"))
#   
#   ## Want to expand the covariates to allow different covariate effects per transition
#   covs <- c("x1", "x2")
#   dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
#   head(dat.mstate.temp.wide)
#   
#   return(list("data.mstate" = dat.mstate.temp.wide,
#               "tmat" = tmat,
#               "data.raw" = data.raw))
#   
# }



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
calc.calib.aj <- function(data.mstate, tmat, t.eval){
  
  #     t.eval <- ceiling(7*365.25)
  #     data.mstate <- data.mstate.obj[["data.mstate"]]
  #     tmat <- data.mstate.obj[["tmat"]]
  # data.mstate <- data.mstate.obj[["data.mstate"]][data.mstate.obj[["data.mstate"]]$patid %in% 1:5000, ]
  
  ### Fit csh's with no predictors
  csh.aj <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data.mstate)
  
  ### Calculate cumulative incidence functions
  msfit.aj <- msfit(csh.aj, trans = tmat)
  
  ### Calculate Aalen-Johansen estimator
  pt.aj <- probtrans(msfit.aj, predt = 0)
  
  ### Extract the closest time in the data to the time we want to evaluate at
  t.eval.dat <- pt.aj[[1]]$time[max(which(pt.aj[[1]]$time <= t.eval))]
  
  ### Extract AJ estimator at this time point
  obs.aj <- pt.aj[[1]][pt.aj[[1]]$time == t.eval.dat, paste("pstate", 1:5, sep = "")]
  
  ### Extract AJ standard error  at this time point
  obs.aj.se <- pt.aj[[1]][pt.aj[[1]]$time == t.eval.dat, paste("se", 1:5, sep = "")]
  
  ### Create output object
  output.object <- list("obs.aj" = obs.aj, "obs.aj.se" = obs.aj.se)
  
  return(output.object)
}


###
### 3.2) Calculate landmark Aalen-Johansen estimator for a cohort at time t.eval
###
calc.calib.lmaj <- function(data.mstate, t.eval){
  
  #   t.eval <- ceiling(7*365.25)
  #     data.mstate <- data.mstate.obj[["data.mstate"]]
  #     tmat <- data.mstate.obj[["tmat"]]
  
  ### Calculate observed risks using landmark Aalen-Johansen
  pt.lmaj <- LMAJ(msdata = data.mstate, s = 0, from = 1)
  
  ### Extract the closest time in the data to the time we want to evaluate at
  t.eval.dat <- pt.lmaj$time[max(which(pt.lmaj$time <= t.eval))]
  
  ### Get the LMAJ estimator at this time point
  obs.lmaj <- pt.lmaj[pt.lmaj$time == t.eval.dat, paste("pstate", 1:5, sep = "")]
  
  return(obs.lmaj)
}


###
### 3.3) Function to extract all individuals in state j, at time t.eval, from a dataset in mstate format
### Used in other functions which assess calibration using blr/mlr at specific time points (3.4 and 3.5)
### 
extract.ids.states <- function(data.mstate, j, t.eval){
  
  ### For non-absorbing states, to be in state j at time t, you must have an observations from state j, where Tstart <= t.eval < Tstop
  if (j < max(data.mstate$to)){
    ## Extract ids
    ids.state.j <- subset(data.mstate, from == j & Tstart <= t.eval & t.eval < Tstop) %>%
      select(patid) %>%
      distinct(patid)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$patid)
  } else if (j == max(data.mstate$to)){
    ### For absorbind state, just have to have moved into it
    ids.state.j <- subset(data.mstate, to == max(data.mstate$to) & t.eval >= Tstop & status == 1) %>%
      select(patid) %>%
      distinct(patid)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$patid)
  }
  
  return(ids.state.j)
}


###
### 3.4A) Assess calibration using binary logistic regression
### Note this does not apply ipcw, and therefore will be incorrect in the presence of censoring
###
calc.calib.blr <- function(data.mstate, data.raw, t.eval, p.est){
  
  ## Let's do state 1 first
  ## Want to know which individual are in state j at time t
#           data.mstate <- data.mstate.reduc
#           data.raw <- data.raw.reduc
#           t.eval <- ceiling(7*365.25)
#           p.est <- p.true
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #     intersect(ids.state.list[[1]], ids.state.list[[3]])
  #   intersect(ids.state.list[[1]], ids.state.list[[4]])
  #   intersect(ids.state.list[[1]], ids.state.list[[5]])
  #   intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5),
           state1.bin = case_when(state.poly == 1 ~ 1,
                                  state.poly != 1 ~ 0),
           state2.bin = case_when(state.poly == 2 ~ 1,
                                  state.poly != 2 ~ 0),
           state3.bin = case_when(state.poly == 3 ~ 1,
                                  state.poly != 3 ~ 0),
           state4.bin = case_when(state.poly == 4 ~ 1,
                                  state.poly != 4 ~ 0),
           state5.bin = case_when(state.poly == 5 ~ 1,
                                  state.poly != 5 ~ 0))
  
  ### Add the predicted risks, and the logit transormation of the predicted risks to the dataset
  p.est.logit <- log(p.est/(1-p.est))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est, p.est.logit)
  
  ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
  ### (to ensure weights/offsets and individuals being modelled are consistent)
  data.raw.noNA <- data.raw %>% subset(!is.na(state.poly))
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1 <- glm(state1.bin ~ p.est.logit1, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm2 <- glm(state2.bin ~ p.est.logit2, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm3 <- glm(state3.bin ~ p.est.logit3, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm4 <- glm(state4.bin ~ p.est.logit4, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm5 <- glm(state5.bin ~ p.est.logit5, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  
  slopes <- c("state1" = as.numeric(coefficients(lrm1)[2]),
              "state2" = as.numeric(coefficients(lrm2)[2]),
              "state3" = as.numeric(coefficients(lrm3)[2]),
              "state4" = as.numeric(coefficients(lrm4)[2]),
              "state5" = as.numeric(coefficients(lrm5)[2]))
  
  slopes.se <- c("state1" = as.numeric(summary(lrm1)$coefficients[2,2]),
                 "state2" = as.numeric(summary(lrm2)$coefficients[2,2]),
                 "state3" = as.numeric(summary(lrm3)$coefficients[2,2]),
                 "state4" = as.numeric(summary(lrm4)$coefficients[2,2]),
                 "state5" = as.numeric(summary(lrm5)$coefficients[2,2]))
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1.offset <- glm(state1.bin ~ offset(p.est.logit1), data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm2.offset <- glm(state2.bin ~ offset(p.est.logit2), data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm3.offset <- glm(state3.bin ~ offset(p.est.logit3), data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm4.offset <- glm(state4.bin ~ offset(p.est.logit4), data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm5.offset <- glm(state5.bin ~ offset(p.est.logit5), data = data.raw.noNA, family = quasibinomial(link = "logit"))
  
  intercepts <- c("state1" = as.numeric(coefficients(lrm1.offset)),
                  "state2" = as.numeric(coefficients(lrm2.offset)),
                  "state3" = as.numeric(coefficients(lrm3.offset)),
                  "state4" = as.numeric(coefficients(lrm4.offset)),
                  "state5" = as.numeric(coefficients(lrm5.offset)))
  
  intercepts.se <- c("state1" = as.numeric(summary(lrm1.offset)$coefficients[2]),
                     "state2" = as.numeric(summary(lrm2.offset)$coefficients[2]),
                     "state3" = as.numeric(summary(lrm3.offset)$coefficients[2]),
                     "state4" = as.numeric(summary(lrm4.offset)$coefficients[2]),
                     "state5" = as.numeric(summary(lrm5.offset)$coefficients[2]))
  
  ###
  ### Create predicted observed risk for each individual
  ###
  pred.obs1 <- predict(lrm1, newdata = data.raw.noNA, type = "response")
  pred.obs2 <- predict(lrm2, newdata = data.raw.noNA, type = "response")
  pred.obs3 <- predict(lrm3, newdata = data.raw.noNA, type = "response")
  pred.obs4 <- predict(lrm4, newdata = data.raw.noNA, type = "response")
  pred.obs5 <- predict(lrm5, newdata = data.raw.noNA, type = "response")
  
  ### Get difference in risk
  diff.pred.obs1 <- mean(pred.obs1 - data.raw.noNA$p.est1)
  diff.pred.obs2 <- mean(pred.obs2 - data.raw.noNA$p.est2)
  diff.pred.obs3 <- mean(pred.obs3 - data.raw.noNA$p.est3)
  diff.pred.obs4 <- mean(pred.obs4 - data.raw.noNA$p.est4)
  diff.pred.obs5 <- mean(pred.obs5 - data.raw.noNA$p.est5)
  diff.pred.obs <- c(diff.pred.obs1, diff.pred.obs2, diff.pred.obs3, diff.pred.obs4, diff.pred.obs5)
  
  ### Create output object
  output.obj <- list("int" = intercepts, "int.se" = intercepts.se,
                     "slopes" = slopes, "slopes.se" = slopes.se,
                     "diff.pred.obs" = diff.pred.obs)
  return(output.obj)
  
}


###
### 3.4B) Assess calibration using binary logistic regression and IPCW weights
###
calc.calib.blr.ipcw <- function(data.mstate, data.raw, t.eval, p.est, max.weight = 10){
  
  ## Let's do state 1 first
  ## Want to know which individual are in state j at time t
#         data.mstate <- data.mstate.reduc
#         data.raw <- data.raw.reduc
#         t.eval <- ceiling(7*365.25)
#         p.est <- p.est.mean[[2]]

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #   Reduce(intersect, ids.state.list)
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5),
           state1.bin = case_when(state.poly == 1 ~ 1,
                                  state.poly != 1 ~ 0),
           state2.bin = case_when(state.poly == 2 ~ 1,
                                  state.poly != 2 ~ 0),
           state3.bin = case_when(state.poly == 3 ~ 1,
                                  state.poly != 3 ~ 0),
           state4.bin = case_when(state.poly == 4 ~ 1,
                                  state.poly != 4 ~ 0),
           state5.bin = case_when(state.poly == 5 ~ 1,
                                  state.poly != 5 ~ 0))
  
  ### Add the predicted risks, and the logit transormation of the predicted risks to the dataset
  p.est.logit <- log(p.est/(1-p.est))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est, p.est.logit)
  
  ### This data could be used to fit invalid logistic regression calibration models
  
  ###
  ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
  ###
  weights <- calc.weights(data.raw, t.eval, max.weight = max.weight)
  weights.DGM <- calc.weights.DGMspec(data.raw, t.eval, max.weight = max.weight, cens_shape = 1, cens_scale, cens_beta_x1, cens_beta_x2)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights, weights.DGM)
  
  ###
  ### Extract slopes and intercepts for miss-specified weighted models
  ###
  
  ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
  ### (to ensure weights/offsets and individuals being modelled are consistent)
  data.raw.noNA <- data.raw %>% subset(!is.na(state.poly))
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1.mspec <- glm(state1.bin ~ p.est.logit1, weights = ipcw.mspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm2.mspec <- glm(state2.bin ~ p.est.logit2, weights = ipcw.mspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm3.mspec <- glm(state3.bin ~ p.est.logit3, weights = ipcw.mspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm4.mspec <- glm(state4.bin ~ p.est.logit4, weights = ipcw.mspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm5.mspec <- glm(state5.bin ~ p.est.logit5, weights = ipcw.mspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  
  slopes.mspec <- c("state1" = as.numeric(coefficients(lrm1.mspec)[2]),
              "state2" = as.numeric(coefficients(lrm2.mspec)[2]),
              "state3" = as.numeric(coefficients(lrm3.mspec)[2]),
              "state4" = as.numeric(coefficients(lrm4.mspec)[2]),
              "state5" = as.numeric(coefficients(lrm5.mspec)[2]))
  
  slopes.mspec.se <- c("state1" = as.numeric(summary(lrm1.mspec)$coefficients[2,2]),
                 "state2" = as.numeric(summary(lrm2.mspec)$coefficients[2,2]),
                 "state3" = as.numeric(summary(lrm3.mspec)$coefficients[2,2]),
                 "state4" = as.numeric(summary(lrm4.mspec)$coefficients[2,2]),
                 "state5" = as.numeric(summary(lrm5.mspec)$coefficients[2,2]))
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1.mspec.offset <- glm(state1.bin ~ offset(p.est.logit1), weights = ipcw.mspec, data = data.raw.noNA, 
                     family = quasibinomial(link = "logit"))
  lrm2.mspec.offset <- glm(state2.bin ~ offset(p.est.logit2), weights = ipcw.mspec, data = data.raw.noNA, 
                     family = quasibinomial(link = "logit"))
  lrm3.mspec.offset <- glm(state3.bin ~ offset(p.est.logit3), weights = ipcw.mspec, data = data.raw.noNA, 
                     family = quasibinomial(link = "logit"))
  lrm4.mspec.offset <- glm(state4.bin ~ offset(p.est.logit4), weights = ipcw.mspec, data = data.raw.noNA, 
                     family = quasibinomial(link = "logit"))
  lrm5.mspec.offset <- glm(state5.bin ~ offset(p.est.logit5), weights = ipcw.mspec, data = data.raw.noNA, 
                     family = quasibinomial(link = "logit"))
  
  intercepts.mspec <- c("state1" = as.numeric(coefficients(lrm1.mspec.offset)),
                  "state2" = as.numeric(coefficients(lrm2.mspec.offset)),
                  "state3" = as.numeric(coefficients(lrm3.mspec.offset)),
                  "state4" = as.numeric(coefficients(lrm4.mspec.offset)),
                  "state5" = as.numeric(coefficients(lrm5.mspec.offset)))
  
  intercepts.mspec.se <- c("state1" = as.numeric(summary(lrm1.mspec.offset)$coefficients[2]),
                     "state2" = as.numeric(summary(lrm2.mspec.offset)$coefficients[2]),
                     "state3" = as.numeric(summary(lrm3.mspec.offset)$coefficients[2]),
                     "state4" = as.numeric(summary(lrm4.mspec.offset)$coefficients[2]),
                     "state5" = as.numeric(summary(lrm5.mspec.offset)$coefficients[2]))
  
  
  ###
  ### Extract slopes and intercepts for perfectly-specified weighted models
  ###
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1.pspec <- glm(state1.bin ~ p.est.logit1, weights = ipcw.pspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm2.pspec <- glm(state2.bin ~ p.est.logit2, weights = ipcw.pspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm3.pspec <- glm(state3.bin ~ p.est.logit3, weights = ipcw.pspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm4.pspec <- glm(state4.bin ~ p.est.logit4, weights = ipcw.pspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm5.pspec <- glm(state5.bin ~ p.est.logit5, weights = ipcw.pspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  
  slopes.pspec <- c("state1" = as.numeric(coefficients(lrm1.pspec)[2]),
                    "state2" = as.numeric(coefficients(lrm2.pspec)[2]),
                    "state3" = as.numeric(coefficients(lrm3.pspec)[2]),
                    "state4" = as.numeric(coefficients(lrm4.pspec)[2]),
                    "state5" = as.numeric(coefficients(lrm5.pspec)[2]))
  
  slopes.pspec.se <- c("state1" = as.numeric(summary(lrm1.pspec)$coefficients[2,2]),
                       "state2" = as.numeric(summary(lrm2.pspec)$coefficients[2,2]),
                       "state3" = as.numeric(summary(lrm3.pspec)$coefficients[2,2]),
                       "state4" = as.numeric(summary(lrm4.pspec)$coefficients[2,2]),
                       "state5" = as.numeric(summary(lrm5.pspec)$coefficients[2,2]))
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1.pspec.offset <- glm(state1.bin ~ offset(p.est.logit1), weights = ipcw.pspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  lrm2.pspec.offset <- glm(state2.bin ~ offset(p.est.logit2), weights = ipcw.pspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  lrm3.pspec.offset <- glm(state3.bin ~ offset(p.est.logit3), weights = ipcw.pspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  lrm4.pspec.offset <- glm(state4.bin ~ offset(p.est.logit4), weights = ipcw.pspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  lrm5.pspec.offset <- glm(state5.bin ~ offset(p.est.logit5), weights = ipcw.pspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  
  intercepts.pspec <- c("state1" = as.numeric(coefficients(lrm1.pspec.offset)),
                        "state2" = as.numeric(coefficients(lrm2.pspec.offset)),
                        "state3" = as.numeric(coefficients(lrm3.pspec.offset)),
                        "state4" = as.numeric(coefficients(lrm4.pspec.offset)),
                        "state5" = as.numeric(coefficients(lrm5.pspec.offset)))
  
  intercepts.pspec.se <- c("state1" = as.numeric(summary(lrm1.pspec.offset)$coefficients[2]),
                           "state2" = as.numeric(summary(lrm2.pspec.offset)$coefficients[2]),
                           "state3" = as.numeric(summary(lrm3.pspec.offset)$coefficients[2]),
                           "state4" = as.numeric(summary(lrm4.pspec.offset)$coefficients[2]),
                           "state5" = as.numeric(summary(lrm5.pspec.offset)$coefficients[2]))
  
  ###
  ### Extract slopes and intercepts for DGM-specified weighted models
  ###
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1.DGMspec <- glm(state1.bin ~ p.est.logit1, weights = ipcw.DGMspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm2.DGMspec <- glm(state2.bin ~ p.est.logit2, weights = ipcw.DGMspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm3.DGMspec <- glm(state3.bin ~ p.est.logit3, weights = ipcw.DGMspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm4.DGMspec <- glm(state4.bin ~ p.est.logit4, weights = ipcw.DGMspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  lrm5.DGMspec <- glm(state5.bin ~ p.est.logit5, weights = ipcw.DGMspec, data = data.raw.noNA, family = quasibinomial(link = "logit"))
  
  slopes.DGMspec <- c("state1" = as.numeric(coefficients(lrm1.DGMspec)[2]),
                    "state2" = as.numeric(coefficients(lrm2.DGMspec)[2]),
                    "state3" = as.numeric(coefficients(lrm3.DGMspec)[2]),
                    "state4" = as.numeric(coefficients(lrm4.DGMspec)[2]),
                    "state5" = as.numeric(coefficients(lrm5.DGMspec)[2]))
  
  slopes.DGMspec.se <- c("state1" = as.numeric(summary(lrm1.DGMspec)$coefficients[2,2]),
                       "state2" = as.numeric(summary(lrm2.DGMspec)$coefficients[2,2]),
                       "state3" = as.numeric(summary(lrm3.DGMspec)$coefficients[2,2]),
                       "state4" = as.numeric(summary(lrm4.DGMspec)$coefficients[2,2]),
                       "state5" = as.numeric(summary(lrm5.DGMspec)$coefficients[2,2]))
  
  ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
  lrm1.DGMspec.offset <- glm(state1.bin ~ offset(p.est.logit1), weights = ipcw.DGMspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  lrm2.DGMspec.offset <- glm(state2.bin ~ offset(p.est.logit2), weights = ipcw.DGMspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  lrm3.DGMspec.offset <- glm(state3.bin ~ offset(p.est.logit3), weights = ipcw.DGMspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  lrm4.DGMspec.offset <- glm(state4.bin ~ offset(p.est.logit4), weights = ipcw.DGMspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  lrm5.DGMspec.offset <- glm(state5.bin ~ offset(p.est.logit5), weights = ipcw.DGMspec, data = data.raw.noNA, 
                           family = quasibinomial(link = "logit"))
  
  intercepts.DGMspec <- c("state1" = as.numeric(coefficients(lrm1.DGMspec.offset)),
                        "state2" = as.numeric(coefficients(lrm2.DGMspec.offset)),
                        "state3" = as.numeric(coefficients(lrm3.DGMspec.offset)),
                        "state4" = as.numeric(coefficients(lrm4.DGMspec.offset)),
                        "state5" = as.numeric(coefficients(lrm5.DGMspec.offset)))
  
  intercepts.DGMspec.se <- c("state1" = as.numeric(summary(lrm1.DGMspec.offset)$coefficients[2]),
                           "state2" = as.numeric(summary(lrm2.DGMspec.offset)$coefficients[2]),
                           "state3" = as.numeric(summary(lrm3.DGMspec.offset)$coefficients[2]),
                           "state4" = as.numeric(summary(lrm4.DGMspec.offset)$coefficients[2]),
                           "state5" = as.numeric(summary(lrm5.DGMspec.offset)$coefficients[2]))
  
  ###
  ### Create predicted observed risk for each individual
  ###
  
  ###
  ### miss-specified model
  pred.obs.mspec1 <- predict(lrm1.mspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.mspec2 <- predict(lrm2.mspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.mspec3 <- predict(lrm3.mspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.mspec4 <- predict(lrm4.mspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.mspec5 <- predict(lrm5.mspec.offset, newdata = data.raw.noNA, type = "response")
  
  ### Get difference in risk
  diff.pred.obs.mspec1 <- mean(pred.obs.mspec1 - data.raw.noNA$p.est1)
  diff.pred.obs.mspec2 <- mean(pred.obs.mspec2 - data.raw.noNA$p.est2)
  diff.pred.obs.mspec3 <- mean(pred.obs.mspec3 - data.raw.noNA$p.est3)
  diff.pred.obs.mspec4 <- mean(pred.obs.mspec4 - data.raw.noNA$p.est4)
  diff.pred.obs.mspec5 <- mean(pred.obs.mspec5 - data.raw.noNA$p.est5)
  diff.pred.obs.mspec <- c(diff.pred.obs.mspec1, diff.pred.obs.mspec2, diff.pred.obs.mspec3, diff.pred.obs.mspec4, diff.pred.obs.mspec5)
  
  ###
  ### perfectly specified model
  pred.obs.pspec1 <- predict(lrm1.pspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.pspec2 <- predict(lrm2.pspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.pspec3 <- predict(lrm3.pspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.pspec4 <- predict(lrm4.pspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.pspec5 <- predict(lrm5.pspec.offset, newdata = data.raw.noNA, type = "response")
  
  ### Get difference in risk
  diff.pred.obs.pspec1 <- mean(pred.obs.pspec1 - data.raw.noNA$p.est1)
  diff.pred.obs.pspec2 <- mean(pred.obs.pspec2 - data.raw.noNA$p.est2)
  diff.pred.obs.pspec3 <- mean(pred.obs.pspec3 - data.raw.noNA$p.est3)
  diff.pred.obs.pspec4 <- mean(pred.obs.pspec4 - data.raw.noNA$p.est4)
  diff.pred.obs.pspec5 <- mean(pred.obs.pspec5 - data.raw.noNA$p.est5)
  diff.pred.obs.pspec <- c(diff.pred.obs.pspec1, diff.pred.obs.pspec2, diff.pred.obs.pspec3, diff.pred.obs.pspec4, diff.pred.obs.pspec5)
  
  ###
  ### DGM spec
  pred.obs.DGMspec1 <- predict(lrm1.DGMspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.DGMspec2 <- predict(lrm2.DGMspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.DGMspec3 <- predict(lrm3.DGMspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.DGMspec4 <- predict(lrm4.DGMspec.offset, newdata = data.raw.noNA, type = "response")
  pred.obs.DGMspec5 <- predict(lrm5.DGMspec.offset, newdata = data.raw.noNA, type = "response")
  
  ### Get difference in risk
  diff.pred.obs.DGMspec1 <- mean(pred.obs.DGMspec1 - data.raw.noNA$p.est1)
  diff.pred.obs.DGMspec2 <- mean(pred.obs.DGMspec2 - data.raw.noNA$p.est2)
  diff.pred.obs.DGMspec3 <- mean(pred.obs.DGMspec3 - data.raw.noNA$p.est3)
  diff.pred.obs.DGMspec4 <- mean(pred.obs.DGMspec4 - data.raw.noNA$p.est4)
  diff.pred.obs.DGMspec5 <- mean(pred.obs.DGMspec5 - data.raw.noNA$p.est5)
  diff.pred.obs.DGMspec <- c(diff.pred.obs.DGMspec1, diff.pred.obs.DGMspec2, diff.pred.obs.DGMspec3, diff.pred.obs.DGMspec4, diff.pred.obs.DGMspec5)
  
  ### Create output object
  output.obj <- list("int.mspec" = intercepts.mspec, "int.mspec.se" = intercepts.mspec.se,
                     "slopes.mspec" = slopes.mspec, "slopes.mspec.se" = slopes.mspec.se,
                     "int.pspec" = intercepts.pspec, "int.pspec.se" = intercepts.pspec.se,
                     "slopes.pspec" = slopes.pspec, "slopes.pspec.se" = slopes.pspec.se,
                     "int.DGMspec" = intercepts.DGMspec, "int.DGMspec.se" = intercepts.DGMspec.se,
                     "slopes.DGMspec" = slopes.DGMspec, "slopes.DGMspec.se" = slopes.DGMspec.se,
                     "diff.pred.obs.mspec" = diff.pred.obs.mspec,
                     "diff.pred.obs.pspec" = diff.pred.obs.pspec,
                     "diff.pred.obs.DGMspec" = diff.pred.obs.DGMspec)
  return(output.obj)
  
}


###
### 3.4C) Calculate standard error of calibration estimate using binary logistic regression and IPCW weights
###
calc.calib.blr.ipcw.boot.se <- function(data.mstate, data.raw, t.eval, p.est, n.boot, max.weight = 10){
  
#     data.mstate <- data.mstate.reduc
#     data.raw <- data.raw.reduc
#     t.eval <- ceiling(7*365.25)
#     p.est <- p.est.mean[[1]]
#     n.boot <- 3
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #   Reduce(intersect, ids.state.list)
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5),
           state1.bin = case_when(state.poly == 1 ~ 1,
                                  state.poly != 1 ~ 0),
           state2.bin = case_when(state.poly == 2 ~ 1,
                                  state.poly != 2 ~ 0),
           state3.bin = case_when(state.poly == 3 ~ 1,
                                  state.poly != 3 ~ 0),
           state4.bin = case_when(state.poly == 4 ~ 1,
                                  state.poly != 4 ~ 0),
           state5.bin = case_when(state.poly == 5 ~ 1,
                                  state.poly != 5 ~ 0))
  
  ### Add the predicted risks, and the logit transormation of the predicted risks to the dataset
  p.est.logit <- log(p.est/(1-p.est))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est, p.est.logit)
  
  ### Calc DGM weights and add to dataset
  weights.DGM <- calc.weights.DGMspec(data.raw, t.eval, max.weight = max.weight, cens_shape = 1, cens_scale, cens_beta_x1, cens_beta_x2)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights.DGM)
  
  ### Calculate and add weights from DGM spec
  ### Now write a function that will be put into a boot function to calculate bootstrapped standard errors for intercept and slopes
  get_int_slope_boot <- function(data, indices){

    ### Create bootstrap sample based on indices
    data.boot <- data[indices, ]
    
    ###
    ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
    ###
    weights <- calc.weights(data.boot, t.eval, max.weight = max.weight)
    
    ### Add these to dataset
    data.boot <- cbind(data.boot, weights)
    
    ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
    ### (to ensure weights/offsets and individuals being modelled are consistent)
    data.boot.noNA <- data.boot %>% subset(!is.na(state.poly))
    
    ###
    ### Extract slopes and intercepts for regular model
    ###
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1 <- glm(state1.bin ~ p.est.logit1, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm2 <- glm(state2.bin ~ p.est.logit2, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm3 <- glm(state3.bin ~ p.est.logit3, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm4 <- glm(state4.bin ~ p.est.logit4, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm5 <- glm(state5.bin ~ p.est.logit5, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    
    slopes <- c("slopes1" = as.numeric(coefficients(lrm1)[2]),
                "slopes2" = as.numeric(coefficients(lrm2)[2]),
                "slopes3" = as.numeric(coefficients(lrm3)[2]),
                "slopes4" = as.numeric(coefficients(lrm4)[2]),
                "slopes5" = as.numeric(coefficients(lrm5)[2]))
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1.offset <- glm(state1.bin ~ offset(p.est.logit1), data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm2.offset <- glm(state2.bin ~ offset(p.est.logit2), data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm3.offset <- glm(state3.bin ~ offset(p.est.logit3), data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm4.offset <- glm(state4.bin ~ offset(p.est.logit4), data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm5.offset <- glm(state5.bin ~ offset(p.est.logit5), data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    
    intercepts <- c("int1" = as.numeric(coefficients(lrm1.offset)),
                          "int2" = as.numeric(coefficients(lrm2.offset)),
                          "int3" = as.numeric(coefficients(lrm3.offset)),
                          "int4" = as.numeric(coefficients(lrm4.offset)),
                          "int5" = as.numeric(coefficients(lrm5.offset)))
    
    ###
    ### Extract slopes and intercepts for miss-specified weighted models
    ###
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1.mspec <- glm(state1.bin ~ p.est.logit1, weights = ipcw.mspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm2.mspec <- glm(state2.bin ~ p.est.logit2, weights = ipcw.mspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm3.mspec <- glm(state3.bin ~ p.est.logit3, weights = ipcw.mspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm4.mspec <- glm(state4.bin ~ p.est.logit4, weights = ipcw.mspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm5.mspec <- glm(state5.bin ~ p.est.logit5, weights = ipcw.mspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    
    slopes.mspec <- c("slopes.mspec1" = as.numeric(coefficients(lrm1.mspec)[2]),
                      "slopes.mspec2" = as.numeric(coefficients(lrm2.mspec)[2]),
                      "slopes.mspec3" = as.numeric(coefficients(lrm3.mspec)[2]),
                      "slopes.mspec4" = as.numeric(coefficients(lrm4.mspec)[2]),
                      "slopes.mspec5" = as.numeric(coefficients(lrm5.mspec)[2]))
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1.mspec.offset <- glm(state1.bin ~ offset(p.est.logit1), weights = ipcw.mspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm2.mspec.offset <- glm(state2.bin ~ offset(p.est.logit2), weights = ipcw.mspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm3.mspec.offset <- glm(state3.bin ~ offset(p.est.logit3), weights = ipcw.mspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm4.mspec.offset <- glm(state4.bin ~ offset(p.est.logit4), weights = ipcw.mspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm5.mspec.offset <- glm(state5.bin ~ offset(p.est.logit5), weights = ipcw.mspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    
    intercepts.mspec <- c("int.mspec1" = as.numeric(coefficients(lrm1.mspec.offset)),
                          "int.mspec2" = as.numeric(coefficients(lrm2.mspec.offset)),
                          "int.mspec3" = as.numeric(coefficients(lrm3.mspec.offset)),
                          "int.mspec4" = as.numeric(coefficients(lrm4.mspec.offset)),
                          "int.mspec5" = as.numeric(coefficients(lrm5.mspec.offset)))
    
    ###
    ### Extract slopes and intercepts for perfectly-specified weighted models
    ###
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1.pspec <- glm(state1.bin ~ p.est.logit1, weights = ipcw.pspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm2.pspec <- glm(state2.bin ~ p.est.logit2, weights = ipcw.pspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm3.pspec <- glm(state3.bin ~ p.est.logit3, weights = ipcw.pspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm4.pspec <- glm(state4.bin ~ p.est.logit4, weights = ipcw.pspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm5.pspec <- glm(state5.bin ~ p.est.logit5, weights = ipcw.pspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    
    slopes.pspec <- c("slopes.pspec1" = as.numeric(coefficients(lrm1.pspec)[2]),
                      "slopes.pspec2" = as.numeric(coefficients(lrm2.pspec)[2]),
                      "slopes.pspec3" = as.numeric(coefficients(lrm3.pspec)[2]),
                      "slopes.pspec4" = as.numeric(coefficients(lrm4.pspec)[2]),
                      "slopes.pspec5" = as.numeric(coefficients(lrm5.pspec)[2]))
    
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1.pspec.offset <- glm(state1.bin ~ offset(p.est.logit1), weights = ipcw.pspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm2.pspec.offset <- glm(state2.bin ~ offset(p.est.logit2), weights = ipcw.pspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm3.pspec.offset <- glm(state3.bin ~ offset(p.est.logit3), weights = ipcw.pspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm4.pspec.offset <- glm(state4.bin ~ offset(p.est.logit4), weights = ipcw.pspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm5.pspec.offset <- glm(state5.bin ~ offset(p.est.logit5), weights = ipcw.pspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    
    intercepts.pspec <- c("int.pspec1" = as.numeric(coefficients(lrm1.pspec.offset)),
                          "int.pspec2" = as.numeric(coefficients(lrm2.pspec.offset)),
                          "int.pspec3" = as.numeric(coefficients(lrm3.pspec.offset)),
                          "int.pspec4" = as.numeric(coefficients(lrm4.pspec.offset)),
                          "int.pspec5" = as.numeric(coefficients(lrm5.pspec.offset)))
    
    ###
    ### Extract slopes and intercepts for DGM-specified weighted models
    ###
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1.DGMspec <- glm(state1.bin ~ p.est.logit1, weights = ipcw.DGMspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm2.DGMspec <- glm(state2.bin ~ p.est.logit2, weights = ipcw.DGMspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm3.DGMspec <- glm(state3.bin ~ p.est.logit3, weights = ipcw.DGMspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm4.DGMspec <- glm(state4.bin ~ p.est.logit4, weights = ipcw.DGMspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    lrm5.DGMspec <- glm(state5.bin ~ p.est.logit5, weights = ipcw.DGMspec, data = data.boot.noNA, family = quasibinomial(link = "logit"))
    
    slopes.DGMspec <- c("slopes.DGMspec1" = as.numeric(coefficients(lrm1.DGMspec)[2]),
                      "slopes.DGMspec2" = as.numeric(coefficients(lrm2.DGMspec)[2]),
                      "slopes.DGMspec3" = as.numeric(coefficients(lrm3.DGMspec)[2]),
                      "slopes.DGMspec4" = as.numeric(coefficients(lrm4.DGMspec)[2]),
                      "slopes.DGMspec5" = as.numeric(coefficients(lrm5.DGMspec)[2]))
    
    
    ### Fit logistic regression recalibration models to the uncensored observations at time t to calculate intercepts
    lrm1.DGMspec.offset <- glm(state1.bin ~ offset(p.est.logit1), weights = ipcw.DGMspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm2.DGMspec.offset <- glm(state2.bin ~ offset(p.est.logit2), weights = ipcw.DGMspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm3.DGMspec.offset <- glm(state3.bin ~ offset(p.est.logit3), weights = ipcw.DGMspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm4.DGMspec.offset <- glm(state4.bin ~ offset(p.est.logit4), weights = ipcw.DGMspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    lrm5.DGMspec.offset <- glm(state5.bin ~ offset(p.est.logit5), weights = ipcw.DGMspec, data = data.boot.noNA, 
                             family = quasibinomial(link = "logit"))
    
    intercepts.DGMspec <- c("int.DGMspec1" = as.numeric(coefficients(lrm1.DGMspec.offset)),
                          "int.DGMspec2" = as.numeric(coefficients(lrm2.DGMspec.offset)),
                          "int.DGMspec3" = as.numeric(coefficients(lrm3.DGMspec.offset)),
                          "int.DGMspec4" = as.numeric(coefficients(lrm4.DGMspec.offset)),
                          "int.DGMspec5" = as.numeric(coefficients(lrm5.DGMspec.offset)))
    
    ###
    ### Create predicted observed risk for each individual
    ###
    
    ###
    ### no weights
    pred.obs1 <- predict(lrm1.offset, newdata = data.boot.noNA, type = "response")
    pred.obs2 <- predict(lrm2.offset, newdata = data.boot.noNA, type = "response")
    pred.obs3 <- predict(lrm3.offset, newdata = data.boot.noNA, type = "response")
    pred.obs4 <- predict(lrm4.offset, newdata = data.boot.noNA, type = "response")
    pred.obs5 <- predict(lrm5.offset, newdata = data.boot.noNA, type = "response")
    
    ### Get difference in risk
    diff.pred.obs1 <- mean(pred.obs1 - data.boot.noNA$p.est1)
    diff.pred.obs2 <- mean(pred.obs2 - data.boot.noNA$p.est2)
    diff.pred.obs3 <- mean(pred.obs3 - data.boot.noNA$p.est3)
    diff.pred.obs4 <- mean(pred.obs4 - data.boot.noNA$p.est4)
    diff.pred.obs5 <- mean(pred.obs5 - data.boot.noNA$p.est5)
    diff.pred.obs <- c(diff.pred.obs1, diff.pred.obs2, diff.pred.obs3, diff.pred.obs4, diff.pred.obs5)
    
    ###
    ### miss-specified model
    pred.obs.mspec1 <- predict(lrm1.mspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.mspec2 <- predict(lrm2.mspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.mspec3 <- predict(lrm3.mspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.mspec4 <- predict(lrm4.mspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.mspec5 <- predict(lrm5.mspec.offset, newdata = data.boot.noNA, type = "response")
    
    ### Get difference in risk
    diff.pred.obs.mspec1 <- mean(pred.obs.mspec1 - data.boot.noNA$p.est1)
    diff.pred.obs.mspec2 <- mean(pred.obs.mspec2 - data.boot.noNA$p.est2)
    diff.pred.obs.mspec3 <- mean(pred.obs.mspec3 - data.boot.noNA$p.est3)
    diff.pred.obs.mspec4 <- mean(pred.obs.mspec4 - data.boot.noNA$p.est4)
    diff.pred.obs.mspec5 <- mean(pred.obs.mspec5 - data.boot.noNA$p.est5)
    diff.pred.obs.mspec <- c(diff.pred.obs.mspec1, diff.pred.obs.mspec2, diff.pred.obs.mspec3, diff.pred.obs.mspec4, diff.pred.obs.mspec5)
    
    ###
    ### perfectly specified model
    pred.obs.pspec1 <- predict(lrm1.pspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.pspec2 <- predict(lrm2.pspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.pspec3 <- predict(lrm3.pspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.pspec4 <- predict(lrm4.pspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.pspec5 <- predict(lrm5.pspec.offset, newdata = data.boot.noNA, type = "response")
    
    ### Get difference in risk
    diff.pred.obs.pspec1 <- mean(pred.obs.pspec1 - data.boot.noNA$p.est1)
    diff.pred.obs.pspec2 <- mean(pred.obs.pspec2 - data.boot.noNA$p.est2)
    diff.pred.obs.pspec3 <- mean(pred.obs.pspec3 - data.boot.noNA$p.est3)
    diff.pred.obs.pspec4 <- mean(pred.obs.pspec4 - data.boot.noNA$p.est4)
    diff.pred.obs.pspec5 <- mean(pred.obs.pspec5 - data.boot.noNA$p.est5)
    diff.pred.obs.pspec <- c(diff.pred.obs.pspec1, diff.pred.obs.pspec2, diff.pred.obs.pspec3, diff.pred.obs.pspec4, diff.pred.obs.pspec5)
    
    ###
    ### DGM spec
    pred.obs.DGMspec1 <- predict(lrm1.DGMspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.DGMspec2 <- predict(lrm2.DGMspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.DGMspec3 <- predict(lrm3.DGMspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.DGMspec4 <- predict(lrm4.DGMspec.offset, newdata = data.boot.noNA, type = "response")
    pred.obs.DGMspec5 <- predict(lrm5.DGMspec.offset, newdata = data.boot.noNA, type = "response")
    
    ### Get difference in risk
    diff.pred.obs.DGMspec1 <- mean(pred.obs.DGMspec1 - data.boot.noNA$p.est1)
    diff.pred.obs.DGMspec2 <- mean(pred.obs.DGMspec2 - data.boot.noNA$p.est2)
    diff.pred.obs.DGMspec3 <- mean(pred.obs.DGMspec3 - data.boot.noNA$p.est3)
    diff.pred.obs.DGMspec4 <- mean(pred.obs.DGMspec4 - data.boot.noNA$p.est4)
    diff.pred.obs.DGMspec5 <- mean(pred.obs.DGMspec5 - data.boot.noNA$p.est5)
    diff.pred.obs.DGMspec <- c(diff.pred.obs.DGMspec1, diff.pred.obs.DGMspec2, diff.pred.obs.DGMspec3, diff.pred.obs.DGMspec4, diff.pred.obs.DGMspec5)
    
    ### create and return output object
    output.obj <- c(intercepts, intercepts.mspec, intercepts.pspec, intercepts.DGMspec, slopes, slopes.mspec, slopes.pspec, slopes.DGMspec,
                    diff.pred.obs, diff.pred.obs.mspec, diff.pred.obs.pspec, diff.pred.obs.DGMspec)
    return(output.obj)
    
  }  
  
  ### Run the bootstrapping
  boot.obj <- boot(data.raw, statistic = get_int_slope_boot, R = n.boot)
  colnames(boot.obj$t) <- c(paste("int", 1:5, sep = ""), paste("int.mspec", 1:5, sep = ""), 
                          paste("int.pspec", 1:5, sep = ""), paste("int.DGMspec", 1:5, sep = ""), 
                          paste("slope", 1:5, sep = ""), paste("slope.mspec", 1:5, sep = ""), 
                          paste("slope.pspec", 1:5, sep = ""), paste("slope.DGMspec", 1:5, sep = ""),
                          paste("diff.pred.obs", 1:5, sep = ""), paste("diff.pred.obs.mspec", 1:5, sep = ""), 
                          paste("diff.pred.obs.pspec", 1:5, sep = ""), paste("diff.pred.obs.DGMspec", 1:5, sep = ""))  
  
  ### Calculate standard errors
  se <- sqrt(apply(boot.obj$t, 2, var))
  names(se) <- c(paste("int", 1:5, sep = ""), paste("int.mspec", 1:5, sep = ""), 
                 paste("int.pspec", 1:5, sep = ""), paste("int.DGMspec", 1:5, sep = ""), 
                 paste("slope", 1:5, sep = ""), paste("slope.mspec", 1:5, sep = ""), 
                 paste("slope.pspec", 1:5, sep = ""), paste("slope.DGMspec", 1:5, sep = ""),
                 paste("diff.pred.obs", 1:5, sep = ""), paste("diff.pred.obs.mspec", 1:5, sep = ""), 
                 paste("diff.pred.obs.pspec", 1:5, sep = ""), paste("diff.pred.obs.DGMspec", 1:5, sep = ""))
  
  ### Return output
  output.object <- list("boot.obj" = boot.obj$t, "se" = se)
  return(output.object)
  
}



###
### 3.5A) Function to calculate calibration after applying the nominal calibration framework of van Hoorde et al,
### Note this does not apply ipcw, and therefore will be incorrect in the presence of censoring
###
calc.calib.mlr <- function(data.mstate, data.raw, t.eval, p.est){
  
  ## Let's do state 1 first
  ## Want to know which individual are in state j at time t
#     data.mstate <- data.mstate.reduc
#     data.raw <- data.raw.reduc
#     t.eval <- ceiling(7*365.25)
#     p.est <- p.est.mean[[2]]
#   str(summary(lrm1.DGMspec.offset))
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #   Reduce(intersect, ids.state.list)
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5))
  
  ### Add p.est to dataset
  data.raw <- data.frame(data.raw, p.est)
  
  ### Add linear predictors from a multinomial framework
  data.raw <- data.raw %>%
    mutate(mlr.lp1 = log(p.est2/p.est1),
           mlr.lp2 = log(p.est3/p.est1),
           mlr.lp3 = log(p.est4/p.est1),
           mlr.lp4 = log(p.est5/p.est1),
           state.poly.fac = as.factor(state.poly))
  
  
  ## Add constraints
  i <- diag(4)
  i1 <- rbind(1, 0, 0, 0)
  i2 <- rbind(0, 1, 0, 0)
  i3 <- rbind(0, 0, 1, 0)
  i4 <- rbind(0, 0, 0, 1)
  clist <- list("(Intercept)" = i, "mlr.lp1" = i1, "mlr.lp2" = i2, "mlr.lp3" = i3, "mlr.lp4" = i4)
  clist
  
  ### Apply nominal recalibration framework
  
  ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
  ### (to ensure weights/offsets and individuals being modelled are consistent)
  data.raw.noNA <- data.raw %>% subset(!is.na(state.poly))
  
  ### Apply models
  calib.model <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4, constraints = clist, 
                      data = data.raw.noNA, family = multinomial(refLevel = "1"))
  
  calib.model.offset <- vgam(data.raw.noNA$state.poly.fac ~ 1, offset = as.matrix(data.raw.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4")]), 
                             family = multinomial(refLevel = "1"))
  
  ###
  ### Create predicted observed risk for each individual
  ###
  pred.obs <- data.frame("mlr.lp1" = coefficients(calib.model.offset)[1] + data.raw.noNA$mlr.lp1,
                         "mlr.lp2" = coefficients(calib.model.offset)[2] + data.raw.noNA$mlr.lp2,
                         "mlr.lp3" = coefficients(calib.model.offset)[3] + data.raw.noNA$mlr.lp3,
                         "mlr.lp4" = coefficients(calib.model.offset)[4] + data.raw.noNA$mlr.lp4)
  pred.obs <- mutate(pred.obs, 
                     p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                     p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                     p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                     p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                     p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)))
  
  ### Get difference in risk
  diff.pred.obs1 <- mean(pred.obs[,"p1"] - data.raw.noNA$p.est1)
  diff.pred.obs2 <- mean(pred.obs[,"p2"] - data.raw.noNA$p.est2)
  diff.pred.obs3 <- mean(pred.obs[,"p3"] - data.raw.noNA$p.est3)
  diff.pred.obs4 <- mean(pred.obs[,"p4"] - data.raw.noNA$p.est4)
  diff.pred.obs5 <- mean(pred.obs[,"p5"] - data.raw.noNA$p.est5)
  diff.pred.obs <- c(diff.pred.obs1, diff.pred.obs2, diff.pred.obs3, diff.pred.obs4, diff.pred.obs5)
  
  ### Create output object
  output.object <- list("int" = calib.model.offset@coefficients, 
                        "int.se" = sqrt(diag(vcov(calib.model.offset))), 
                        "slopes" = calib.model@coefficients[paste("mlr.lp", 1:4, sep = "")], 
                        "slopes.se" = sqrt(diag(vcov(calib.model))[paste("mlr.lp", 1:4, sep = "")]),
                        "diff.pred.obs" = diff.pred.obs)
  return(output.object)
  
}

###
### 3.5B) Function to calculate calibration after applying the nominal calibration framework of van Hoorde et al,
### weighted using IPCW weights
###
calc.calib.mlr.ipcw <- function(data.mstate, data.raw, t.eval, p.est, max.weight = 10){
  
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
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #   Reduce(intersect, ids.state.list)
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5))
  
  ### Add p.est to dataset
  data.raw <- data.frame(data.raw, p.est)

  ### Add linear predictors from a multinomial framework
  data.raw <- data.raw %>%
    mutate(mlr.lp1 = log(p.est2/p.est1),
           mlr.lp2 = log(p.est3/p.est1),
           mlr.lp3 = log(p.est4/p.est1),
           mlr.lp4 = log(p.est5/p.est1),
           state.poly.fac = as.factor(state.poly))
  
  ###
  ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
  ###
  weights <- calc.weights(data.raw, t.eval, max.weight = max.weight)
  weights.DGM <- calc.weights.DGMspec(data.raw, t.eval, max.weight = max.weight, cens_shape = 1, cens_scale, cens_beta_x1, cens_beta_x2)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights, weights.DGM)
  
  ###
  ### Now fit the nominal calibration framework using the weights
  ###
  
  ## Add constraints
  i <- diag(4)
  i1 <- rbind(1, 0, 0, 0)
  i2 <- rbind(0, 1, 0, 0)
  i3 <- rbind(0, 0, 1, 0)
  i4 <- rbind(0, 0, 0, 1)
  clist <- list("(Intercept)" = i, "mlr.lp1" = i1, "mlr.lp2" = i2, "mlr.lp3" = i3, "mlr.lp4" = i4)
  clist
  
  ### Apply nominal recalibration framework
  ###
  
  ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
  ### (to ensure weights and individuals being modelled are consistent)
  data.raw.noNA <- data.raw %>% subset(!is.na(state.poly))
  
  ### Miss-specified weights
  calib.model.mspec <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4, constraints = clist, weights = ipcw.mspec,
                      data = data.raw.noNA, family = multinomial(refLevel = "1"))
  
  calib.model.mspec.offset <- vgam(data.raw.noNA$state.poly.fac ~ 1, offset = as.matrix(data.raw.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4")]),
                             weights = data.raw.noNA$ipcw.mspec,
                             family = multinomial(refLevel = "1"))
  
  ### Perfectly specified weights
  calib.model.pspec <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4, constraints = clist, weights = ipcw.pspec,
                            data = data.raw.noNA, family = multinomial(refLevel = "1"))
  
  calib.model.pspec.offset <- vgam(data.raw.noNA$state.poly.fac ~ 1, offset = as.matrix(data.raw.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4")]),
                                   weights = data.raw.noNA$ipcw.pspec,
                                   family = multinomial(refLevel = "1"))
  
  
  ### DGM specified weights
  calib.model.DGMspec <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4, constraints = clist, weights = ipcw.DGMspec,
                            data = data.raw.noNA, family = multinomial(refLevel = "1"))
  
  calib.model.DGMspec.offset <- vgam(data.raw.noNA$state.poly.fac ~ 1, offset = as.matrix(data.raw.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4")]),
                                   weights = data.raw.noNA$ipcw.DGMspec,
                                   family = multinomial(refLevel = "1"))
  
  ###
  ### Create predicted observed risk for each individual
  ###
  
  ###
  ### miss-specified model spec
  pred.obs.mspec <- data.frame("mlr.lp1" = coefficients(calib.model.mspec.offset)[1] + data.raw.noNA$mlr.lp1,
                                 "mlr.lp2" = coefficients(calib.model.mspec.offset)[2] + data.raw.noNA$mlr.lp2,
                                 "mlr.lp3" = coefficients(calib.model.mspec.offset)[3] + data.raw.noNA$mlr.lp3,
                                 "mlr.lp4" = coefficients(calib.model.mspec.offset)[4] + data.raw.noNA$mlr.lp4)
  pred.obs.mspec <- mutate(pred.obs.mspec, 
                             p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)))
  
  ### Get difference in risk
  diff.pred.obs.mspec1 <- mean(pred.obs.mspec[,"p1"] - data.raw.noNA$p.est1)
  diff.pred.obs.mspec2 <- mean(pred.obs.mspec[,"p2"] - data.raw.noNA$p.est2)
  diff.pred.obs.mspec3 <- mean(pred.obs.mspec[,"p3"] - data.raw.noNA$p.est3)
  diff.pred.obs.mspec4 <- mean(pred.obs.mspec[,"p4"] - data.raw.noNA$p.est4)
  diff.pred.obs.mspec5 <- mean(pred.obs.mspec[,"p5"] - data.raw.noNA$p.est5)
  diff.pred.obs.mspec <- c(diff.pred.obs.mspec1, diff.pred.obs.mspec2, diff.pred.obs.mspec3, diff.pred.obs.mspec4, diff.pred.obs.mspec5)
  
  ###
  ### Perfectly specified model 
  pred.obs.pspec <- data.frame("mlr.lp1" = coefficients(calib.model.pspec.offset)[1] + data.raw.noNA$mlr.lp1,
                                 "mlr.lp2" = coefficients(calib.model.pspec.offset)[2] + data.raw.noNA$mlr.lp2,
                                 "mlr.lp3" = coefficients(calib.model.pspec.offset)[3] + data.raw.noNA$mlr.lp3,
                                 "mlr.lp4" = coefficients(calib.model.pspec.offset)[4] + data.raw.noNA$mlr.lp4)
  pred.obs.pspec <- mutate(pred.obs.pspec, 
                             p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)))
  
  ### Get difference in risk
  diff.pred.obs.pspec1 <- mean(pred.obs.pspec[,"p1"] - data.raw.noNA$p.est1)
  diff.pred.obs.pspec2 <- mean(pred.obs.pspec[,"p2"] - data.raw.noNA$p.est2)
  diff.pred.obs.pspec3 <- mean(pred.obs.pspec[,"p3"] - data.raw.noNA$p.est3)
  diff.pred.obs.pspec4 <- mean(pred.obs.pspec[,"p4"] - data.raw.noNA$p.est4)
  diff.pred.obs.pspec5 <- mean(pred.obs.pspec[,"p5"] - data.raw.noNA$p.est5)
  diff.pred.obs.pspec <- c(diff.pred.obs.pspec1, diff.pred.obs.pspec2, diff.pred.obs.pspec3, diff.pred.obs.pspec4, diff.pred.obs.pspec5)
  
  ###
  ### DGM spec
  pred.obs.DGMspec <- data.frame("mlr.lp1" = coefficients(calib.model.DGMspec.offset)[1] + data.raw.noNA$mlr.lp1,
                         "mlr.lp2" = coefficients(calib.model.DGMspec.offset)[2] + data.raw.noNA$mlr.lp2,
                         "mlr.lp3" = coefficients(calib.model.DGMspec.offset)[3] + data.raw.noNA$mlr.lp3,
                         "mlr.lp4" = coefficients(calib.model.DGMspec.offset)[4] + data.raw.noNA$mlr.lp4)
  pred.obs.DGMspec <- mutate(pred.obs.DGMspec, 
                     p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                     p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                     p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                     p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                     p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)))

  ### Get difference in risk
  diff.pred.obs.DGMspec1 <- mean(pred.obs.DGMspec[,"p1"] - data.raw.noNA$p.est1)
  diff.pred.obs.DGMspec2 <- mean(pred.obs.DGMspec[,"p2"] - data.raw.noNA$p.est2)
  diff.pred.obs.DGMspec3 <- mean(pred.obs.DGMspec[,"p3"] - data.raw.noNA$p.est3)
  diff.pred.obs.DGMspec4 <- mean(pred.obs.DGMspec[,"p4"] - data.raw.noNA$p.est4)
  diff.pred.obs.DGMspec5 <- mean(pred.obs.DGMspec[,"p5"] - data.raw.noNA$p.est5)
  diff.pred.obs.DGMspec <- c(diff.pred.obs.DGMspec1, diff.pred.obs.DGMspec2, diff.pred.obs.DGMspec3, diff.pred.obs.DGMspec4, diff.pred.obs.DGMspec5)
  
  ### Create output object
  output.object <- list("int.mspec" = calib.model.mspec.offset@coefficients, 
                        "int.mspec.se" = sqrt(diag(vcov(calib.model.mspec.offset))), 
                        "slopes.mspec" = calib.model.mspec@coefficients[paste("mlr.lp", 1:4, sep = "")], 
                        "slopes.mspec.se" = sqrt(diag(vcov(calib.model.mspec))[paste("mlr.lp", 1:4, sep = "")]),
                        "int.pspec" = calib.model.pspec.offset@coefficients, 
                        "int.pspec.se" = sqrt(diag(vcov(calib.model.pspec.offset))), 
                        "slopes.pspec" = calib.model.pspec@coefficients[paste("mlr.lp", 1:4, sep = "")], 
                        "slopes.pspec.se" = sqrt(diag(vcov(calib.model.pspec))[paste("mlr.lp", 1:4, sep = "")]),
                        "int.DGMspec" = calib.model.DGMspec.offset@coefficients, 
                        "int.DGMspec.se" = sqrt(diag(vcov(calib.model.DGMspec.offset))), 
                        "slopes.DGMspec" = calib.model.DGMspec@coefficients[paste("mlr.lp", 1:4, sep = "")], 
                        "slopes.DGMspec.se" = sqrt(diag(vcov(calib.model.DGMspec))[paste("mlr.lp", 1:4, sep = "")]),
                        "diff.pred.obs.mspec" = diff.pred.obs.mspec,
                        "diff.pred.obs.pspec" = diff.pred.obs.pspec,
                        "diff.pred.obs.DGMspec" = diff.pred.obs.DGMspec)
  
  return(output.object)
  
}


###
### 3.5C) Calculate standard error of calibration estimate using multinomial nominal recalibration framework and IPCW weights
###
calc.calib.mlr.ipcw.boot.se <- function(data.mstate, data.raw, t.eval, p.est, n.boot, max.weight = 10){
  
  #     data.mstate <- data.mstate.reduc
  #     data.raw <- data.raw.reduc
  #     t.eval <- ceiling(7*365.25)
  #     p.est <- p.true
  #     n.boot <- 3
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  #   ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #   Reduce(intersect, ids.state.list)
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5))
  
  ### Add p.est to dataset
  data.raw <- data.frame(data.raw, p.est)
  
  ### Add linear predictors from a multinomial framework
  data.raw <- data.raw %>%
    mutate(mlr.lp1 = log(p.est2/p.est1),
           mlr.lp2 = log(p.est3/p.est1),
           mlr.lp3 = log(p.est4/p.est1),
           mlr.lp4 = log(p.est5/p.est1),
           state.poly.fac = as.factor(state.poly))
  
  ### Calc DGM weights and add to dataset
  weights.DGM <- calc.weights.DGMspec(data.raw, t.eval, max.weight = max.weight, cens_shape = 1, cens_scale, cens_beta_x1, cens_beta_x2)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights.DGM)
  
  ### Now write a function that will be put into a boot function to calculate bootstrapped standard errors for intercept and slopes
  get_int_slope_boot <- function(data, indices){
    
    ### Create bootstrap sample based on indices
    data.boot <- data[indices, ]
    
    ###
    ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
    ###
    weights <- calc.weights(data.boot, t.eval, max.weight = max.weight)
    
    ### Add these to dataset
    data.boot <- cbind(data.boot, weights)
    
    ###
    ### Now fit the nominal calibration framework using the weights
    ###
    
    ## Add constraints
    i <- diag(4)
    i1 <- rbind(1, 0, 0, 0)
    i2 <- rbind(0, 1, 0, 0)
    i3 <- rbind(0, 0, 1, 0)
    i4 <- rbind(0, 0, 0, 1)
    clist <- list("(Intercept)" = i, "mlr.lp1" = i1, "mlr.lp2" = i2, "mlr.lp3" = i3, "mlr.lp4" = i4)
    clist
    
    ### Apply nominal recalibration framework
    ###
    
    ### Create a dataset which is only individuals who don't have NA values for state.poly.fac 
    ### (to ensure weights and individuals being modelled are consistent)
    data.boot.noNA <- data.boot %>% subset(!is.na(state.poly))
    
    ### No weights
    calib.model <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4, constraints = clist,
                              data = data.boot.noNA, family = multinomial(refLevel = "1"))
    
    calib.model.offset <- vgam(data.boot.noNA$state.poly.fac ~ 1, offset = as.matrix(data.boot.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4")]),
                                     family = multinomial(refLevel = "1"))
    
    ### Miss-specified weights
    calib.model.mspec <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4, constraints = clist, weights = ipcw.mspec,
                              data = data.boot.noNA, family = multinomial(refLevel = "1"))
    
    calib.model.mspec.offset <- vgam(data.boot.noNA$state.poly.fac ~ 1, offset = as.matrix(data.boot.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4")]),
                                     weights = data.boot.noNA$ipcw.mspec,
                                     family = multinomial(refLevel = "1"))
    
    ### Perfectly specified weights
    calib.model.pspec <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4, constraints = clist, weights = ipcw.pspec,
                              data = data.boot.noNA, family = multinomial(refLevel = "1"))
    
    calib.model.pspec.offset <- vgam(data.boot.noNA$state.poly.fac ~ 1, offset = as.matrix(data.boot.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4")]),
                                     weights = data.boot.noNA$ipcw.pspec,
                                     family = multinomial(refLevel = "1"))
    
    ### DGM specified weights
    calib.model.DGMspec <- vgam(state.poly.fac ~ mlr.lp1 + mlr.lp2 + mlr.lp3 + mlr.lp4, constraints = clist, weights = ipcw.DGMspec,
                              data = data.boot.noNA, family = multinomial(refLevel = "1"))
    
    calib.model.DGMspec.offset <- vgam(data.boot.noNA$state.poly.fac ~ 1, offset = as.matrix(data.boot.noNA[, c("mlr.lp1", "mlr.lp2", "mlr.lp3", "mlr.lp4")]),
                                     weights = data.boot.noNA$ipcw.DGMspec,
                                     family = multinomial(refLevel = "1"))
    
    ###
    ### Create predicted observed risk for each individual
    ###
    
    ###
    ### no weights
    pred.obs <- data.frame("mlr.lp1" = coefficients(calib.model.offset)[1] + data.boot.noNA$mlr.lp1,
                                 "mlr.lp2" = coefficients(calib.model.offset)[2] + data.boot.noNA$mlr.lp2,
                                 "mlr.lp3" = coefficients(calib.model.offset)[3] + data.boot.noNA$mlr.lp3,
                                 "mlr.lp4" = coefficients(calib.model.offset)[4] + data.boot.noNA$mlr.lp4)
    pred.obs <- mutate(pred.obs, 
                             p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)))
    
    ### Get difference in risk
    diff.pred.obs1 <- mean(pred.obs[,"p1"] - data.boot.noNA$p.est1)
    diff.pred.obs2 <- mean(pred.obs[,"p2"] - data.boot.noNA$p.est2)
    diff.pred.obs3 <- mean(pred.obs[,"p3"] - data.boot.noNA$p.est3)
    diff.pred.obs4 <- mean(pred.obs[,"p4"] - data.boot.noNA$p.est4)
    diff.pred.obs5 <- mean(pred.obs[,"p5"] - data.boot.noNA$p.est5)
    diff.pred.obs <- c(diff.pred.obs1, diff.pred.obs2, diff.pred.obs3, diff.pred.obs4, diff.pred.obs5)
    
    ###
    ### miss-specified model spec
    pred.obs.mspec <- data.frame("mlr.lp1" = coefficients(calib.model.mspec.offset)[1] + data.boot.noNA$mlr.lp1,
                                 "mlr.lp2" = coefficients(calib.model.mspec.offset)[2] + data.boot.noNA$mlr.lp2,
                                 "mlr.lp3" = coefficients(calib.model.mspec.offset)[3] + data.boot.noNA$mlr.lp3,
                                 "mlr.lp4" = coefficients(calib.model.mspec.offset)[4] + data.boot.noNA$mlr.lp4)
    pred.obs.mspec <- mutate(pred.obs.mspec, 
                             p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)))
    
    ### Get difference in risk
    diff.pred.obs.mspec1 <- mean(pred.obs.mspec[,"p1"] - data.boot.noNA$p.est1)
    diff.pred.obs.mspec2 <- mean(pred.obs.mspec[,"p2"] - data.boot.noNA$p.est2)
    diff.pred.obs.mspec3 <- mean(pred.obs.mspec[,"p3"] - data.boot.noNA$p.est3)
    diff.pred.obs.mspec4 <- mean(pred.obs.mspec[,"p4"] - data.boot.noNA$p.est4)
    diff.pred.obs.mspec5 <- mean(pred.obs.mspec[,"p5"] - data.boot.noNA$p.est5)
    diff.pred.obs.mspec <- c(diff.pred.obs.mspec1, diff.pred.obs.mspec2, diff.pred.obs.mspec3, diff.pred.obs.mspec4, diff.pred.obs.mspec5)
    
    ###
    ### Perfectly specified model 
    pred.obs.pspec <- data.frame("mlr.lp1" = coefficients(calib.model.pspec.offset)[1] + data.boot.noNA$mlr.lp1,
                                 "mlr.lp2" = coefficients(calib.model.pspec.offset)[2] + data.boot.noNA$mlr.lp2,
                                 "mlr.lp3" = coefficients(calib.model.pspec.offset)[3] + data.boot.noNA$mlr.lp3,
                                 "mlr.lp4" = coefficients(calib.model.pspec.offset)[4] + data.boot.noNA$mlr.lp4)
    pred.obs.pspec <- mutate(pred.obs.pspec, 
                             p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                             p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)))
    
    ### Get difference in risk
    diff.pred.obs.pspec1 <- mean(pred.obs.pspec[,"p1"] - data.boot.noNA$p.est1)
    diff.pred.obs.pspec2 <- mean(pred.obs.pspec[,"p2"] - data.boot.noNA$p.est2)
    diff.pred.obs.pspec3 <- mean(pred.obs.pspec[,"p3"] - data.boot.noNA$p.est3)
    diff.pred.obs.pspec4 <- mean(pred.obs.pspec[,"p4"] - data.boot.noNA$p.est4)
    diff.pred.obs.pspec5 <- mean(pred.obs.pspec[,"p5"] - data.boot.noNA$p.est5)
    diff.pred.obs.pspec <- c(diff.pred.obs.pspec1, diff.pred.obs.pspec2, diff.pred.obs.pspec3, diff.pred.obs.pspec4, diff.pred.obs.pspec5)
    
    ###
    ### DGM spec
    pred.obs.DGMspec <- data.frame("mlr.lp1" = coefficients(calib.model.DGMspec.offset)[1] + data.boot.noNA$mlr.lp1,
                                   "mlr.lp2" = coefficients(calib.model.DGMspec.offset)[2] + data.boot.noNA$mlr.lp2,
                                   "mlr.lp3" = coefficients(calib.model.DGMspec.offset)[3] + data.boot.noNA$mlr.lp3,
                                   "mlr.lp4" = coefficients(calib.model.DGMspec.offset)[4] + data.boot.noNA$mlr.lp4)
    pred.obs.DGMspec <- mutate(pred.obs.DGMspec, 
                               p1 = 1/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                               p2 = exp(mlr.lp1)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                               p3 = exp(mlr.lp2)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                               p4 = exp(mlr.lp3)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)),
                               p5 = exp(mlr.lp4)/(1 + exp(mlr.lp1) + exp(mlr.lp2) + exp(mlr.lp3) + exp(mlr.lp4)))
    
    ### Get difference in risk
    diff.pred.obs.DGMspec1 <- mean(pred.obs.DGMspec[,"p1"] - data.boot.noNA$p.est1)
    diff.pred.obs.DGMspec2 <- mean(pred.obs.DGMspec[,"p2"] - data.boot.noNA$p.est2)
    diff.pred.obs.DGMspec3 <- mean(pred.obs.DGMspec[,"p3"] - data.boot.noNA$p.est3)
    diff.pred.obs.DGMspec4 <- mean(pred.obs.DGMspec[,"p4"] - data.boot.noNA$p.est4)
    diff.pred.obs.DGMspec5 <- mean(pred.obs.DGMspec[,"p5"] - data.boot.noNA$p.est5)
    diff.pred.obs.DGMspec <- c(diff.pred.obs.DGMspec1, diff.pred.obs.DGMspec2, diff.pred.obs.DGMspec3, diff.pred.obs.DGMspec4, diff.pred.obs.DGMspec5)
    
    
    ### Create and return output object
    output.obj <- c(calib.model.offset@coefficients,
                    calib.model.mspec.offset@coefficients, 
                    calib.model.pspec.offset@coefficients,
                    calib.model.DGMspec.offset@coefficients,
                    calib.model@coefficients[paste("mlr.lp", 1:4, sep = "")],
                    calib.model.mspec@coefficients[paste("mlr.lp", 1:4, sep = "")],
                    calib.model.pspec@coefficients[paste("mlr.lp", 1:4, sep = "")],
                    calib.model.DGMspec@coefficients[paste("mlr.lp", 1:4, sep = "")],
                    diff.pred.obs,
                    diff.pred.obs.mspec,
                    diff.pred.obs.pspec,
                    diff.pred.obs.DGMspec)
    
    return(output.obj)
    
  }  
  
  
  ### Run the bootstrapping
  boot.obj <- boot(data.raw, statistic = get_int_slope_boot, R = n.boot)
  colnames(boot.obj$t) <- c(paste("int", 2:5, sep = ""), paste("int.mspec", 2:5, sep = ""), 
                            paste("int.pspec", 2:5, sep = ""), paste("int.DGMspec", 2:5, sep = ""),
                            paste("slope", 2:5, sep = ""), paste("slope.mspec", 2:5, sep = ""), 
                            paste("slope.pspec", 2:5, sep = ""), paste("slope.DGMspec", 2:5, sep = ""),
                            paste("diff.pred.obs", 1:5, sep = ""), paste("diff.pred.obs.mspec", 1:5, sep = ""), 
                            paste("diff.pred.obs.pspec", 1:5, sep = ""), paste("diff.pred.obs.DGMspec", 1:5, sep = ""))
    
    ### Calculate standard errors
    se <- sqrt(apply(boot.obj$t, 2, var))
  names(se) <- c(paste("int", 2:5, sep = ""), paste("int.mspec", 2:5, sep = ""), 
                 paste("int.pspec", 2:5, sep = ""), paste("int.DGMspec", 2:5, sep = ""),
                 paste("slope", 2:5, sep = ""), paste("slope.mspec", 2:5, sep = ""), 
                 paste("slope.pspec", 2:5, sep = ""), paste("slope.DGMspec", 2:5, sep = ""),
                 paste("diff.pred.obs", 1:5, sep = ""), paste("diff.pred.obs.mspec", 1:5, sep = ""), 
                 paste("diff.pred.obs.pspec", 1:5, sep = ""), paste("diff.pred.obs.DGMspec", 1:5, sep = ""))
  
  ### Return output
  output.object <- list("boot.obj" = boot.obj$t, "se" = se)
  return(output.object)
  
}


###
### 3.6A) Define a function to calculate pseudo-value for an individual, using the Aalen-Johansen estimator
### 
### obs.aj is pre-calculated Aalen-Johansen estimator in entire cohort data.mstate
func.calc.pv.aj <- function(patid.eval, data.mstate, obs.aj, tmat, n.cohort, t.eval){
  
  ### Calculate AJ estimate without patient in dataset
  est.drop.pat <- calc.calib.aj(subset(data.mstate, patid != patid.eval), 
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
### pv.same is the pre-calculated pseudo-value
func.calc.pv.aj.efficient <- function(patid.eval, data.mstate, obs.aj, tmat, n.cohort, t.eval, pv.same){
  
  ### Create temporary dataset of the first row of the data of the individuals of interest
  data.mstate.temp <- subset(data.mstate, patid == patid.eval & from == 1 & to == 2)
  
  if (data.mstate.temp$Tstop > t.eval){
    ### Assign pv.pat to be the pre-calculated value of pseudo value
    pv.pat <- pv.same
    
  } else if (data.mstate.temp$Tstop <= t.eval){
    ### Calculate AJ estimate without patient in dataset
    pv.pat <- func.calc.pv.aj(patid.eval = patid.eval, 
                              data.mstate = data.mstate, 
                              obs.aj = obs.aj, 
                              tmat = tmat, 
                              n.cohort = n.cohort, 
                              t.eval = t.eval)
  }
  
  return(pv.pat)
}


###############################################################################
###############################################################################
###
### This section contains functions used for assessing moderate calibration ###
###
###############################################################################
###############################################################################


###
### 4.1A) Assess moderate calibration using binary logistic regression
### Note this does not apply ipcw, and therefore will be incorrect in the presence of censoring
###
calc.calib.blr.mod <- function(data.mstate, data.raw, t.eval, p.est){
  
  # Let's do state 1 first
  # Want to know which individual are in state j at time t
#           data.mstate <- data.mstate.reduc
#           data.raw <- data.raw.reduc
#           t.eval <- ceiling(7*365.25)
#           p.est <- p.true
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #     intersect(ids.state.list[[1]], ids.state.list[[3]])
  #   intersect(ids.state.list[[1]], ids.state.list[[4]])
  #   intersect(ids.state.list[[1]], ids.state.list[[5]])
  #   intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5),
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
           state.poly.fac = as.factor(state.poly))
  
  ### Add the predicted risks, and the logit transormation of the predicted risks to the dataset
  p.est.logit <- log(p.est/(1-p.est))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est, p.est.logit)
  
  ### Fit loess recalibration models to the uncensored observations at time t to calculate intercepts
  loess1 <- loess(state1.bin ~ p.est1, data = data.raw)
  loess2 <- loess(state2.bin ~ p.est2, data = data.raw)
  loess3 <- loess(state3.bin ~ p.est3, data = data.raw)
  loess4 <- loess(state4.bin ~ p.est4, data = data.raw)
  loess5 <- loess(state5.bin ~ p.est5, data = data.raw)
  
  ### Created 'predicted observed' probabilities for each individual
  data.raw$loess.pred.obs1 <- predict(loess1, newdata = data.raw)
  data.raw$loess.pred.obs2 <- predict(loess2, newdata = data.raw)
  data.raw$loess.pred.obs3 <- predict(loess3, newdata = data.raw)
  data.raw$loess.pred.obs4 <- predict(loess4, newdata = data.raw)
  data.raw$loess.pred.obs5 <- predict(loess5, newdata = data.raw)
  
  ### Calculate the ECI
  ECI <- calc.ECI(data.raw, p.est, data.raw[, paste("loess.pred.obs", 1:5, sep = "")])
  
  ### Produce plots for each and store in a list
  plots.list <- vector("list", 5)
  for (i in 1:5){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("loess.pred.obs", i, sep = "")]
    data.raw$obs.true <- data.raw[, paste("p.true", i, sep = "")]
    
    ### Create the plots
    plots.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>% 
                                arrange(pred) %>% select(patid, pred, obs, obs.true)) +
      geom_line(aes(x = pred, y = obs), color = "red") +
      geom_line(aes(x = pred, y = obs.true), color = "blue", linetype = "dotted") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, sep = ""))
  }
  
  ### Return output object
  output.object <- list("plots.list" = plots.list, "ECI" = ECI)
  return(output.object)
  
}



###
### 4.1B) Assess moderate calibration using binary logistic regression and IPCW
###
calc.calib.blr.ipcw.mod <- function(data.mstate, data.raw, t.eval, p.est, max.weight = 10){
  
  # Let's do state 1 first
  # Want to know which individual are in state j at time t
#   data.mstate <- data.mstate.reduc
#   data.raw <- data.raw.reduc
#   t.eval <- ceiling(7*365.25)
#   p.est <- p.true
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #     intersect(ids.state.list[[1]], ids.state.list[[3]])
  #   intersect(ids.state.list[[1]], ids.state.list[[4]])
  #   intersect(ids.state.list[[1]], ids.state.list[[5]])
  #   intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5),
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
           state.poly.fac = as.factor(state.poly))
  
  ### Add the predicted risks, and the logit transormation of the predicted risks to the dataset
  p.est.logit <- log(p.est/(1-p.est))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est, p.est.logit)
  
  ###
  ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
  ###
  weights <- calc.weights(data.raw, t.eval, max.weight = max.weight)
  weights.DGM <- calc.weights.DGMspec(data.raw, t.eval, max.weight = max.weight, cens_shape = 1, cens_scale, cens_beta_x1, cens_beta_x2)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights, weights.DGM)
  
  ###
  ### Fit the miss-specified and perfectly specified calibration models
  ###

  ### Fit loess recalibration models to the uncensored observations at time t to calculate intercepts
  loess.mspec1 <- loess(state1.bin ~ p.est1, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.mspec"])
  loess.mspec2 <- loess(state2.bin ~ p.est2, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.mspec"])
  loess.mspec3 <- loess(state3.bin ~ p.est3, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.mspec"])
  loess.mspec4 <- loess(state4.bin ~ p.est4, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.mspec"])
  loess.mspec5 <- loess(state5.bin ~ p.est5, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.mspec"])
  
  ### Created 'predicted observed' probabilities for each individual
  data.raw$loess.mspec.pred.obs1 <- predict(loess.mspec1, newdata = data.raw)
  data.raw$loess.mspec.pred.obs2 <- predict(loess.mspec2, newdata = data.raw)
  data.raw$loess.mspec.pred.obs3 <- predict(loess.mspec3, newdata = data.raw)
  data.raw$loess.mspec.pred.obs4 <- predict(loess.mspec4, newdata = data.raw)
  data.raw$loess.mspec.pred.obs5 <- predict(loess.mspec5, newdata = data.raw)
  
  ### Produce plots for each and store in a list
  plots.mspec.list <- vector("list", 5)
  for (i in 1:5){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("loess.mspec.pred.obs", i, sep = "")]
    data.raw$obs.true <- data.raw[, paste("p.true", i, sep = "")]
    
    ### Create the plots
    plots.mspec.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>% 
                                      arrange(pred) %>% select(patid, pred, obs, obs.true)) +
      geom_line(aes(x = pred, y = obs), color = "red") +
      geom_line(aes(x = pred, y = obs.true), color = "blue", linetype = "dotted") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, sep = ""))
  }

  ### Fit loess recalibration models to the uncensored observations at time t to calculate intercepts
  loess.pspec1 <- loess(state1.bin ~ p.est1, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.pspec"])
  loess.pspec2 <- loess(state2.bin ~ p.est2, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.pspec"])
  loess.pspec3 <- loess(state3.bin ~ p.est3, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.pspec"])
  loess.pspec4 <- loess(state4.bin ~ p.est4, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.pspec"])
  loess.pspec5 <- loess(state5.bin ~ p.est5, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.pspec"])
  
  ### Created 'predicted observed' probabilities for each individual
  data.raw$loess.pspec.pred.obs1 <- predict(loess.pspec1, newdata = data.raw)
  data.raw$loess.pspec.pred.obs2 <- predict(loess.pspec2, newdata = data.raw)
  data.raw$loess.pspec.pred.obs3 <- predict(loess.pspec3, newdata = data.raw)
  data.raw$loess.pspec.pred.obs4 <- predict(loess.pspec4, newdata = data.raw)
  data.raw$loess.pspec.pred.obs5 <- predict(loess.pspec5, newdata = data.raw)
  
  ### Produce plots for each and store in a list
  plots.pspec.list <- vector("list", 5)
  for (i in 1:5){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("loess.pspec.pred.obs", i, sep = "")]
    data.raw$obs.true <- data.raw[, paste("p.true", i, sep = "")]
    
    ### Create the plots
    plots.pspec.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>% 
                                arrange(pred) %>% select(patid, pred, obs, obs.true)) +
      geom_line(aes(x = pred, y = obs), color = "red") +
      geom_line(aes(x = pred, y = obs.true), color = "blue", linetype = "dotted") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, sep = ""))
  }
  
  ### Fit loess recalibration models to the uncensored observations at time t to calculate intercepts
  loess.DGMspec1 <- loess(state1.bin ~ p.est1, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.DGMspec"])
  loess.DGMspec2 <- loess(state2.bin ~ p.est2, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.DGMspec"])
  loess.DGMspec3 <- loess(state3.bin ~ p.est3, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.DGMspec"])
  loess.DGMspec4 <- loess(state4.bin ~ p.est4, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.DGMspec"])
  loess.DGMspec5 <- loess(state5.bin ~ p.est5, data = data.raw[!is.na(data.raw$state.poly), ], weights = data.raw[!is.na(data.raw$state.poly), "ipcw.DGMspec"])
  
  ### Created 'predicted observed' probabilities for each individual
  data.raw$loess.DGMspec.pred.obs1 <- predict(loess.DGMspec1, newdata = data.raw)
  data.raw$loess.DGMspec.pred.obs2 <- predict(loess.DGMspec2, newdata = data.raw)
  data.raw$loess.DGMspec.pred.obs3 <- predict(loess.DGMspec3, newdata = data.raw)
  data.raw$loess.DGMspec.pred.obs4 <- predict(loess.DGMspec4, newdata = data.raw)
  data.raw$loess.DGMspec.pred.obs5 <- predict(loess.DGMspec5, newdata = data.raw)
  
  ### Produce plots for each and store in a list
  plots.DGMspec.list <- vector("list", 5)
  for (i in 1:5){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("loess.DGMspec.pred.obs", i, sep = "")]
    data.raw$obs.true <- data.raw[, paste("p.true", i, sep = "")]
    
    ### Create the plots
    plots.DGMspec.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>% 
                                        arrange(pred) %>% select(patid, pred, obs, obs.true)) +
      geom_line(aes(x = pred, y = obs), color = "red") +
      geom_line(aes(x = pred, y = obs.true), color = "blue", linetype = "dotted") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, sep = ""))
  }
  
  ### Calculate the ECI
  ECI.mspec <- calc.ECI(data.raw, p.est, data.raw[, paste("loess.mspec.pred.obs", 1:5, sep = "")])
  ECI.pspec <- calc.ECI(data.raw, p.est, data.raw[, paste("loess.pspec.pred.obs", 1:5, sep = "")])
  ECI.DGMspec <- calc.ECI(data.raw, p.est, data.raw[, paste("loess.DGMspec.pred.obs", 1:5, sep = "")])
  
  ### Return output object
  output.object <- list("plots.mspec.list" = plots.mspec.list, "plots.pspec.list" = plots.pspec.list, "plots.DGMspec.list" = plots.DGMspec.list, 
                        "ECI.mspec" = ECI.mspec, "ECI.pspec" = ECI.pspec, "ECI.DGMspec" = ECI.DGMspec)
  return(output.object)
  
}


###
### 4.2A) Assess moderate calibration using multinomial nominal recalibration framework
### Note this does not apply ipcw, and therefore will be incorrect in the presence of censoring
###
calc.calib.mlr.mod <- function(data.mstate, data.raw, t.eval, p.est, ps.int = 4, degree = 3){
  
  # Let's do state 1 first
  # Want to know which individual are in state j at time t
#             data.mstate <- data.mstate.reduc
#             data.raw <- data.raw.reduc
#             t.eval <- ceiling(7*365.25)
#             p.est <- p.true
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #     intersect(ids.state.list[[1]], ids.state.list[[3]])
  #   intersect(ids.state.list[[1]], ids.state.list[[4]])
  #   intersect(ids.state.list[[1]], ids.state.list[[5]])
  #   intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5))
  
  ### Add p.est to dataset
  data.raw <- data.frame(data.raw, p.est)
  
  ### Add linear predictors from a multinomial framework
  data.raw <- data.raw %>%
    mutate(mlr.lp1 = log(p.est2/p.est1),
           mlr.lp2 = log(p.est3/p.est1),
           mlr.lp3 = log(p.est4/p.est1),
           mlr.lp4 = log(p.est5/p.est1),
           state.poly.fac = as.factor(state.poly))
  
  ### Apply nominal recalibration framework with vector spline smoothers
  calib.model <- vgam(state.poly.fac ~ sm.ps(mlr.lp1, ps.int = ps.int, degree = degree) + sm.ps(mlr.lp2, ps.int = ps.int, degree = degree) + sm.ps(mlr.lp3, ps.int = ps.int, degree = degree) + 
                        sm.ps(mlr.lp4, ps.int = ps.int, degree = degree), 
                      data = data.raw, family = multinomial(refLevel = "1"))
  
  ### Now generate predicted risks for each individual either A), in the validation cohort, or B) fo rsome vector of predicted risks
  mlr.pred.obs <- predict(calib.model, newdata = data.raw, type = "response")
  colnames(mlr.pred.obs) <- paste("mlr.pred.obs", 1:ncol(mlr.pred.obs), sep = "")

  ### Add plot data to data.raw
  data.raw <- cbind(data.raw, mlr.pred.obs)
  
  ### Calculate the ECI
  ECI <- calc.ECI(data.raw, p.est, mlr.pred.obs)
  
  ### Create plots
  ### Produce plots for each and store in a list
  plots.list <- vector("list", 5)
  for (i in 1:5){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("mlr.pred.obs", i, sep = "")]
    data.raw$obs.true <- data.raw[, paste("p.true", i, sep = "")]
    
    ### Create the plots
    plots.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>%
                                arrange(pred) %>% select(patid, pred, obs, obs.true)) +
      geom_point(aes(x = pred, y = obs), color = "red", size = 0.5, alpha = 0.05) +
      geom_line(aes(x = pred, y = obs.true), color = "blue", linetype = "dotted") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, sep = ""))
  }
  
  ### Return output object
  output.object <- list("plots.list" = plots.list, "ECI" = ECI)
  return(output.object)
  
}


###
### 4.2B) Assess moderate calibration using multinomial nominal recalibration framework and IPCW
###
calc.calib.mlr.ipcw.mod <- function(data.mstate, data.raw, t.eval, p.est, ps.int = 4, degree = 3, max.weight = 10){
  
  # Let's do state 1 first
  # Want to know which individual are in state j at time t
#               data.mstate <- data.mstate.reduc
#               data.raw <- data.raw.reduc
#               t.eval <- ceiling(7*365.25)
#               p.est <- p.true
#   ps.int <- 4
#   degree <- 3
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Etract ids for each state at this time
  ids.state.list <- vector("list", max(data.mstate$to))
  for (j in 1:max(data.mstate$to)){
    ids.state.list[[j]] <- extract.ids.states(data.mstate, j, t.eval)
  }
  
  ### Check there is no intersection (i.e. any individual is only in one state, and my code isn't wrong)
  #     intersect(ids.state.list[[1]], ids.state.list[[3]])
  #   intersect(ids.state.list[[1]], ids.state.list[[4]])
  #   intersect(ids.state.list[[1]], ids.state.list[[5]])
  #   intersect(ids.state.list[[4]], ids.state.list[[5]])
  
  ###
  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>% 
    mutate(state.poly = case_when(patid %in% ids.state.list[[1]] ~ 1,
                                  patid %in% ids.state.list[[2]] ~ 2,
                                  patid %in% ids.state.list[[3]] ~ 3,
                                  patid %in% ids.state.list[[4]] ~ 4,
                                  patid %in% ids.state.list[[5]] ~ 5))
  
  ### Add p.est to dataset
  data.raw <- data.frame(data.raw, p.est)
  
  ### Add linear predictors from a multinomial framework
  data.raw <- data.raw %>%
    mutate(mlr.lp1 = log(p.est2/p.est1),
           mlr.lp2 = log(p.est3/p.est1),
           mlr.lp3 = log(p.est4/p.est1),
           mlr.lp4 = log(p.est5/p.est1),
           state.poly.fac = as.factor(state.poly))
  
  ###
  ### Next step is to calculate the ipcw, to apply a weighted (and valid) calibration model
  ###
  weights <- calc.weights(data.raw, t.eval, max.weight = max.weight)
  weights.DGM <- calc.weights.DGMspec(data.raw, t.eval, max.weight = max.weight, cens_shape = 1, cens_scale, cens_beta_x1, cens_beta_x2)
  
  ### Add these to dataset
  data.raw <- cbind(data.raw, weights, weights.DGM)
  
  ###
  ### Miss-specified weighted model
  ###
  
  ### Apply nominal recalibration framework with vector spline smoothers
  calib.model.mspec <- vgam(state.poly.fac ~ sm.ps(mlr.lp1, ps.int = ps.int, degree = degree) + sm.ps(mlr.lp2, ps.int = ps.int, degree = degree) + sm.ps(mlr.lp3, ps.int = ps.int, degree = degree) + 
                        sm.ps(mlr.lp4, ps.int = ps.int, degree = degree), weights = data.raw[!is.na(data.raw$state.poly), "ipcw.mspec"],
                      data = data.raw[!is.na(data.raw$state.poly), ], family = multinomial(refLevel = "1"))
  
  ###
  ### Perfectly specified weighted model
  ###
  
  ### Apply nominal recalibration framework with vector spline smoothers
  calib.model.pspec <- vgam(state.poly.fac ~ sm.ps(mlr.lp1, ps.int = ps.int, degree = degree) + sm.ps(mlr.lp2, ps.int = ps.int, degree = degree) + sm.ps(mlr.lp3, ps.int = ps.int, degree = degree) + 
                              sm.ps(mlr.lp4, ps.int = ps.int, degree = degree), weights = data.raw[!is.na(data.raw$state.poly), "ipcw.pspec"],
                            data = data.raw[!is.na(data.raw$state.poly), ], family = multinomial(refLevel = "1"))
  
  
  ###
  ### DGM specified weighted model
  ###
  
  ### Apply nominal recalibration framework with vector spline smoothers
  calib.model.DGMspec <- vgam(state.poly.fac ~ sm.ps(mlr.lp1, ps.int = ps.int, degree = degree) + sm.ps(mlr.lp2, ps.int = ps.int, degree = degree) + sm.ps(mlr.lp3, ps.int = ps.int, degree = degree) + 
                              sm.ps(mlr.lp4, ps.int = ps.int, degree = degree), weights = data.raw[!is.na(data.raw$state.poly), "ipcw.DGMspec"],
                            data = data.raw[!is.na(data.raw$state.poly), ], family = multinomial(refLevel = "1"))
  
  ###
  ### Generate predicted-observed risks and add to data.raw
  ###
  
  ### For all other functions, I just generate predicted risks for all individuals, and then just plot for those who were uncensored
  ### However, some of the censored individuals are causing an error, therefore I must generate NA vectors, then assign the
  ### predicted observed probabilities to the crrect individuals
  
  ### mspec
  ## Create dataframe to store
  dat.mlr.mspec.pred.obs <- data.frame(matrix(NA, ncol = 5, nrow = nrow(data.raw)))
  ## Assign colnames
  colnames(dat.mlr.mspec.pred.obs) <- paste("mlr.mspec.pred.obs", 1:ncol(dat.mlr.mspec.pred.obs), sep = "")
  ## Calc pred.obs for those who are uncesored
  mlr.mspec.pred.obs <- predict(calib.model.mspec, newdata = data.raw[!is.na(data.raw$state.poly), ], type = "response")
  ## Assign to appropriate individuals
  dat.mlr.mspec.pred.obs[!is.na(data.raw$state.poly), ] <- mlr.mspec.pred.obs
  
  ### pspec
  ## Create dataframe to store
  dat.mlr.pspec.pred.obs <- data.frame(matrix(NA, ncol = 5, nrow = nrow(data.raw)))
  ## Assign colnames
  colnames(dat.mlr.pspec.pred.obs) <- paste("mlr.pspec.pred.obs", 1:ncol(dat.mlr.pspec.pred.obs), sep = "")
  ## Calc pred.obs for those who are uncesored
  mlr.pspec.pred.obs <- predict(calib.model.pspec, newdata = data.raw[!is.na(data.raw$state.poly), ], type = "response")
  ## Assign to appropriate individuals
  dat.mlr.pspec.pred.obs[!is.na(data.raw$state.poly), ] <- mlr.pspec.pred.obs
  
  ### DGMspec
  ## Create dataframe to store
  dat.mlr.DGMspec.pred.obs <- data.frame(matrix(NA, ncol = 5, nrow = nrow(data.raw)))
  ## Assign colnames
  colnames(dat.mlr.DGMspec.pred.obs) <- paste("mlr.DGMspec.pred.obs", 1:ncol(dat.mlr.DGMspec.pred.obs), sep = "")
  ## Calc pred.obs for those who are uncesored
  mlr.DGMspec.pred.obs <- predict(calib.model.DGMspec, newdata = data.raw[!is.na(data.raw$state.poly), ], type = "response")
  ## Assign to appropriate individuals
  dat.mlr.DGMspec.pred.obs[!is.na(data.raw$state.poly), ] <- mlr.DGMspec.pred.obs
  
  ### Then add it to data.raw
  data.raw <- cbind(data.raw, dat.mlr.mspec.pred.obs, dat.mlr.pspec.pred.obs, dat.mlr.DGMspec.pred.obs)
  
  ### Create plots
  
  ### mspec
  ### Produce plots for each and store in a list
  plots.mspec.list <- vector("list", 5)
  for (i in 1:5){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("mlr.mspec.pred.obs", i, sep = "")]
    data.raw$obs.true <- data.raw[, paste("p.true", i, sep = "")]
    
    ### Create the plots
    plots.mspec.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>%
                                      arrange(pred) %>%  select(patid, pred, obs, obs.true)) +
      geom_point(aes(x = pred, y = obs), color = "red", size = 0.5, alpha = 0.05) +
      geom_line(aes(x = pred, y = obs.true), color = "blue", linetype = "dotted") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, sep = ""))
  }
  
  ### pspec
  ### Produce plots for each and store in a list
  plots.pspec.list <- vector("list", 5)
  for (i in 1:5){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("mlr.pspec.pred.obs", i, sep = "")]
    data.raw$obs.true <- data.raw[, paste("p.true", i, sep = "")]
    
    ### Create the plots
    plots.pspec.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>%
                                      arrange(pred) %>%  select(patid, pred, obs, obs.true)) +
      geom_point(aes(x = pred, y = obs), color = "red", size = 0.5, alpha = 0.05) +
      geom_line(aes(x = pred, y = obs.true), color = "blue", linetype = "dotted") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, sep = ""))
  }
  
  ### DGMspec
  ### Produce plots for each and store in a list
  plots.DGMspec.list <- vector("list", 5)
  for (i in 1:5){
    
    ### Creaet variables to plot 
    data.raw$pred <- data.raw[, paste("p.est", i, sep = "")]
    data.raw$obs <- data.raw[, paste("mlr.DGMspec.pred.obs", i, sep = "")]
    data.raw$obs.true <- data.raw[, paste("p.true", i, sep = "")]
    
    ### Create the plots
    plots.DGMspec.list[[i]] <- ggplot(data = data.raw %>% subset(!is.na(data.raw$state.poly)) %>% slice(1:10000) %>%
                                        arrange(pred) %>%  select(patid, pred, obs, obs.true)) +
      geom_point(aes(x = pred, y = obs), color = "red", size = 0.5, alpha = 0.05) +
      geom_line(aes(x = pred, y = obs.true), color = "blue", linetype = "dotted") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      ylim(c(0, max(data.raw$pred[!is.na(data.raw$state.poly)]))) + 
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .005)) + 
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, sep = ""))
  }
  
  ### Calculate the ECI
  ECI.mspec <- calc.ECI(data.raw, p.est, dat.mlr.mspec.pred.obs)
  ECI.pspec <- calc.ECI(data.raw, p.est, dat.mlr.pspec.pred.obs)
  ECI.DGMspec <- calc.ECI(data.raw, p.est, dat.mlr.DGMspec.pred.obs)
  
  ### Return output object
  output.object <- list("plots.mspec.list" = plots.mspec.list, "plots.pspec.list" = plots.pspec.list, "plots.DGMspec.list" = plots.DGMspec.list, 
                        "ECI.mspec" = ECI.mspec, "ECI.pspec" = ECI.pspec, "ECI.DGMspec" = ECI.DGMspec)
  return(output.object)
  
}

###
### 4.3) Calculate Aalen-Johansen estimator seperately for each transition
###
calc.calib.aj.moderate <- function(data.mstate, data.raw, tmat, t.eval, p.est, num.groups){
#     data.mstate <- data.mstate.reduc
#     data.raw <- data.raw.reduc
#     p.est <- p.est.perf
#     tmat
#     t.eval <- ceiling(7*365.25)
#     num.groups <- 10
  
  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")
  
  ### Combine data.raw and predicted probabilities
  data.raw <- cbind(data.raw, p.est)
  
  ### Transition to state2
  aj.decile.pred <- vector("list", max(data.mstate$to))
  aj.decile.obs <- vector("list", max(data.mstate$to))
  
  ### Create object to store sub-cohorts
  data.mstate.deciles <- vector("list", max(data.mstate$to))
  data.raw.deciles <- vector("list", max(data.mstate$to))
  
  ## Create object to store plots
  data.ggplot <- vector("list", max(data.mstate$to))
  plots.ggplot <- vector("list", max(data.mstate$to))
  
  for (i in 1:max(data.mstate$to)){
    
    print(paste("i = ", i, Sys.time()))
    
    ### Create output vector for observed and predicted
    aj.decile.pred[[i]] <- vector(mode = "numeric", num.groups)
    aj.decile.obs[[i]] <- vector(mode = "numeric", num.groups)
    
    ### Arrange data by the risk for each transition
    data.temp <- arrange(data.raw, !! rlang::sym(paste("p.est", i, sep = "")))
    
    ### Create num.groups subcohorts
    data.mstate.deciles[[i]] <- vector("list", num.groups)
    data.raw.deciles[[i]] <- data.temp %>%
      group_by((row_number() - 1) %/% (n()/num.groups)) %>%
      nest %>% pull(data)
    
    ### Calculate AJ within each cohort, and the mean predicted risk
    for (j in 1:num.groups){
      
      print(paste("j = ", j, Sys.time()))
      ### Turn into dataframe from tibble
      data.raw.deciles[[i]][[j]] <- data.frame(data.raw.deciles[[i]][[j]])
      
      ### Create mstate dataframes, based on who is in the data.raw.deciles dataframes
      data.mstate.deciles[[i]][[j]] <- subset(data.mstate, patid %in% data.raw.deciles[[i]][[j]]$patid)
      
      ### Calc AJ
      obs.aj.temp <- calc.calib.aj(data.mstate = data.mstate.deciles[[i]][[j]],
                                   tmat = tmat, 
                                   t.eval = t.eval)[["obs.aj"]]
      
      ### Only interested in the observed risk estimate for state 2
      aj.decile.obs[[i]][j] <- as.numeric(obs.aj.temp[paste("pstate", i, sep = "")])
      aj.decile.pred[[i]][j] <- mean(select(data.raw.deciles[[i]][[j]], paste("p.est", i, sep = ""))[,1])
      
    }
    
    ### Create plots from each of these
    
    ## First create dataset
    data.ggplot[[i]] <- data.frame("obs" = aj.decile.obs[[i]], "pred" = aj.decile.pred[[i]])
    
    ## Now create plot
    plots.ggplot[[i]] <- ggplot(data = data.ggplot[[i]], aes(x = pred, y = obs, color = "red")) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      xlab("Predicted risk") + ylab("Observed risk") + 
      xlim(c(0, max(c(data.ggplot[[i]]$pred, data.ggplot[[i]]$obs)))) + 
      ylim(c(0, max(c(data.ggplot[[i]]$pred, data.ggplot[[i]]$obs)))) +
      theme(legend.position = "none") +
      ggtitle(paste("State ", i, sep = ""))
  }
  
  ### Create output object and return it
  output.obj <- list("obs" = aj.decile.obs, "pred" = aj.decile.pred, "plots.list" = plots.ggplot)
  return(output.obj)
  
}



###
### 4.4) Calculate calibration using pseudo values
###
calc.calib.pv.moderate <- function(data.raw, pv.comb, p.est){

#   str(pv.comb)
#   p.est <- p.est.perf[1:20000, ]
#   data.raw <- data.raw.reduc[1:20000, ]
  
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
  
  ### Created 'predicted observed' probabilities for each individual
  data.pv$loess.pred.obs1 <- predict(loess1, newdata = data.pv)
  data.pv$loess.pred.obs2 <- predict(loess2, newdata = data.pv)
  data.pv$loess.pred.obs3 <- predict(loess3, newdata = data.pv)
  data.pv$loess.pred.obs4 <- predict(loess4, newdata = data.pv)
  data.pv$loess.pred.obs5 <- predict(loess5, newdata = data.pv)
  
  ### Calculate the ECI
  ECI <- calc.ECI(data.raw, p.est, data.pv[, paste("loess.pred.obs", 1:5, sep = "")])
  
  ### Produce plots for each and store in a list
  plots.list <- vector("list", 5)
  for (i in 1:5){
    
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
      ggtitle(paste("State ", i, sep = ""))
  }
  
  
  ### Create output object and return it
  output.object <- list("plots.list" = plots.list, "ECI" = ECI)
  return(output.object)
  
}




#########################################################################
###
### This section contains misc functions used for a variety of things ###
###
#########################################################################

###
### Misc 1) This function will calculate the required scale, to result in a survival probability of p, after 7 years, 
### assuming a exponential baseline hazard
###
calc.scale <- function(p){
  return(7*365.25/(-log(p)))
}

###
### Misc 2) A program to calculate ECI
###
calc.ECI <- function(data.raw.in, p.est.in, p.obs.in){
  
  stopifnot(ncol(p.est.in) == 5)
  stopifnot(ncol(p.obs.in) == 5)
  
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
  
  # Calc numerator
  ECI.numer <- (sum((p.est.in[,1] - p.obs.in[,1])^2) + sum((p.est.in[,2] - p.obs.in[,2])^2) +
                  sum((p.est.in[,3] - p.obs.in[,3])^2) + sum((p.est.in[,4] - p.obs.in[,4])^2) +
                  sum((p.est.in[,5] - p.obs.in[,5])^2))
  # Calc denominator
  ECI.denom <- (sum((p.est.in[,1] - K1)^2) + sum((p.est.in[,2] - K2)^2) +
                  sum((p.est.in[,3] - K3)^2) + sum((p.est.in[,4] - K4)^2) +
                  sum((p.est.in[,5] - K5)^2))
  # Calc ECI
  ECI <- ECI.numer/ECI.denom
  
  return(ECI)
}


###
### Misc 3) Function to calculate IPC weights estimating the model from the data
###
calc.weights <- function(data.raw, t.eval, max.weight = 10){
  
  ### First need to estimate the censoring mechanism from the data

  ### We add a censoring time that would be observed from the data (this would be the min of censoring time, and death)
  ### We take the min of these as the event time
  ### It's an "event" if censoring happens, and "censored" if death happens (as we want to estimate the rate of censoring)
  if (!("State.6" %in% colnames(data.raw))){
    data.raw <- data.raw %>%
      mutate(State.5.noNA = case_when(is.na(State.5) ~ t.eval + 1,
                                      TRUE ~ State.5),
             model.cens.times = pmin(State.5.noNA, cens.times),
             model.cens.s = case_when(cens.times <= State.5.noNA ~ 1,
                                      cens.times > State.5.noNA ~ 0))
  } else if ("State.6" %in% colnames(data.raw)){
    data.raw <- data.raw %>%
      mutate(State.6.noNA = case_when(is.na(State.6) ~ t.eval + 1,
                                      TRUE ~ State.6),
             model.cens.times = pmin(State.6.noNA, cens.times),
             model.cens.s = case_when(cens.times <= State.6.noNA ~ 1,
                                      cens.times > State.6.noNA ~ 0))
  }
  
  ###
  ### Create models for censoring in order to calculate the IPCW weights
  ###
  
  ### A model where we do not adjust for predictors covariates (i.e. miss-specified ipcw weights)
  cens.model.mspec <- coxph(Surv(model.cens.times, model.cens.s) ~ 1, data = data.raw)
  
  ### And a model with perfectly specified ipcw weights
  cens.model.pspec <- coxph(Surv(model.cens.times, model.cens.s) ~ x1 + x2, data = data.raw)
  
  ### Calculate a data frame containing probability of censored and uncenosred at each time point
  ### The weights will be the probability of being uncensored, at the time of the event for each individual
  
  ### For miss-specified model
  ###
  data.weights.mspec <- basehaz(cens.model.mspec)
  data.weights.mspec$p.uncens <- exp(-data.weights.mspec$hazard)
  data.weights.mspec$p.cens <- 1 - data.weights.mspec$p.uncens
  
  
  ### Perfectly specified model
  ###
  
  ## Extract baseline hazard
  data.weights.pspec <- basehaz(cens.model.pspec, centered = FALSE)
  ## Extract linear predictors for each individual
  lp.pspec <- as.matrix(data.raw[, c("x1", "x2")]) %*% as.numeric(coefficients(cens.model.pspec))
  ## Add lp.pspec to data,raw
  data.raw$lp.pspec <- lp.pspec[1:nrow(data.raw),]
  
  ### We do not convert hazards to probabilities yet (like for mspec), as weights are different for each individual, so we do this
  ### within prob.uncens.func.pspec below
  
  
  ### Create weights for the cohort at time t.eval
  ### Note for individuals who died, we take the probability of them being uncensored at the time of death
  ### For individuals still alive, we take the probability of being uncensored at time t.eval
  
  ### If individual has death prior to the time we are evaluating, assign weight at time of death
  ### Get location of individuals who had death prior to evaluation time (or censoring)
  if (!("State.6" %in% colnames(data.raw))){
    obs.death.prior <- !is.na(data.raw$state.poly) & data.raw$State.5.noNA < t.eval
  } else if ("State.6" %in% colnames(data.raw)){
    obs.death.prior <- !is.na(data.raw$state.poly) & data.raw$State.6.noNA < t.eval
  }
  
  ###
  ### Now create miss-specified probability of censoring weights
  ###
  
  ### First assign all individuals a weight of the probability of being uncesored at time t.eval
  data.raw$pcw.mspec <- data.weights.mspec$p.uncens[max(which(data.weights.mspec$time <= t.eval))]
  
  ## Write a function which will extract the uncensored probability for an individual at a given time
  prob.uncens.func.mspec <- function(t){
    return(data.weights.mspec$p.uncens[max(which(data.weights.mspec$time <= t))])
  }
  
  ## Apply this function to all the times at which individuals have died prior to censoring, and assign to the appropriate individuals
  if (!("State.6" %in% colnames(data.raw))){
    data.raw$pcw.mspec[obs.death.prior] <- sapply(data.raw$State.5.noNA[obs.death.prior], prob.uncens.func.mspec)
  } else if ("State.6" %in% colnames(data.raw)){
    data.raw$pcw.mspec[obs.death.prior] <- sapply(data.raw$State.6.noNA[obs.death.prior], prob.uncens.func.mspec)
  }
  
  ### Finally, assign an NA to individuals who are censored prior to time of interest
  data.raw$pcw.mspec[is.na(data.raw$state.poly)] <- NA
  
  ### Invert these
  data.raw$ipcw.mspec <- 1/data.raw$pcw.mspec
  
  ### Finally cap these at 10
  data.raw$ipcw.mspec <- pmin(data.raw$ipcw.mspec, max.weight)
  
  ###
  ### Now create perfectly specified probability of censoring weights
  ###
  
  ### First assign all individuals a weight of the probability of being uncesored at time t.eval
  ### This is the linear predictor times the cumulative hazard at time t.eval, and appropriate transformation to get a risk
  data.raw$pcw.pspec <- as.numeric(exp(-exp(data.raw$lp.pspec)*data.weights.pspec$hazard[max(which(data.weights.pspec$time <= t.eval))]))
  
  ## Write a function which will extract the uncensored probability for an individual at a given time
  prob.uncens.func.pspec <- function(input){
    ## Assign t and patid
    t <- input[1]
    lp <- input[2]
  
    ## Get hazard at appropriate time
    bhaz.t <- data.weights.pspec$hazard[max(which(data.weights.pspec$time <= t))]
    
    ## Return risk
    return(exp(-exp(lp)*bhaz.t))
  }
  
  ## Apply this function to all the times at which individuals have died prior to censoring, and assign to the appropriate individuals
  if (!("State.6" %in% colnames(data.raw))){
    data.raw$pcw.pspec[obs.death.prior] <- apply(data.raw[obs.death.prior, c("State.5.noNA", "lp.pspec")], 1, FUN = prob.uncens.func.pspec)
#     mapply(prob.uncens.func.pspec, t = data.raw$State.5.noNA[obs.death.prior], 
#                                                   patid = data.raw$patid[obs.death.prior])
    
  } else if ("State.6" %in% colnames(data.raw)){
    data.raw$pcw.pspec[obs.death.prior] <- apply(data.raw[obs.death.prior, c("State.6.noNA", "lp.pspec")], 1, FUN = prob.uncens.func.pspec)
#       mapply(prob.uncens.func.pspec, t = data.raw$State.6.noNA[obs.death.prior], 
#                                                   patid = data.raw$patid[obs.death.prior])
  }

  ### Invert these
  data.raw$ipcw.pspec <- 1/data.raw$pcw.pspec
  
  ### Finally cap these at 10
  data.raw$ipcw.pspec <- pmin(data.raw$ipcw.pspec, max.weight)
  
  ### Create output object
  output.object <- data.frame("ipcw.mspec" = data.raw$ipcw.mspec, "ipcw.pspec" = data.raw$ipcw.pspec)
  
  return(output.object)
  
}


###
### Misc 4) Function to calculate IPC weights using the DGM to assign them (i.e. perfect weights)
###
calc.weights.DGMspec <- function(data.raw, t.eval, max.weight = 10, cens_shape, cens_scale, cens_beta_x1, cens_beta_x2){
  
  ### First need to estimate the censoring mechanism from the data
  
  ### We add a censoring time that would be observed from the data (this would be the min of censoring time, and death)
  ### We take the min of these as the event time
  ### It's an "event" if censoring happens, and "censored" if death happens (as we want to estimate the rate of censoring)
  if (!("State.6" %in% colnames(data.raw))){
    data.raw <- data.raw %>%
      mutate(State.5.noNA = case_when(is.na(State.5) ~ t.eval + 1,
                                      TRUE ~ State.5),
             model.cens.times = pmin(State.5.noNA, cens.times),
             model.cens.s = case_when(cens.times <= State.5.noNA ~ 1,
                                      cens.times > State.5.noNA ~ 0))
  } else if ("State.6" %in% colnames(data.raw)){
    data.raw <- data.raw %>%
      mutate(State.6.noNA = case_when(is.na(State.6) ~ t.eval + 1,
                                      TRUE ~ State.6),
             model.cens.times = pmin(State.6.noNA, cens.times),
             model.cens.s = case_when(cens.times <= State.6.noNA ~ 1,
                                      cens.times > State.6.noNA ~ 0))
  }
  
  ### Create weights for the cohort at time t.eval
  ### Note for individuals who died, we take the probability of them being uncensored at the time of death
  ### For individuals still alive, we take the probability of being uncensored at time t.eval
  
  ### If individual has death prior to the time we are evaluating, assign weight at time of death
  ### Get location of individuals who had death prior to evaluation time (or censoring)
  if (!("State.6" %in% colnames(data.raw))){
    obs.death.prior <- !is.na(data.raw$state.poly) & data.raw$State.5.noNA < t.eval
  } else if ("State.6" %in% colnames(data.raw)){
    obs.death.prior <- !is.na(data.raw$state.poly) & data.raw$State.6.noNA < t.eval
  }
  
  ###
  ### Now assign the correct weights, as identified from the DGM
  ###
  prob_uncens_DGM <- function(input){
    
    ## Input is a vector of x1, x2 and t.eval
    x1 <- input[1]
    x2 <- input[2]
    t <- input[3]
  
    ## Return probability of being uncensored at time t
    prob <- exp(-((t/cens_scale)^(cens_shape))*exp(x1*cens_beta_x1 + x2*cens_beta_x2))
    return(prob)
  }
  
  ### First assign all individuals the weights at time t.eval
  temp.t.eval.data <- cbind(data.raw[, c("x1", "x2")], rep(t.eval, nrow(data.raw)))
  data.raw$pcw.DGMspec <- apply(temp.t.eval.data, 1, FUN = prob_uncens_DGM)
  rm(temp.t.eval.data)
  
  ## Now apply this function to all the times at which individuals have died prior to censoring, and assign to the appropriate individuals
  Sys.time()
  if (!("State.6" %in% colnames(data.raw))){
    data.raw$pcw.DGMspec[obs.death.prior] <- apply(data.raw[obs.death.prior, c("x1", "x2", "State.5.noNA")], 1, FUN = prob_uncens_DGM)
  } else if ("State.6" %in% colnames(data.raw)){
    data.raw$pcw.DGMspec[obs.death.prior] <- apply(data.raw[obs.death.prior, c("x1", "x2", "State.6.noNA")], 1, FUN = prob_uncens_DGM)
  }

  ### Invert these
  data.raw$ipcw.DGMspec <- 1/data.raw$pcw.DGMspec
  
  ### Finally cap these at 10
  data.raw$ipcw.DGMspec <- pmin(data.raw$ipcw.DGMspec, max.weight)
  
  ### Create output object
  output.object <- data.frame("ipcw.DGMspec" = data.raw$ipcw.DGMspec)
  
  return(output.object)
  
}


