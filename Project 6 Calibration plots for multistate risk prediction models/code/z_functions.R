#############################################################
### Program with all functions to be used in this project ###
### This will be sourced at the start of each .R file     ###
#############################################################

### Section 1: Functions for generating data.
### Section 2: Functions for calculating pseudo-values in simulation.
### Section 3: Functions to plot results.
### Section 4: Functions for clinical example.

### NOTES
### - Functions for estimating the calibration plots themselves, are contained within the calibmsm package.
### - Functions for calculating the 'true' transition probabilities, are contained in z_functions_true_transition_probs.R

################################################
### Section 1: Functions for generating data ###
################################################

###
### 1.1) Generate data according to DGM
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
                         shape45, scale45, #shape and scale for weibull baseline hazard for transition 4 -> 5
                         beta.x12,  #covariate effects for transiion 12
                         beta.x13,  #covariate effects for transiion 13
                         beta.x15,  #covariate effects for transiion 15
                         beta.x24,  #covariate effects for transiion 24
                         beta.x25,  #covariate effects for transiion 25
                         beta.x34,  #covariate effects for transiion 34
                         beta.x35,  #covariate effects for transiion 35
                         beta.x45,  #covariate effects for transiion 45
                         x.in, #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
                         numsteps) #number of sampler steps in gems data generation process
{
  
  ## Generate a baseline covariate data frame
  bl <- x.in
  
  ## Generate an empty hazard matrix
  hf <- generateHazardMatrix(5)
  
  ## Change the entries of the transitions we want to allow
  ## Define the transitions as weibull
  hf[[1, 2]] <- function(t, shape, scale, beta) {
    exp(bl["x12"]*beta)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[1, 3]] <- function(t, shape, scale, beta) {
    exp(bl["x13"]*beta)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[1, 5]] <- function(t, shape, scale, beta) {
    exp(bl["x15"]*beta)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[2, 4]] <- function(t, shape, scale, beta) {
    exp(bl["x24"]*beta)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[2, 5]] <- function(t, shape, scale, beta) {
    exp(bl["x25"]*beta)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[3, 4]] <- function(t, shape, scale, beta) {
    exp(bl["x34"]*beta)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[3, 5]] <- function(t, shape, scale, beta) {
    exp(bl["x35"]*beta)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  hf[[4, 5]] <- function(t, shape, scale, beta) {
    exp(bl["x45"]*beta)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
  
  #print(hf)
  
  ## Note that by using (t + sum(history)) we are implementing a clock forward approach
  
  ## Generate an empty parameter matrix
  par <- generateParameterMatrix(hf)
  
  ## Use the vector of scales in each transition hazard
  par[[1, 2]] <- list(shape = shape12, scale = scale12, 
                      beta = beta.x12)
  par[[1, 3]] <- list(shape = shape13, scale = scale13, 
                      beta = beta.x13)
  par[[1, 5]] <- list(shape = shape15, scale = scale15, 
                      beta = beta.x15)
  par[[2, 4]] <- list(shape = shape24, scale = scale24, 
                      beta = beta.x24)
  par[[2, 5]] <- list(shape = shape25, scale = scale25, 
                      beta = beta.x25)
  par[[3, 4]] <- list(shape = shape34, scale = scale34, 
                      beta = beta.x34)
  par[[3, 5]] <- list(shape = shape35, scale = scale35, 
                      beta = beta.x35)
  par[[4, 5]] <- list(shape = shape45, scale = scale45, 
                      beta = beta.x45)
  
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
### 1.2) Convert data generated from gen.dat.DGM1 into a format that can be analysed using the mstate package.
### This is an 'msdata' class dataset.
### The simulated continuous values are rounded to integers (days).
### Censoring mechanism is also applied to the data.
###
convert.mstate.integer.DGM1.cens <- function(cohort.in,
                                             max.follow,
                                             cens_shape,
                                             cens_scale,
                                             cens_beta_x12,
                                             cens_beta_x13,
                                             cens_beta_x15,
                                             cens_beta_x24,
                                             cens_beta_x25,
                                             cens_beta_x34,
                                             cens_beta_x35,
                                             cens_beta_x45){
  
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
  cens.times <- simsurv("weibull", lambdas = 1/cens_scale, gammas = cens_shape, 
                        x = dplyr::select(cohort.in, x12, x13, x15, x24, x25, x34, x35, x45), 
                        betas = c("x12" = cens_beta_x12,
                                  "x13" = cens_beta_x13, 
                                  "x15" = cens_beta_x15, 
                                  "x24" = cens_beta_x24, 
                                  "x25" = cens_beta_x25, 
                                  "x34" = cens_beta_x34, 
                                  "x35" = cens_beta_x35, 
                                  "x45" = cens_beta_x45))
  
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
  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
  colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
    paste0("state", 1:ncol(dat.mstate.temp), ".s")
  
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat.mstate.temp))){
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
  }
  
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
  dat.mstate.temp$x12 <- cohort.in$x12
  dat.mstate.temp$x13 <- cohort.in$x13
  dat.mstate.temp$x15 <- cohort.in$x15
  dat.mstate.temp$x24 <- cohort.in$x24
  dat.mstate.temp$x25 <- cohort.in$x25
  dat.mstate.temp$x34 <- cohort.in$x34
  dat.mstate.temp$x35 <- cohort.in$x35
  dat.mstate.temp$x45 <- cohort.in$x45
  dat.mstate.temp$patid <- cohort.in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3,5), c(4,5), c(4,5), c(5), c()),
                   names = paste0("state", 1:5))
  
  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:5)),
                                 status = c(NA, paste0("state", 2:5, ".s")), 
                                 keep = c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45", "cens.time", "patid"))
  
  return(list("data.mstate" = dat.mstate.temp.wide,
              "tmat" = tmat,
              "data.raw" = data.raw))
  
}


###
### 1.3) Convert data generated from gen.dat.DGM1 into a format that can be analysed using the mstate package.
### This is an 'msdata' class dataset.
### Censoring is not applied. This is used in 3_validate_true_transition_probabilities.R, where we want an uncensored
### cohort, to compare event rates with the transition probabilities calculated from the formulae.
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
  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2 )] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
  colnames(dat.mstate.temp.noNA)[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- 
    paste0("state", 1:ncol(dat.mstate.temp), ".s")
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat.mstate.temp))){
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,i]),(i+ncol(dat.mstate.temp))] <- 1
  }
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat.mstate.temp <- dat.mstate.temp.noNA
  
  ## Now need to add baseline data
  dat.mstate.temp$x12 <- cohort.in$x12
  dat.mstate.temp$x13 <- cohort.in$x13
  dat.mstate.temp$x15 <- cohort.in$x15
  dat.mstate.temp$x24 <- cohort.in$x24
  dat.mstate.temp$x25 <- cohort.in$x25
  dat.mstate.temp$x34 <- cohort.in$x34
  dat.mstate.temp$x35 <- cohort.in$x35
  dat.mstate.temp$x45 <- cohort.in$x45
  dat.mstate.temp$patid <- cohort.in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3,5), c(4,5), c(4,5), c(5), c()),
                   names = paste0("state", 1:5))
  
  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:5)),
                                 status = c(NA, paste0("state", 2:5, ".s")), 
                                 keep = c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45", "patid"))
  
  return(list("data.mstate" = dat.mstate.temp.wide,
              "tmat" = tmat,
              "data.raw" = cohort.in))
  
}


###
### 1.4) Function to calculate true IPC weights based on the DGM (i.e. perfect weights).
### This is used to calculate perfectly specified weights in the sensitivity analyses.
###
calc_weights_DGM <- function(data.raw, 
                             t.eval, 
                             max.weight = 10, 
                             cens_shape, 
                             cens_scale, 
                             cens_beta_x12, 
                             cens_beta_x13, 
                             cens_beta_x15, 
                             cens_beta_x24, 
                             cens_beta_x25, 
                             cens_beta_x34, 
                             cens_beta_x35, 
                             cens_beta_x45){
  
  ### First need to estimate the censoring mechanism from the data
  
  ### Create weights for the cohort at time t.eval
  ### Note for individuals who died, we take the probability of them being uncensored at the time of death
  ### For individuals still alive, we take the probability of being uncensored at time t.eval
  
  ### If individual has death prior to the time we are evaluating, assign weight at time of death
  ### Get location of individuals who had death prior to evaluation time (or censoring)
  obs.death.prior <- (data.raw$dtcens.s == 0) & data.raw$dtcens < t.eval
  
  ###
  ### Now assign the correct weights, as identified from the DGM
  ###
  prob_uncens_DGM <- function(input){
    
    ## Input is a vector of the x's and t.eval
    x12 <- input[1]
    x13 <- input[2]
    x15 <- input[3]
    x24 <- input[4]
    x25 <- input[5]
    x34 <- input[6]
    x35 <- input[7]
    x45 <- input[8]
    t <- input[9]
    
    ## Return probability of being uncensored at time t
    prob <- exp(-((t/cens_scale)^(cens_shape))*exp(x12*cens_beta_x12 + 
                                                     x13*cens_beta_x13 + 
                                                     x15*cens_beta_x15 + 
                                                     x24*cens_beta_x24 + 
                                                     x25*cens_beta_x25 + 
                                                     x34*cens_beta_x34 + 
                                                     x35*cens_beta_x35 + 
                                                     x45*cens_beta_x45))
    return(prob)
  }
  
  ### First assign all individuals the weights at time t.eval
  temp.t.eval.data <- cbind(data.raw[, c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45")], 
                            rep(t.eval, nrow(data.raw)))
  data.raw$pcw.DGMspec <- apply(temp.t.eval.data, 1, FUN = prob_uncens_DGM)
  rm(temp.t.eval.data)
  
  ## Now apply this function to all the times at which individuals have died prior to censoring, and assign to the appropriate individuals
  data.raw$pcw.DGMspec[obs.death.prior] <- 
    apply(data.raw[obs.death.prior, c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45", "dtcens")], 
          1, 
          FUN = prob_uncens_DGM)
  
  ### Invert these
  data.raw$ipcw.DGMspec <- 1/data.raw$pcw.DGMspec
  
  ### Finally cap these at 10
  data.raw$ipcw.DGMspec <- pmin(data.raw$ipcw.DGMspec, max.weight)
  
  ### Create output object
  output.object <- data.frame("ipcw.DGM" = data.raw$ipcw.DGMspec)
  
  return(output.object)
  
}


###
### 1.5) This function will calculate the required scale, to result in a survival probability of p, after 7 years, 
### assuming an exponential hazard.
###
calc.scale <- function(p){
  return(7*365.25/(-log(p)))
}


#########################################################################
### Section 2: Functions for calculating pseudo-values in simulation. ###
#########################################################################

###
### 2.1) Calculate Aalen-Johansen estimator for a cohort of class 'msdata', at time t.eval.
###
calc.calib.aj <- function(data.mstate, tmat, t.eval)
  
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
### 2.2) Define a function to calculate pseudo-value for an individual, using the Aalen-Johansen estimator.
###
### obs.aj is pre-calculated Aalen-Johansen estimator in entire cohort data.mstate
###
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
### 2.3) Define an efficient function to calculate pseudo-value for an individual, using the Aalen-Johansen estimator.
### Efficiency is gained through the fact that the pseudo-value is the same for all individuals still in state 1 and
### uncensored at time t.eval. There is therefore an extra input argument containing the pseudo-value for someone who 
### this is the case. It then only has to be calculated once, prior to calculating pseudo-values for the entire cohort.
###
### obs.aj is pre-calculated Aalen-Johansen estimator in entire cohort data.mstate
### pv.same is the pre-calculated pseudo-value
###
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


#############################################
### Section 3: Functions to plot results. ###
#############################################

###
### 3.1) Create plots for calibration data estimated with BLR-IPCW and PV
### Plots estimated calibration curves alongside the "true" calibration curve.
###
plot_calib_msm <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL,
                           marg.density = FALSE, marg.density.size = 5, marg.density.type = "density",
                           marg.rug = FALSE, marg.rug.transparency = 0.1, 
                           inclu.legend = TRUE, legend.seperate = TRUE, legend.title = NULL, # Note legend.seperate only impacts when marg.density = TRUE
                           inclu.title = TRUE,
                           axis.titles.x = NULL, axis.titles.text.x = "Predicted risk",
                           axis.titles.y = NULL, axis.titles.text.y = "Observed risk",
                           bottom = NULL, left = NULL, size = 12, CI = TRUE){
  
  ### Extract plot data and relevant metadata
  plot.data <- x
  
  ### Create list to store plots
  plots.list <- vector("list", 5)
  
  for (state.k in 1:5){
    
    ### Assign plot data
    plot.data.k <- plot.data[[state.k]]
    
    ### Transform data if CI == TRUE
    if (CI == TRUE){
      ### Pivot longer to create data for ggplot and assign appropriate labels
      plot.data.k.longer <- tidyr::pivot_longer(plot.data.k, cols = c(obs, obs.upper, obs.lower, true), names_to = "line.group")
      plot.data.k.longer <- dplyr::mutate(plot.data.k.longer,
                                          line.group = base::factor(line.group),
                                          mapping = dplyr::case_when(line.group == "obs" ~ 1,
                                                                     line.group %in% c("obs.upper", "obs.lower") ~ 2,
                                                                     line.group == "true" ~ 3),
                                          mapping = base::factor(mapping))
      
      levels(plot.data.k.longer$line.group) <- c("Calibration", "Upper", "Lower", "True")
      levels(plot.data.k.longer$mapping) <- c("Calibration", "95% CI", "True")
    } else {
      ### Transform data without CI's
      plot.data.k.longer <- tidyr::pivot_longer(plot.data.k, cols = c(obs, true), names_to = "line.group")
      plot.data.k.longer <- dplyr::mutate(plot.data.k.longer,
                                          line.group = base::factor(line.group),
                                          mapping = dplyr::case_when(line.group == "obs" ~ 1,
                                                                     line.group == "true" ~ 2),
                                          mapping = base::factor(mapping))
      
      levels(plot.data.k.longer$line.group) <- c("Calibration", "True")
      levels(plot.data.k.longer$mapping) <- c("Calibration", "True")
      
    }
    
    ### Create the plots
    plots.list[[state.k]] <- ggplot2::ggplot(data = plot.data.k.longer |> 
                                               dplyr::arrange(pred) |> 
                                               dplyr::select(pred, line.group, value, mapping)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = value, group = line.group, color = mapping, lty = mapping)) +
      ggplot2::xlim(c(min(min(plot.data.k.longer$value), min(plot.data.k.longer$pred)),
                      max(max(plot.data.k.longer$value), max(plot.data.k.longer$pred)))) +
      ggplot2::ylim(c(min(min(plot.data.k.longer$value), min(plot.data.k.longer$pred)),
                      max(max(plot.data.k.longer$value), max(plot.data.k.longer$pred)))) +
      ggplot2::theme(text = ggplot2::element_text(size = size), 
                     legend.text = ggplot2::element_text(size = size)) +
      ggplot2::labs(x = NULL, y = NULL)
    
    ### Specify colours manually to be color blind friendly. Number of colours depends on if CI is true or FALSE
    if (CI == TRUE){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::scale_color_manual(values = c("red", "red", "blue")) + 
        ggplot2::scale_linetype_manual(values = c("solid", "longdash", "solid"))
    } else if (CI == FALSE){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::scale_color_manual(values = c("red", "blue")) + 
        ggplot2::scale_linetype_manual(values = c("solid", "solid"))
    }
    
    ### Add legend title if specified
    if (is.null(legend.title)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::theme(legend.title = ggplot2::element_blank())
    } else if (!is.null(legend.title)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::theme(legend.title = ggplot2::element_text(size = size, face = "bold")) + 
        ggplot2::guides(color = ggplot2::guide_legend(title = legend.title), lty = ggplot2::guide_legend(title = legend.title))
    }
    
    ### Add ggtitles if specified
    if (inclu.title == TRUE){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::ggtitle(paste("State ", state.k, sep = ""))
    }
    
    ### Add axis titles
    if (is.null(axis.titles.x)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::xlab(axis.titles.text.x)
    } else if (!is.null(axis.titles.x)){
      if (state.k %in% axis.titles.x){
        plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::xlab(axis.titles.text.x)
      } 
    }
    
    if (is.null(axis.titles.y)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::ylab(axis.titles.text.y)
    } else if (!is.null(axis.titles.y)){
      if (state.k %in% axis.titles.y){
        plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::ylab(axis.titles.text.y)
      }
    }
    
    ### If marginal density plot has been requested add density plot
    if (marg.density == TRUE){
      ## Save legend
      if (state.k == 1){
        legend.save <- ggpubr::get_legend(plots.list[[state.k]], position = "bottom")
      }
      plots.list[[state.k]] <- plots.list[[state.k]] +
        ## Add a geom_point object of the line and set to invisible (scatter plot required for marginal density using ggMarginal)
        ## Subset to ignore the confidence intervals when doing the density plots
        ggplot2::geom_point(data = plot.data.k.longer |> 
                              dplyr::arrange(pred) |> 
                              dplyr::select(pred, line.group, value, mapping) |> 
                              subset(line.group == "Calibration"),
                            ggplot2::aes(x = pred, y = value), col = grDevices::rgb(0, 0, 0, alpha = 0)) + 
        ## Remove legend
        ggplot2::theme(legend.position = "none")
      
      ## Add ggMarginal
      plots.list[[state.k]] <- ggExtra::ggMarginal(plots.list[[state.k]], 
                                                   margins = "x", 
                                                   x = pred, 
                                                   size = marg.density.size, 
                                                   type = marg.density.type, 
                                                   colour = "red")
      
      ### If marginal density plot was not requested
    } else {
      ## Remove legend
      if (inclu.legend == FALSE){
        plots.list[[state.k]] <- plots.list[[state.k]] + 
          ggplot2::theme(legend.position = "none")
      }
      ## Add marginal rug plot if requested
      if (marg.rug == TRUE){
        plots.list[[state.k]] <- plots.list[[state.k]] +
          ggplot2::geom_rug(data = plot.data.k.longer |> dplyr::arrange(pred) |> dplyr::select(pred, line.group, value, mapping) |> subset(line.group == "Calibration"),
                            ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = marg.rug.transparency))
      }
    }
  }
  
  ### Assign nrow and ncol if not provided by user
  if (is.null(nrow)){
    nrow <- 2
  }
  if (is.null(ncol)){
    ncol <- base::ceiling(length(plots.list)/2)
  }
  
  ### Combine plots into single ggplot
  if (combine == TRUE){
    if (marg.density == FALSE){
      if (inclu.legend == TRUE){
        plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, common.legend = TRUE, legend = "bottom")
      } else {
        plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol)
      }
      
    } else if (marg.density == TRUE){
      plots.list <- gridExtra::arrangeGrob(grobs = plots.list,
                                           layout_matrix = base::matrix(base::seq_len(nrow*ncol),
                                                                        nrow = nrow,
                                                                        ncol = ncol,
                                                                        byrow = TRUE),
                                           top = NULL, bottom = bottom, left = left)
      
      if (inclu.legend == TRUE){
        if (legend.seperate == FALSE){
          plots.list <- gridExtra::arrangeGrob(plots.list, legend.save, nrow = 2, heights = c(15, 1))
        } else {
          plots.list <- list("plots" = plots.list, "legend" = legend.save)
        }
        
      }
    }
  }
  
  ### Return output object
  return(plots.list)
  
}


###
### 3.2) Create plots for calibration data estimated with MLR-IPCW.
### Superimposes estimated calibration scatter plot alongside the "true" calibration scatter plot.
###
plot_calib_mlr <- function(x, ..., combine = TRUE, ncol = NULL, nrow = NULL, point.size = 0.5, transparency.plot = 0.25,
                           marg.density = FALSE, marg.density.size = 5, marg.density.type = "density",
                           marg.rug = FALSE, marg.rug.transparency = 0.1, 
                           inclu.legend = TRUE, legend.seperate = TRUE, legend.title = NULL, # Note legend.seperate only impacts when marg.density = TRUE
                           inclu.title = TRUE,
                           axis.titles.x = NULL, axis.titles.text.x = "Predicted risk",
                           axis.titles.y = NULL, axis.titles.text.y = "Observed risk",
                           bottom = NULL, left = NULL, size = 12){
  
  ### Extract plot data and relevant metadata
  plot.data <- x
  
  ### Create list to store plots
  plots.list <- vector("list", 5)
  
  ### Cycle through states
  for (state.k in 1:5){
    
    ### Assign plot data
    plot.data.k <- plot.data[[state.k]]
    
    ### Pivot longer to create data for ggplot and assign appropriate labels
    plot.data.k.longer <- tidyr::pivot_longer(plot.data.k, cols = c(obs, true), names_to = "group")
    plot.data.k.longer <- dplyr::mutate(plot.data.k.longer,
                                        group = base::factor(group),
                                        mapping = dplyr::case_when(group == "obs" ~ 1,
                                                                   group == "true" ~ 2),
                                        mapping = base::factor(mapping))
    
    levels(plot.data.k.longer$group) <- c("Calibration", "True")
    levels(plot.data.k.longer$mapping) <- c("Calibration", "True")
    
    ### Create the plots
    plots.list[[state.k]] <- ggplot2::ggplot(data = plot.data.k.longer |> 
                                               dplyr::arrange(dplyr::desc(group)) |>  
                                               dplyr::select(id, pred, value, group)) +
      ggplot2::geom_point(ggplot2::aes(x = pred, y = value, color = group), alpha = transparency.plot, size = point.size) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      ggplot2::xlim(c(min(plot.data.k.longer$value), max(plot.data.k.longer$pred))) +
      ggplot2::ylim(c(min(plot.data.k.longer$value), max(plot.data.k.longer$value))) +
      ggplot2::theme(text = ggplot2::element_text(size = size), 
                     legend.text = ggplot2::element_text(size = size)) +
      ggplot2::labs(x = NULL, y = NULL) + 
      ggplot2::scale_color_manual(values = c("red", "blue"))
    
    ### Add legend title if specified
    if (is.null(legend.title)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::theme(legend.title = ggplot2::element_blank())
    } else if (!is.null(legend.title)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::theme(legend.title = ggplot2::element_text(size = size, face = "bold")) + 
        ggplot2::guides(color = ggplot2::guide_legend(title = legend.title))
    }
    
    ### Add ggtitles if specified
    if (inclu.title == TRUE){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::ggtitle(paste("State ", state.k, sep = ""))
    }
    
    ### Add axis titles
    if (is.null(axis.titles.x)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::xlab(axis.titles.text.x)
    } else if (!is.null(axis.titles.x)){
      if (state.k %in% axis.titles.x){
        plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::xlab(axis.titles.text.x)
      } 
    }
    
    if (is.null(axis.titles.y)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::ylab(axis.titles.text.y)
    } else if (!is.null(axis.titles.y)){
      if (state.k %in% axis.titles.y){
        plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::ylab(axis.titles.text.y)
      }
    }
    
    ### If marginal density plot has been requested add density plot
    if (marg.density == TRUE){
      if (state.k == 1){
        ### Create a plot for the legend (note we do not want transparency in the legend so have to create it seperately)
        temp.plot <- 
          ggplot2::ggplot(data = plot.data.k.longer 
                          |> dplyr::arrange(pred) 
                          |>  dplyr::select(id, pred, value, group)) +
          ggplot2::geom_point(ggplot2::aes(x = pred, y = value, color = group), size = point.size) + 
          ggplot2::scale_color_manual(values = c("red", "blue")) +
          ggplot2::theme(legend.text = ggplot2::element_text(size = size))
        
        ### Add legend title if specified
        if (is.null(legend.title)){
          temp.plot <- temp.plot + ggplot2::theme(legend.title = ggplot2::element_blank())
        } else if (!is.null(legend.title)){
          temp.plot <- temp.plot + ggplot2::theme(legend.title = ggplot2::element_text(size = size, face = "bold")) + 
            ggplot2::guides(color = ggplot2::guide_legend(title = legend.title))
        }
        
        ### Extract and save legend
        legend.save <- ggpubr::get_legend(temp.plot, position = "bottom") 
      }
      
      ### Remove legend from main plots
      plots.list[[state.k]] <- plots.list[[state.k]]+ 
        ggplot2::theme(legend.position = "none")
      
      ## Add ggMarginal
      plots.list[[state.k]] <- ggExtra::ggMarginal(plots.list[[state.k]], 
                                                   margins = "x", 
                                                   size = marg.density.size, 
                                                   type = marg.density.type, 
                                                   colour = "red")
      ### If marginal rug plot has been requested
    } else {
      ## Remove legend
      if (inclu.legend == FALSE){
        plots.list[[state.k]] <- plots.list[[state.k]] + 
          ggplot2::theme(legend.position = "none")
      }
      ## Add marginal rug plot if requested
      if (marg.rug == TRUE){
        plots.list[[state.k]] <- plots.list[[state.k]] +
          ggplot2::geom_rug(ggplot2::aes(x = pred, y = value), col = grDevices::rgb(1, 0, 0, alpha = marg.rug.transparency))
      }
    }
    
  }
  
  ### Assign nrow and ncol if not provided by user
  if (is.null(nrow)){
    nrow <- 2
  }
  if (is.null(ncol)){
    ncol <- base::ceiling(length(plots.list)/2)
  }
  
  ### Combine plots into single ggplot
  if (combine == TRUE){
    if (marg.density == FALSE){
      if (inclu.legend == TRUE){
        plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, common.legend = TRUE, legend = "bottom")
      } else {
        plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol)
      }
      
    } else if (marg.density == TRUE){
      plots.list <- gridExtra::arrangeGrob(grobs = plots.list,
                                           layout_matrix = base::matrix(base::seq_len(nrow*ncol),
                                                                        nrow = nrow,
                                                                        ncol = ncol,
                                                                        byrow = TRUE),
                                           top = NULL, bottom = bottom, left = left)
      
      if (inclu.legend == TRUE){
        if (legend.seperate == FALSE){
          plots.list <- gridExtra::arrangeGrob(plots.list, legend.save, nrow = 2, heights = c(15, 1))
        } else {
          plots.list <- list("plots" = plots.list, "legend" = legend.save)
        }
        
      }
    }
  }
  
  ### Return output object
  return(plots.list)
  
}

###
### 3.3) Create plots for calibration data estimated with BLR-IPCW and PV in the small sample simulation.
### Superimposes estimated calibration curves from each simulation alongside the "true" calibration curve.
###
plot_calib_msm_small_sample <- function(x, ..., calib.type, combine = TRUE, ncol = NULL, nrow = NULL, plot.transparency = 0.1,
                                        marg.density = FALSE, marg.density.size = 5, marg.density.type = "density",
                                        marg.rug = FALSE, marg.rug.transparency = 0.1, 
                                        inclu.legend = TRUE, legend.seperate = TRUE, legend.title = NULL, # Note legend.seperate only impacts when marg.density = TRUE
                                        inclu.title = TRUE,
                                        axis.titles.x = NULL, axis.titles.text.x = "Predicted risk",
                                        axis.titles.y = NULL, axis.titles.text.y = "Observed risk",
                                        bottom = NULL, left = NULL, size = 12){

  ### Extract plot data and relevant metadata
  plot.data <- x
  
  ### Create list to store plots
  plots.list <- vector("list", 5)
  
  for (state.k in 1:5){
    
    ### Assign plot data
    plot.data.k <- plot.data[[state.k]]
    
    ### Transform data
    plot.data.k <- dplyr::mutate(plot.data.k,
                                 mapping = dplyr::case_when(plot != 0 ~ 1,
                                                            plot == 0 ~ 2),
                                 mapping = base::factor(mapping),
                                 plot = base::factor(plot),
                                 plot = stats::relevel(plot, ref = max(levels(plot))),
                                 plot = factor(plot, levels = rev(levels(plot))))
    levels(plot.data.k$mapping) <- c("Calibration", "True")
    
    ### Create the plots
    plots.list[[state.k]] <- ggplot2::ggplot(data = plot.data.k |> 
                                               dplyr::arrange(desc(plot), pred) |> 
                                               dplyr::select(pred, obs, mapping, plot))
    
    ### Different plot if scatter or lineplot
    if (calib.type == "line"){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
        ggplot2::geom_line(ggplot2::aes(x = pred, y = obs, group = plot, color = mapping, alpha = mapping))
    } else if (calib.type == "scatter"){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
        ggplot2::geom_point(ggplot2::aes(x = pred, y = obs, group = mapping, color = mapping, alpha = mapping))
    }
    
    ### Add rest of stuff to plots
    plots.list[[state.k]] <- plots.list[[state.k]] +
      ggplot2::xlim(c(min(min(plot.data.k$obs), min(plot.data.k$pred)),
                      min(max(plot.data.k$obs), max(plot.data.k$pred)) + 0.1)) +
      ggplot2::ylim(c(min(min(plot.data.k$obs), min(plot.data.k$pred)),
                      min(max(plot.data.k$obs), max(plot.data.k$pred)) + 0.1)) +
      ggplot2::theme(text = ggplot2::element_text(size = size), 
                     legend.text = ggplot2::element_text(size = size)) +
      ggplot2::labs(x = NULL, y = NULL)
    
    ### Specify colours manually to be color blind friendly. Number of colours depends on if CI is true or FALSE
    plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::scale_color_manual(values = c("red", "deepskyblue")) + 
      ggplot2::scale_alpha_manual(values = c(plot.transparency, 1), guide = "none")
    
    ### Add legend title if specified
    if (is.null(legend.title)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::theme(legend.title = ggplot2::element_blank())
    } else if (!is.null(legend.title)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::theme(legend.title = ggplot2::element_text(size = size, face = "bold")) + 
        ggplot2::guides(color = ggplot2::guide_legend(title = legend.title))
    }
    
    ### Add ggtitles if specified
    if (inclu.title == TRUE){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::ggtitle(paste("State ", state.k, sep = ""))
    }
    
    ### Add axis titles
    if (is.null(axis.titles.x)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::xlab(axis.titles.text.x)
    } else if (!is.null(axis.titles.x)){
      if (state.k %in% axis.titles.x){
        plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::xlab(axis.titles.text.x)
      } 
    }
    
    if (is.null(axis.titles.y)){
      plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::ylab(axis.titles.text.y)
    } else if (!is.null(axis.titles.y)){
      if (state.k %in% axis.titles.y){
        plots.list[[state.k]] <- plots.list[[state.k]] + ggplot2::ylab(axis.titles.text.y)
      }
    }
    
    ### If marginal density plot has been requested add density plot
    if (marg.density == TRUE){
      ## Save legend
      if (state.k == 1){
        legend.save <- ggpubr::get_legend(plots.list[[state.k]], position = "bottom")
      }
      plots.list[[state.k]] <- plots.list[[state.k]] +
        ## Add a geom_point object of the line and set to invisible (scatter plot required for marginal density using ggMarginal)
        ## Subset to ignore the confidence intervals when doing the density plots
        ggplot2::geom_point(data = plot.data.k |> 
                              dplyr::arrange(pred) |> 
                              dplyr::select(pred, obs, mapping) |> 
                              subset(mapping == "Calibration"),
                            ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(0, 0, 0, alpha = 0)) + 
        ## Remove legend
        ggplot2::theme(legend.position = "none")
      
      ## Add ggMarginal
      plots.list[[state.k]] <- ggExtra::ggMarginal(plots.list[[state.k]], 
                                                   margins = "x", 
                                                   x = pred, 
                                                   size = marg.density.size, 
                                                   type = marg.density.type, 
                                                   colour = "red")
      
      ### If marginal density plot was not requested
    } else {
      ## Remove legend
      if (inclu.legend == FALSE){
        plots.list[[state.k]] <- plots.list[[state.k]] + 
          ggplot2::theme(legend.position = "none")
      }
      ## Add marginal rug plot if requested
      if (marg.rug == TRUE){
        plots.list[[state.k]] <- plots.list[[state.k]] +
          ggplot2::geom_rug(data = plot.data.k |> dplyr::arrange(pred) |> dplyr::select(pred, obs, mapping) |> subset(mapping == "Calibration"),
                            ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(1, 0, 0, alpha = marg.rug.transparency))
      }
    }
  }
  
  ### Assign nrow and ncol if not provided by user
  if (is.null(nrow)){
    nrow <- 2
  }
  if (is.null(ncol)){
    ncol <- base::ceiling(length(plots.list)/2)
  }
  
  ### Combine plots into single ggplot
  if (combine == TRUE){
    if (marg.density == FALSE){
      if (inclu.legend == TRUE){
        plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, common.legend = TRUE, legend = "bottom")
      } else {
        plots.list <- ggpubr::ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol)
      }
      
    } else if (marg.density == TRUE){
      plots.list <- gridExtra::arrangeGrob(grobs = plots.list,
                                           layout_matrix = base::matrix(base::seq_len(nrow*ncol),
                                                                        nrow = nrow,
                                                                        ncol = ncol,
                                                                        byrow = TRUE),
                                           top = NULL, bottom = bottom, left = left)
      
      if (inclu.legend == TRUE){
        if (legend.seperate == FALSE){
          plots.list <- gridExtra::arrangeGrob(plots.list, legend.save, nrow = 2, heights = c(15, 1))
        } else {
          plots.list <- list("plots" = plots.list, "legend" = legend.save)
        }
        
      }
    }
  }
  
  ### Return output object
  return(plots.list)
  
}


##################################################
### Section 4: Functions for clinical example. ###
##################################################

###
### 4.1) Calculate Aalen-Johansen estimator for a cohort (dataset class 'msdata') from the clinical example at time t.eval
###
calc.calib.aj.ce <- function(data.mstate, tmat, t.eval){

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
### 4.2) Define a function to calculate pseudo-value for an individual, using the Aalen-Johansen estimator,
### for data in the clinical eample.
### 
### obs.aj is pre-calculated Aalen-Johansen estimator in entire cohort
###
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
### 4.3) Define an efficient function to calculate pseudo-value for an individual, using the Aalen-Johansen estimator,
### for data in the clinical example.
### Efficiency is gained through the fact that the pseudo-value is the same for all individuals still in state 1 and
### uncensored at time t.eval. There is therefore an extra input argument containing the pseudo-value for someone who 
### this is the case. It then only has to be calculated once, prior to calculating pseudo-values for the entire cohort.
###
### obs.aj is pre-calculated Aalen-Johansen estimator in entire cohort data.mstate
### pv.same is the pre-calculated pseudo-value
###
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


############################
### Misc other functions ###
### Currently unused     ###
############################

###
### Misc 1) A program to calculate ECI.
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
### Misc 2) A program to calculate ECI in the clinical example.
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


