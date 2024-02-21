
###
### Define DGM
###
DGM1 <- function(n, #number of patients to simulate
                 max.follow, #maximum follow up
                 shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                 shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                 shape23, scale23, #shape and scale for weibull baseline hazard for transition 2 -> 3
                 beta12.x1, beta12.x2, #covariate effects for transiion 12
                 beta13.x1, beta13.x2, #covariate effects for transiion 13
                 beta23.x1, beta23.x2, #covariate effects for transiion 23
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
  #   shape23 <- 1
  #   scale23 <- 5*1588.598
  #
  #   #qweibull(0.8, 1, 1588.598)
  #
  #   ## Covariate effects
  #   beta12.x1 <- 1
  #   beta12.x2 <- 1
  #   beta13.x1 <- 0.5
  #   beta13.x2 <- 0.5
  #   beta23.x1 <- 1
  #   beta23.x2 <- 0.5
  #
  #   x.in <- x.baseline
  #   numsteps <- max.follow

  ## Generate a baseline covariate data frame
  bl <- x.in

  ## Generate an empty hazard matrix
  hf <- generateHazardMatrix(3)
  #hf

  ## Change the entries of the transitions we want to allow
  ## Define the transitions as weibull
  hf[[1, 2]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*(t/scale)^(shape - 1)}

  hf[[1, 3]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*(t/scale)^(shape - 1)}

  hf[[2, 3]] <- function(t, shape, scale, beta.x1, beta.x2) {
    exp(bl["x1"]*beta.x1 + bl["x2"]*beta.x2)*(shape/scale)*(t/scale)^(shape - 1)}

  print(hf)


  ## We are using a clock reset approach to generate data
  ## Replace t with (t + sum(history)) to implement a clock forward approach

  ## Generate an empty parameter matrix
  par <- generateParameterMatrix(hf)

  ## Use the vector of scales in each transition hazard
  par[[1, 2]] <- list(shape = shape12, scale = scale12,
                      beta.x1 = beta12.x1, beta.x2 = beta12.x2)
  par[[1, 3]] <- list(shape = shape13, scale = scale13,
                      beta.x1 = beta13.x1, beta.x2 = beta13.x2)
  par[[2, 3]] <- list(shape = shape23, scale = scale23,
                      beta.x1 = beta23.x1, beta.x2 = beta23.x2)

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


# ###
# ### 1.2B) Convert data into a format that can be analysed using the mstate package
# ### This is what will be used in the actual simulation, is it includes censoring
# ###
# ###
# convert.mstate.cens <- function(cohort.in,
#                                 max.follow,
#                                 cens_shape,
#                                 cens_scale,
#                                 cens_beta_x1,
#                                 cens_beta_x2){
#
#   #   cohort.in <- cohort[["cohort"]]
#   #   max.follow <- cohort[["max.follow"]]
#   #   cens_shape <- 1
#   #   cens_scale <- 4000
#   #   cens_beta_x1 <- 0
#   #   cens_beta_x2 <- 0
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
#   dat.mstate.temp <- select(cohort.in, paste("State.", 1:3, sep = ""))
#   colnames(dat.mstate.temp) <- paste0("state", 1:3)
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
#                               state3 >= cens.time ~ cens.time)) %>%
#     mutate(state2.s = case_when(state2 < cens.time ~ state2.s,
#                                 state2 == cens.time ~ 0),
#            state3.s = case_when(state3 < cens.time ~ state3.s,
#                                 state3 == cens.time ~ 0))
#
#   ## Now need to add baseline data
#   dat.mstate.temp$x1 <- cohort.in$x1
#   dat.mstate.temp$x2 <- cohort.in$x2
#   dat.mstate.temp$patid <- cohort.in$patid
#
#   ### Now we can use msprep from the mstate package to turn into wide format
#   ## First create a transition matrix corresponding to the columns
#   tmat <- transMat(x = list(c(2,3), c(3), c()),
#                    names = paste0("state", 1:3))
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
#                                  time = c(NA, paste0("state", 2:3)),
#                                  status = c(NA, paste0("state", 2:3, ".s")),
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




convert.mstate.cens <- function(cohort.in,
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
           State.3 = round(State.3)) %>%
    mutate(State.2 = case_when(!is.na(State.2) & !is.na(State.1) & State.1 == State.2 ~ State.2 + 1,
                               TRUE ~ State.2),
           State.3 = case_when(!is.na(State.3) & !is.na(State.1) & State.3 == State.1 ~ State.3 + 1,
                               !is.na(State.3) & !is.na(State.2) & State.3 == State.2 ~ State.3 + 1,
                               TRUE ~ State.3)
           )

  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat.mstate.temp <- select(cohort.in, paste("State.", 1:3, sep = ""))
  colnames(dat.mstate.temp) <- paste0("state", 1:3)

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
                              state3 >= cens.time ~ cens.time)) %>%
    mutate(state2.s = case_when(state2 < cens.time ~ state2.s,
                                state2 == cens.time ~ 0),
           state3.s = case_when(state3 < cens.time ~ state3.s,
                                state3 == cens.time ~ 0))

  ## Now need to add baseline data
  dat.mstate.temp$x1 <- cohort.in$x1
  dat.mstate.temp$x2 <- cohort.in$x2
  dat.mstate.temp$patid <- cohort.in$patid

  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- transMat(x = list(c(2,3), c(3), c()),
                   names = paste0("state", 1:3))


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
                                 time = c(NA, paste0("state", 2:3)),
                                 status = c(NA, paste0("state", 2:3, ".s")),
                                 keep = c("x1","x2","cens.time","patid"))

  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("x1", "x2")
  dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
  head(dat.mstate.temp.wide)

  return(list("data.mstate" = dat.mstate.temp.wide,
              "tmat" = tmat,
              "data.raw" = data.raw))

}





