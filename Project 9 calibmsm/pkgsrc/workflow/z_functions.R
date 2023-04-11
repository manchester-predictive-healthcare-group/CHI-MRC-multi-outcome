#############################################################
###                                                       ###
### Section 1) Misc functions used within other functions ###
###                                                       ###
#############################################################

###
### 1.1) Function to extract all individuals in state j, at time t.eval, from a dataset in mstate format
### Used in other functions which assess calibration using blr/mlr at specific time points (3.4 and 3.5)
###
extract.ids.states <- function(data.mstate, tmat, j, t.eval){

  ### Define maximum state number
  max.state <- max(data.mstate$to)

  ### Identify which states are absorbing states
  absorbing.states <- which(apply(tmat, 1, function(x) {sum(!is.na(x))}) == 0)

  ### For non-absorbing states, to be in state j at time t, you must have an observations from state j, where Tstart <= t.eval < Tstop
  if (!(j %in% absorbing.states)){
    ## Extract ids
    ids.state.j <- subset(data.mstate, from == j & Tstart <= t.eval & t.eval < Tstop) %>%
      select(id) %>%
      distinct(id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$id)
  } else if (j %in% absorbing.states){
    ### For absorbing state, just have to have moved into it
    ids.state.j <- subset(data.mstate, to == j & t.eval >= Tstop & status == 1) %>%
      select(id) %>%
      distinct(id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$id)
  }

  return(ids.state.j)
}


###
### 1.2) Function to calculate IPC weights estimating the model from the data
###

### The dataset data.raw must have a time until censoring variable (dtcens) and a an event indicator (dtcens.s)
### - dtcens.s = 1 if censoring is observed at time dtcens, dtcens.s = 0 otherwise.
### - The variable dtcens is the time until an individual has been censored, or the time until they have reached an absorbing state. If
###   and individual reached an absorbing state, dtcens.s = 0 (event indicator = 0, as censoring has not been observed)

calc.weights <- function(data.mstate, data.raw, covs, j = NULL, landmark.type = NULL, s, t.eval, max.weight = 10, stabilised = FALSE){

#   data.mstate <- msebmt
#   data.raw <- ebmt
#   covs <- NULL
#   j
#   landmark.type <- "state"
#   s
#   t.eval <- t.eval
#   max.weight <- 10
#
  ### Create a new outcome, which is the time until censored from s
  data.raw$dtcens.modified <- data.raw$dtcens - s

  ### Save a copy of data.raw
  data.raw.save <- data.raw

  ### If landmark.type = "state", calculate weights only in individuals in state j at time s
  ### If landmark.type = "all", calculate weights in all uncensored individuals at time s (note that this excludes individuals
  ### who have reached absorbing states, who have been 'censored' from the survival distribution is censoring)
  if (landmark.type == "state"){
    ### Identify individuals who are uncensored in state j at time s
    ids.uncens <- subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
      select(id) %>%
      distinct(id) %>%
      pull(id)

  } else if (landmark.type == "all"){
    ### Identify individuals who are uncensored time s
    ids.uncens <- subset(data.mstate, Tstart <= s & s < Tstop) %>%
      select(id) %>%
      distinct(id) %>%
      pull(id)

  }

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate %>% subset(id %in% ids.uncens)
  data.raw <- data.raw %>% subset(id %in% ids.uncens)

  ###
  ### Create models for censoring in order to calculate the IPCW weights
  ### Seperate models for estimating the weights, and stabilising the weights (intercept only model)
  ###
  if (!is.null(covs)){
    ### A model where we adjust for predictor variables
    cens.model <- coxph(as.formula(paste("Surv(dtcens.modified, dtcens.s) ~ ",
                                         paste(covs, collapse = "+"),
                                         sep = "")),
                        data = data.raw)

    ### Intercept only model (numerator for stabilised weights)
    cens.model.int <- coxph(as.formula(paste("Surv(dtcens.modified, dtcens.s) ~ 1",
                                             sep = "")),
                            data = data.raw)
  } else if (is.null(covs)){
    ### If user has not input any predictors for estimating weights, the model for estimating the weights is the intercept only model (i.e. Kaplan Meier estimator)

    ### Intercept only model (numerator for stabilised weights)
    cens.model.int <- coxph(as.formula(paste("Surv(dtcens.modified, dtcens.s) ~ 1",
                                             sep = "")),
                            data = data.raw)
    ### Assign cens.model to be the same
    cens.model <- cens.model.int


  }

  ### Calculate a data frame containing probability of censored and uncenosred at each time point
  ### The weights will be the probability of being uncensored, at the time of the event for each individual

  ## Extract baseline hazard
  data.weights <- basehaz(cens.model, centered = FALSE)
  ## Add lp to data.raw.save
  data.raw.save$lp <- predict(cens.model, newdata = data.raw.save, type = "lp", reference = "zero")

  ### Create weights for the cohort at time t.eval - s
  ### Note for individuals who reached an absorbing state, we take the probability of them being uncensored at the time of reached the
  ### abosrbing state. For individuals still alive, we take the probability of being uncensored at time t.eval - s.

  ### Get location of individuals who entered absorbing states or were censored prior to evaluation time
  obs.absorbed.prior <- which(data.raw.save$dtcens < t.eval & data.raw.save$dtcens.s == 0)
  obs.censored.prior <- which(data.raw.save$dtcens < t.eval & data.raw.save$dtcens.s == 1)

  ###
  ### Now create unstabilised probability of (un)censoring weights
  ### Note that weights are the probability of being uncensored, so if an individual has low probability of being uncesored,
  ### the inervse of this will be big, weighting them strongly
  ###

  ### First assign all individuals a weight of the probability of being uncesored at time t.eval
  ### This is the linear predictor times the cumulative hazard at time t.eval, and appropriate transformation to get a risk
  data.raw.save$pcw <- as.numeric(exp(-exp(data.raw.save$lp)*data.weights$hazard[max(which(data.weights$time <= t.eval - s))]))

  ## Write a function which will extract the uncensored probability for an individual with linear predictor lp at a given time t
  prob.uncens.func <- function(input){

    ## Assign t and person_id
    t <- input[1]
    lp <- input[2]

    if (t <= 0){
      return(NA)
    } else if (t > 0){
      ## Get hazard at appropriate time
      if (t < min(data.weights$time)){
        bhaz.t <- 0
      } else if (t >= min(data.weights$time)){
        bhaz.t <- data.weights$hazard[max(which(data.weights$time <= t))]
      }

      ## Return risk
      return(exp(-exp(lp)*bhaz.t))
    }
  }

  ### Apply this function to all the times at which individuals have entered an absorbing state prior to censoring
  data.raw.save$pcw[obs.absorbed.prior] <- apply(data.raw.save[obs.absorbed.prior, c("dtcens.modified", "lp")], 1, FUN = prob.uncens.func)

  ### For individuals who were censored prior to t.eval, assign the weight as NA
  data.raw.save$pcw[obs.censored.prior] <- NA

  ### Invert these
  data.raw.save$ipcw <- 1/data.raw.save$pcw

  ###
  ### Stabilise these weights dependent on user-input
  ###
  if (stabilised == TRUE){

    ## Extract baseline hazard
    data.weights.numer <- basehaz(cens.model.int, centered = TRUE)

    ### Assign all individuals a weight of the probability of being uncesored at time t.eval
    data.raw.save$pcw.numer <- as.numeric(exp(-data.weights.numer$hazard[max(which(data.weights.numer$time <= t.eval - s))]))

#     ### Write a function which will extract the uncensored probability for an individual with linear predictor lp at a given time t
#     prob.uncens.func.numer <- function(input){
#       ## Assign t
#       t <- input
#
#       ## Get hazard at appropriate time
#       bhaz.t <- data.weights.numer$hazard[max(which(data.weights.numer$time <= t))]
#
#       ## Return risk
#       return(exp(-bhaz.t))
#     }
#
#     ### Apply this function to all the times at which individuals have died prior to censoring, and assign to the appropriate individuals
#     data.raw.save$pcw.numer[obs.absorbed.prior] <- sapply(data.raw.save[obs.absorbed.prior, c("dtcens.modified")], FUN = prob.uncens.func.numer)
#
    ### Create stabilised weight
    data.raw.save$ipcw.stab <- data.raw.save$pcw.numer*data.raw.save$ipcw
  }

  ### Finally cap these at 10 and create output object

  ### Create output object
  if (stabilised == FALSE){
    data.raw.save$ipcw <- pmin(data.raw.save$ipcw, max.weight)
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw, "pcw" = data.raw.save$pcw)
  } else if (stabilised == TRUE){
    data.raw.save$ipcw <- pmin(data.raw.save$ipcw, max.weight)
    data.raw.save$ipcw.stab <- pmin(data.raw.save$ipcw.stab, max.weight)
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw, "ipcw.stab" = data.raw.save$ipcw.stab, "pcw" = data.raw.save$pcw)
  }

  return(output.weights)

}


##########################################################################
###                                                                    ###
### Section 2) Functions to create calibration curves or scatter plots ###
###                                                                    ###
##########################################################################

###
### 2.1) Calculate calibration using binary logistic regression framework with inverse probability of censoring weights.
### Calibration model uses loess smoothers.
###
calc.calib.blr.loess <- function(data.mstate, data.raw, tmat, j, s, t.eval, p.est, span, degree, CI = FALSE, R.boot, stabilised = FALSE){

#   ### Calculate ipcw weights
#   span <- 1
#   degree.blr <- 2
#   weights.fixed <- calc.weights(data.mstate = msebmt,
#                                 data.raw = ebmt,
#                                 covs = covs,
#                                 j = j,
#                                 landmark.type = "all",
#                                 s = s,
#                                 t.eval = t.eval)
#
#   ### Combineweights with ebmt dataset
#   ebmt.weights <- cbind(ebmt, weights.fixed)
#   data.mstate = msebmt
#   data.raw = ebmt.weights
#   tmat
#   j=j
#   s=s
#   t.eval = t.eval
#   p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = "")))
#   span = span
#   degree = degree.blr
#   CI = 95
#   R.boot = 250

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(p.est) != 0)

  ### Add the predicted risks, and the logit transormation of the predicted risks to data.raw
  p.est.logit <- log(p.est[,valid.transitions]/(1-p.est[,valid.transitions]))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est[,valid.transitions], p.est.logit)

  ### Identify individuals who are in state j at time s
  ids.state.j <- subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    select(id) %>%
    distinct(id) %>%
    pull(id)

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate %>% subset(id %in% ids.state.j)
  data.raw <- data.raw %>% subset(id %in% ids.state.j)

  ### Extract which state individuals are in at time t.eval
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract.ids.states(data.mstate, tmat, k, t.eval)
  }

  #   ### Create a function which identifies which states an individual is in
  #   identify.state <- function(id){
  #     if (sum(lapply(ids.state.list, function(x){id %in% x}) == TRUE) == 1){
  #       group <- which(lapply(ids.state.list, function(x){id %in% x}) == TRUE)
  #     } else if (sum(lapply(ids.state.list, function(x){id %in% x}) == TRUE) == 0){
  #       group <- NA
  #     }
  #     return(group)
  #   }
  #   names(ids.state.list) <- 1:6
  #   id <- 1
  #   time.in1 <- Sys.time()
  #   test <- sapply(1:100000, function(id) {as.numeric(names(ids.state.list)[vapply(ids.state.list, is.element, el = id, FUN.VALUE = FALSE)])})
  #   time.out1 <- Sys.time()
  #   time.out1 - time.in1

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- 1:nrow(data.raw)
  m1 <- outer(v1, ids.state.list, FUN = Vectorize('%in%'))
  state.poly <- lapply(split(m1, row(m1)), function(x) (1:max.state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state.poly, length)
  state.poly[idx] <- NA

  ## Add to data.raw
  data.raw <- mutate(data.raw, state.poly = unlist(state.poly),
                     state.poly.fac = factor(state.poly))

  ### Create binary variables for each possible state that can be transitioned to
  temp.dummy <- to_dummy(data.raw, state.poly.fac)
  colnames(temp.dummy) <- paste("state", valid.transitions, ".bin", sep = "")

  ### Add to dataset
  data.raw <- cbind(data.raw, temp.dummy)
  rm(temp.dummy)

  ### Reduce data.raw to individuals who are uncensored at time t.eval
  data.raw.uncens <- subset(data.raw, !is.na(state.poly))

  ###
  ### Calculate observed risks
  ###

  ### Create function which we will apply boot to if calculating confidence intervals
  calc.obs.loess.func <- function(data.in, indices, state.k, data.in.uncens){

    ### Create bootstrapped dataset
    data.boot <- data.in[indices, c("id", paste("p.est", state.k, sep = ""), paste("state", state.k, ".bin", sep = ""), "state.poly", "ipcw")]

    ### Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- subset(data.boot, !is.na(state.poly))

    ### Define equation
    eq.loess <- formula(paste("state", state.k, ".bin ~ p.est", state.k, sep = ""))

    ### Fit model
    if (stabilised == FALSE){
      loess.model <- loess(eq.loess,
                           data = data.boot.uncens,
                           weights = data.boot.uncens[, "ipcw"],
                           span = span,
                           degree = degree)
    } else if (stabilised == TRUE){
      loess.model <- loess(eq.loess,
                           data = data.boot.uncens,
                           weights = data.boot.uncens[, "ipcw.stab"],
                           span = span,
                           degree = degree)
    }

    ### Create predictions for the vector of predicted probabilities for people included in original the calibration curve (this is the individuals
    ### who are in the unbootstrapped dataset, and are uncensored at time t.eval. This would be the vector of predicted probabilities for the
    ### calibration curve if no bootstrapping was done)

    ### Create predicted observed risks for these individuals
    loess.pred.obs <- predict(loess.model, newdata = data.in.uncens)

    return(loess.pred.obs)
  }

  ### Create object to store output
  output.object <- vector("list", length(valid.transitions))

  ### Loop through and fit models
  for (k in 1:length(valid.transitions)){

    ### Assign state of interest
    state.k <- valid.transitions[k]

    ### Calculate predicted observed probabilities
    loess.pred.obs <- calc.obs.loess.func(data.raw, 1:nrow(data.raw), state.k, data.raw.uncens)

    ### Assign output
    output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                     "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                     "obs" = loess.pred.obs)

    ### Calculate confidence intervals
    if (CI != FALSE){

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Run bootstrapping
      boot.obs <- boot(data.raw, calc.obs.loess.func, R = R.boot, state.k = valid.transitions[k], data.in.uncens = data.raw.uncens)$t

      ### Extract confidence bands
      lower <- apply(boot.obs, 2, quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.obs, 2, quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot.obs)) > 0){
        print(paste("WARNING, SOME BOOTSTRAPPED OBSERVED RISKS WERE NA FOR STATE", valid.transitions[k]))
        print(paste("THERE ARE ", sum(apply(boot.obs, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's FOR OBSERVED RISKS"))
        print(paste("THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs, 1, function(x) {sum(is.na(x))}))))
      }

      ### Assign output
      output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                       "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                       "obs" = loess.pred.obs,
                                       "obs.lower" = lower,
                                       "obs.upper" = upper)
    }
  }

  ### Create metadata object
  metadata <- list("valid.transitions" = valid.transitions, "CI" = CI)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plot.data" = output.object, "metadata" = metadata)
  return(output.object.comb)
}


###
### 2.2) Calculate calibration using binary logistic regression framework with inverse probability of censoring weights.
### Calibration model uses restricted cubic splines.
###
calc.calib.blr.rcs <- function(data.mstate, data.raw, tmat, j, s, t.eval, p.est, nk, CI = FALSE, R.boot, stabilised = FALSE){

#       ### Calculate ipcw weights
#       weights <- calc.weights(data.mstate = msebmt,
#                               data.raw = ebmt,
#                               covs = covs,
#                               j = j,
#                               landmark.type = "all",
#                               s = s,
#                               t.eval = t.eval)
#
#       ### Combineweights with ebmt dataset
#       ebmt.weights <- cbind(ebmt, weights)
#
#       data.mstate <- msebmt
#       data.raw <- ebmt.weights
#       tmat
#       j<-1
#       s<-s
#       t.eval <- t.eval
#       p.est <- tp.all[[1]] %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
#       span <- span
#       degree <- degree.blr
#       CI <- 95
#       R.boot <- 500
#       nk <- 3

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(p.est) != 0)

  ### Add the predicted risks, and the logit transormation of the predicted risks to data.raw
  p.est.logit <- log(p.est[,valid.transitions]/(1-p.est[,valid.transitions]))
  colnames(p.est.logit) <- paste("p.est.logit", valid.transitions, sep = "")
  data.raw <- data.frame(data.raw, p.est[,valid.transitions], p.est.logit)

  ### Identify individuals who are in state j at time s
  ids.state.j <- subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    select(id) %>%
    distinct(id) %>%
    pull(id)

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate %>% subset(id %in% ids.state.j)
  data.raw <- data.raw %>% subset(id %in% ids.state.j)

  ### Extract which state individuals are in at time t.eval
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract.ids.states(data.mstate, tmat, k, t.eval)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- 1:nrow(data.raw)
  m1 <- outer(v1, ids.state.list, FUN = Vectorize('%in%'))
  state.poly <- lapply(split(m1, row(m1)), function(x) (1:max.state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state.poly, length)
  state.poly[idx] <- NA

  ## Add to data.raw
  data.raw <- mutate(data.raw, state.poly = unlist(state.poly),
                     state.poly.fac = factor(state.poly))

  ### Create binary variables for each possible state that can be transitioned to
  temp.dummy <- to_dummy(data.raw, state.poly.fac)
  colnames(temp.dummy) <- paste("state", valid.transitions, ".bin", sep = "")

  ### Add to dataset
  data.raw <- cbind(data.raw, temp.dummy)
  rm(temp.dummy)

  ### Reduce data.raw to individuals who are uncensored at time t.eval
  data.raw.uncens <- subset(data.raw, !is.na(state.poly))

  ###
  ### Calculate predicted-observed risks
  ###

  ### Create function which we will apply boot to if calculating confidence intervals
  calc.obs.rcs.func <- function(data.in, indices, state.k, data.in.uncens){

    ### Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ### Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- subset(data.boot, !is.na(state.poly))

    ### Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- rcspline.eval(data.boot.uncens[,paste("p.est.logit", state.k, sep = "")],nk=nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ### Add the cubic splines for logit of the predicted probability to data.boot.uncens
    data.boot.uncens <- data.frame(cbind(data.boot.uncens, rcs.mat))

    ### Define equation
    eq.LHS <- paste("state", state.k, ".bin ~ ", sep = "")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq.rcs <- formula(paste(eq.LHS, eq.RHS, sep = ""))

    ### Fit model
    if (stabilised == FALSE){
      rcs.model <- lrm(eq.rcs,
                       data = data.boot.uncens,
                       weights = data.boot.uncens[, "ipcw"])
    } else if (stabilised == TRUE){
      rcs.model <- lrm(eq.rcs,
                       data = data.boot.uncens,
                       weights = data.boot.uncens[, "ipcw.stab"])
    }


    ### Create predicted observed probabilities. Create predictions for the vector of predicted probabilities form the original model
    ### So that we always plot over the same range.

    ### For this, need to calculate the correct splines for the original dataset.
    rcs.mat.data.in.uncens <- rcspline.eval(data.in.uncens[,paste("p.est.logit", state.k, sep = "")],knots = knots,inclx=T)
    colnames(rcs.mat.data.in.uncens) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    #attr(rcs.mat.data.in.uncens,"knots")

    ### Add the cubic splines for logit of the predicted probability to data.in.uncens
    data.in.uncens <- data.frame(cbind(data.in.uncens, rcs.mat.data.in.uncens))

    ### Calculate predicted observed probabilities
    rcs.pred.obs <- predict(rcs.model, newdata = data.in.uncens, type = "fitted")

    return(rcs.pred.obs)
  }

  ### Create object to store output
  output.object <- vector("list", length(valid.transitions))

  ### Loop through and fit models
  for (k in 1:length(valid.transitions)){

    ### Assign state of interest
    state.k <- valid.transitions[k]

    ### Calculate predicted observed probabilities
    rcs.pred.obs <- calc.obs.rcs.func(data.raw, 1:nrow(data.raw), state.k, data.raw.uncens)

    ### Assign output
    output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                     "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                     "obs" = rcs.pred.obs)

    ### Calculate confidence intervals
    if (CI != FALSE){

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Run bootstrapping
      boot.obs <- boot(data.raw, calc.obs.rcs.func, R = R.boot, state.k = valid.transitions[k], data.in.uncens = data.raw.uncens)$t

      ### Extract confidence bands
      lower <- apply(boot.obs, 2, quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.obs, 2, quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot.obs)) > 0){
        print(paste("WARNING, SOME BOOTSTRAPPED OBSERVED RISKS WERE NA FOR STATE", valid.transitions[k]))
        print(paste("THERE ARE ", sum(apply(boot.obs, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's FOR OBSERVED RISKS"))
        print(paste("THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs, 1, function(x) {sum(is.na(x))}))))
      }

      ### Assign output
      output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                       "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                       "obs" = rcs.pred.obs,
                                       "obs.lower" = lower,
                                       "obs.upper" = upper)
    }
  }

  ### Create metadata object
  metadata <- list("valid.transitions" = valid.transitions, "CI" = CI)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plot.data" = output.object, "metadata" = metadata)
  return(output.object.comb)
}


###
### 2.3) Calculate calibration using multinomial logistic regression framework with inverse probability of censoring weights.
### Calibration model uses nominal recalibration framework of van Hoorde.
###
calc.calib.mlr <- function(data.mstate, data.raw, tmat, j, s, t.eval, p.est, ps.int = 4, degree = 3){

      ### Calculate ipcw weights
      weights <- calc.weights(data.mstate = msebmtcal,
                              data.raw = ebmtcal,
                              covs = covs,
                              j = j,
                              landmark.type = "all",
                              s = s,
                              t.eval = t.eval)

      ### Combineweights with ebmt dataset
      ebmt.weights <- cbind(ebmtcal, weights)

      data.mstate <- msebmt
      data.raw <- ebmt.weights
      tmat
      j<-j
      s<-s
      t.eval <- t.eval
      p.est <- tp.all[[1]] %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
      ps.int <- 4
      degree <- 3

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(p.est) != 0)

  ### Add linear predictors from a multinomial framework
  ### Start by reducing to non-zero columns
  p.est.mlr <- p.est[,valid.transitions]

  ### Calculate linear predictors
  p.est.mlr <- log(p.est.mlr[,2:ncol(p.est.mlr)]/p.est.mlr[,1])
  colnames(p.est.mlr) <- paste("mlr.lp", 1:(ncol(p.est.mlr)), sep = "")

  ### Add to data frame
  data.raw <- data.frame(data.raw, p.est, p.est.mlr)

  ### Identify individuals who are in state j at time s
  ids.state.j <- subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    select(id) %>%
    distinct(id) %>%
    pull(id)

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate %>% subset(id %in% ids.state.j)
  data.raw <- data.raw %>% subset(id %in% ids.state.j)

  ### Extract which state individuals are in at time t.eval
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract.ids.states(data.mstate, tmat, k, t.eval)
  }

  ### Create a function which identifies which states an individual is in
  identify.state <- function(id){
    if (sum(lapply(ids.state.list, function(x){id %in% x}) == TRUE) == 1){
      group <- which(lapply(ids.state.list, function(x){id %in% x}) == TRUE)
    } else if (sum(lapply(ids.state.list, function(x){id %in% x}) == TRUE) == 0){
      group <- NA
    }
    return(group)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  data.raw <- data.raw %>%
    mutate(state.poly = sapply(id, identify.state),
           state.poly.fac = factor(state.poly))

  ### Reduce data.raw to individuals who are uncensored at time t.eval
  data.raw.uncens <- subset(data.raw, !is.na(state.poly))

  ### Define equation
  eq.LHS <- paste("state.poly.fac ~ ")
  eq.RHS <- paste("sm.ps(mlr.lp", 1:ncol(p.est.mlr), ", ps.int = ps.int, degree = degree)", sep = "", collapse = "+")
  eq.mlr <- as.formula(paste(eq.LHS, eq.RHS, sep =""))

  ### Assign reference category
  ref.cat <- paste(valid.transitions[1])

  ### Apply nominal recalibration framework with vector spline smoothers
  calib.model <- vgam(eq.mlr, weights = data.raw.uncens[, "ipcw"],
                      data = data.raw.uncens, family = multinomial(refLevel = ref.cat))

  ###
  ### Generate predicted-observed risks and add to data.raw.uncens
  ###

  ### For all other functions, I just generate predicted risks for all individuals, and then just plot for those who were uncensored
  ### However, some of the censored individuals are causing an error, therefore I must generate NA vectors, then assign the
  ### predicted observed probabilities to the crrect individuals

  ## Create dataframe to store
  dat.mlr.pred.obs <- data.frame(matrix(NA, ncol = length(valid.transitions), nrow = nrow(data.raw.uncens)))
  ## Assign colnames
  colnames(dat.mlr.pred.obs) <- paste("mlr.pred.obs", valid.transitions, sep = "")
  ## Calc pred.obs for those who are uncesored
  mlr.pred.obs <- predict(calib.model, newdata = data.raw.uncens, type = "response")
  ## Assign to appropriate individuals
  dat.mlr.pred.obs[, ] <- mlr.pred.obs

  ### Then add it to data.raw.uncens
  data.raw.uncens <- cbind(data.raw.uncens, dat.mlr.pred.obs)

  ### Assign output
  output.object <- select(data.raw.uncens, id, paste("p.est", valid.transitions, sep = ""), paste("mlr.pred.obs", valid.transitions, sep = ""))

  ### Create metadata object
  metadata <- list("valid.transitions" = valid.transitions)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plot.data" = output.object, "metadata" = metadata)
  return(output.object.comb)

}


###
### 2.4) Calculate calibration curves for a competing risks sub-model using graphical calibration curves (Austin et al.)
###
calc.calib.cmprsk.mod <- function(data.mstate, data.raw, j, s, t.eval, p.est, nk = 3){

  #       data.mstate <- msebmt
  #       data.raw <- ebmt
  #       p.est <- tp.all[,paste("pstate", 1:6, sep = "")]
  #       j <- 3
  #       s <- 30
  #       t.eval <- 1826
  #       nk <- 3
  #
  ###
  ### Create a plot for each possible transition seperately
  ### Note there it is direct transitions (i.e. from state 1 you cannot transition to state 4 directly,
  ### so there will not be a calibration plot for this, even though there was a plot for 1 -> 4 for the transitio probabilities.)
  ###

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  ### Also drop the staet that an individual is already in, because there is no cmprsk model for staying in the same state
  valid.transitions <- which(colSums(p.est) != 0)
  valid.transitions <- valid.transitions[-(valid.transitions == j)]

  ### Add the predicted risks, and the complementary log log transormation of the predicted risks to data.raw
  p.est.cll <- log(-log(1 - p.est[,valid.transitions]))
  colnames(p.est.cll) <- paste("p.est.cll", valid.transitions, sep = "")

  ### Identify individuals who are in state j at time s
  ids.state.j <- subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    select(id) %>%
    distinct(id) %>%
    pull(id)

  ### Subset data.mstate and data.raw to these individuals.
  data.mstate <- data.mstate %>% subset(id %in% ids.state.j)
  data.raw <- data.raw %>% subset(id %in% ids.state.j)

  ### Add the cloglog risks and predicted risks to the landmark dataset
  data.raw <- cbind(data.raw, p.est[,valid.transitions], p.est.cll)

  ### Finally, identify individuals which are censored before experiencing any events (used to maniuplate data for Fine-Gray regression later)
  ids.cens  <- data.mstate %>% subset(from == j) %>% group_by(id) %>% dplyr::summarize(sum = sum(status)) %>% subset(sum == 0) %>% pull(id)

  ###
  ### Produce calibration plots for each possible transition
  ###

  ### Start by creating a list to store the plots
  plots.list <- vector("list", length(valid.transitions))

  for (k in 1:length(valid.transitions)){

    ### Assign state.k
    state.k <- as.numeric(valid.transitions[k])

    ### Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- rcspline.eval(data.raw[,paste("p.est.cll", state.k, sep = "")],nk=nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ### Create a new dataframe for the validation, to avoid recurison with data.raw
    ### Add the cubic splines for thecomplementary loglog of the predicted probability, and the predicted probability itself
    valid.df <- data.frame(data.raw$id, data.raw[,paste("p.est", state.k, sep = "")], rcs.mat)
    colnames(valid.df) <- c("id", "pred", colnames(rcs.mat))

    ### Want to validate the competing risks model out of state j at time s, into state k, so remove individuals not in state k at time s,
    ### and only retain transitions into state k. Also deduct immortal time from time variable
    data.mstate.j.k.s <- subset(data.mstate, from == j & to == state.k & Tstart <= s & s < Tstop) %>%
      mutate(time = Tstop - s) %>%
      select(c(time, status))

    ### Add to valid.df
    valid.df <- cbind(valid.df, data.mstate.j.k.s)

    ### For individuals who do not have the event of interest, and also are not censored (i.e. they have a different competing event),
    ### set the follow up time to the maximum
    valid.df <- mutate(valid.df, time = case_when(status == 0 & !(id %in% ids.cens) ~ max(time),
                                                  TRUE ~ time))

    ### Create dataset to fit the recalibration model of Austin et al (Graphical calibration curves, BMC Diagnostic and Prognostic, DOI10.1186/s41512-021-00114-6)
    valid.df.crprep <- crprep(Tstop='time',status='status',trans=1,
                              keep=colnames(rcs.mat),valid.df)

    ### Create formula and fit the Fine-Gray recalibration model
    eq.LHS <- paste("Surv(Tstart,Tstop,status==1)~")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq <- formula(paste(eq.LHS, eq.RHS, sep = ""))
    model.calibrate.fg <- cph(eq,weights=weight.cens,x=T,y=T,surv=T,data=valid.df.crprep)

    ### Generate predicted probabilities and standard errors
    valid.df$obs.fg <- 1-survest(model.calibrate.fg,newdata=valid.df.crprep,time=t.eval-s)$surv
    valid.df$obs.fg.upper<-1-survest(model.calibrate.fg,newdata=valid.df.crprep,time=t.eval-s)$lower
    valid.df$obs.fg.lower<-1-survest(model.calibrate.fg,newdata=valid.df.crprep,time=t.eval-s)$upper

    ### Produce plots for each and store in a list

    ### Pivot longer to create data for ggplot and assign appropriate labels
    valid.df.longer <- pivot_longer(valid.df, cols = c(obs.fg, obs.fg.upper, obs.fg.lower), names_to = "line.group")
    valid.df.longer <- mutate(valid.df.longer,
                              line.group = factor(line.group),
                              mapping = case_when(line.group == "obs.fg" ~ 1,
                                                  line.group %in% c("obs.fg.upper", "obs.fg.lower") ~ 2),
                              mapping = factor(mapping))

    levels(valid.df.longer$line.group) <- c("Calibration", "Upper", "Lower")
    levels(valid.df.longer$mapping) <- c("Calibration", "95% CI")

    ### Create the plot
    plots.list[[k]] <- ggplot(data = valid.df.longer %>% arrange(pred) %>% select(id, pred, line.group, value, mapping)) +
      geom_line(aes(x = pred, y = value, group = line.group, color = mapping)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      xlab("Predicted risk") + ylab("Observed risk") +
      xlim(c(0, max(valid.df.longer$pred,
                    valid.df.longer$value))) +
      ylim(c(0, max(valid.df.longer$pred,
                    valid.df.longer$value))) +
      geom_rug(data = valid.df.longer %>% arrange(pred) %>% select(id, pred, line.group, value, mapping) %>% subset(line.group == "Calibration"),
               aes(x = pred, y = value), col = rgb(1, 0, 0, alpha = .1)) +
      theme(legend.position = "none") +
      ggtitle(paste("State ", state.k, sep = ""))

  }

  ### Return plots
  return(plots.list)
}


##########################################################################
###                                                                    ###
### Section 3) Functions to create plots from calibration data         ###
###                                                                    ###
##########################################################################

###
### 3.1) Create plots for calibration curves from binary logistic regression framework
###
plot.calib.blr <- function(object.in, combine = FALSE, ncol, nrow){

  ### Extract plot dat aand relevant metadata
  plot.data <- object.in[["plot.data"]]
  valid.transitions <- object.in[["metadata"]][["valid.transitions"]]
  CI <- object.in[["metadata"]][["CI"]]

  if (CI != FALSE){
    ### Create list to store plots
    plots.list <- vector("list", length(valid.transitions))

    for (k in 1:length(valid.transitions)){
      ### Assign plot data
      plot.data.k <- plot.data[[k]]

      ### Assign state of interest
      state.k <- valid.transitions[k]

      ### Pivot longer to create data for ggplot and assign appropriate labels
      plot.data.k.longer <- pivot_longer(plot.data.k, cols = c(obs, obs.upper, obs.lower), names_to = "line.group")
      plot.data.k.longer <- mutate(plot.data.k.longer,
                                   line.group = factor(line.group),
                                   mapping = case_when(line.group == "obs" ~ 1,
                                                       line.group %in% c("obs.upper", "obs.lower") ~ 2),
                                   mapping = factor(mapping))

      levels(plot.data.k.longer$line.group) <- c("Calibration", "Upper", "Lower")
      levels(plot.data.k.longer$mapping) <- c("Calibration", "95% CI")

      ### Create the plots
      plots.list[[k]] <- ggplot(data = plot.data.k.longer %>% arrange(pred) %>% select(id, pred, line.group, value, mapping)) +
        geom_line(aes(x = pred, y = value, group = line.group, color = mapping)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        xlab("Predicted risk") + ylab("Observed risk") +
#                 xlim(c(0, max(plot.data.k.longer$pred,
#                               plot.data.k.longer$value))) +
#                 ylim(c(0, max(plot.data.k.longer$pred,
#                               plot.data.k.longer$value))) +
  xlim(c(0, 0.5)) +
  ylim(c(0, 0.5)) +
#         geom_rug(data = plot.data.k.longer %>% arrange(pred) %>% select(id, pred, line.group, value, mapping) %>% subset(line.group == "Calibration"),
#                  aes(x = pred, y = value), col = rgb(1, 0, 0, alpha = .1)) +
        geom_rug(data = plot.data.k.longer %>% arrange(pred) %>% select(id, pred, line.group, value, mapping) %>% subset(line.group == "Calibration"),
                 aes(x = pred, y = value), col = rgb(1, 0, 0)) +
        ggtitle(paste("State ", state.k, sep = ""))
    }
  } else if (CI == FALSE){
    ### Create list to store plots
    plots.list <- vector("list", length(valid.transitions))
    for (k in 1:length(valid.transitions)){

      ### Assign plot data
      plot.data.k <- plot.data[[k]]

      ### Assign state of interest
      state.k <- valid.transitions[k]

      ### Create the plots
      plots.list[[k]] <- ggplot(data = plot.data.k %>% arrange(pred) %>% select(id, pred, obs)) +
        geom_line(aes(x = pred, y = obs)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        xlab("Predicted risk") + ylab("Observed risk") +
#         xlim(c(0, max(plot.data.k$pred,
#                       plot.data.k$obs))) +
#         ylim(c(0, max(plot.data.k$pred,
#                       plot.data.k$obs))) +
        xlim(c(0, 0.5)) +
        ylim(c(0, 0.5)) +
        #geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .1)) +
        geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0)) +
        theme(legend.position = "none") +
        ggtitle(paste("State ", state.k, sep = ""))
    }
  }

  ### Combine plots into single ggplot
  if (combine == TRUE){
    plots.list <- ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol, common.legend = TRUE)
  }

  ### Return output object
  return(plots.list)

}


###
### 3.2) Create plots for calibration curves from binary logistic regression framework
###
plot.calib.mlr <- function(object.in, combine = FALSE, ncol, nrow){

  ### Extract plot dat aand relevant metadata
  plot.data <- object.in[["plot.data"]]
  valid.transitions <- object.in[["metadata"]][["valid.transitions"]]

  ### Create list to store plots
  plots.list <- vector("list", length(valid.transitions))
  for (k in 1:length(valid.transitions)){

    ### Assign state of interest
    state.k <- valid.transitions[k]

    ### Create variables to plot
    plot.data$pred <- plot.data[, paste("p.est", state.k, sep = "")]
    plot.data$obs <- plot.data[, paste("mlr.pred.obs", state.k, sep = "")]

    ### Create the plots
    plots.list[[k]] <- ggplot(data = plot.data  %>% arrange(pred) %>%  select(id, pred, obs)) +
      #geom_point(aes(x = pred, y = obs), color = "red", size = 0.5, alpha = .25) +
      geom_point(aes(x = pred, y = obs), color = "red", size = 0.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      xlab("Predicted risk") + ylab("Observed risk") +
      xlim(c(0, max(plot.data$pred,
                    plot.data$obs))) +
      ylim(c(0, max(plot.data$pred,
                    plot.data$obs))) +
      #geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0, alpha = .1)) +
      geom_rug(aes(x = pred, y = obs), col = rgb(1, 0, 0)) +
      theme(legend.position = "none") +
      ggtitle(paste("State ", state.k, sep = ""))

  }

  ### Combine plots into single ggplot
  if (combine == TRUE){
    plots.list <- ggarrange(plotlist = plots.list, nrow = nrow, ncol = ncol)
  }

  ### Return output object
  return(plots.list)

}


##########################################################################
###                                                                    ###
### Section 2=TEMP) Functions to create calibration curves or scatter plots ###
### LIVE EDITS OF EXISTING FUNCTIONS ABOVE
###                                                                    ###
##########################################################################

###
### 2.1) Calculate calibration using binary logistic regression framework with inverse probability of censoring weights.
### Calibration model uses loess smoothers.
###
calc.calib.blr.loess.ihw <- function(data.mstate, data.raw, tmat, j, s, t.eval, p.est, span, degree, CI = FALSE, R.boot, stabilised = FALSE){

#     ### Choose span and degree for loess plots
#     span <- 1
#     degree.blr <- 2
#
#     ### Choose ps.int and degree for mlr plots
#     ps.int <- 4
#     degree.mlr <- 3
#
#     ### Choose number of knots for rcs
#     nk <- 3
#
#     data.mstate = msebmt
#     data.raw = ebmt
#     tmat
#     j=1
#     s=s
#     t.eval = t.eval
#     p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = "")))
#     span = span
#     degree = degree.blr
#     CI = 95
#     R.boot = 25
#     stabilised = FALSE
#     head(p.est)

  ### Extract transition matrix
  tmat <- attributes(data.mstate)[["trans"]]

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(p.est) != 0)

  ### Add the predicted risks, and the logit transormation of the predicted risks to data.raw
  p.est.logit <- log(p.est[,valid.transitions]/(1-p.est[,valid.transitions]))
  colnames(p.est.logit) <- paste("p.est.logit", 1:ncol(p.est.logit), sep = "")
  data.raw <- data.frame(data.raw, p.est[,valid.transitions], p.est.logit)

  ### Identify individuals who are in state j at time s
  ids.state.j <- subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    select(id) %>%
    distinct(id) %>%
    pull(id)

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate %>% subset(id %in% ids.state.j)
  data.raw <- data.raw %>% subset(id %in% ids.state.j)

  ### Extract which state individuals are in at time t.eval
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract.ids.states(data.mstate, tmat, k, t.eval)
  }

  #   ### Create a function which identifies which states an individual is in
  #   identify.state <- function(id){
  #     if (sum(lapply(ids.state.list, function(x){id %in% x}) == TRUE) == 1){
  #       group <- which(lapply(ids.state.list, function(x){id %in% x}) == TRUE)
  #     } else if (sum(lapply(ids.state.list, function(x){id %in% x}) == TRUE) == 0){
  #       group <- NA
  #     }
  #     return(group)
  #   }
  #   names(ids.state.list) <- 1:6
  #   id <- 1
  #   time.in1 <- Sys.time()
  #   test <- sapply(1:100000, function(id) {as.numeric(names(ids.state.list)[vapply(ids.state.list, is.element, el = id, FUN.VALUE = FALSE)])})
  #   time.out1 <- Sys.time()
  #   time.out1 - time.in1

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- 1:nrow(data.raw)
  m1 <- outer(v1, ids.state.list, FUN = Vectorize('%in%'))
  state.poly <- lapply(split(m1, row(m1)), function(x) (1:max.state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state.poly, length)
  state.poly[idx] <- NA

  ## Add to data.raw
  data.raw <- mutate(data.raw, state.poly = unlist(state.poly),
                     state.poly.fac = factor(state.poly))

  ### Create binary variables for each possible state that can be transitioned to
  temp.dummy <- to_dummy(data.raw, state.poly.fac)
  colnames(temp.dummy) <- paste("state", valid.transitions, ".bin", sep = "")

  ### Add to dataset
  data.raw <- cbind(data.raw, temp.dummy)
  rm(temp.dummy)

  ### Create dataset of individuals who are uncensored at time t.eval
  data.raw.uncens <- subset(data.raw, !is.na(state.poly))

  ###
  ### Calculate observed risks
  ###

  ### Create function which we will apply boot to if calculating confidence intervals
  calc.obs.loess.func <- function(data.in, indices, state.k, data.in.uncens){

    ### Create bootstrapped dataset
    data.boot <- data.in[indices, c("id", paste("p.est", state.k, sep = ""), paste("state", state.k, ".bin", sep = ""), "state.poly", covs, "dtcens", "dtcens.s")]

    ### Calculate weights
    weights <- calc.weights(data.mstate, data.boot, covs, j = j, landmark.type = "all", s, t.eval, max.weight = 10, stabilised = stabilised)

    ### Add to data.boot
    data.boot$ipcw <- weights$ipcw
    if (stabilised == TRUE){
      data.boot$ipcw.stab <- weights$ipcw.stab
    }

    ### Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- subset(data.boot, !is.na(state.poly))

    ### Define equation
    eq.loess <- formula(paste("state", state.k, ".bin ~ p.est", state.k, sep = ""))

    ### Fit model
    if (stabilised == FALSE){
      loess.model <- loess(eq.loess,
                           data = data.boot.uncens,
                           weights = data.boot.uncens[, "ipcw"],
                           span = span,
                           degree = degree)
    } else if (stabilised == TRUE){
      loess.model <- loess(eq.loess,
                           data = data.boot.uncens,
                           weights = data.boot.uncens[, "ipcw.stab"],
                           span = span,
                           degree = degree)
    }

    ### Create predictions for the vector of predicted probabilities for people included in original the calibration curve (this is the individuals
    ### who are in the unbootstrapped dataset, and are uncensored at time t.eval. This would be the vector of predicted probabilities for the
    ### calibration curve if no bootstrapping was done)

    ### Create predicted observed risks for these individuals
    loess.pred.obs <- predict(loess.model, newdata = data.in.uncens)

    return(loess.pred.obs)
  }

  ### Create object to store output
  output.object <- vector("list", length(valid.transitions))

  ### Loop through and fit models
  for (k in 1:length(valid.transitions)){

    ### Assign state of interest
    state.k <- valid.transitions[k]

    ### Calculate predicted observed probabilities
    loess.pred.obs <- calc.obs.loess.func(data.raw, 1:nrow(data.raw), state.k, data.raw.uncens)

    ### Assign output
    output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                     "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                     "obs" = loess.pred.obs)

    ### Calculate confidence intervals
    if (CI != FALSE){

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Run bootstrapping
      boot.obs <- boot(data.raw, calc.obs.loess.func, R = 20, state.k = valid.transitions[k], data.in.uncens = data.raw.uncens)$t

      ### Extract confidence bands
      lower <- apply(boot.obs, 2, quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.obs, 2, quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot.obs)) > 0){
        print(paste("WARNING, SOME BOOTSTRAPPED OBSERVED RISKS WERE NA FOR STATE", valid.transitions[k]))
        print(paste("THERE ARE ", sum(apply(boot.obs, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's FOR OBSERVED RISKS"))
        print(paste("THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs, 1, function(x) {sum(is.na(x))}))))
      }

      ### Assign output
      output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                       "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                       "obs" = loess.pred.obs,
                                       "obs.lower" = lower,
                                       "obs.upper" = upper)
    }
  }

  ### Create metadata object
  metadata <- list("valid.transitions" = valid.transitions, "CI" = CI)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plot.data" = output.object, "metadata" = metadata)
  return(output.object.comb)
}






###
### 2.2) Calculate calibration using binary logistic regression framework with inverse probability of censoring weights.
### Calibration model uses restricted cubic splines.
###
calc.calib.blr.rcs.ihw <- function(data.mstate, data.raw, tmat, j, s, t.eval, p.est, curve.nk, CI = FALSE, R.boot, stabilised = FALSE){

  #       ### Calculate ipcw weights
  #       weights <- calc.weights(data.mstate = msebmt,
  #                               data.raw = ebmt,
  #                               covs = covs,
  #                               j = j,
  #                               landmark.type = "all",
  #                               s = s,
  #                               t.eval = t.eval)
  #
  #       ### Combineweights with ebmt dataset
  #       ebmt.weights <- cbind(ebmt, weights)
  #
  #       data.mstate <- msebmt
  #       data.raw <- ebmt.weights
  #       tmat
  #       j<-1
  #       s<-s
  #       t.eval <- t.eval
  #       p.est <- tp.all[[1]] %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  #       span <- span
  #       degree <- degree.blr
  #       CI <- 95
  #       R.boot <- 500
  #       curve.nk <- 3

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(p.est) != 0)

  ### Add the predicted risks, and the logit transormation of the predicted risks to data.raw
  p.est.logit <- log(p.est[,valid.transitions]/(1-p.est[,valid.transitions]))
  colnames(p.est.logit) <- paste("p.est.logit", valid.transitions, sep = "")
  data.raw <- data.frame(data.raw, p.est[,valid.transitions], p.est.logit)

  ### Identify individuals who are in state j at time s
  ids.state.j <- subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    select(id) %>%
    distinct(id) %>%
    pull(id)

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate %>% subset(id %in% ids.state.j)
  data.raw <- data.raw %>% subset(id %in% ids.state.j)

  ### Extract which state individuals are in at time t.eval
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract.ids.states(data.mstate, tmat, k, t.eval)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- 1:nrow(data.raw)
  m1 <- outer(v1, ids.state.list, FUN = Vectorize('%in%'))
  state.poly <- lapply(split(m1, row(m1)), function(x) (1:max.state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state.poly, length)
  state.poly[idx] <- NA

  ## Add to data.raw
  data.raw <- mutate(data.raw, state.poly = unlist(state.poly),
                     state.poly.fac = factor(state.poly))

  ### Create binary variables for each possible state that can be transitioned to
  temp.dummy <- to_dummy(data.raw, state.poly.fac)
  colnames(temp.dummy) <- paste("state", valid.transitions, ".bin", sep = "")

  ### Add to dataset
  data.raw <- cbind(data.raw, temp.dummy)
  rm(temp.dummy)

  ### Reduce data.raw to individuals who are uncensored at time t.eval
  data.raw.uncens <- subset(data.raw, !is.na(state.poly))

  ###
  ### Calculate predicted-observed risks
  ###

  ### Create function which we will apply boot to if calculating confidence intervals
  calc.obs.rcs.func <- function(data.in, indices, state.k, data.in.uncens){

    ### Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ### Calculate weights
    weights <- calc.weights(data.mstate, data.boot, covs, j = j, landmark.type = "all", s, t.eval, max.weight = 10, stabilised = stabilised)

    ### Add to data.boot
    data.boot$ipcw <- weights$ipcw
    if (stabilised == TRUE){
      data.boot$ipcw.stab <- weights$ipcw.stab
    }

    ### Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- subset(data.boot, !is.na(state.poly))

    ### Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- rcspline.eval(data.boot.uncens[,paste("p.est.logit", state.k, sep = "")],nk=curve.nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ### Add the cubic splines for logit of the predicted probability to data.boot.uncens
    data.boot.uncens <- data.frame(cbind(data.boot.uncens, rcs.mat))

    ### Define equation
    eq.LHS <- paste("state", state.k, ".bin ~ ", sep = "")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq.rcs <- formula(paste(eq.LHS, eq.RHS, sep = ""))

    ### Fit model
    if (stabilised == FALSE){
      rcs.model <- lrm(eq.rcs,
                       data = data.boot.uncens,
                       weights = data.boot.uncens[, "ipcw"])
    } else if (stabilised == TRUE){
      rcs.model <- lrm(eq.rcs,
                       data = data.boot.uncens,
                       weights = data.boot.uncens[, "ipcw.stab"])
    }


    ### Create predicted observed probabilities. Create predictions for the vector of predicted probabilities form the original model
    ### So that we always plot over the same range.

    ### For this, need to calculate the correct splines for the original dataset.
    rcs.mat.data.in.uncens <- rcspline.eval(data.in.uncens[,paste("p.est.logit", state.k, sep = "")],knots = knots,inclx=T)
    colnames(rcs.mat.data.in.uncens) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    #attr(rcs.mat.data.in.uncens,"knots")

    ### Add the cubic splines for logit of the predicted probability to data.in.uncens
    data.in.uncens <- data.frame(cbind(data.in.uncens, rcs.mat.data.in.uncens))

    ### Calculate predicted observed probabilities
    rcs.pred.obs <- predict(rcs.model, newdata = data.in.uncens, type = "fitted")

    return(rcs.pred.obs)
  }

  ### Create object to store output
  output.object <- vector("list", length(valid.transitions))

  ### Loop through and fit models
  for (k in 1:length(valid.transitions)){

    ### Assign state of interest
    state.k <- valid.transitions[k]

    ### Calculate predicted observed probabilities
    rcs.pred.obs <- calc.obs.rcs.func(data.raw, 1:nrow(data.raw), state.k, data.raw.uncens)

    ### Assign output
    output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                     "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                     "obs" = rcs.pred.obs)

    ### Calculate confidence intervals
    if (CI != FALSE){

      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Run bootstrapping
      boot.obs <- boot(data.raw, calc.obs.rcs.func, R = R.boot, state.k = valid.transitions[k], data.in.uncens = data.raw.uncens)$t

      ### Extract confidence bands
      lower <- apply(boot.obs, 2, quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.obs, 2, quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot.obs)) > 0){
        print(paste("WARNING, SOME BOOTSTRAPPED OBSERVED RISKS WERE NA FOR STATE", valid.transitions[k]))
        print(paste("THERE ARE ", sum(apply(boot.obs, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's FOR OBSERVED RISKS"))
        print(paste("THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs, 1, function(x) {sum(is.na(x))}))))
      }

      ### Assign output
      output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                       "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                       "obs" = rcs.pred.obs,
                                       "obs.lower" = lower,
                                       "obs.upper" = upper)
    }
  }

  ### Create metadata object
  metadata <- list("valid.transitions" = valid.transitions, "CI" = CI)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plot.data" = output.object, "metadata" = metadata)
  return(output.object.comb)
}











###
### 2.2) Calculate calibration using binary logistic regression framework with inverse probability of censoring weights.
### Calibration model uses restricted cubic splines.
###
calc.calib.blr <- function(data.mstate, data.raw, tmat, j, s, t.eval, p.est, curve.type = "rcs", rcs.nk, loess.span, loess.degree,
                               weights = NULL, w.covs, w.landmark.type, w.max, w.stabilised, CI = FALSE, CI.R.boot){

        ### Calculate ipcw weights
        weights <- calc.weights(data.mstate = msebmt,
                                data.raw = ebmt,
                                covs = covs,
                                j = j,
                                landmark.type = "all",
                                s = s,
                                t.eval = t.eval)

        ### Combineweights with ebmt dataset
        ebmt.weights <- cbind(ebmt, weights)

        data.mstate <- msebmt
        data.raw <- ebmt.weights
        tmat
        j<-1
        s<-s
        t.eval <- t.eval
        p.est <- tp.all[[1]] %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
        span <- span
        degree <- degree.blr
        CI <- 95
        R.boot <- 500
        curve.nk <- 3

  ### If a vector of weights has been provided, add it to the dataset
  if (!is.null(weights)){
    ### First check whether it is the correct length (NA's should be present)
    if (length(weights) != nrow(data.raw)){
      print("Weights vector not same length as data.raw")
      stop()
    } else {
      data.raw$ipcw <- weights
    }
  }

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(p.est) != 0)

  ### Add the predicted risks, and the logit transormation of the predicted risks to data.raw
  p.est.logit <- log(p.est[,valid.transitions]/(1-p.est[,valid.transitions]))
  colnames(p.est.logit) <- paste("p.est.logit", valid.transitions, sep = "")
  data.raw <- data.frame(data.raw, p.est[,valid.transitions], p.est.logit)

  ### Identify individuals who are in state j at time s
  ids.state.j <- subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    select(id) %>%
    distinct(id) %>%
    pull(id)

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate %>% subset(id %in% ids.state.j)
  data.raw <- data.raw %>% subset(id %in% ids.state.j)

  ### Extract which state individuals are in at time t.eval
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract.ids.states(data.mstate, tmat, k, t.eval)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- 1:nrow(data.raw)
  m1 <- outer(v1, ids.state.list, FUN = Vectorize('%in%'))
  state.poly <- lapply(split(m1, row(m1)), function(x) (1:max.state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state.poly, length)
  state.poly[idx] <- NA

  ## Add to data.raw
  data.raw <- mutate(data.raw, state.poly = unlist(state.poly),
                     state.poly.fac = factor(state.poly))
  test <- model.matrix(~state.poly.fac - 1, data.raw)
  ### Create binary variables for each possible state that can be transitioned to
  ## Start by creating NA data.frame
  temp.dummy <- data.frame(matrix(NA, ncol = length(valid.transitions), nrow = nrow(data.raw)))

  ## Create dummy variables
  temp.dummy.calc <- model.matrix(~state.poly.fac - 1, data.raw)

  ## Assign to temp.dummy, for rows where data.raw is not NA
  temp.dummy[!is.na(data.raw$state.poly), ] <- temp.dummy.calc

  ## Assign colnames
  colnames(temp.dummy) <- paste("state", valid.transitions, ".bin", sep = "")

  ### Add to dataset
  data.raw <- cbind(data.raw, temp.dummy)
  rm(temp.dummy)

  ### Reduce data.raw to individuals who are uncensored at time t.eval
  data.raw.uncens <- subset(data.raw, !is.na(state.poly))

  ############################################################
  ### Write functions to estimate predicted observed risks ###
  ############################################################

  ### Functions are written to allow bootstrapping and have the following arguments:
  ### Arg1: data.in = dataset which we want to assess calibration in
  ### Arg2: indices = vector of indices to sample dataset
  ### Arg3: state we want to assess calibration for
  ### Arg4: data.in.uncens = dataset of uncensored individuals at time t.eval, it is these we will generate the predicted observed risks for

  ###
  ### Function 1 uses loess smoothers
  ### Weights are calculated within the function to allow for proper bootstrapping
  ###
  calc.obs.loess.weights.func <- function(data.in, indices, state.k, data.in.uncens){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Calculate weights
    weights <- calc.weights(data.mstate, data.boot, w.covs, j, w.landmark.type, s, t.eval, w.max, w.stabilised)

    ## Add to data.boot
    data.boot$ipcw <- weights$ipcw
    if (w.stabilised == TRUE){
      data.boot$ipcw.stab <- weights$ipcw.stab
    }

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- subset(data.boot, !is.na(state.poly))

    ## Define equation
    eq.loess <- formula(paste("state", state.k, ".bin ~ p.est", state.k, sep = ""))

    ## Fit model
    if (w.stabilised == FALSE){
      loess.model <- loess(eq.loess,
                           data = data.boot.uncens,
                           weights = data.boot.uncens[, "ipcw"],
                           span = loess.span,
                           degree = loess.degree)
    } else if (w.stabilised == TRUE){
      loess.model <- loess(eq.loess,
                           data = data.boot.uncens,
                           weights = data.boot.uncens[, "ipcw.stab"],
                           span = loess.span,
                           degree = loess.degree)
    }

    ## Create predicted observed probabilities. Create predictions for the vector of predicted probabilities for uncensored individuals from original dataset.
    loess.pred.obs <- predict(loess.model, newdata = data.in.uncens)

    return(loess.pred.obs)
  }

  ###
  ### Function 2 uses restricted cubic splines
  ### Weights are calculated within the function to allow for proper bootstrapping
  ###
  calc.obs.rcs.weights.func <- function(data.in, indices, state.k, data.in.uncens){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Calculate weights
    weights <- calc.weights(data.mstate, data.boot, w.covs, j, w.landmark.type, s, t.eval, w.max, w.stabilised)

    ## Add to data.boot
    data.boot$ipcw <- weights$ipcw
    if (w.stabilised == TRUE){
      data.boot$ipcw.stab <- weights$ipcw.stab
    }

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- subset(data.boot, !is.na(state.poly))

    ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- rcspline.eval(data.boot.uncens[,paste("p.est.logit", state.k, sep = "")],nk=rcs.nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.boot.uncens
    data.boot.uncens <- data.frame(cbind(data.boot.uncens, rcs.mat))

    ## Define equation
    eq.LHS <- paste("state", state.k, ".bin ~ ", sep = "")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq.rcs <- formula(paste(eq.LHS, eq.RHS, sep = ""))

    ## Fit model
    if (w.stabilised == FALSE){
      rcs.model <- lrm(eq.rcs,
                       data = data.boot.uncens,
                       weights = data.boot.uncens[, "ipcw"])
    } else if (w.stabilised == TRUE){
      rcs.model <- lrm(eq.rcs,
                       data = data.boot.uncens,
                       weights = data.boot.uncens[, "ipcw.stab"])
    }

    ## Create predicted observed probabilities. Create predictions for the vector of predicted probabilities for uncensored individuals from original dataset.

    ## For this, need to calculate the correct splines for the original dataset.
    rcs.mat.data.in.uncens <- rcspline.eval(data.in.uncens[,paste("p.est.logit", state.k, sep = "")],knots = knots,inclx=T)
    colnames(rcs.mat.data.in.uncens) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    #attr(rcs.mat.data.in.uncens,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.in.uncens
    data.in.uncens <- data.frame(cbind(data.in.uncens, rcs.mat.data.in.uncens))

    ## Calculate predicted observed probabilities
    rcs.pred.obs <- predict(rcs.model, newdata = data.in.uncens, type = "fitted")

    return(rcs.pred.obs)
  }


  ###
  ### Function 3 uses loess smoothers
  ### Weights are not calculated within the function (assumes they have been inputted by user)
  ###
  calc.obs.loess.func <- function(data.in, indices, state.k, data.in.uncens){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, c("id", paste("p.est", state.k, sep = ""), paste("state", state.k, ".bin", sep = ""), "state.poly", "ipcw")]

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- subset(data.boot, !is.na(state.poly))

    ## Define equation
    eq.loess <- formula(paste("state", state.k, ".bin ~ p.est", state.k, sep = ""))

    ## Fit model
    loess.model <- loess(eq.loess,
                         data = data.boot.uncens,
                         weights = data.boot.uncens[, "ipcw"],
                         span = loess.span,
                         degree = loess.degree)

    ## Create predictions for the vector of predicted probabilities for people included in original the calibration curve (this is the individuals
    ## who are in the unbootstrapped dataset, and are uncensored at time t.eval. This would be the vector of predicted probabilities for the
    ## calibration curve if no bootstrapping was done)

    ## Create predicted observed risks for these individuals
    loess.pred.obs <- predict(loess.model, newdata = data.in.uncens)

    return(loess.pred.obs)
  }


  ###
  ### Function 4 uses restricted cubic splines
  ### Weights are not calculated within the function (assumes they have been inputted by user)
  ###
  calc.obs.rcs.func <- function(data.in, indices, state.k, data.in.uncens){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- subset(data.boot, !is.na(state.poly))

    ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- rcspline.eval(data.boot.uncens[,paste("p.est.logit", state.k, sep = "")],nk=rcs.nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.boot.uncens
    data.boot.uncens <- data.frame(cbind(data.boot.uncens, rcs.mat))

    ## Define equation
    eq.LHS <- paste("state", state.k, ".bin ~ ", sep = "")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq.rcs <- formula(paste(eq.LHS, eq.RHS, sep = ""))

    ## Fit model
    rcs.model <- lrm(eq.rcs,
                     data = data.boot.uncens,
                     weights = data.boot.uncens[, "ipcw"])


    ## Create predicted observed probabilities. Create predictions for the vector of predicted probabilities form the original model
    ## So that we always plot over the same range.

    ## For this, need to calculate the correct splines for the original dataset.
    rcs.mat.data.in.uncens <- rcspline.eval(data.in.uncens[,paste("p.est.logit", state.k, sep = "")],knots = knots,inclx=T)
    colnames(rcs.mat.data.in.uncens) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    #attr(rcs.mat.data.in.uncens,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.in.uncens
    data.in.uncens <- data.frame(cbind(data.in.uncens, rcs.mat.data.in.uncens))

    ## Calculate predicted observed probabilities
    rcs.pred.obs <- predict(rcs.model, newdata = data.in.uncens, type = "fitted")

    return(rcs.pred.obs)
  }


  ### Assign function to calculate predicted observed risks based on user input
  if (is.null(weights)){
    if (curve.type == "loess"){
      calc.obs.func <- calc.obs.loess.weights.func
    } else if (curve.type == "rcs"){
      calc.obs.func <- calc.obs.rcs.weights.func
    }
  } else {
    if (curve.type == "loess"){
      calc.obs.func <- calc.obs.loess.func
    } else if (curve.type == "rcs"){
      calc.obs.func <- calc.obs.rcs.func
    }
  }


  ### Create object to store output
  output.object <- vector("list", length(valid.transitions))

  ### Loop through and fit models
  for (k in 1:length(valid.transitions)){

    ### Assign state of interest
    state.k <- valid.transitions[k]

    ### Calculate predicted observed probabilities for calibration curve (note the chosen set of indices samples every patient once)
    rcs.pred.obs <- calc.obs.func(data.raw, 1:nrow(data.raw), state.k, data.raw.uncens)

    ### Assign output
    output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                     "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                     "obs" = rcs.pred.obs)

    ### Calculate confidence intervals (only if user-specifies them, and user didn't input their own weights, which would lead to incorrect confidence intervals)
#     if (CI != FALSE & (is.null(weights))){
  if (CI != FALSE ){
      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Run bootstrapping
      boot.obs <- boot(data.raw, calc.obs.func, R = CI.R.boot, state.k = valid.transitions[k], data.in.uncens = data.raw.uncens)$t

      ### Extract confidence bands
      lower <- apply(boot.obs, 2, quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.obs, 2, quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot.obs)) > 0){
        print(paste("WARNING, SOME BOOTSTRAPPED OBSERVED RISKS WERE NA FOR STATE", valid.transitions[k]))
        print(paste("THERE ARE ", sum(apply(boot.obs, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's FOR OBSERVED RISKS"))
        print(paste("THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs, 1, function(x) {sum(is.na(x))}))))
      }

      ### Assign output
      output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                       "pred" = data.raw.uncens[, paste("p.est", valid.transitions[k], sep = "")],
                                       "obs" = rcs.pred.obs,
                                       "obs.lower" = lower,
                                       "obs.upper" = upper)
    }
  }

  ### Create metadata object
  metadata <- list("valid.transitions" = valid.transitions, "CI" = CI, "curve.type" = curve.type)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plot.data" = output.object, "metadata" = metadata)
  return(output.object.comb)
}



###
###
###
calc.weights.nonlinear <- function(data.mstate, data.raw, covs, j = NULL, landmark.type = NULL, s, t.eval, max.weight = 10, stabilised = FALSE){

    data.mstate <- msebmt
    data.raw <- ebmt
    covs <- covs
    j
    landmark.type <- "state"
    s
    t.eval <- t.eval
    max.weight <- 10
    w.nk <- 3
  #
  ### Create a new outcome, which is the time until censored from s
  data.raw$dtcens.modified <- data.raw$dtcens - s

  ### Save a copy of data.raw
  data.raw.save <- data.raw

  ### If landmark.type = "state", calculate weights only in individuals in state j at time s
  ### If landmark.type = "all", calculate weights in all uncensored individuals at time s (note that this excludes individuals
  ### who have reached absorbing states, who have been 'censored' from the survival distribution is censoring)
  if (landmark.type == "state"){
    ### Identify individuals who are uncensored in state j at time s
    ids.uncens <- subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
      select(id) %>%
      distinct(id) %>%
      pull(id)

  } else if (landmark.type == "all"){
    ### Identify individuals who are uncensored time s
    ids.uncens <- subset(data.mstate, Tstart <= s & s < Tstop) %>%
      select(id) %>%
      distinct(id) %>%
      pull(id)

  }

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate %>% subset(id %in% ids.uncens)
  data.raw <- data.raw %>% subset(id %in% ids.uncens)

  ### Create appropriate predictor variables if we have splines/nonlinear components in model for weights
  covs.spline <- c("x1", "x2")
  covs.spline.nk <- c(3,4)

  data.raw$x1 <- rnorm(nrow(data.raw), 0, 1)
  data.raw$x2 <- rnorm(nrow(data.raw), 0, 1)

  ### Create object to store expanded covariates and knot positions
  covar.expand <- vector("list", length(covs.spline))
  knots <- vector("list", length(covs.spline))
  for (i in 1:length(covs.spline)){
    covar.expand[[i]] <- rcspline.eval(data.raw[,covs.spline[i]],nk=covs.spline.nk[i],inclx=T)
    colnames(covar.expand[[i]]) <- paste("rcs.", covs.spline[i], ".", 1:ncol(covar.expand[[i]]), sep = "")
    knots[[i]] <- attr(covar.expand[[i]], "knots")
  }

  head(covar.expand[[1]])
  as.formula(paste("Surv(dtcens.modified, dtcens.s) ~ ",
                   paste(covs, collapse = "+"),

                   sep = ""))

  paste(covs.spline, collapse = "+")
  ###
  ### Create models for censoring in order to calculate the IPCW weights
  ### Seperate models for estimating the weights, and stabilising the weights (intercept only model)
  ###
  if (!is.null(covs)){
    ### A model where we adjust for predictor variables
    cens.model <- coxph(as.formula(paste("Surv(dtcens.modified, dtcens.s) ~ ",
                                         paste(covs, collapse = "+"),
                                         sep = "")),
                        data = data.raw)

    ### Intercept only model (numerator for stabilised weights)
    cens.model.int <- coxph(as.formula(paste("Surv(dtcens.modified, dtcens.s) ~ 1",
                                             sep = "")),
                            data = data.raw)
  } else if (is.null(covs)){
    ### If user has not input any predictors for estimating weights, the model for estimating the weights is the intercept only model (i.e. Kaplan Meier estimator)

    ### Intercept only model (numerator for stabilised weights)
    cens.model.int <- coxph(as.formula(paste("Surv(dtcens.modified, dtcens.s) ~ 1",
                                             sep = "")),
                            data = data.raw)
    ### Assign cens.model to be the same
    cens.model <- cens.model.int


  }

  ### Calculate a data frame containing probability of censored and uncenosred at each time point
  ### The weights will be the probability of being uncensored, at the time of the event for each individual

  ## Extract baseline hazard
  data.weights <- basehaz(cens.model, centered = FALSE)
  ## Add lp to data.raw.save
  data.raw.save$lp <- predict(cens.model, newdata = data.raw.save, type = "lp", reference = "zero")

  ### Create weights for the cohort at time t.eval - s
  ### Note for individuals who reached an absorbing state, we take the probability of them being uncensored at the time of reached the
  ### abosrbing state. For individuals still alive, we take the probability of being uncensored at time t.eval - s.

  ### Get location of individuals who entered absorbing states or were censored prior to evaluation time
  obs.absorbed.prior <- which(data.raw.save$dtcens < t.eval & data.raw.save$dtcens.s == 0)
  obs.censored.prior <- which(data.raw.save$dtcens < t.eval & data.raw.save$dtcens.s == 1)

  ###
  ### Now create unstabilised probability of (un)censoring weights
  ### Note that weights are the probability of being uncensored, so if an individual has low probability of being uncesored,
  ### the inervse of this will be big, weighting them strongly
  ###

  ### First assign all individuals a weight of the probability of being uncesored at time t.eval
  ### This is the linear predictor times the cumulative hazard at time t.eval, and appropriate transformation to get a risk
  data.raw.save$pcw <- as.numeric(exp(-exp(data.raw.save$lp)*data.weights$hazard[max(which(data.weights$time <= t.eval - s))]))

  ## Write a function which will extract the uncensored probability for an individual with linear predictor lp at a given time t
  prob.uncens.func <- function(input){

    ## Assign t and person_id
    t <- input[1]
    lp <- input[2]

    if (t <= 0){
      return(NA)
    } else if (t > 0){
      ## Get hazard at appropriate time
      if (t < min(data.weights$time)){
        bhaz.t <- 0
      } else if (t >= min(data.weights$time)){
        bhaz.t <- data.weights$hazard[max(which(data.weights$time <= t))]
      }

      ## Return risk
      return(exp(-exp(lp)*bhaz.t))
    }
  }

  ### Apply this function to all the times at which individuals have entered an absorbing state prior to censoring
  data.raw.save$pcw[obs.absorbed.prior] <- apply(data.raw.save[obs.absorbed.prior, c("dtcens.modified", "lp")], 1, FUN = prob.uncens.func)

  ### For individuals who were censored prior to t.eval, assign the weight as NA
  data.raw.save$pcw[obs.censored.prior] <- NA

  ### Invert these
  data.raw.save$ipcw <- 1/data.raw.save$pcw

  ###
  ### Stabilise these weights dependent on user-input
  ###
  if (stabilised == TRUE){

    ## Extract baseline hazard
    data.weights.numer <- basehaz(cens.model.int, centered = TRUE)

    ### Assign all individuals a weight of the probability of being uncesored at time t.eval
    data.raw.save$pcw.numer <- as.numeric(exp(-data.weights.numer$hazard[max(which(data.weights.numer$time <= t.eval - s))]))

    #     ### Write a function which will extract the uncensored probability for an individual with linear predictor lp at a given time t
    #     prob.uncens.func.numer <- function(input){
    #       ## Assign t
    #       t <- input
    #
    #       ## Get hazard at appropriate time
    #       bhaz.t <- data.weights.numer$hazard[max(which(data.weights.numer$time <= t))]
    #
    #       ## Return risk
    #       return(exp(-bhaz.t))
    #     }
    #
    #     ### Apply this function to all the times at which individuals have died prior to censoring, and assign to the appropriate individuals
    #     data.raw.save$pcw.numer[obs.absorbed.prior] <- sapply(data.raw.save[obs.absorbed.prior, c("dtcens.modified")], FUN = prob.uncens.func.numer)
    #
    ### Create stabilised weight
    data.raw.save$ipcw.stab <- data.raw.save$pcw.numer*data.raw.save$ipcw
  }

  ### Finally cap these at 10 and create output object

  ### Create output object
  if (stabilised == FALSE){
    data.raw.save$ipcw <- pmin(data.raw.save$ipcw, max.weight)
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw, "pcw" = data.raw.save$pcw)
  } else if (stabilised == TRUE){
    data.raw.save$ipcw <- pmin(data.raw.save$ipcw, max.weight)
    data.raw.save$ipcw.stab <- pmin(data.raw.save$ipcw.stab, max.weight)
    output.weights <- data.frame("id" = data.raw.save$id, "ipcw" = data.raw.save$ipcw, "ipcw.stab" = data.raw.save$ipcw.stab, "pcw" = data.raw.save$pcw)
  }

  return(output.weights)

}




