#' Create data for calibration curves using a binary logistic regression
#' framework with inverse probability of censoring weights
#'
#' @description
#' Creates the underlying data for the calibration plots. Observed event
#' probabilities at time `t.eval` are estimated for inputted predicted
#' transition probabilities `tp.pred` out of state `j` at time `s`.
#' `calc_calib_blr` estimates calibration curves using a binary logistic
#' framework in combination with landmarking and inverse probability of
#' censoring weights. A choice between restricted cubic splines and loess
#' smoothers for estimating the calibration curve can be made using `curve.type`.
#'
#' Two datasets for the same cohort of inidividuals must be provided. A `msdata`
#' format dataset generated using the `mstate` package. A `data.frame` with one
#' row per individual, relevant variables for estimating the weights, and a time
#' until censoring varaible (`dtcens`) and indicator (`dtcens.s`). Weights are
#' estimated using a cox-proportional hazard model and assuming linear
#' functional form of the variables defined in `w.covs`. We urge users to
#' specify their own model for estimating the weights. Confidence intervals for
#' the calibration curves can be estimated using bootstrapping. This procedure
#' uses the internal method for estimating weights, we therefore encourage
#' users to specify their own bootstrapping procedure, which incorporates their
#' own model for estimating the weights. Details on how to do this are provided
#' in the vignette.
#'
#' @param data.mstate Validation data in `msdata` format.
#' @param data.raw Validation data in data.frame (one row per individual).
#' @param t.eval Follow up time at which to calculate weights.
#' @param j Landmark state at which predictions were made
#' @param s Landmark time at which predictions were made
#' @param t.eval Follow up time at which calibration is to be assessed
#' @param tp.pred Matrix of predicted transition probabilities at time t.eval, if in state j at time s. There must be a seperate column for the predicted transition probabilities into every state, even if these predicted transition probabilities are 0.
#' @param curve.type Whether calibration curves are estimated using restricted cubic splines ('rcs') or loess smoothers ('loess')
#' @param rcs.nk Number of knots when curves are estimated using restricted cubic splines
#' @param loess.span Span when curves are estimated using loess smoothers
#' @param loess.degree Degree when curves are estimated using loess smoothers
#' @param weights Vector of inverse probability of censoring weights
#' @param w.covs Character vector of variable names to adjust for when calculating inverse probability of censoring weights
#' @param w.landmark.type Whether weights are estimated in all individuals uncensored at time s ('all') or only in individuals uncensored and in state j at time s ('state')
#' @param w.max Maximum bound for inverse probability of censoring weights
#' @param w.stabilised Indicates whether inverse probability of censoring weights should be stabilised or not
#' @param w.max.follow Maximum follow up for model calculating inverse probability of censoring weights. Reducing this to `t.eval` + 1 may aid in the proportional hazards assumption being met in this model.
#' @param CI Size of confidence intervals as a %
#' @param CI.R.boot Number of bootstrap replicates when estimating the confidence interval for the calibration curve
#' @param data.pred.plot Data frame or matrix of predicted risks for each possible transition over which to plot the calibration curves. Must have one column for every possible transition.
#' @param transitions.out Transitions for which to calculate calibration curves. Will do all possible transitions if left as NULL.
#'
#' @export
calc_calib_blr <- function(data.mstate, data.raw, j, s, t.eval, tp.pred, curve.type = "rcs", rcs.nk = 3, loess.span, loess.degree,
                           weights = NULL, w.covs, w.landmark.type = "state", w.max = 10, w.stabilised = FALSE, w.max.follow = NULL, CI = FALSE, CI.R.boot,
                           data.pred.plot = NULL, transitions.out = NULL){


# data.mstate <- msebmtcal
# data.raw <- data.boot
# j <- 1
# j.in <- 1
# s<-0
# t.eval <- 1826
# tp.pred = tps0 %>% dplyr::filter(j == j.in) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
# curve.type = "rcs"
# rcs.nk = 3
# weights <- weights.manual
# w.covs = c("year", "agecl", "proph", "match")
# w.landmark.type = "all"
# w.max = 10
# w.stabilised = FALSE
# CI = FALSE
# str(data.pred.plot)
# data.pred.plot <- tps0 %>% dplyr::filter(j == j.in) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
# data.pred.plot$id <- 1:nrow(data.pred.plot)

  ###
  ### FIRST DO WARNINGS
  ###

  ### Warning if confidence interval requested internally, and a set of predicted risks to plot over has also been inputted (data.pred.plot ignored)

  ### Warning if weights inputted manually, and confidence interval requested internally

  ### Check if people have specified internal bootstrap and only for specific weights...

  ### Check if transitions.out is only specified for non-zero columns

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

  ### Extract transition matrix from msdata object
  tmat <- attributes(data.mstate)$trans

  ### Assign the maximum state an individual may enter
  max.state <- max(data.mstate$to)

  ### Assign colnames to predicted transition probabilities (and in data.pred.plot)
  colnames(tp.pred) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  if (!is.null(data.pred.plot)){
    colnames(data.pred.plot) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")
  }

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(tp.pred) != 0)

  ### Add the predicted risks, and the logit transformation of the predicted risks to data.raw and data.pred.plot
  tp.pred.logit <- log(tp.pred[,valid.transitions]/(1-tp.pred[,valid.transitions]))
  colnames(tp.pred.logit) <- paste("tp.pred.logit", valid.transitions, sep = "")
  data.raw <- data.frame(data.raw, tp.pred[,valid.transitions], tp.pred.logit)

  if (!is.null(data.pred.plot)){
    data.pred.plot.logit <- log(data.pred.plot[,valid.transitions]/(1-data.pred.plot[,valid.transitions]))
    colnames(data.pred.plot.logit) <- paste("tp.pred.logit", valid.transitions, sep = "")
    data.pred.plot <- data.frame(data.pred.plot[,valid.transitions], data.pred.plot.logit)

  }

  ### Extract which state individuals are in at time t.eval
  ids.state.list <- vector("list", max.state)
  for (k in valid.transitions){
    ids.state.list[[k]] <- extract_ids_states(data.mstate, tmat, k, t.eval)
  }

  ### Create a variable to say which state an individual was in at the time of interest
  ## Create list containing the relevant data
  v1 <- data.raw$id
  m1 <- outer(v1, ids.state.list, FUN = Vectorize('%in%'))
  state.poly <- lapply(split(m1, row(m1)), function(x) (1:max.state)[x])

  ## Change integer(0) values to NA's
  idx <- !sapply(state.poly, length)
  state.poly[idx] <- NA

  ## Add to data.raw
  data.raw <- dplyr::mutate(data.raw, state.poly = unlist(state.poly),
                     state.poly.fac = factor(state.poly))

  ### Create binary variables for each possible state that can be transitioned to
  ## Start by creating NA data.frame
  temp.dummy <- data.frame(matrix(NA, ncol = length(valid.transitions), nrow = nrow(data.raw)))

  ## Create dummy variables
  temp.dummy.calc <- stats::model.matrix(~state.poly.fac - 1, data.raw)

  ## Assign to temp.dummy, for rows where data.raw is not NA
  temp.dummy[!is.na(data.raw$state.poly), ] <- temp.dummy.calc

  ## Assign colnames
  colnames(temp.dummy) <- paste("state", valid.transitions, ".bin", sep = "")

  ### Add to dataset
  data.raw <- cbind(data.raw, temp.dummy)
  rm(temp.dummy)

  ### Identify individuals who are in state j at time s (will be used for landmarking)
  ids.state.js <- base::subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    dplyr::select(id) %>%
    dplyr::distinct(id) %>%
    dplyr::pull(id)

  ### Reduce data.raw to landmarked dataset of individuals who are uncensored at time t.eval,
  ### this is the set of predicted risks over which we plot calibration curves
  data.raw.lmk.js.uncens <- data.raw %>% base::subset(id %in% ids.state.js) %>% base::subset(!is.na(state.poly))


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
  calc_obs_loess_weights_func <- function(data.in, indices, state.k, data.in.uncens, data.pred.plot){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Create landmarked dataset
    data.boot.lmk.js <-  data.boot %>% base::subset(id %in% ids.state.js)

    ## Calculate weights
    ## Note this is done in the entire dataset data.boot, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
    ## calculating the weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t.eval = t.eval,
                            s = s,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised,
                            max.follow = w.max.follow)

    ## Add to data.boot
    data.boot.lmk.js <- dplyr::left_join(data.boot.lmk.js, dplyr::distinct(weights), by = dplyr::join_by(id))

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.lmk.js.uncens <- base::subset(data.boot.lmk.js, !is.na(state.poly))

    ## Define equation
    eq.loess <- stats::formula(paste("state", state.k, ".bin ~ tp.pred", state.k, sep = ""))

    ## Fit model
    if (w.stabilised == FALSE){
      loess.model <- stats::loess(eq.loess,
                                  data = data.boot.lmk.js.uncens,
                                  weights = data.boot.lmk.js.uncens[, "ipcw"],
                                  span = loess.span,
                                  degree = loess.degree)
    } else if (w.stabilised == TRUE){
      loess.model <- stats::loess(eq.loess,
                                  data = data.boot.lmk.js.uncens,
                                  weights = data.boot.lmk.js.uncens[, "ipcw.stab"],
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
  calc_obs_rcs_weights_func <- function(data.in, indices, state.k, data.in.uncens, data.pred.plot){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Create landmarked dataset
    data.boot.lmk.js <-  data.boot %>% base::subset(id %in% ids.state.js)

    ## Calculate weights
    ## Note this is done in the entire dataset data.boot, which has its own functionality (w.landmark.type) to landmark on j and s, or just s, before
    ## calculating the weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t.eval = t.eval,
                            s = s,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised,
                            max.follow = w.max.follow)

    ## Add to data.boot
    data.boot.lmk.js <- dplyr::left_join(data.boot.lmk.js, dplyr::distinct(weights), by = dplyr::join_by(id))

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.lmk.js.uncens <- base::subset(data.boot.lmk.js, !is.na(state.poly))

    ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- Hmisc::rcspline.eval(data.boot.lmk.js.uncens[,paste("tp.pred.logit", state.k, sep = "")],nk=rcs.nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ## Add the cubic splines for logit of the predicted probability to  data.boot.lmk.js.uncens
    data.boot.lmk.js.uncens <- data.frame(cbind(data.boot.lmk.js.uncens, rcs.mat))

    ## Define equation
    eq.LHS <- paste("state", state.k, ".bin ~ ", sep = "")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq.rcs <- stats::formula(paste(eq.LHS, eq.RHS, sep = ""))

    ## Fit model
    ## NB: Warnings are suppressed beceause rms gives the following warning:
      ## In rms::lrm(eq.rcs, data = data.boot.lmk.js.uncens, weights = data.boot.lmk.js.uncens[,:
      ## currently weights are ignored in model validation and bootstrapping lrm fits
    ## We are not using the model validation or bootstrapping aspects of rms::lrm (we are applying bootstrapping ourselves),
    ## meaning the warning is not neccesary.
    suppressWarnings(
      if (w.stabilised == FALSE){
        rcs.model <- rms::lrm(eq.rcs,
                              data = data.boot.lmk.js.uncens,
                              weights = data.boot.lmk.js.uncens[, "ipcw"])
      } else if (w.stabilised == TRUE){
        rcs.model <- rms::lrm(eq.rcs,
                              data = data.boot.lmk.js.uncens,
                              weights = data.boot.lmk.js.uncens[, "ipcw.stab"])
      }
    )

    ## Create predicted observed probabilities. Create predictions for the vector of predicted probabilities for uncensored individuals from original dataset.

    ## For this, need to calculate the correct splines for the original dataset.
    rcs.mat.data.in.uncens <- Hmisc::rcspline.eval(data.in.uncens[,paste("tp.pred.logit", state.k, sep = "")],knots = knots,inclx=T)
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
  calc_obs_loess_func <- function(data.in, indices, state.k, data.in.uncens, data.pred.plot){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Create landmarked dataset
    data.boot.lmk.js <-  data.boot %>% base::subset(id %in% ids.state.js)

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.lmk.js.uncens <- base::subset(data.boot.lmk.js, !is.na(state.poly))

    ## Define equation
    eq.loess <- stats::formula(paste("state", state.k, ".bin ~ tp.pred", state.k, sep = ""))

    ## Fit model
    loess.model <- stats::loess(eq.loess,
                                data = data.boot.lmk.js.uncens,
                                weights = data.boot.lmk.js.uncens[, "ipcw"],
                                span = loess.span,
                                degree = loess.degree)

    ## Create predictions for the vector of predicted probabilities for people included in original the calibration curve (this is the individuals
    ## who are in the unbootstrapped dataset, and are uncensored at time t.eval. This would be the vector of predicted probabilities for the
    ## calibration curve if no bootstrapping was done)

    ## Create predicted observed risks
    if (is.null(data.pred.plot)){
      loess.pred.obs <- predict(loess.model, newdata = data.in.uncens)
    } else {
      loess.pred.obs <- predict(loess.model, newdata = data.pred.plot)
    }


    return(loess.pred.obs)
  }


  ###
  ### Function 4 uses restricted cubic splines
  ### Weights are not calculated within the function (assumes they have been inputted by user)
  ###
  calc_obs_rcs_func <- function(data.in, indices, state.k, data.in.uncens, data.pred.plot){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Create landmarked dataset
    data.boot.lmk.js <-  data.boot %>% base::subset(id %in% ids.state.js)

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.lmk.js.uncens <- base::subset(data.boot.lmk.js, !is.na(state.poly))

    ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- Hmisc::rcspline.eval(data.boot.lmk.js.uncens[,paste("tp.pred.logit", state.k, sep = "")],nk=rcs.nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.boot.lmk.js.uncens
    data.boot.lmk.js.uncens <- data.frame(cbind(data.boot.lmk.js.uncens, rcs.mat))

    ## Define equation
    eq.LHS <- paste("state", state.k, ".bin ~ ", sep = "")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq.rcs <- stats::formula(paste(eq.LHS, eq.RHS, sep = ""))

    ## Fit model
    suppressWarnings(
    rcs.model <- rms::lrm(eq.rcs,
                          data = data.boot.lmk.js.uncens,
                          weights = data.boot.lmk.js.uncens[, "ipcw"])
    )


    ## Create predicted observed probabilities. Create predictions for the vector of predicted probabilities form the original model
    ## So that we always plot over the same range.

    ## For this, need to calculate the correct splines for the original dataset.
    rcs.mat.data.in.uncens <- Hmisc::rcspline.eval(data.in.uncens[,paste("tp.pred.logit", state.k, sep = "")],knots = knots,inclx=T)
    colnames(rcs.mat.data.in.uncens) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    #attr(rcs.mat.data.in.uncens,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.in.uncens
    data.in.uncens <- data.frame(cbind(data.in.uncens, rcs.mat.data.in.uncens))

    ## Calculate predicted observed probabilities
    if (is.null(data.pred.plot)){
      rcs.pred.obs <- predict(rcs.model, newdata = data.in.uncens, type = "fitted")
    } else {
      ## Create spline terms for data.pred.plot
      rcs.mat.data.pred.plot <- Hmisc::rcspline.eval(data.pred.plot[,paste("tp.pred.logit", state.k, sep = "")],knots = knots,inclx=T)
      colnames(rcs.mat.data.pred.plot) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")

      ## Add to dataset
      data.pred.plot <- data.frame(cbind(data.pred.plot, rcs.mat.data.pred.plot))

      ## Calculate predicted observed probabilities
      rcs.pred.obs <- predict(rcs.model, newdata = data.pred.plot, type = "fitted")
    }

    return(rcs.pred.obs)
  }

  ###################################
  ### FUNCTIONS HAVE BEEN DEFINED ###
  ###################################

  ### Assign function to calculate predicted observed risks based on user input
  if (is.null(weights)){
    if (curve.type == "loess"){
      calc_obs_func <- calc_obs_loess_weights_func
    } else if (curve.type == "rcs"){
      calc_obs_func <- calc_obs_rcs_weights_func
    }
  } else {
    if (curve.type == "loess"){
      calc_obs_func <- calc_obs_loess_func
    } else if (curve.type == "rcs"){
      calc_obs_func <- calc_obs_rcs_func
    }
  }


  ### Create object to store output
  if (is.null(transitions.out)){
    transitions.out <- valid.transitions
  }

  output.object <- vector("list", length(transitions.out))
  names(output.object) <- paste("state", transitions.out, sep = "")

  ### Loop through and fit models
  for (k in 1:length(transitions.out)){

    ### Assign state of interest
    state.k <- transitions.out[k]

    ### Calculate predicted observed probabilities for calibration curve (note the chosen set of indices samples every patient once)
    rcs.pred.obs <- calc_obs_func(data.in = data.raw, indices = 1:nrow(data.raw), state.k = state.k,
                                  data.in.uncens = data.raw.lmk.js.uncens, data.pred.plot = data.pred.plot)

    ### Assign output
    if (is.null(data.pred.plot)){
      output.object[[k]] <- data.frame("id" = data.raw.lmk.js.uncens[, "id"],
                                       "pred" = data.raw.lmk.js.uncens[, paste("tp.pred", transitions.out[k], sep = "")],
                                       "obs" = rcs.pred.obs)
    } else {
      output.object[[k]] <- data.frame(
                                       "pred" = data.pred.plot[, paste("tp.pred", transitions.out[k], sep = "")],
                                       "obs" = rcs.pred.obs)
    }



    ### Calculate confidence intervals (only if user-specifies them, and user didn't input their own weights, which would lead to incorrect confidence intervals)
    #     if (CI != FALSE & (is.null(weights))){
    if (CI != FALSE ){
      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Run bootstrapping
      boot.obs <- boot::boot(data.raw, calc_obs_func, R = CI.R.boot, state.k = transitions.out[k], data.in.uncens = data.raw.lmk.js.uncens)$t

      ### Extract confidence bands
      lower <- apply(boot.obs, 2, stats::quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.obs, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot.obs)) > 0){
        print(paste("WARNING, SOME BOOTSTRAPPED OBSERVED RISKS WERE NA FOR STATE", transitions.out[k]))
        print(paste("THERE ARE ", sum(apply(boot.obs, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's FOR OBSERVED RISKS"))
        print(paste("THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs, 1, function(x) {sum(is.na(x))}))))
      }

      ### Assign output
      output.object[[k]] <- data.frame("id" = data.raw.lmk.js.uncens[, "id"],
                                       "pred" = data.raw.lmk.js.uncens[, paste("tp.pred", transitions.out[k], sep = "")],
                                       "obs" = rcs.pred.obs,
                                       "obs.lower" = lower,
                                       "obs.upper" = upper)
    }
  }

  ### Create metadata object
  metadata <- list("valid.transitions" = valid.transitions, "CI" = CI, "curve.type" = curve.type)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plotdata" = output.object, "metadata" = metadata)

  ### Assign calib_blr class
  attr(output.object.comb, "class") <- "calib_blr"

  return(output.object.comb)
}
