###
### 1.1) Function to extract all individuals in state j, at time t.eval, from a dataset in mstate format
### Used in other functions which assess calibration using blr/mlr at specific time points (3.4 and 3.5)
###
extract_ids_states <- function(data.mstate, tmat, j, t.eval){

  ### Define maximum state number
  max.state <- max(data.mstate$to)

  ### Identify which states are absorbing states
  absorbing.states <- which(apply(tmat, 1, function(x) {sum(!is.na(x))}) == 0)

  ### For non-absorbing states, to be in state j at time t, you must have an observations from state j, where Tstart <= t.eval < Tstop
  if (!(j %in% absorbing.states)){
    ## Extract ids
    ids.state.j <- base::subset(data.mstate, from == j & Tstart <= t.eval & t.eval < Tstop) %>%
      dplyr::select(id) %>%
      dplyr::distinct(id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$id)
  } else if (j %in% absorbing.states){
    ### For absorbing state, just have to have moved into it
    ids.state.j <- base::subset(data.mstate, to == j & t.eval >= Tstop & status == 1) %>%
      dplyr::select(id) %>%
      dplyr::distinct(id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$id)
  }

  return(ids.state.j)
}




###
### 2.2) Calculate calibration using binary logistic regression framework with inverse probability of censoring weights.
### Calibration model uses restricted cubic splines.
###
#' Create data for calibration curves using a binary logistic regression framework with inverse probability of censoring weights
#'
#' @param data.mstate Validation data in msdata format.
#' @param data.raw Validation data in data.frame (one row per individual).
#' @param t.eval Follow up time at which to calculate weights.
#' @param j Landmark state at which predictions were made
#' @param s Landmark time at which predictions were made
#' @param t.eval Follow up time at which calibration is to be assessed
#' @param tp.pred Vector of predicted transition probabilities at time t.eval
#' @param curve.type Whether calibration curves are estimated using restricted cubic splines ('rcs') or loess smoothers ('loess')
#' @param rcs.nk Number of knots when curves are estimated using restricted cubic splines
#' @param loess.span Span when curves are estimated using loess smoothers
#' @param loess.degree Degree when curves are estimated using loess smoothers
#' @param weights Vector of inverse probability of censoring weights
#' @param w.covs Character vector of variable names to adjust for when calculating inverse probability of censoring weights
#' @param w.landmark.type Whether weights are estimated in all individuals uncensored at time s ('all') or only in individuals uncensored and in state j at time s ('state')
#' @param w.max Maximum bound forinverse probability of censoring weights
#' @param w.stabilised Indicates whether inverse probability of censoring weights should be stabilised or not
#' @param CI Size of confidence intervals as a %
#' @param CI.R.boot Number of bootstrap replicates when estimating the confidence interval for the calibration curve
#'
#' @export
calc_calib_blr <- function(data.mstate, data.raw, j, s, t.eval, tp.pred, curve.type = "rcs", rcs.nk = 3, loess.span, loess.degree,
                           weights = NULL, w.covs, w.landmark.type = "state", w.max = 10, w.stabilised = FALSE, CI = FALSE, CI.R.boot){


        # data.mstate <- msebmtcal
        # data.raw <- ebmtcal
        # j<-1
        # s<-100
        # t.eval <- 1826
        # tp.pred = tps100 %>% dplyr::filter(j == 1) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
        # curve.type = "rcs"
        # rcs.nk = 3
        # weights <- NULL
        # w.covs = c("year", "agecl", "proph", "match")
        # w.landmark.type = "all"
        # w.max = 10
        # w.stabilised = FALSE
        # CI = FALSE

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

  ### Assign colnames to predicted transition probabilities
  colnames(tp.pred) <- paste("tp.pred", 1:ncol(tp.pred), sep = "")

  ### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
  valid.transitions <- which(colSums(tp.pred) != 0)

  ### Add the predicted risks, and the logit transormation of the predicted risks to data.raw
  tp.pred.logit <- log(tp.pred[,valid.transitions]/(1-tp.pred[,valid.transitions]))
  colnames(tp.pred.logit) <- paste("tp.pred.logit", valid.transitions, sep = "")
  data.raw <- data.frame(data.raw, tp.pred[,valid.transitions], tp.pred.logit)

  ### Identify individuals who are in state j at time s
  ids.state.j <- base::subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
    dplyr::select(id) %>%
    dplyr::distinct(id) %>%
    dplyr::pull(id)

  ### Subset data.mstate and data.raw to these individuals
  data.mstate <- data.mstate %>% base::subset(id %in% ids.state.j)
  data.raw <- data.raw %>% base::subset(id %in% ids.state.j)

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

  ### Reduce data.raw to individuals who are uncensored at time t.eval
  data.raw.uncens <- base::subset(data.raw, !is.na(state.poly))

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
  calc_obs_loess_weights_func <- function(data.in, indices, state.k, data.in.uncens){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Calculate weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t.eval = t.eval,
                            s = s,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised)

    ## Add to data.boot
    data.boot$ipcw <- weights$ipcw
    if (w.stabilised == TRUE){
      data.boot$ipcw.stab <- weights$ipcw.stab
    }

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- base::subset(data.boot, !is.na(state.poly))

    ## Define equation
    eq.loess <- stats::formula(paste("state", state.k, ".bin ~ tp.pred", state.k, sep = ""))

    ## Fit model
    if (w.stabilised == FALSE){
      loess.model <- stats::loess(eq.loess,
                                  data = data.boot.uncens,
                                  weights = data.boot.uncens[, "ipcw"],
                                  span = loess.span,
                                  degree = loess.degree)
    } else if (w.stabilised == TRUE){
      loess.model <- stats::loess(eq.loess,
                                  data = data.boot.uncens,
                                  weights = data.boot.uncens[, "ipcw.stab"],
                                  span = loess.span,
                                  degree = loess.degree)
    }

    ## Create predicted observed probabilities. Create predictions for the vector of predicted probabilities for uncensored individuals from original dataset.
    loess.pred.obs <- stats::predict(loess.model, newdata = data.in.uncens)

    return(loess.pred.obs)
  }

  ###
  ### Function 2 uses restricted cubic splines
  ### Weights are calculated within the function to allow for proper bootstrapping
  ###
  calc_obs_rcs_weights_func <- function(data.in, indices, state.k, data.in.uncens){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Calculate weights
    weights <- calc_weights(data.mstate = data.mstate,
                            data.raw = data.boot,
                            covs = w.covs,
                            t.eval = t.eval,
                            s = s,
                            landmark.type = w.landmark.type,
                            j = j,
                            max.weight = w.max,
                            stabilised = w.stabilised)

    ## Add to data.boot
    data.boot$ipcw <- weights$ipcw
    if (w.stabilised == TRUE){
      data.boot$ipcw.stab <- weights$ipcw.stab
    }

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- base::subset(data.boot, !is.na(state.poly))

    ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- Hmisc::rcspline.eval(data.boot.uncens[,paste("tp.pred.logit", state.k, sep = "")],nk=rcs.nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.boot.uncens
    data.boot.uncens <- data.frame(cbind(data.boot.uncens, rcs.mat))

    ## Define equation
    eq.LHS <- paste("state", state.k, ".bin ~ ", sep = "")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq.rcs <- stats::formula(paste(eq.LHS, eq.RHS, sep = ""))

    ## Fit model
    ## NB: Warnings are suppressed beceause rms gives the following warning:
      ## In rms::lrm(eq.rcs, data = data.boot.uncens, weights = data.boot.uncens[,:
      ## currently weights are ignored in model validation and bootstrapping lrm fits
    ## We are not using the model validation or bootstrapping aspects of rms::lrm (we are applying bootstrapping ourselves),
    ## meaning the warning is not neccesary.
    suppressWarnings(
      if (w.stabilised == FALSE){
        rcs.model <- rms::lrm(eq.rcs,
                              data = data.boot.uncens,
                              weights = data.boot.uncens[, "ipcw"])
      } else if (w.stabilised == TRUE){
        rcs.model <- rms::lrm(eq.rcs,
                              data = data.boot.uncens,
                              weights = data.boot.uncens[, "ipcw.stab"])
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
    rcs.pred.obs <- stats::predict(rcs.model, newdata = data.in.uncens, type = "fitted")

    return(rcs.pred.obs)
  }


  ###
  ### Function 3 uses loess smoothers
  ### Weights are not calculated within the function (assumes they have been inputted by user)
  ###
  calc_obs_loess_func <- function(data.in, indices, state.k, data.in.uncens){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- base::subset(data.boot, !is.na(state.poly))

    ## Define equation
    eq.loess <- stats::formula(paste("state", state.k, ".bin ~ tp.pred", state.k, sep = ""))

    ## Fit model
    loess.model <- stats::loess(eq.loess,
                                data = data.boot.uncens,
                                weights = data.boot.uncens[, "ipcw"],
                                span = loess.span,
                                degree = loess.degree)

    ## Create predictions for the vector of predicted probabilities for people included in original the calibration curve (this is the individuals
    ## who are in the unbootstrapped dataset, and are uncensored at time t.eval. This would be the vector of predicted probabilities for the
    ## calibration curve if no bootstrapping was done)

    ## Create predicted observed risks for these individuals
    loess.pred.obs <- stats::predict(loess.model, newdata = data.in.uncens)

    return(loess.pred.obs)
  }


  ###
  ### Function 4 uses restricted cubic splines
  ### Weights are not calculated within the function (assumes they have been inputted by user)
  ###
  calc_obs_rcs_func <- function(data.in, indices, state.k, data.in.uncens){

    ## Create bootstrapped dataset
    data.boot <- data.in[indices, ]

    ## Reduce data.boot to individuals who are uncensored at time t.eval
    data.boot.uncens <- base::subset(data.boot, !is.na(state.poly))

    ## Create restricted cubic splines for the cloglog of the linear predictor for the state of interst
    rcs.mat <- Hmisc::rcspline.eval(data.boot.uncens[,paste("tp.pred.logit", state.k, sep = "")],nk=rcs.nk,inclx=T)
    colnames(rcs.mat) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    knots <- attr(rcs.mat,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.boot.uncens
    data.boot.uncens <- data.frame(cbind(data.boot.uncens, rcs.mat))

    ## Define equation
    eq.LHS <- paste("state", state.k, ".bin ~ ", sep = "")
    eq.RHS <- paste("rcs.x", 1:ncol(rcs.mat), sep = "", collapse = "+")
    eq.rcs <- stats::formula(paste(eq.LHS, eq.RHS, sep = ""))

    ## Fit model
    rcs.model <- rms::lrm(eq.rcs,
                          data = data.boot.uncens,
                          weights = data.boot.uncens[, "ipcw"])


    ## Create predicted observed probabilities. Create predictions for the vector of predicted probabilities form the original model
    ## So that we always plot over the same range.

    ## For this, need to calculate the correct splines for the original dataset.
    rcs.mat.data.in.uncens <- Hmisc::rcspline.eval(data.in.uncens[,paste("tp.pred.logit", state.k, sep = "")],knots = knots,inclx=T)
    colnames(rcs.mat.data.in.uncens) <- paste("rcs.x", 1:ncol(rcs.mat), sep = "")
    #attr(rcs.mat.data.in.uncens,"knots")

    ## Add the cubic splines for logit of the predicted probability to data.in.uncens
    data.in.uncens <- data.frame(cbind(data.in.uncens, rcs.mat.data.in.uncens))

    ## Calculate predicted observed probabilities
    rcs.pred.obs <- stats::predict(rcs.model, newdata = data.in.uncens, type = "fitted")

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
  output.object <- vector("list", length(valid.transitions))

  ### Loop through and fit models
  for (k in 1:length(valid.transitions)){

    ### Assign state of interest
    state.k <- valid.transitions[k]

    ### Calculate predicted observed probabilities for calibration curve (note the chosen set of indices samples every patient once)
    rcs.pred.obs <- calc_obs_func(data.in = data.raw, indices = 1:nrow(data.raw), state.k = state.k, data.in.uncens = data.raw.uncens)

    ### Assign output
    output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                     "pred" = data.raw.uncens[, paste("tp.pred", valid.transitions[k], sep = "")],
                                     "obs" = rcs.pred.obs)

    ### Calculate confidence intervals (only if user-specifies them, and user didn't input their own weights, which would lead to incorrect confidence intervals)
    #     if (CI != FALSE & (is.null(weights))){
    if (CI != FALSE ){
      ### Define alpha for CI's
      alpha <- (1-CI/100)/2

      ### Run bootstrapping
      boot.obs <- boot::boot(data.raw, calc_obs_func, R = CI.R.boot, state.k = valid.transitions[k], data.in.uncens = data.raw.uncens)$t

      ### Extract confidence bands
      lower <- apply(boot.obs, 2, stats::quantile, probs = alpha, na.rm = TRUE)
      upper <- apply(boot.obs, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

      ### Produce a warning if any NA values
      if(sum(is.na(boot.obs)) > 0){
        print(paste("WARNING, SOME BOOTSTRAPPED OBSERVED RISKS WERE NA FOR STATE", valid.transitions[k]))
        print(paste("THERE ARE ", sum(apply(boot.obs, 1, function(x) {sum(is.na(x)) > 0})), " ITERATIONS WITH NA's FOR OBSERVED RISKS"))
        print(paste("THE MEAN NUMBER OF NA's IN EACH ITERATION IS", mean(apply(boot.obs, 1, function(x) {sum(is.na(x))}))))
      }

      ### Assign output
      output.object[[k]] <- data.frame("id" = data.raw.uncens[, "id"],
                                       "pred" = data.raw.uncens[, paste("tp.pred", valid.transitions[k], sep = "")],
                                       "obs" = rcs.pred.obs,
                                       "obs.lower" = lower,
                                       "obs.upper" = upper)
    }
  }

  ### Create metadata object
  metadata <- list("valid.transitions" = valid.transitions, "CI" = CI, "curve.type" = curve.type)

  ### Crate a combined output object with metadata, as well as plot data
  output.object.comb <- list("plotdata" = output.object, "metadata" = metadata)

  ### Assign calibmsm class
  attr(output.object.comb, "class") <- "calib_blr"

  return(output.object.comb)
}



###
### 2.3) Calculate calibration using multinomial logistic regression framework with inverse probability of censoring weights.
### Calibration model uses nominal recalibration framework of van Hoorde.
###
calc_calib_mlr <- function(data.mstate, data.raw, tmat, j, s, t.eval, p.est, ps.int = 4, degree = 3,
                           weights = NULL, w.covs, w.landmark.type = "state", w.max = 10, w.stabilised = FALSE){

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
  #       j<-j
  #       s<-s
  #       t.eval <- t.eval
  #       p.est <- tp.all[[1]] %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
  #       ps.int <- 4
  #       degree <- 3

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
  ### predicted observed probabilities to the correct individuals

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
