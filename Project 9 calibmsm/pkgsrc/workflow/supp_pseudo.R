###
### Supplementary material 1 pseudo-value analysis
### The aim of this analysis is to repeat the analysis but using the pseudo-value approach
### Focus on j  1, s = 0
###

### Clear workspace
rm(list=ls())

### Load calibmsm and required data
library("calibmsm")
library("dplyr")
library("mstate")
data("ebmtcal")
data("msebmtcal")
data("tps0")
data("tps100")

########################################
### Define functions for the program ###
########################################

###
### Define a function to calculate Aalen-Johansen estimator of transition probabilities for a group of patients
calc.calib.aj.ce <- function(data.mstate, tmat, t.eval){

  ### Assign max state number
  max.state <- max(data.mstate$to)

  ### Fit csh's with no predictors
  csh.aj <- survival::coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data.mstate)

  ### Calculate cumulative incidence functions
  msfit.aj <- mstate::msfit(csh.aj, trans = tmat)

  ### Calculate Aalen-Johansen estimator
  pt.aj <- mstate::probtrans(msfit.aj, predt = 0)

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
### Define a function to calculate pseudo-value for an individual, using the Aalen-Johansen estimator
func.calc.pv.aj.ce <- function(person_id.eval, data.mstate, obs.aj, tmat, n.cohort, t.eval){

  ### Calculate AJ estimate without patient in dataset
  est.drop.pat <- calc.calib.aj.ce(subset(data.mstate, id != person_id.eval),
                                   tmat = tmat, t.eval = t.eval)

  ### Retain just the estimate (not the standard error)
  est.drop.pat <- est.drop.pat[["obs.aj"]]

  ### Calculate the pseudo-value
  pv.pat <- n.cohort*obs.aj - (n.cohort-1)*est.drop.pat
  pv.pat <- as.numeric(pv.pat[1, ])

  return(pv.pat)
}


########################################################
### Do the prep before calculating the pseudo-values ###
########################################################

### Specifically want to create datasets of individuals grouped by predicted risk into each state k, and calculate
### the Aalen-Johansen estimator within each of these groups

### Assign t.eval
t.eval <- 1826

### Assign number of pctls
n.pctls <- 10

### Create objects to store data sorted by predicted risk for each state k, and the Aalen-Johansen estimator for each subgroup
data.sort.pctls <- vector("list", 6)
obs.aj.pctls <- vector("list", 6)

### Assign data.raw by merging tps0 for j = 1 with embtcal
data.raw <- dplyr::inner_join(ebmtcal, tps0 %>% subset(j == 1), by = join_by("id"))

### Define data.mstate
data.mstate <- msebmtcal

### For each states k, order patients by predicted risk and split into percentiles
### Loop through the states of interest
for (state in 1:6){

  print(paste("state = ", state, Sys.time()))

  ###
  ### Create n.pctls datasets, grouped by predicted risk
  ###

  ### Arrange data by risk of category of interest
  data.raw.sort.temp <- dplyr::arrange(data.raw, !!sym(paste("pstate", state, sep = "")))

  ### Create a list of length n.pctls
  data.sort.pctls[[state]] <- vector("list", n.pctls)
  group.size <- nrow(data.raw.sort.temp)/n.pctls
  for (i in 1:n.pctls){
    data.sort.pctls[[state]][[i]] <- dplyr::slice(data.raw.sort.temp, ((round((i-1)*group.size))+1):(round(i*group.size)))
  }

  ###
  ### Calculate the Aalen-Johansen estimator within each group
  ###

  ### Create a list to store them
  obs.aj.pctls[[state]] <- vector("list", n.pctls)

  ### Calculate obs.aj for each group
  for (i in 1:n.pctls){
    print(paste("Calc obs.aj = ", state, "pctl = ", i, Sys.time()))
    obs.aj.pctls[[state]][[i]] <- calc.calib.aj.ce(data.mstate = base::subset(data.mstate, id %in% data.sort.pctls[[state]][[i]]$id),
                                                   tmat = attributes(data.mstate)$trans, t.eval = t.eval)$obs.aj
  }

}


##########################################
### Now to calculate the pseudo-values ###
##########################################

### Create list to store
pv.out.pctls <- vector("list", 6)

### Run through and calculate pseudo-values
for (state in 1:6){
  ### Create output list
  pv.out.pctls[[state]] <- vector("list", n.pctls)
  for(pctl in 1:n.pctls){
    print(paste("state = ", state, "pctl = ", pctl, Sys.time()))
    ### Get vector of patient ids from the state and pctl combo
    ids.vec <- data.sort.pctls[[state]][[pctl]]$id
    ### Calculate pseudo-values and store
    pv.temp <- lapply(ids.vec, func.calc.pv.aj.ce,
                      data.mstate = subset(data.mstate, id %in% data.sort.pctls[[state]][[pctl]]$id),
                      obs.aj = obs.aj.pctls[[state]][[pctl]],
                      tmat = attributes(data.mstate)$trans,
                      n.cohort = nrow(data.sort.pctls[[state]][[pctl]]),
                      t.eval = t.eval)

    ### Only interested in pseudo-values for the state of interest (which we sorted by)
    pv.temp <- do.call("rbind", pv.temp)[, state]

    ### Combine with person_ids, so we have a record of who these pseudo-values are for
    pv.out.pctls[[state]][[pctl]] <- data.frame(cbind(ids.vec, pv.temp))
    colnames(pv.out.pctls[[state]][[pctl]]) <- c("id", paste("pv.", state, sep = ""))

  }

  ### Reduce the output to one dataset for each state and sort by id
  pv.out.pctls[[state]] <- do.call("rbind", pv.out.pctls[[state]])
  pv.out.pctls[[state]] <- dplyr::arrange(pv.out.pctls[[state]], id)

}

### combine the pseudo-values into one dataset
pv.comb <- inner_join(pv.out.pctls[[1]], pv.out.pctls[[2]], by = join_by("id")) %>%
  inner_join(pv.out.pctls[[3]], by = join_by("id")) %>%
  inner_join(pv.out.pctls[[4]], by = join_by("id")) %>%
  inner_join(pv.out.pctls[[5]], by = join_by("id")) %>%
  inner_join(pv.out.pctls[[6]], by = join_by("id"))

save.image("workflow/supp_pseudo.RData")

################################################
### Now calculate moderate calibration plots ###
################################################

###
### 4.3) Calculate calibration using pseudo values
###
calc.calib.pv.moderate.ce <- function(pv.comb, p.est){

  ### Assign colnames to p.est
  colnames(p.est) <- paste("p.est", 1:ncol(p.est), sep = "")

  ### Combine data.raw and predicted probabilities
  data.pv <- cbind(pv.comb, p.est)

  ### Fit loess recalibration models to the uncensored observations at time t to calculate intercepts
  loess1 <- loess(pv.1 ~ p.est1, data = data.pv)
  loess2 <- loess(pv.2 ~ p.est2, data = data.pv)
  loess3 <- loess(pv.3 ~ p.est3, data = data.pv)
  loess4 <- loess(pv.4 ~ p.est4, data = data.pv)
  loess5 <- loess(pv.5 ~ p.est5, data = data.pv)
  loess6 <- loess(pv.6 ~ p.est6, data = data.pv)

  ### Created 'predicted observed' probabilities for each individual
  data.pv$loess.pred.obs1 <- predict(loess1, newdata = data.pv)
  data.pv$loess.pred.obs2 <- predict(loess2, newdata = data.pv)
  data.pv$loess.pred.obs3 <- predict(loess3, newdata = data.pv)
  data.pv$loess.pred.obs4 <- predict(loess4, newdata = data.pv)
  data.pv$loess.pred.obs5 <- predict(loess5, newdata = data.pv)
  data.pv$loess.pred.obs6 <- predict(loess6, newdata = data.pv)

  ### Create plot titles that will be used
  plot.titles <- c("Transplant", "Recovery", "AE", "Recovery + AE", "Relapse", "Death")
  plots.list <- vector("list", 6)
  for (i in 1:6){

    ### Creaet variables to plot
    data.pv$pred <- data.pv[, paste("p.est", i, sep = "")]
    data.pv$obs <- data.pv[, paste("loess.pred.obs", i, sep = "")]

    ### Create the plots
    plots.list[[i]] <- ggplot2::ggplot(data = data.pv %>%
                                dplyr::arrange(pred) %>% dplyr::select(pred, obs),
                                ggplot2::aes(x = pred, y = obs, color = "red")) +
      ggplot2::geom_line() +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk") +
      ggplot2::xlim(c(0, max(data.pv$pred))) +
      ggplot2::ylim(c(0, max(data.pv$pred))) +
      ggplot2::geom_rug(ggplot2::aes(x = pred, y = obs), col = grDevices::rgb(1, 0, 0, alpha = .3)) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(paste("State ", i, " (", plot.titles[i], ")", sep = ""))
  }


  ### Create output object and return it
  output.object <- list("plots.list" = plots.list)
  return(output.object)

}

### Create plots
plots.pv.j1s0 <- calc.calib.pv.moderate.ce(pv.comb, dplyr::select(data.raw, paste("pstate", 1:6, sep = "")))

plots.ggcombine <- ggpubr::ggarrange(plotlist = plots.pv.j1s0[["plots.list"]], nrow = 2, ncol = 3, common.legend = TRUE)

Cairo::CairoPNG(paste("workflow/figures/supp_pseudo_j1s0_pseudo.png", sep = ""),
         dpi = 300, width = 15, height = 10, unit = "in")
print(plots.ggcombine)
dev.off()

save.image("workflow/supp_pseudo.RData")


##################################################
### Calibration plots using our normal methods ###
##################################################


###
### With vars, no CI
###

### Create calibration plots data
dat.calib.blr <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = 1826,
                 tp.pred = tps0 %>% subset(j == 1) %>%
                   dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                 curve.type = "loess",
                 loess.span = 0.75,
                 loess.degree = 2,
                 w.covs = c("agecl", "year", "proph", "match"))

### Create calibration plot
plot.calibmsm.loess <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)

### Save it
Cairo::CairoPNG(paste("workflow/figures/supp_pseudo_j1s0_calibmsm.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plot.calibmsm.loess)
dev.off()


### Create calibration plots data
dat.calib.blr <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = 1826,
                 tp.pred = tps0 %>% subset(j == 1) %>%
                   dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                 curve.type = "loess",
                 loess.span = 0.75,
                 loess.degree = 2,
                 w.covs = c("agecl", "proph", "match"))

### Create calibration plot
plot.calibmsm.loess <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)

### Save it
Cairo::CairoPNG(paste("workflow/figures/supp_pseudo_j1s0_calibmsm_noyear.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plot.calibmsm.loess)
dev.off()



###
### With vars, no CI, max follow
###

### Create calibration plots data
dat.calib.blr <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = 1826,
                 tp.pred = tps0 %>% subset(j == 1) %>%
                   dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                 curve.type = "loess",
                 loess.span = 0.75,
                 loess.degree = 2,
                 w.covs = c("agecl", "year", "proph", "match"),
                 w.max.follow = "t.eval")


### Create calibration plot
plot.calibmsm.loess <- plot(dat.calib.blr2, combine = TRUE, nrow = 2, ncol = 3)

### Save it
Cairo::CairoPNG(paste("workflow/figures/supp_pseudo_j1s0_calibmsm_maxfollow.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plot.calibmsm.loess)
dev.off()


###
### With vars, with CI
###

### Create calibration plots data
dat.calib.blr <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = 1826,
                 tp.pred = tps0 %>% subset(j == 1) %>%
                   dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                 curve.type = "loess",
                 loess.span = 0.75,
                 loess.degree = 2,
                 w.covs = c("agecl", "year", "proph", "match"),
                 CI = 95,
                 CI.R.boot = 200)

### Create calibration plot
plot.calibmsm.loess <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)

### Save it
Cairo::CairoPNG(paste("workflow/figures/supp_pseudo_j1s0_calibmsm_ci.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plot.calibmsm.loess)
dev.off()


###
### null vars, no CI
###

### Create calibration plots data
dat.calib.blr <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = 1826,
                 tp.pred = tps0 %>% subset(j == 1) %>%
                   dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                 curve.type = "loess",
                 loess.span = 0.75,
                 loess.degree = 2,
                 w.covs = NULL)

### Create calibration plot
plot.calibmsm.loess <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)

### Save it
Cairo::CairoPNG(paste("workflow/figures/supp_pseudo_j1s0_calibmsm_null.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plot.calibmsm.loess)
dev.off()


###
### null vars, no CI
###

### Create calibration plots data
dat.calib.blr <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = 1826,
                 tp.pred = tps0 %>% subset(j == 1) %>%
                   dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                 curve.type = "loess",
                 loess.span = 0.75,
                 loess.degree = 2,
                 w.covs = NULL,
                 CI = 95,
                 CI.R.boot = 200)

### Create calibration plot
plot.calibmsm.loess <- plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)

### Save it
Cairo::CairoPNG(paste("workflow/figures/supp_pseudo_j1s0_calibmsm_null_ci.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plot.calibmsm.loess)
dev.off()

###############################################################################################################
### Lets explore the mechanism for estimating the weights, as the calibration plot for state 3 seems weird! ###
###############################################################################################################

# ### Define parameters
# t.eval <- 1826
# j <- 1
# s <- 0
# landmark.type <- "state"
# covs <- c("agecl", "year", "proph", "match")
# max.follow <- "t.eval"
# data.raw <- ebmtcal
# data.mstate <- msebmtcal
#
# ### Modify everybody to be censored after time t.eval, if a max.follow has been specified
# if(!is.null(max.follow)){
#   if (max.follow == "t.eval"){
#     data.raw <- dplyr::mutate(data.raw,
#                               dtcens.s = dplyr::case_when(dtcens < t.eval + 2 ~ dtcens.s,
#                                                           dtcens >= t.eval + 2 ~ 0),
#                               dtcens = dplyr::case_when(dtcens < t.eval + 2 ~ dtcens,
#                                                         dtcens >= t.eval + 2 ~ t.eval + 2))
#   } else {
#     data.raw <- dplyr::mutate(data.raw,
#                               dtcens.s = dplyr::case_when(dtcens < max.follow + 2 ~ dtcens.s,
#                                                           dtcens >= max.follow + 2 ~ 0),
#                               dtcens = dplyr::case_when(dtcens < max.follow + 2 ~ dtcens,
#                                                         dtcens >= max.follow + 2 ~ max.follow))
#   }
# }
#
# exp(-0.03)
#
# sum(ebmtcal$dtcens.s == 1 & ebmtcal$dtcens < 1826)
# ebmtcal %>% group_by(year, dtcens.s) %>% summarize(mean = mean(dtcens),
#                                                     median = median(dtcens),
#                                                     max = max(dtcens))
#
# data.raw %>% group_by(year, dtcens.s) %>% summarize(mean = mean(dtcens),
#                                                    median = median(dtcens),
#                                                    max = max(dtcens))
#
# ebmtcal %>% group_by(year) %>% summarize(mean = mean(dtcens.s))
# data.raw %>% group_by(year) %>% summarize(mean = mean(dtcens.s))
#
# data.raw.test <- ebmtcal %>% subset(year == "1995-1998")
#
# data.raw %>% group_by(year, rel.s) %>% summarize(mean = mean(rel),
#                                                     median = median(rel),
#                                                     max = max(rel))
#
# data.raw %>% group_by(year, srv.s) %>% summarize(mean = mean(srv),
#                                                     median = median(srv),
#                                                     max = max(srv))
#
#
# hist(ebmtcal %>% subset(year == "1995-1998" & dtcens.s == 0) %>% pull(dtcens))
# hist(ebmtcal %>% subset(year == "1985-1989" & dtcens.s == 0) %>% pull(dtcens))
#
# ### Create a new outcome, which is the time until censored from s
# data.raw$dtcens.modified <- data.raw$dtcens - s
#
# ### Save a copy of data.raw
# data.raw.save <- data.raw
#
# ### If landmark.type = "state", calculate weights only in individuals in state j at time s
# ### If landmark.type = "all", calculate weights in all uncensored individuals at time s (note that this excludes individuals
# ### who have reached absorbing states, who have been 'censored' from the survival distribution is censoring)
# if (landmark.type == "state"){
#   ### Identify individuals who are uncensored in state j at time s
#   ids.uncens <- base::subset(data.mstate, from == j & Tstart <= s & s < Tstop) %>%
#     dplyr::select(id) %>%
#     dplyr::distinct(id) %>%
#     dplyr::pull(id)
#
# } else if (landmark.type == "all"){
#   ### Identify individuals who are uncensored time s
#   ids.uncens <- base::subset(data.mstate, Tstart <= s & s < Tstop) %>%
#     dplyr::select(id) %>%
#     dplyr::distinct(id) %>%
#     dplyr::pull(id)
#
# }
#
# ### Subset data.mstate and data.raw to these individuals
# data.mstate <- data.mstate %>% base::subset(id %in% ids.uncens)
# data.raw <- data.raw %>% base::subset(id %in% ids.uncens)
#
# ###
# ### Create models for censoring in order to calculate the IPCW weights
# ### Seperate models for estimating the weights, and stabilising the weights (intercept only model)
# ###
# if (!is.null(covs)){
#   ### A model where we adjust for predictor variables
#   cens.model <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens.modified, dtcens.s) ~ ",
#                                                         paste(covs, collapse = "+"),
#                                                         sep = "")),
#                                 data = data.raw)
#
#   ### Intercept only model (numerator for stabilised weights)
#   cens.model.int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens.modified, dtcens.s) ~ 1",
#                                                             sep = "")),
#                                     data = data.raw)
# } else if (is.null(covs)){
#   ### If user has not input any predictors for estimating weights, the model for estimating the weights is the intercept only model (i.e. Kaplan Meier estimator)
#
#   ### Intercept only model (numerator for stabilised weights)
#   cens.model.int <- survival::coxph(stats::as.formula(paste("survival::Surv(dtcens.modified, dtcens.s) ~ 1",
#                                                             sep = "")),
#                                     data = data.raw)
#   ### Assign cens.model to be the same
#   cens.model <- cens.model.int
#
#
# }
#
# ### The coefficients for the "year" weight is really big
# max(data.mstate$time[data.mstate$to == 3 & data.mstate$status == 1])
#
# ### Get id's of people who have transitions to state 3
# ids.2 <- data.mstate$id[data.mstate$from == 1 & data.mstate$to == 2 & data.mstate$status == 1]
# ids.3 <- data.mstate$id[data.mstate$from == 1 & data.mstate$to == 3 & data.mstate$status == 1]
# ids.5 <- data.mstate$id[data.mstate$from == 1 & data.mstate$to == 5 & data.mstate$status == 1]
# ids.6 <- data.mstate$id[data.mstate$from == 1 & data.mstate$to == 6 & data.mstate$status == 1]
#
# ### Calculate a data frame containing probability of censored and uncenosred at each time point
# ### The weights will be the probability of being uncensored, at the time of the event for each individual
#
# ## Extract baseline hazard
# data.weights <- survival::basehaz(cens.model, centered = TRUE)
# ## Add lp to data.raw.save
# data.raw.save$lp <- stats::predict(cens.model, newdata = data.raw.save, type = "lp", reference = "zero")
#
#
# data.weights$hazard[max(which(data.weights$time <= t.eval - s))]
#
# tempsurvfit <- survfit(cens.model, data = data.raw)
# #str(tempsurvfit)
# tempsurvfit$cumhaz[max(which(tempsurvfit$time <= t.eval - s))]
# exp(-0.05)
#
# tempsurvfit <- survfit(Surv(dtcens, dtcens.s) ~ 1, data = data.raw)
# tempsurvfit$cumhaz[max(which(tempsurvfit$time <= t.eval - s))]
#
#
# ### Create weights for the cohort at time t.eval - s
# ### Note for individuals who reached an absorbing state, we take the probability of them being uncensored at the time of reached the
# ### abosrbing state. For individuals still alive, we take the probability of being uncensored at time t.eval - s.
#
# ### Get location of individuals who entered absorbing states or were censored prior to evaluation time
# obs.absorbed.prior <- which(data.raw.save$dtcens <= t.eval & data.raw.save$dtcens.s == 0)
# obs.censored.prior <- which(data.raw.save$dtcens <= t.eval & data.raw.save$dtcens.s == 1)
#
# ###
# ### Now create unstabilised probability of (un)censoring weights
# ### Note that weights are the probability of being uncensored, so if an individual has low probability of being uncesored,
# ### the inervse of this will be big, weighting them strongly
# ###
#
# ### First assign all individuals a weight of the probability of being uncensored at time t.eval
# ### This is the linear predictor times the cumulative hazard at time t.eval, and appropriate transformation to get a risk
# data.raw.save$pcw <- as.numeric(exp(-exp(data.raw.save$lp)*data.weights$hazard[max(which(data.weights$time <= t.eval - s))]))
#
# ## Write a function which will extract the uncensored probability for an individual with linear predictor lp at a given time t
# prob.uncens.func <- function(input){
#
#   ## Assign t and person_id
#   t <- input[1]
#   lp <- input[2]
#
#   if (t <= 0){
#     return(NA)
#   } else if (t > 0){
#     ## Get hazard at appropriate time
#     if (t < min(data.weights$time)){
#       bhaz.t <- 0
#     } else if (t >= min(data.weights$time)){
#       bhaz.t <- data.weights$hazard[max(which(data.weights$time <= t))]
#     }
#
#     ## Return risk
#     return(exp(-exp(lp)*bhaz.t))
#   }
# }
#
# ### Apply this function to all the times at which individuals have entered an absorbing state prior to censoring
# data.raw.save$pcw[obs.absorbed.prior] <- apply(data.raw.save[obs.absorbed.prior, c("dtcens.modified", "lp")], 1, FUN = prob.uncens.func)
#
# ### For individuals who were censored prior to t.eval, assign the weight as NA
# data.raw.save$pcw[obs.censored.prior] <- NA
#
# ### Invert these
# data.raw.save$ipcw <- 1/data.raw.save$pcw
#
# save.image("workflow/supp_pseudo_explore_weights.RData")
