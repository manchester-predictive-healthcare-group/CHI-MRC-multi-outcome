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

### Assign n.pctls
n.pctls <- 3

### Create objects to store data sorted by predicted risk for each state k, and the Aalen-Johansen estimator for each subgroup
data.sort.groups <- vector("list", 6)
obs.aj.groups <- vector("list", 6)

### Assign data.raw by merging tps0 for j = 1 with embtcal
data.raw <- dplyr::inner_join(ebmtcal, tps0 %>% subset(j == 1), by = join_by("id"))

### Split ebmtcal into three datasets defined by year
data.raw.list <- vector("list", 3)
data.raw.list[[1]] <- subset(data.raw, year == "1985-1989")
data.raw.list[[2]] <- subset(data.raw, year == "1990-1994")
data.raw.list[[3]] <- subset(data.raw, year == "1995-1998")

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
  data.raw.list.sort.temp <- vector("list", 3)
  for (group in 1:3){
    data.raw.list.sort.temp[[group]] <- dplyr::arrange(data.raw.list[[group]], !!sym(paste("pstate", state, sep = "")))
  }

  ### Create a list of length n.pctls
  data.sort.groups[[state]] <- vector("list", 3)
  group.size <- vector("numeric", 3)
  for (group in 1:3){
    data.sort.groups[[state]][[group]] <- vector("list", n.pctls)
    group.size[group] <- nrow(data.raw.list.sort.temp[[group]])/n.pctls
    for (pctl in 1:n.pctls){
      data.sort.groups[[state]][[group]][[pctl]] <-
        dplyr::slice(data.raw.list.sort.temp[[group]], ((round((pctl-1)*group.size[group]))+1):(round(pctl*group.size[group])))
    }
  }

  ###
  ### Calculate the Aalen-Johansen estimator within each group
  ###

  ### Create a list to store them
  obs.aj.groups[[state]] <- vector("list", 3)
  for (group in 1:3){
    obs.aj.groups[[state]][[group]] <- vector("list", n.pctls)
  }

  ### Calculate obs.aj for each group
  for (group in 1:3){
    for (pctl in 1:n.pctls){
      print(paste("Calc obs.aj = ", state, "group = ", group, "pctl = ", pctl, Sys.time()))
      obs.aj.groups[[state]][[group]][[pctl]] <-
        calc.calib.aj.ce(data.mstate = base::subset(data.mstate, id %in% data.sort.groups[[state]][[group]][[pctl]]$id),
                         tmat = attributes(data.mstate)$trans,
                         t.eval = t.eval)$obs.aj

    }
  }
}

##########################################
### Now to calculate the pseudo-values ###
##########################################

### Create list to store
pv.out.groups <- vector("list", 6)

### Run through and calculate pseudo-values
for (state in 2:6){
  ### Create output list
  pv.out.groups[[state]] <- vector("list", 3)
  for(group in 1:3){
    pv.out.groups[[state]][[group]] <- vector("list", 3)
    for (pctl in 1:3){
      print(paste("state = ", state, "group = ", group, Sys.time()))
      ### Get vector of patient ids from the state and group combo
      ids.vec <- data.sort.groups[[state]][[group]][[pctl]]$id
      ### Calculate pseudo-values and store
      pv.temp <- lapply(ids.vec, func.calc.pv.aj.ce,
                        data.mstate = subset(data.mstate, id %in% data.sort.groups[[state]][[group]][[pctl]]$id),
                        obs.aj = obs.aj.groups[[state]][[group]][[pctl]],
                        tmat = attributes(data.mstate)$trans,
                        n.cohort = nrow(data.sort.groups[[state]][[group]][[pctl]]),
                        t.eval = t.eval)

      ### Only interested in pseudo-values for the state of interest (which we sorted by)
      pv.temp <- do.call("rbind", pv.temp)[, state]

      ### Combine with person_ids, so we have a record of who these pseudo-values are for
      pv.out.groups[[state]][[group]][[pctl]] <- data.frame(cbind(ids.vec, pv.temp))
      colnames(pv.out.groups[[state]][[group]][[pctl]]) <- c("id", paste("pv.", state, sep = ""))

    }

    ### Rbind pseudo-values for each percentile within each group
    pv.out.groups[[state]][[group]] <- do.call("rbind", pv.out.groups[[state]][[group]])

  }

  ### Rbind pseudo-values for each each group and sort by id
  pv.out.groups[[state]] <- do.call("rbind", pv.out.groups[[state]])
  pv.out.groups[[state]] <- dplyr::arrange(pv.out.groups[[state]], id)

}

### combine the pseudo-values into one dataset
pv.comb <- inner_join(pv.out.groups[[1]], pv.out.groups[[2]], by = join_by("id")) %>%
  inner_join(pv.out.groups[[3]], by = join_by("id")) %>%
  inner_join(pv.out.groups[[4]], by = join_by("id")) %>%
  inner_join(pv.out.groups[[5]], by = join_by("id")) %>%
  inner_join(pv.out.groups[[6]], by = join_by("id"))

save.image("workflow/supp_pseudo_year_pctl.RData")

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

Cairo::CairoPNG(paste("workflow/figures/supp_pseudo_j1s0_pseudo_year_pctl.png", sep = ""),
                dpi = 300, width = 15, height = 10, unit = "in")
print(plots.ggcombine)
dev.off()

save.image("workflow/supp_pseudo_year.RData")

