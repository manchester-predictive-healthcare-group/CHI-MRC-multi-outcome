###
### This file will test the bootstrap manually
###
rm(list=ls())
load_all()

data("ebmtcal")
data("msebmtcal")
data("tps0")
data("tps100")

### Start by defining key variables
data.mstate <- msebmtcal
data.raw <- ebmtcal
j<-1
s<-0
t.eval <- 1826
tp.pred = tps0 %>% dplyr::filter(j == !!j) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))
curve.type = "rcs"
rcs.nk = 3
weights <- NULL
w.covs = c("year", "agecl", "proph", "match")
w.landmark.type = "all"
w.max = 10
w.stabilised = FALSE
CI = FALSE

###
### Start by producing data using built in bootstrap
###
data.calib.boot.internal <- calc_calib_blr(data.mstate = msebmtcal,
                                           data.raw = ebmtcal,
                                           j=j,
                                           s=s,
                                           t.eval = t.eval,
                                           tp.pred = tp.pred,
                                           curve.type = curve.type,
                                           rcs.nk = rcs.nk,
                                           w.covs = w.covs,
                                           w.landmark.type = "state",
                                           w.max = 10,
                                           w.stabilised = FALSE,
                                           CI = 95,
                                           CI.R.boot = 200)


#######################################################################################################
### Now produce the data using the manual bootstrap approach and user defined estimation of weights ###
#######################################################################################################

###
### Step 1: Define the set of predicted risks over which we want to plot the calibration plots
###

### This should be the individuals uncensored at time t.eval
ids.uncens <- data.raw %>% subset(dtcens > t.eval | (dtcens < t.eval & dtcens.s == 0)) %>% dplyr::pull(id)
### Extract the predicted risks  for these individuals
data.pred.plot <- tps0 %>% dplyr::filter(j == !!j & id %in% ids.uncens) %>% dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

###
### Step 2: Create the calibration curve
###

### Calculate inverse probability of censoring weights
### We use the same function that is used internally in calc_calib_blr, but the point of this exercise is to define a function to calculate the weights yourself
### based on knowledge of the clinical settings
weights.manual <- calc_weights(data.mstate = msebmtcal,
                                     data.raw = ebmtcal,
                                     covs = c("year", "agecl", "proph", "match"),
                                     t.eval = t.eval,
                                     s = s,
                                     landmark.type = "state",
                                     j = j,
                                     max.weight = 10,
                                     stabilised = FALSE)$ipcw

### The structure of the weights object should be a numeric vector of length nrow(data.raw).
str(weights.manual)

### We have assigned weights of NA to individuals who are censored at t.eval, but this is irrelevant, as these individuals will not be included
### in the assessment of calibration at time t.eval.

### Calculate calibration curve using the manually calculated weights and extract the observed event rates
dat.calib.boot.manual <- calc_calib_blr(data.mstate = msebmtcal,
                                        data.raw = ebmtcal,
                                        j=j,
                                        s=s,
                                        t.eval = t.eval,
                                        tp.pred = tp.pred,
                                        curve.type = curve.type,
                                        rcs.nk = rcs.nk,
                                        weights = weights.manual,
                                        data.pred.plot = data.pred.plot
                                        )

###
### Step 3: Bootstrap this process manually
###

### Repeat write a function to repeat this process using boot for each state
calc_obs_boot <- function(data, indices, tp.pred, state.k){

  ### Bootstrap dataset
  data.boot <- data[indices,]
  tp.pred.boot <- tp.pred[indices, ]

  ### Calculate weights
  weights.manual <- calc_weights(data.mstate = msebmtcal,
                                 data.raw = data.boot,
                                 covs = c("year", "agecl", "proph", "match"),
                                 t.eval = t.eval,
                                 s = s,
                                 landmark.type = "state",
                                 j = j,
                                 max.weight = 10,
                                 stabilised = FALSE)$ipcw

  ### Calculate calibration curve and extract observed event rates
  curve.obs <- calc_calib_blr(data.mstate = msebmtcal,
                              data.raw = data.boot,
                              j=j,
                              s=s,
                              t.eval = t.eval,
                              tp.pred = tp.pred.boot,
                              curve.type = curve.type,
                              rcs.nk = rcs.nk,
                              weights = weights.manual,
                              data.pred.plot = data.pred.plot,
                              states.out = state.k)[["plotdata"]][[paste("state", state.k, sep = "")]]$obs

  return(curve.obs)

}

### Define alpha for CI's
alpha <- (1-95/100)/2

### Extract what states an individual can move into from state j (states with a non-zero predicted risk)
valid.transitions <- which(colSums(tp.pred) != 0)

### Create dataset for output
plot.data.list <- vector("list", length(valid.transitions))

### Run bootstrapping and create data for plots
for (k in 1:length(valid.transitions)){

  ### Assign state k
  state.k <- valid.transitions[k]

  ### Run bootstrapping
  boot.obs <- boot::boot(ebmtcal, calc_obs_boot, R = 100, tp.pred = tp.pred, state.k = state.k)$t

  ### Extract confidence bands
  lower <- apply(boot.obs, 2, stats::quantile, probs = alpha, na.rm = TRUE)
  upper <- apply(boot.obs, 2, stats::quantile, probs = 1-alpha, na.rm = TRUE)

  ### Assign output
  plot.data.list[[k]] <- data.frame(
    "pred" = dat.calib.boot.manual[["plotdata"]][[k]]$pred,
    "obs" = dat.calib.boot.manual[["plotdata"]][[k]]$obs,
    "obs.lower" = lower,
    "obs.upper" = upper)

}

### Want to put into same format as data.calib.boot.internal
### Create metadata
metadata <- list("valid.transitions"= valid.transitions,
                 "CI" = 95,
                 "curve.type" = "rcs")

data.calib.boot.external <- list("plotdata" = plot.data.list,
                                 "metadata" = metadata)


### Pivot longer to create data for ggplot and assign appropriate labels
plot1 <- plot.calib_blr(data.calib.boot.internal, combine = TRUE)
plot2 <- plot.calib_blr(data.calib.boot.external, combine = TRUE)

Cairo::CairoPNG("workflow/plot_internal.png", dpi = 300, width = 15, height = 10, unit = "in")
print(plot1)
dev.off()
Cairo::CairoPNG("workflow/plot_external.png", dpi = 300, width = 15, height = 10, unit = "in")
print(plot2)
dev.off()

