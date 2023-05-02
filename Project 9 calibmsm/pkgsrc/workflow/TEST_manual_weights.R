### Load packages
library(calibmsm)
library(dplyr)

########################################################
### Section 3: Clinical setting and data preparation ###
########################################################

### Load and view data
data("ebmtcal")
head(ebmtcal)

data("msebmtcal")
head(msebmtcal)

data("tps0")
head(tps0)

data("tps100")
head(tps100)

#######################################################
### Section 4: Calibration curves and scatter plots ###
#######################################################

###
### Section 4.1: Calibration plots from state j = 1 at time s = 0
###

### Assign follow up time to assess calibration (5 years)
t.eval <- 1826

### Fit model to generate observed event probabilities using BLR-IPCW
dat.calib.blr <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=0,
                 t.eval = t.eval,
                 tp.pred = tps0 %>%
                   dplyr::filter(j == 1) %>%
                   dplyr::select(any_of(paste("pstate", 1:6, sep = ""))),
                 curve.type = "rcs",
                 rcs.nk = 3,
                 w.covs = c("year", "agecl", "proph", "match"))

### View structure of output
str(dat.calib.blr[["plotdata"]])
str(dat.calib.blr[["metadata"]])

### Plot Figure 1
plot(dat.calib.blr, combine = TRUE, nrow = 2, ncol = 3)


###
### Section 5.2: Estimating confidence interval manually using bootstrapping
###

### Store predicted transition probabilities in an dataframe so they can be bootstrapped alongside ebmtcal
tp.pred.s0 <- tps0 %>%
  dplyr::filter(j == 1) %>%
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

tp.pred.s100 <- tps100 %>%
  dplyr::filter(j == 1) %>%
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))

### Identify individuals who are uncensored at time t.eval
ids.uncens <- ebmtcal %>%
  subset(dtcens > t.eval | (dtcens < t.eval & dtcens.s == 0)) %>%
  dplyr::pull(id)

### Extract the predicted risks out of state 1 for these individuals
### All bootstrapped calibration curves will be plotted over these values
data.pred.plot <- tps0 %>%
  dplyr::filter(j == 1 & id %in% ids.uncens) %>%
  dplyr::select(any_of(paste("pstate", 1:6, sep = "")))


###
### Estimate calibration curve using manually estimated weights
###

### First estimate the weights
### We are using the internal function from calibmsm, but in practice your own function should be defined here
s.in <- 0
weights.manual <-
  calc_weights(data.mstate = msebmtcal,
               data.raw = ebmtcal,
               covs = c("year", "agecl", "proph", "match"),
               t.eval = t.eval,
               s = s.in,
               landmark.type = "state",
               j = 1,
               max.weight = 10,
               stabilised = FALSE)$ipcw
str(weights.manual)
sum(!is.na(weights.manual))

test.data <- ebmtcal[!is.na(weights.manual),]
### Estimate the calibration curve, manually specifying the weights
dat.calib.boot.manual <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=s.in,
                 t.eval = t.eval,
                 tp.pred = tp.pred.s0,
                 curve.type = "rcs",
                 rcs.nk = 3,
                 weights = weights.manual)
str(dat.calib.boot.manual[[1]][[1]])

str(ids.plot)
ids.plot <- dat.calib.boot.manual[[1]][[1]]$id
test.data.ids <-  test.data %>% subset(!(id %in% ids.plot))





s.in <- 100
weights.manual <-
  calc_weights(data.mstate = msebmtcal,
               data.raw = ebmtcal,
               covs = c("year", "agecl", "proph", "match"),
               t.eval = t.eval,
               s = s.in,
               landmark.type = "state",
               j = 1,
               max.weight = 10,
               stabilised = FALSE)$ipcw
str(weights.manual)
sum(!is.na(weights.manual))
### Estimate the calibration curve, manually specifying the weights
dat.calib.boot.manual <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=s.in,
                 t.eval = t.eval,
                 tp.pred = tp.pred.s100,
                 curve.type = "rcs",
                 rcs.nk = 3,
                 weights = weights.manual)
str(dat.calib.boot.manual[[1]][[1]])

dat.calib.boot.manual <-
  calc_calib_blr(data.mstate = msebmtcal,
                 data.raw = ebmtcal,
                 j=1,
                 s=s.in,
                 t.eval = t.eval,
                 tp.pred = tp.pred.s100,
                 curve.type = "rcs",
                 rcs.nk = 3,
                 w.covs = c("year", "agecl", "match"))
str(dat.calib.boot.manual[[1]][[1]])


