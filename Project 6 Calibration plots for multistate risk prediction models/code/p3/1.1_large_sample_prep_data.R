###
### This program will prep the data ready to calculate calibration for the cohort size and scenario of interest
### Aalen-Johansen estimator of the entire cohort will be calculated
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
print(paste("scen = ", scen, sep = ""))

### Load packages
source("code/z_functions.R")
source("code/z_load_packages.R")
source("code/p2/0_input_parameters.R")

### Load simulated data
data.pop.mstate <- readRDS(paste("data/sim/popdata/pop.mstate.dgm1.", scen, ".RData", sep = ""))
data.pop.raw <- readRDS(paste("data/sim/popdata/pop.raw.dgm1.", scen, ".RData", sep = ""))
tmat <- readRDS(paste("data/sim/popdata/tmat.dgm1.RData", sep = ""))

### Set seed
set.seed(101)

### Extract 200,000 individuals, which is what we're doing the large sample analysis with
n.cohort <- 200000
data.mstate <- data.pop.mstate[data.pop.mstate$patid %in% 1:n.cohort, ]
data.raw <- data.pop.raw[data.pop.raw$patid %in% 1:n.cohort, ]
str(data.mstate)
str(data.raw)

### Extract true transition probabilities
tp.true <- data.raw %>% select(paste("p.true", 1:5, sep = ""))

### Now defined the predicted transition probabilities based off these
tp.pred <- vector("list", 3)

### First element is perfectly prediction probabilities
tp.pred[[1]] <- dplyr::mutate(tp.true,
                              tp1 = p.true1,
                              tp2 = p.true2,
                              tp3 = p.true3,
                              tp4 = p.true4,
                              tp5 = p.true5) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Second set of predicted transition probabilities we add vector (0.5, 0.25, 0, -0.25, -0.5) to the log-odds, convert back and normalise
tp.pred[[2]] <- dplyr::mutate(tp.true,
                              tp1 = exp(log(p.true1/(1 - p.true1)) + 0.5)/(1 + exp(log(p.true1/(1 - p.true1)) + 0.5)),
                              tp2 = exp(log(p.true2/(1 - p.true2)) + 0.25)/(1 + exp(log(p.true2/(1 - p.true2)) + 0.25)),
                              tp3 = p.true3,
                              tp4 = exp(log(p.true4/(1 - p.true4)) - 0.25)/(1 + exp(log(p.true4/(1 - p.true4)) - 0.25)),
                              tp5 = exp(log(p.true5/(1 - p.true5)) - 0.5)/(1 + exp(log(p.true5/(1 - p.true5)) - 0.5)),
                              tp1 = tp1/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp2 = tp2/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp3 = tp3/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp4 = tp4/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp5 = tp5/(tp1 + tp2 + tp3 + tp4 + tp5)) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)


### Third set of predicted transition probabilities  (0.5, -0.5, -0.5, -0.5, 0.5) to the log-odds, convert back and normalise
tp.pred[[3]] <- dplyr::mutate(tp.true,
                              tp1 = exp(log(p.true1/(1 - p.true1)) + 0.5)/(1 + exp(log(p.true1/(1 - p.true1)) + 0.5)),
                              tp2 = exp(log(p.true2/(1 - p.true2)) - 0.5)/(1 + exp(log(p.true2/(1 - p.true2)) - 0.5)),
                              tp3 = exp(log(p.true3/(1 - p.true3)) - 0.5)/(1 + exp(log(p.true3/(1 - p.true3)) - 0.5)),
                              tp4 = exp(log(p.true4/(1 - p.true4)) - 0.5)/(1 + exp(log(p.true4/(1 - p.true4)) - 0.5)),
                              tp5 = exp(log(p.true5/(1 - p.true5)) + 0.5)/(1 + exp(log(p.true5/(1 - p.true5)) + 0.5)),
                              tp1 = tp1/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp2 = tp2/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp3 = tp3/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp4 = tp4/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp5 = tp5/(tp1 + tp2 + tp3 + tp4 + tp5)) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)


### Need to create the censoring variable which will be used to estimate th weights when assessing calibration
### dtcens is time until either censored, or entering an absorbing state
### dtcens.s == 1 if individual is censored, dtcens.s == 0 if absorbing state is entered
### Also add an "id" variable which is required for calibmsm functionality
data.raw <- mutate(data.raw,
                   dtcens = pmin(State.5, cens.times, na.rm = TRUE),
                   dtcens.s = case_when(dtcens == cens.times ~ 1,
                                        dtcens == State.5 ~ 0),
                   id = patid)

data.mstate$id <- data.mstate$patid

### 
### Calculate Aalen-Johansen estimator in entire cohort
print(paste("AJ ", Sys.time()))
obs.aj.object <- calc.calib.aj(data.mstate, 
                               tmat, 
                               t.eval)
print(paste("AJ ", Sys.time()))

### Extract standard errors and estimates seperately
obs.aj <- obs.aj.object[["obs.aj"]]
obs.aj.se <- obs.aj.object[["obs.aj.se"]]


###
### Save image
rm(data.pop.raw, data.pop.mstate)
save.image(paste("data/sim/large_sample_prep_data_", scen, ".RData", sep = ""))
print("IMAGE SAVED")