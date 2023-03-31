### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages and functions
source("code/z_functions.R")
source("code/z_load_packages.R")
source("code/z_functions_true_transition_probs.R")

### Read in seed/iter
args <- commandArgs(trailingOnly = T)
set <- as.numeric(args[1])
print(paste("set = ", set, sep = ""))

### Generate the x1 and x2 values we want to generate for, according to the set argument
covars <- rbind(c(0, 0),
                c(1, 0),
                c(0, 1),
                c(1, 1),
                c(-1, -1))
x1.eval <- covars[set, 1]
x2.eval <- covars[set, 2]
print(paste("x1 = ", x1.eval, ", x2 = ", x2.eval, sep = ""))

### Generate some baseline data
n.cohort <- 500000
x.baseline <- data.frame("x1" = rep(x1.eval, n.cohort), "x2" = rep(x2.eval, n.cohort))

### Generate multistate data for these individuals
time.in <- Sys.time()
cohort <- gen.dat.DGM1(n = n.cohort, #number of patients to simulate
                         max.follow = ceiling(7*365.25), #maximum follow up
                         shape12 = 1, scale12 = 1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 2
                         shape13 = 1, scale13 = 0.5*1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 3
                         shape15 = 1, scale15 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 5
                         shape24 = 1, scale24 = 1588.598, #shape and scale for weibull baseline hazard for transition 2 -> 4
                         shape25 = 1, scale25 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 2 -> 5
                         shape34 = 1, scale34 = 0.5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 4
                         shape35 = 1, scale35 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 5
                         shape45 = 1, scale45 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 5
                         beta12.x1 = 1, beta12.x2 = 1, #covariate effects for transiion 12
                         beta13.x1 = 0.5, beta13.x2 = 0.5, #covariate effects for transiion 13
                         beta15.x1 = 1, beta15.x2 = 0.5, #covariate effects for transiion 15
                         beta24.x1 = 0.5, beta24.x2 = 1, #covariate effects for transiion 24
                         beta25.x1 = 1, beta25.x2 = 1, #covariate effects for transiion 25
                         beta34.x1 = 0.5, beta34.x2 = 0.5, #covariate effects for transiion 34
                         beta35.x1 = 1, beta35.x2 = 0.5, #covariate effects for transiion 35
                         beta45.x1 = 0.5, beta45.x2 = 1, #covariate effects for transiion 45
                         x.in = x.baseline, #baseline predictors, dataframe with two columns (x1 continuous, x2 continuous)
                         numsteps = ceiling(7*365.25))
time.out <- Sys.time()
time.diff <- time.out - time.in
print(paste("cohort generation took ", time.diff, sep = ""))
time.diff

### Convert into msate format
cohort.mstate <- convert.mstate.DGM1.nocens(cohort[["cohort"]], cohort[["max.follow"]])

### Create true trans probs
true.trans.probs <- calc.true.transition.probs.DGM1(u.eval = 0, t.eval = ceiling(7*365.25), x1.eval = x1.eval, x2.eval = x2.eval,
                                                    shape12 = 1, scale12 = 1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 2
                                                    shape13 = 1, scale13 = 0.5*1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 3
                                                    shape15 = 1, scale15 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 5
                                                    shape24 = 1, scale24 = 1588.598, #shape and scale for weibull baseline hazard for transition 2 -> 4
                                                    shape25 = 1, scale25 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 2 -> 5
                                                    shape34 = 1, scale34 = 0.5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 4
                                                    shape35 = 1, scale35 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 5
                                                    shape45 = 1, scale45 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 5
                                                    beta12.x1 = 1, beta12.x2 = 1, #covariate effects for transiion 12
                                                    beta13.x1 = 0.5, beta13.x2 = 0.5, #covariate effects for transiion 13
                                                    beta15.x1 = 1, beta15.x2 = 0.5, #covariate effects for transiion 15
                                                    beta24.x1 = 0.5, beta24.x2 = 1, #covariate effects for transiion 24
                                                    beta25.x1 = 1, beta25.x2 = 1, #covariate effects for transiion 25
                                                    beta34.x1 = 0.5, beta34.x2 = 0.5, #covariate effects for transiion 34
                                                    beta35.x1 = 1, beta35.x2 = 0.5, #covariate effects for transiion 35
                                                    beta45.x1 = 0.5, beta45.x2 = 1 #covariate effects for transiion 45
)

### Save image
save.image(paste("data/validate_true_transition_probabilities_DGM1_set", set, ".RData", sep = ""))
print("IMAGE SAVED")

# for (set in 1:5){
#   print(paste("set =", set))
#   load(paste("data/validate_true_transition_probabilities_DGM2_set", set, ".RData", sep = ""))
#   print(events(cohort.mstate[["data.wide"]])[["Frequencies"]][,"no event"]/n.cohort)
#   print(true.trans.probs)
#   print(events(cohort.mstate[["data.wide"]])[["Frequencies"]][,"no event"]/n.cohort - true.trans.probs)
# }