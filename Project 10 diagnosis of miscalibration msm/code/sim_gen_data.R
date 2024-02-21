### Data generating mechanism for illness-death model

### Set working directory
rm(list=ls())

### Load packages
source("z_functions.R")
source("z_load_packages.R")

### Extract arguments
args <- commandArgs(trailingOnly = T)
scenario <- as.numeric(args[1])

### Set parameters
n.cohort <- 10000
max.follow <- 1827
numsteps <- max.follow

shape12 <- 1
scale12 <- 5000
shape13 <- 1
scale13 <- 5000
shape23 <- 1
scale23 <- 5000

beta12.x1 <- 0.5
beta12.x2 <- -0.5
beta13.x1 <- 0.5
beta13.x2 <- -0.5
beta23.x1 <- 0.5
beta23.x2 <- -0.5

### Change input paramters depending on scenario
if (scenario == 2){
  ## Scale of transition 12 reduced
  scale12 <- 3700
} else if (scenario == 3){
  ## Scale of transition 13 reduced
  scale13 <- 3700
} else if (scenario == 4){
  ## Scale of transition 23 reduced
  scale23 <- 3700
} else if (scenario == 5){
  ## Effect of beta12.x2 changed
  beta12.x2 <- 0
} else if (scenario == 6){
  ## Effect of beta13.x2 changed
  beta13.x2 <- 0
} else if (scenario == 7){
  ## Effect of beta23.x2 changed
  beta23.x2 <- 0
}

### Set seed for censoring which happens at random
set.seed(101)

### Generate baseline data
x.baseline <- data.frame("x1" = rnorm(n.cohort, 0, 1), "x2" = rnorm(n.cohort, 0, 1))

### Generate transition data
dat1 <- DGM1(n = n.cohort, #number of patients to simulate
             max.follow = max.follow, #maximum follow up
             shape12 = shape12, scale12 = scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
             shape13 = shape13, scale13 = scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
             shape23 = shape23, scale23 = scale23, #shape and scale for weibull baseline hazard for transition 2 -> 3
             beta12.x1 = beta12.x1, beta12.x2 = beta12.x2, #covariate effects for transiion 12
             beta13.x1 = beta13.x1, beta13.x2 = beta13.x2, #covariate effects for transiion 13
             beta23.x1 = beta23.x1, beta23.x2 = beta23.x2, #covariate effects for transiion 23
             x.in = x.baseline, #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
             numsteps = numsteps)

### Structure of data
str(dat1$cohort)
head(dat1$cohort, n = 100)

### Define beta for censoring mechanism
cens_beta_x1 <- 0
cens_beta_x2 <- 0

### Covert to mstate format and apply censoring
data.convert.cens <- convert.mstate.cens(cohort.in = dat1[["cohort"]],
                                   max.follow = 1827,
                                   cens_shape = 1,
                                   cens_scale = 25000,
                                   cens_beta_x1 = cens_beta_x1,
                                   cens_beta_x2 = cens_beta_x2)
head(data.convert.cens[["data.mstate"]])
str(data.convert.cens[["data.mstate"]])
head(data.convert.cens[["data.raw"]])

print(paste("FINISH CONVERT MSTATE ", Sys.time(), sep = ""))

### Extract data.mstate and data.raw
data.mstate <- data.convert.cens[["data.mstate"]]
data.raw <- data.convert.cens[["data.raw"]]

### Add patid to data.raw
data.raw$id <- data.raw$patid

### Define formula for multistate model
adj.msm <- paste(paste("x1.", (1:3), sep = "", collapse = "+"),
                 paste("x2.", (1:3), sep = "", collapse = "+"),
                 sep = "+")
formula.msm <- paste("Surv(Tstart, Tstop, status) ~ ", adj.msm, " + strata(trans)", sep = "")

### Fit multistsate model
msm.cox.fit <- coxph(as.formula(formula.msm), data = data.mstate)

### Save this
save.image(paste("sim_gen_data.RData")


