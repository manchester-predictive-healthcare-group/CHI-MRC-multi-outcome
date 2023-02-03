##### This program is to come up with the input parameters so the scenarios match Table 1 #####

### Clear workspace
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Source packages
source("code/sim_function_load_packages.R")

### Load functions for the simulation

source("code/sim_function_fit_models_all.R")
source("code/sim_function_generate_data_all.R")
source("code/sim_function_calc_true_risk_all.R")

## Followup time at which to evaluate
t.eval.sim <- 10

### The following paramters are specific to the data generating mechanism
## Shape and scale paramters for baseline hazard of each transition
  shape12.sim <- 1
  scale12.sim <- 15
  shape13.sim <- 1
  scale13.sim <- 10
  shape24.sim <- 1
  scale24.sim <- 10
  shape34.sim <- 1 
  scale34.sim <- 15
  
  ## Covariate effects for each transition
  beta12.cont.sim <- 0.25
  beta12.cat.sim <- 0.5
  beta13.cont.sim <- 0.3
  beta13.cat.sim <- 0.25
  beta24.cont.sim <- 0.1
  beta24.cat.sim <- 0.1
  beta34.cont.sim <- 0.25
  beta34.cat.sim <- 0.5

# beta12.cont.sim <- 0
# beta12.cat.sim <- 0
# beta13.cont.sim <- 0
# beta13.cat.sim <- 0
# beta24.cont.sim <- 0
# beta24.cat.sim <- 0
# beta34.cont.sim <- 0
# beta34.cat.sim <- 0
  
  ## Number of sampling steps when simulating data using multistate model
  numsteps.sim <- 50000
  

  ## Baseline hazard parameters for marginal distributions
  baselineA.sim <- c(1,5)
  baselineB.sim <- c(1,2.5)
  ## Covariate effects
  COV_betaA.sim <- c(0.25, 0.75)
  COV_betaB.sim <- c(0.5, 0.5)
  ## Copula association parameter
  copula.param.sim <- 0.
  ## Name of copula and rotation

    copula.sim <- "clayton"
    rotate.cop.sim <- 0

    copula.sim <- "clayton"
    rotate.cop.sim <- 90

    copula.sim <- "gumbel"
    rotate.cop.sim <- 0

    copula.sim <- "gumbel"
    rotate.cop.sim <- 90

    copula.sim <- "joe"
    rotate.cop.sim <- 0

    copula.sim <- "joe"
    rotate.cop.sim <- 90

    copula.sim <- "fgm"
    rotate.cop.sim <- 0


  ## Baseline hazard parameters for marginal distributions
  baselineA.sim <- c(1,5)
  baselineB.sim <- c(1,2.5)
  ## Covariate effects
  COV_betaA.sim <- c(0.5, 0.5)
  COV_betaB.sim <- c(0.25, 0.75)
  ## Frailty association parameter
  frail.var.sim <- 0.0001
  ## Assigm distribution

    frail.dist.sim <- "normal"
    int.method.sim <- "pcubature"

    frail.dist.sim <- "gamma"
    int.method.sim <- "hcubature"

#####################
### Suppose condition A the risk is 13%
### Suppose condition B the risk is 6%
t.eval.sim <- 10
x1.eval.sim <- 0
x2.eval.sim <- 0

## Baseline hazard parameters for marginal distributions
baselineA.sim <- c(1,70)
baselineB.sim <- c(1,140)
## Covariate effects
COV_betaA.sim <- c(0.25, 0.75)
COV_betaB.sim <- c(0.5, 0.5)


calc.true.risk.copula(
  t.eval = t.eval.sim, x1.eval = x1.eval.sim, x2.eval = x2.eval.sim,
  copula = "clayton", copula.param = 3, rotate.cop = rotate.cop.sim,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)


##################### Calc the true risks

true.risk.frailty.gamma <- calc.true.risk.frailty(
  t.eval = t.eval.sim, x1.eval = x1.eval.sim, x2.eval = x2.eval.sim,
  frail.dist = "gamma", frail.var = frail.var.sim,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
  int.method = "hcubature")

true.risk.frailty.normal <- calc.true.risk.frailty(
  t.eval = t.eval.sim, x1.eval = x1.eval.sim, x2.eval = x2.eval.sim,
  frail.dist = "normal", frail.var = frail.var.sim,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim,
  int.method = "pcubature")

true.risk.msm <- calc.true.risk.msm(t.eval = t.eval.sim, x1.eval = x1.eval.sim, x2.eval = x2.eval.sim,
                                    shape12 = shape12.sim, scale12 = scale12.sim,
                                    shape13 = shape13.sim, scale13 = scale13.sim,
                                    shape24 = shape24.sim, scale24 = scale24.sim,
                                    shape34 = shape34.sim, scale34 = scale34.sim,
                                    beta12.cont = beta12.cont.sim, beta12.cat = beta12.cat.sim,
                                    beta13.cont = beta13.cont.sim, beta13.cat = beta13.cat.sim,
                                    beta24.cont = beta24.cont.sim, beta24.cat = beta24.cat.sim,
                                    beta34.cont = beta34.cont.sim, beta34.cat = beta34.cat.sim)

true.risk.copula.clayton <- calc.true.risk.copula(
  t.eval = t.eval.sim, x1.eval = x1.eval.sim, x2.eval = x2.eval.sim,
  copula = "clayton", copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)

true.risk.copula.gumb <- calc.true.risk.copula(
  t.eval = t.eval.sim, x1.eval = x1.eval.sim, x2.eval = x2.eval.sim,
  copula = "gumbel", copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)

true.risk.copula.fgm <- calc.true.risk.copula(
  t.eval = t.eval.sim, x1.eval = x1.eval.sim, x2.eval = x2.eval.sim,
  copula = "fgm", copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)


########### Suppose we want to calculate some input parameters for copula approach ###############

## Baseline hazard parameters for marginal distributions
## Baseline hazard parameters for marginal distributions
baselineA.sim <- c(1,60)
baselineB.sim <- c(1,130)
## Covariate effects
COV_betaA.sim <- c(0.25, 0.75)
COV_betaB.sim <- c(0.5, 0.5)

## Copula association parameter
copula.param.sim <- 3
## Name of copula and rotation

copula.sim <- "clayton"
rotate.cop.sim <- 0


calc.true.risk.copula(
  t.eval = t.eval.sim, x1.eval = 0, x2.eval = 0,
  copula = "clayton", copula.param = copula.param.sim, rotate.cop = 0,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)

calc.true.risk.copula(
  t.eval = t.eval.sim, x1.eval = 0, x2.eval = 0,
  copula = "clayton", copula.param = copula.param.sim, rotate.cop = 0,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.A"]*calc.true.risk.copula(
    t.eval = t.eval.sim, x1.eval = 0, x2.eval = 0,
    copula = "clayton", copula.param = copula.param.sim, rotate.cop = 0,
    baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
    baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.B"]

calc.true.risk.copula(
  t.eval = t.eval.sim, x1.eval = 0, x2.eval = 1,
  copula = "clayton", copula.param = copula.param.sim, rotate.cop = 0,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.A"]*calc.true.risk.copula(
    t.eval = t.eval.sim, x1.eval = 0, x2.eval = 1,
    copula = "clayton", copula.param = copula.param.sim, rotate.cop = 0,
    baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
    baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.B"]

### Calculate the true overall joint risk in the cohort


### Want to make this into a function, so I can easily compare stuff
est.input.func <- function(baselineA.sim,
                           baselineB.sim,
                           ## Covariate effects
                           COV_betaA.sim,
                           COV_betaB.sim,
                           
                           ## Copula association parameter
                           copula.param.sim,
                           ## Name of copula and rotation
                           
                           copula.sim,
                           rotate.cop.sim){
  true.J.risk.c.clay.int.x20 <- function(x1){
    return(
      dnorm(x1, 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.sim, x1.eval = x1, x2.eval = 0,
        copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
        baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
        baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.joint.true"])
  }
  
  true.J.risk.c.clay.int.x21 <- function(x1){
    return(
      dnorm(x1, 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.sim, x1.eval = x1, x2.eval = 1,
        copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
        baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
        baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.joint.true"])
  }
  
  true.A.risk.c.clay.int.x20 <- function(x1){
    return(
      dnorm(x1, 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.sim, x1.eval = x1, x2.eval = 0,
        copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
        baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
        baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.A"])
  }
  
  true.A.risk.c.clay.int.x21 <- function(x1){
    return(
      dnorm(x1, 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.sim, x1.eval = x1, x2.eval = 1,
        copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
        baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
        baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.A"])
  }
  
  true.B.risk.c.clay.int.x20 <- function(x1){
    return(
      dnorm(x1, 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.sim, x1.eval = x1, x2.eval = 0,
        copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
        baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
        baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.B"])
  }
  
  true.B.risk.c.clay.int.x21 <- function(x1){
    return(
      dnorm(x1, 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.sim, x1.eval = x1, x2.eval = 1,
        copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
        baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
        baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.B"])
  }
  
  true.J.N.risk.c.clay.int.x20 <- function(x1){
    return(
      dnorm(x1, 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.sim, x1.eval = x1, x2.eval = 0,
        copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
        baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
        baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.A"]*calc.true.risk.copula(
          t.eval = t.eval.sim, x1.eval = x1, x2.eval = 0,
          copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
          baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
          baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.B"])
  }
  
  true.J.N.risk.c.clay.int.x21 <- function(x1){
    return(
      dnorm(x1, 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.sim, x1.eval = x1, x2.eval = 1,
        copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
        baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
        baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.A"]*calc.true.risk.copula(
          t.eval = t.eval.sim, x1.eval = x1, x2.eval = 1,
          copula = copula.sim, copula.param = copula.param.sim, rotate.cop = rotate.cop.sim,
          baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
          baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)["risk.marg.true.B"])
  }
  
  ### True joint risk over entire cohort
  true.joint.risk <- (cubintegrate(true.J.risk.c.clay.int.x20, lower = qnorm(0.00001, 0, 1), upper = qnorm(0.99999, 0, 1))$integral +
     cubintegrate(true.J.risk.c.clay.int.x21, lower = qnorm(0.00001, 0, 1), upper = qnorm(0.99999, 0, 1))$integral)/2
  
  ### Naive joint risk over entire cohort
  naive.joint.risk <- (cubintegrate(true.J.N.risk.c.clay.int.x20, lower = qnorm(0.00001, 0, 1), upper = qnorm(0.99999, 0, 1))$integral +
     cubintegrate(true.J.N.risk.c.clay.int.x21, lower = qnorm(0.00001, 0, 1), upper = qnorm(0.99999, 0, 1))$integral)/2
  
  ### Marginal risk of A
  true.marg.risk.A <- (cubintegrate(true.A.risk.c.clay.int.x20, lower = qnorm(0.00001, 0, 1), upper = qnorm(0.99999, 0, 1))$integral +
                    cubintegrate(true.A.risk.c.clay.int.x21, lower = qnorm(0.00001, 0, 1), upper = qnorm(0.99999, 0, 1))$integral)/2
  
  ### Marginal risk of B
  true.marg.risk.B <- (cubintegrate(true.B.risk.c.clay.int.x20, lower = qnorm(0.00001, 0, 1), upper = qnorm(0.99999, 0, 1))$integral +
                    cubintegrate(true.B.risk.c.clay.int.x21, lower = qnorm(0.00001, 0, 1), upper = qnorm(0.99999, 0, 1))$integral)/2
  
  return(c("true.joint.risk" = true.joint.risk,
         "naive.joint.risk" = naive.joint.risk,
         "true.marg.risk.A" = true.marg.risk.A,
         "true.marg.risk.B" = true.marg.risk.B))
}



### Want to make this into a function, so I can easily compare stuff
est.input.func(baselineA.sim = c(1,70),
               baselineB.sim = c(1,108.5),
               ## Covariate effects
               COV_betaA.sim = c(0.25, 0.75),
               COV_betaB.sim = c(0.5, 0.5),
                           
               ## Copula association parameter
               copula.param.sim = 1.79,
               ## Name of copula and rotation
               
               copula.sim = "clayton",
               rotate.cop.sim = 0)

### NEED TO REACH A FINAL DECISION ON COVARIATE EFFECTS FIRST AND LOCK THESE IN
### I THINK BIGGER EFFECTS IS A GOOD IDEA, MAYBE 1 AND 2?

est.input.func(baselineA.sim = c(1,70),
               baselineB.sim = c(1,108.5),
               ## Covariate effects
               COV_betaA.sim = c(0.25, 0.75),
               COV_betaB.sim = c(0.5, 0.5),
               
               ## Copula association parameter
               copula.param.sim = 1.79,
               ## Name of copula and rotation
               
               copula.sim = "clayton",
               rotate.cop.sim = 0)


t.eval.sim <- 10
x1.eval.sim <- -2
x2.eval.sim <- 1

## Baseline hazard parameters for marginal distributions
baselineA.sim <- c(1,70)
baselineB.sim <- c(1,108.5)
## Covariate effects ## I THINK THESE VALUES ARE PRETTY GOOD!! MEANS THE RISKS RANGE BETWEEN LIKE 3% to 50%!!!!
COV_betaA.sim <- c(0.5, 0.75)
COV_betaB.sim <- c(0.75, 0.5)

calc.true.risk.copula(
  t.eval = t.eval.sim, x1.eval = x1.eval.sim, x2.eval = x2.eval.sim,
  copula = "clayton", copula.param = 1.79, rotate.cop = 0,
  baselineA = baselineA.sim, COV_betaA = COV_betaA.sim, 
  baselineB = baselineB.sim, COV_betaB = COV_betaB.sim)