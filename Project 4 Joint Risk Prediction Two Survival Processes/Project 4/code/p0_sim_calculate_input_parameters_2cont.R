##### This program is to come up with some input parameters #####

### This is a test run of a what a simulation will look like ###

### Source packages
source("code/sim_function_load_packages.R")

### Clear workspace
rm(list=ls())

### Assign t.eval.sim
t.eval.sim <- 10

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Load functions for the simulation
source("code/sim_function_fit_models_all_2cont.R")
source("code/sim_function_generate_data_all_2cont.R")
source("code/sim_function_calc_true_risk_all_2cont.R")


### Write functions to return the marginal, joint, naive joint adjusted for X, and naive unadjusted joint risks


##################
### Copula DGM ###
##################
est.pop.risks.copula.func <- function(t.eval.in, 
                                  baselineA.in, ## Baseline shape and scale
                                  baselineB.in,
                                  
                                  COV_betaA.in, ## Covariate effects
                                  COV_betaB.in,
                                  
                                  copula.in, ## copula type
                                  copula.param.in, ## copula association parameter
                                  rotate.cop.in ## copula rotation
                                  ){
  # Function to integrate over to get joint risk
  true.J.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.joint.true"])
  }
  
  # Function to integrate over to get marginal risk of A
  true.A.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.marg.true.A"])
  }
  
  # Function to integrate over to get marginal risk of B
  true.B.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.marg.true.B"])
  }
  
  # Function to integrate over to get joint risk adjusted for X, but not taking into account residual correlation
  true.J.N.ADJ.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.copula(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.marg.true.A"]*calc.true.risk.copula(
          t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
          copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
          baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
          baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.marg.true.B"])
  }
  
  ### Create lower and upper limits
  lower.lim <- qnorm(0.00001,0,1)
  upper.lim <- qnorm(0.99999,0,1)
  
  ### True joint risk over entire cohort
  true.joint.risk <- cubintegrate(true.J.risk.int, 
                                  lower = c(lower.lim, lower.lim), 
                                  upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  ### Joint risk adjusted for X over entire cohort
  joint.adj.X.risk <- cubintegrate(true.J.N.ADJ.risk.int, 
                                   lower = c(lower.lim, lower.lim), 
                                   upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  ### Marginal risk of A
  true.marg.risk.A <- cubintegrate(true.A.risk.int, 
                                   lower = c(lower.lim, lower.lim), 
                                   upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  ### Marginal risk of B
  true.marg.risk.B <- cubintegrate(true.B.risk.int, 
                                   lower = c(lower.lim, lower.lim), 
                                   upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  return(c("true.joint.risk" = true.joint.risk,
         "joint.adj.X.risk" = joint.adj.X.risk,
         "naive.joint.risk" = true.marg.risk.A*true.marg.risk.B,
         "true.marg.risk.A" = true.marg.risk.A,
         "true.marg.risk.B" = true.marg.risk.B))
}

#############################
### Copula DGM vectorized ###
#############################
est.pop.risks.copula.func.vec <- function(t.eval.in, 
                                          baselineA.in, ## Baseline shape and scale
                                          baselineB.in,
                                          
                                          COV_betaA.in, ## Covariate effects
                                          COV_betaB.in,
                                          
                                          copula.in, ## copula type
                                          copula.param.in, ## copula association parameter
                                          rotate.cop.in ## copula rotation
){
  # Function to integrate over to get joint risk
  true.J.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.copula(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
      baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
      baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.joint.true"]}), ncol = ncol(x))
  }
  
  # Function to integrate over to get marginal risk of A
  true.A.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.copula(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
      baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
      baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.marg.true.A"]}), ncol = ncol(x))
  }
  
  # Function to integrate over to get marginal risk of B
  true.B.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.copula(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
      baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
      baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.marg.true.B"]}), ncol = ncol(x))
  }
  
  # Function to integrate over to get joint risk adjusted for X, but not taking into account residual correlation
  true.J.N.ADJ.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.copula(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
      baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
      baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.marg.true.A"]*calc.true.risk.copula(
        t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
        copula = copula.in, copula.param = copula.param.in, rotate.cop = rotate.cop.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in)["risk.marg.true.B"]}), ncol = ncol(x))
  }
  
  ### Create lower and upper limits
  lower.lim <- qnorm(0.00001,0,1)
  upper.lim <- qnorm(0.99999,0,1)
  
  ### True joint risk over entire cohort
  true.joint.risk <- cubature::pcubature(true.J.risk.int, 
                                         lower = c(lower.lim, lower.lim), 
                                         upper = c(upper.lim, upper.lim), 
                                         vectorInterface = TRUE)$integral
  
  ### Joint risk adjusted for X over entire cohort
  joint.adj.X.risk <- cubature::pcubature(true.J.N.ADJ.risk.int, 
                                              lower = c(lower.lim, lower.lim), 
                                              upper = c(upper.lim, upper.lim),
                                              vectorInterface = TRUE)$integral
  
  ### Marginal risk of A
  true.marg.risk.A <- cubature::pcubature(true.A.risk.int, 
                                          lower = c(lower.lim, lower.lim),
                                          upper = c(upper.lim, upper.lim),
                                          vectorInterface = TRUE)$integral
  
  ### Marginal risk of B
  true.marg.risk.B <- cubature::pcubature(true.B.risk.int, 
                                          lower = c(lower.lim, lower.lim),
                                          upper = c(upper.lim, upper.lim),
                                          vectorInterface = TRUE)$integral
  
  return(c("true.joint.risk" = true.joint.risk,
           "joint.adj.X.risk" = joint.adj.X.risk,
           "naive.joint.risk" = true.marg.risk.A*true.marg.risk.B,
           "true.marg.risk.A" = true.marg.risk.A,
           "true.marg.risk.B" = true.marg.risk.B))
}


###################
### Frailty DGM ###
###################  
est.pop.risks.frailty.func <- function(t.eval.in, 
                                       baselineA.in, ## Baseline shape and scale
                                       baselineB.in,
                                       
                                       COV_betaA.in, ## Covariate effects
                                       COV_betaB.in,
                                       
                                       frail.dist.in, ## frailty distribution
                                       frail.var.in ## frailty association parameter
){
  # Function to integrate over to get joint risk
  true.J.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.frailty(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        frail.dist = frail.dist.in, frail.var = frail.var.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.joint.true"])
  }
  
  # Function to integrate over to get marginal risk of A
  true.A.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.frailty(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        frail.dist = frail.dist.in, frail.var = frail.var.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.marg.true.A"])
  }
  
  # Function to integrate over to get marginal risk of B
  true.B.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.frailty(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        frail.dist = frail.dist.in, frail.var = frail.var.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.marg.true.B"])
  }
  
  # Function to integrate over to get joint risk adjusted for X, but not taking into account residual correlation
  true.J.N.ADJ.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.frailty(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        frail.dist = frail.dist.in, frail.var = frail.var.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.marg.true.A"]*calc.true.risk.frailty(
          t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
          frail.dist = frail.dist.in, frail.var = frail.var.in,
          baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
          baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.marg.true.B"])
  }
  
  ### Create lower and upper limits
  lower.lim <- qnorm(0.00001,0,1)
  upper.lim <- qnorm(0.99999,0,1)
  
  ### True joint risk over entire cohort
  true.joint.risk <- cubintegrate(true.J.risk.int, 
                                  lower = c(lower.lim, lower.lim), 
                                  upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  ### Joint risk adjusted for X over entire cohort
  joint.adj.X.risk <- cubintegrate(true.J.N.ADJ.risk.int, 
                                       lower = c(lower.lim, lower.lim), 
                                       upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  ### Marginal risk of A
  true.marg.risk.A <- cubintegrate(true.A.risk.int, 
                                   lower = c(lower.lim, lower.lim), 
                                   upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  ### Marginal risk of B
  true.marg.risk.B <- cubintegrate(true.B.risk.int, 
                                   lower = c(lower.lim, lower.lim), 
                                   upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  return(c("true.joint.risk" = true.joint.risk,
           "joint.adj.X.risk" = joint.adj.X.risk,
           "naive.joint.risk" = true.marg.risk.A*true.marg.risk.B,
           "true.marg.risk.A" = true.marg.risk.A,
           "true.marg.risk.B" = true.marg.risk.B))
}


##############################
### Frailty DGM vectorized ###
##############################  
est.pop.risks.frailty.func.vec <- function(t.eval.in, 
                                           baselineA.in, ## Baseline shape and scale
                                           baselineB.in,
                                           
                                           COV_betaA.in, ## Covariate effects
                                           COV_betaB.in,
                                           
                                           frail.dist.in, ## frailty distribution
                                           frail.var.in ## frailty association parameter
){
  # Function to integrate over to get joint risk
  true.J.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.frailty(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      frail.dist = frail.dist.in, frail.var = frail.var.in,
      baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
      baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.joint.true"]}), ncol = ncol(x))
  }
  
  # Function to integrate over to get marginal risk of A
  true.A.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.frailty(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      frail.dist = frail.dist.in, frail.var = frail.var.in,
      baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
      baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.marg.true.A"]}), ncol = ncol(x))
  }
  
  # Function to integrate over to get marginal risk of B
  true.B.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.frailty(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      frail.dist = frail.dist.in, frail.var = frail.var.in,
      baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
      baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.marg.true.B"]}), ncol = ncol(x))
  }
  
  # Function to integrate over to get joint risk adjusted for X, but not taking into account residual correlation
  true.J.N.ADJ.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.frailty(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      frail.dist = frail.dist.in, frail.var = frail.var.in,
      baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
      baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.marg.true.A"]*calc.true.risk.frailty(
        t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
        frail.dist = frail.dist.in, frail.var = frail.var.in,
        baselineA = baselineA.in, COV_betaA = COV_betaA.in, 
        baselineB = baselineB.in, COV_betaB = COV_betaB.in, int.method = "hcubature")["risk.marg.true.B"]}), ncol = ncol(x))
  }
  
  ### Create lower and upper limits
  lower.lim <- qnorm(0.001,0,1)
  upper.lim <- qnorm(0.999,0,1)
  
  ### True joint risk over entire cohort
  true.joint.risk <- cubature::pcubature(true.J.risk.int, 
                                         lower = c(lower.lim, lower.lim), 
                                         upper = c(upper.lim, upper.lim), 
                                         vectorInterface = TRUE)$integral
  
  ### Joint risk adjusted for X over entire cohort
  joint.adj.X.risk <- cubature::pcubature(true.J.N.ADJ.risk.int, 
                                              lower = c(lower.lim, lower.lim), 
                                              upper = c(upper.lim, upper.lim),
                                              vectorInterface = TRUE)$integral
  
  ### Marginal risk of A
  true.marg.risk.A <- cubature::pcubature(true.A.risk.int, 
                                          lower = c(lower.lim, lower.lim),
                                          upper = c(upper.lim, upper.lim),
                                          vectorInterface = TRUE)$integral
  
  ### Marginal risk of B
  true.marg.risk.B <- cubature::pcubature(true.B.risk.int, 
                                          lower = c(lower.lim, lower.lim),
                                          upper = c(upper.lim, upper.lim),
                                          vectorInterface = TRUE)$integral
  
  return(c("true.joint.risk" = true.joint.risk,
           "joint.adj.X.risk" = joint.adj.X.risk,
           "naive.joint.risk" = true.marg.risk.A*true.marg.risk.B,
           "true.marg.risk.A" = true.marg.risk.A,
           "true.marg.risk.B" = true.marg.risk.B))
}


###############
### MSM DGM ###
###############                              
est.pop.risks.msm.func <- function(t.eval.in, 
                                   shape12.in, scale12.in, #shape and scale for weibull baseline hazard for transition 1 -> 2
                                   shape13.in, scale13.in, #shape and scale for weibull baseline hazard for transition 1 -> 3
                                   shape24.in, scale24.in, #shape and scale for weibull baseline hazard for transition 2 -> 4
                                   shape34.in, scale34.in, #shape and scale for weibull baseline hazard for transition 3 -> 4
                                   beta12.cont.in, beta12.cont2.in, #covariate effects for transiion 12
                                   beta13.cont.in, beta13.cont2.in, #covariate effects for transiion 13
                                   beta24.cont.in, beta24.cont2.in, #covariate effects for transiion 24
                                   beta34.cont.in, beta34.cont2.in #covariate effects for transiion 34
){
  # Function to integrate over to get joint risk
  true.J.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.msm(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        shape12.in, scale12.in, 
        shape13.in, scale13.in, 
        shape24.in, scale24.in, 
        shape34.in, scale34.in,
        beta12.cont.in, beta12.cont2.in, 
        beta13.cont.in, beta13.cont2.in, 
        beta24.cont.in, beta24.cont2.in, 
        beta34.cont.in, beta34.cont2.in)["risk.joint.true"])
  }
  
  # Function to integrate over to get marginal risk of A
  true.A.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.msm(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        shape12.in, scale12.in, 
        shape13.in, scale13.in, 
        shape24.in, scale24.in, 
        shape34.in, scale34.in,
        beta12.cont.in, beta12.cont2.in, 
        beta13.cont.in, beta13.cont2.in, 
        beta24.cont.in, beta24.cont2.in, 
        beta34.cont.in, beta34.cont2.in)["risk.marg.true.A"])
  }
  
  # Function to integrate over to get marginal risk of B
  true.B.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.msm(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        shape12.in, scale12.in, 
        shape13.in, scale13.in, 
        shape24.in, scale24.in, 
        shape34.in, scale34.in,
        beta12.cont.in, beta12.cont2.in, 
        beta13.cont.in, beta13.cont2.in, 
        beta24.cont.in, beta24.cont2.in, 
        beta34.cont.in, beta34.cont2.in)["risk.marg.true.B"])
  }
  
  # Function to integrate over to get joint risk adjusted for X, but not taking into account residual correlation
  true.J.N.ADJ.risk.int <- function(x){
    return(
      dnorm(x[1], 0, 1)*dnorm(x[2], 0, 1)*calc.true.risk.msm(
        t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
        shape12.in, scale12.in, 
        shape13.in, scale13.in, 
        shape24.in, scale24.in, 
        shape34.in, scale34.in,
        beta12.cont.in, beta12.cont2.in, 
        beta13.cont.in, beta13.cont2.in, 
        beta24.cont.in, beta24.cont2.in, 
        beta34.cont.in, beta34.cont2.in)["risk.marg.true.A"]*calc.true.risk.msm(
          t.eval = t.eval.in, x1.eval = x[1], x2.eval = x[2],
          shape12.in, scale12.in, 
          shape13.in, scale13.in, 
          shape24.in, scale24.in, 
          shape34.in, scale34.in,
          beta12.cont.in, beta12.cont2.in, 
          beta13.cont.in, beta13.cont2.in, 
          beta24.cont.in, beta24.cont2.in, 
          beta34.cont.in, beta34.cont2.in)["risk.marg.true.B"])
  }
  
  ### Create lower and upper limits
  lower.lim <- qnorm(0.00001,0,1)
  upper.lim <- qnorm(0.99999,0,1)
  
  ### True joint risk over entire cohort
  true.joint.risk <- cubintegrate(true.J.risk.int, 
                                  lower = c(lower.lim, lower.lim), 
                                  upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  ### Joint risk adjusted for X over entire cohort
  joint.adj.X.risk <- cubintegrate(true.J.N.ADJ.risk.int, 
                                       lower = c(lower.lim, lower.lim), 
                                       upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  ### Marginal risk of A
  true.marg.risk.A <- cubintegrate(true.A.risk.int, 
                                   lower = c(lower.lim, lower.lim), 
                                   upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  ### Marginal risk of B
  true.marg.risk.B <- cubintegrate(true.B.risk.int, 
                                   lower = c(lower.lim, lower.lim), 
                                   upper = c(upper.lim, upper.lim), method = "hcubature")$integral
  
  return(c("true.joint.risk" = true.joint.risk,
           "joint.adj.X.risk" = joint.adj.X.risk,
           "naive.joint.risk" = true.marg.risk.A*true.marg.risk.B,
           "true.marg.risk.A" = true.marg.risk.A,
           "true.marg.risk.B" = true.marg.risk.B))
}


##########################
### msm DGM vectorized ###
##########################
est.pop.risks.msm.func.vec <- function(t.eval.in,                                    
                                       shape12.in, scale12.in, #shape and scale for weibull baseline hazard for transition 1 -> 2
                                       shape13.in, scale13.in, #shape and scale for weibull baseline hazard for transition 1 -> 3
                                       shape24.in, scale24.in, #shape and scale for weibull baseline hazard for transition 2 -> 4
                                       shape34.in, scale34.in, #shape and scale for weibull baseline hazard for transition 3 -> 4
                                       beta12.cont.in, beta12.cont2.in, #covariate effects for transiion 12
                                       beta13.cont.in, beta13.cont2.in, #covariate effects for transiion 13
                                       beta24.cont.in, beta24.cont2.in, #covariate effects for transiion 24
                                       beta34.cont.in, beta34.cont2.in #covariate effects for transiion 34
){
  # Function to integrate over to get joint risk
  true.J.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.msm(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      shape12.in, scale12.in, 
      shape13.in, scale13.in, 
      shape24.in, scale24.in, 
      shape34.in, scale34.in,
      beta12.cont.in, beta12.cont2.in, 
      beta13.cont.in, beta13.cont2.in, 
      beta24.cont.in, beta24.cont2.in, 
      beta34.cont.in, beta34.cont2.in)["risk.joint.true"]}), ncol = ncol(x))
  }
  
  # Function to integrate over to get marginal risk of A
  true.A.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.msm(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      shape12.in, scale12.in, 
      shape13.in, scale13.in, 
      shape24.in, scale24.in, 
      shape34.in, scale34.in,
      beta12.cont.in, beta12.cont2.in, 
      beta13.cont.in, beta13.cont2.in, 
      beta24.cont.in, beta24.cont2.in, 
      beta34.cont.in, beta34.cont2.in)["risk.marg.true.A"]}), ncol = ncol(x))
  }
  
  # Function to integrate over to get marginal risk of B
  true.B.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.msm(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      shape12.in, scale12.in, 
      shape13.in, scale13.in, 
      shape24.in, scale24.in, 
      shape34.in, scale34.in,
      beta12.cont.in, beta12.cont2.in, 
      beta13.cont.in, beta13.cont2.in, 
      beta24.cont.in, beta24.cont2.in, 
      beta34.cont.in, beta34.cont2.in)["risk.marg.true.B"]}), ncol = ncol(x))
  }
  
  # Function to integrate over to get joint risk adjusted for X, but not taking into account residual correlation
  true.J.N.ADJ.risk.int <- function(x){
    matrix(apply(x, 2, function(z) {dnorm(z[1], 0, 1)*dnorm(z[2], 0, 1)*calc.true.risk.msm(
      t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
      shape12.in, scale12.in, 
      shape13.in, scale13.in, 
      shape24.in, scale24.in, 
      shape34.in, scale34.in,
      beta12.cont.in, beta12.cont2.in, 
      beta13.cont.in, beta13.cont2.in, 
      beta24.cont.in, beta24.cont2.in, 
      beta34.cont.in, beta34.cont2.in)["risk.marg.true.A"]*calc.true.risk.msm(
        t.eval = t.eval.in, x1.eval = z[1], x2.eval = z[2],
        shape12.in, scale12.in, 
        shape13.in, scale13.in, 
        shape24.in, scale24.in, 
        shape34.in, scale34.in,
        beta12.cont.in, beta12.cont2.in, 
        beta13.cont.in, beta13.cont2.in, 
        beta24.cont.in, beta24.cont2.in, 
        beta34.cont.in, beta34.cont2.in)["risk.marg.true.B"]}), ncol = ncol(x))
  }
  
  ### Create lower and upper limits
  lower.lim <- qnorm(0.00001,0,1)
  upper.lim <- qnorm(0.99999,0,1)
  
  ### True joint risk over entire cohort
  true.joint.risk <- cubature::pcubature(true.J.risk.int, 
                                         lower = c(lower.lim, lower.lim), 
                                         upper = c(upper.lim, upper.lim), 
                                         vectorInterface = TRUE)$integral
  
  ### Joint risk adjusted for X over entire cohort
  joint.adj.X.risk <- cubature::pcubature(true.J.N.ADJ.risk.int, 
                                              lower = c(lower.lim, lower.lim), 
                                              upper = c(upper.lim, upper.lim),
                                              vectorInterface = TRUE)$integral
  
  ### Marginal risk of A
  true.marg.risk.A <- cubature::pcubature(true.A.risk.int, 
                                          lower = c(lower.lim, lower.lim),
                                          upper = c(upper.lim, upper.lim),
                                          vectorInterface = TRUE)$integral
  
  ### Marginal risk of B
  true.marg.risk.B <- cubature::pcubature(true.B.risk.int, 
                                          lower = c(lower.lim, lower.lim),
                                          upper = c(upper.lim, upper.lim),
                                          vectorInterface = TRUE)$integral
  
  return(c("true.joint.risk" = true.joint.risk,
           "joint.adj.X.risk" = joint.adj.X.risk,
           "naive.joint.risk" = true.marg.risk.A*true.marg.risk.B,
           "true.marg.risk.A" = true.marg.risk.A,
           "true.marg.risk.B" = true.marg.risk.B))
}

###########################
### Censoring mechanism ###

calc.prop.cens <- function(n.in, baseline_cens.in, COV_beta_cens.in){
  
  ### Create the dataset of baseline characteristics
  x.baseline.in <- data.frame("x1" = rnorm(n.in, 0, 1), "x2" = rnorm(n.in, 0, 1))
  
  ## Generate the censoring times for each individual
  cens.times <- simsurv("weibull", lambdas = 1/baseline_cens.in[2], gammas = baseline_cens.in[1], x = x.baseline.in, 
                        betas = c("x1" = COV_beta_cens.in[1], "x2" = COV_beta_cens.in[2]))
  
  ## Return proportion of censoring events prior to 10 years
  return(sum(cens.times$eventtime < 10)/nrow(cens.times))
}


### Results in approx 10% censoring by 10 years
baseline_cens.sim <- c(1, 95)
COV_beta_cens.sim <- c(0.1, 0.1)

prop.cens.10 <- calc.prop.cens(n.in = 1000000, baseline_cens.in = baseline_cens.sim, COV_beta_cens.in = COV_beta_cens.sim)
prop.cens.10



####################################################
####################################################
### Calculate input parameters for main scenario ###
####################################################
####################################################

### Denote this scenario s1.1

### When we vary the residual dependence, we will consider the ratio between the mean joint risk adjusted for X, and the
### mean joint risk adjusted for X and residual dependence.

### In the CPRD dataset we have found mean marginal risks of approx 10% for each outcome, a mean naive unadjusted joint risk of 1%, and a mean 
### joint risk adjusted for X + residual dependence of 1.5%. This is a ratio of 1.5 between the mean joint risk adjusted for X + residual 
### dependence, and the naive unadjusted mean joint risk. We will therefore target this ratio 1.5 in this scenario. However, it's important to note,
### it's unclear from the CPRD data how much of this difference can be explained by adjusting for X alone, until we have explored the clinical example
### in detail. We are attempting to have concordence between simulation and clinical example, but they do not need to be identical. I therefore propose
### targeting a mean joint risk adjusted for X of 1.25%. This will mean:
## Naive unadjusted mean joint risk: 1%
## Mean joint risk adjusted for X: 1.25%
## Mean joint risk adjusted for X and residual dependence: 1.5%

### This means the residual correlation ratio is 1.5/1.25 = 20% increase in the joint risk induced by the residual dependence.
### And The ratio of naive unadjusted mean joint risk/mean joint risk adjusted for X and residual dependence is 1.5, as in the real data.


######################
### Clayton copula ###

t.eval.sim <- 10
baselineA.s1.1.clay <- c(1,115)
baselineB.s1.1.clay <- c(1,120)
COV_betaA.s1.1.clay <- c(0.25, 0.6)
COV_betaB.s1.1.clay <- c(0.73, 0.15)
copula.param.s1.1.clay <- 0.28
copula.s1.1.clay <- "clayton"
rotate.cop.s1.1.clay <- 0

pop.risks.s1.1.clay <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                   baselineA.in = baselineA.s1.1.clay,
                                                   baselineB.in = baselineB.s1.1.clay,
                                                   COV_betaA.in = COV_betaA.s1.1.clay,
                                                   COV_betaB.in = COV_betaB.s1.1.clay,
                                                   copula.param.in = copula.param.s1.1.clay,
                                                   copula.in = copula.s1.1.clay,
                                                   rotate.cop.in = rotate.cop.s1.1.clay)

pop.risks.s1.1.clay


#####################
### Gumbel copula ###

baselineA.s1.1.gumb <- c(1,115)
baselineB.s1.1.gumb <- c(1,120)
COV_betaA.s1.1.gumb <- c(0.25, 0.6)
COV_betaB.s1.1.gumb <- c(0.73, 0.15)
copula.param.s1.1.gumb <- 1.03
copula.s1.1.gumb <- "gumbel"
rotate.cop.s1.1.gumb <- 0

pop.risks.s1.1.gumb <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                   baselineA.in = baselineA.s1.1.gumb,
                                                   baselineB.in = baselineB.s1.1.gumb,
                                                   COV_betaA.in = COV_betaA.s1.1.gumb,
                                                   COV_betaB.in = COV_betaB.s1.1.gumb,
                                                   copula.param.in = copula.param.s1.1.gumb,
                                                   copula.in = copula.s1.1.gumb,
                                                   rotate.cop.in = rotate.cop.s1.1.gumb)

pop.risks.s1.1.gumb



####################
### Frank copula ###

baselineA.s1.1.frank <- c(1,115)
baselineB.s1.1.frank <- c(1,120)
COV_betaA.s1.1.frank <- c(0.25, 0.6)
COV_betaB.s1.1.frank <- c(0.73, 0.15)
copula.param.s1.1.frank <- 0.65
copula.s1.1.frank <- "frank"
rotate.cop.s1.1.frank <- 0

pop.risks.s1.1.frank <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                    baselineA.in = baselineA.s1.1.frank,
                                                    baselineB.in = baselineB.s1.1.frank,
                                                    COV_betaA.in = COV_betaA.s1.1.frank,
                                                    COV_betaB.in = COV_betaB.s1.1.frank,
                                                    copula.param.in = copula.param.s1.1.frank,
                                                    copula.in = copula.s1.1.frank,
                                                    rotate.cop.in = rotate.cop.s1.1.frank)

pop.risks.s1.1.frank


######################
### Normal frailty ###

baselineA.s1.1.norm <- c(1,126)
baselineB.s1.1.norm <- c(1,131)
COV_betaA.s1.1.norm <- c(0.25, 0.6)
COV_betaB.s1.1.norm <- c(0.73, 0.15)
frail.var.s1.1.norm <- 0.5
frail.dist.s1.1.norm <- "normal"

pop.risks.s1.1.norm <- est.pop.risks.frailty.func.vec(t.eval.in = t.eval.sim,
                                                    baselineA.in = baselineA.s1.1.norm,
                                                    baselineB.in = baselineB.s1.1.norm,
                                                    COV_betaA.in = COV_betaA.s1.1.norm,
                                                    COV_betaB.in = COV_betaB.s1.1.norm,
                                                    frail.dist.in = frail.dist.s1.1.norm,
                                                    frail.var.in = frail.var.s1.1.norm)

pop.risks.s1.1.norm


#####################
### Gamma frailty ###

baselineA.s1.1.gamma <- c(1,112)
baselineB.s1.1.gamma <- c(1,117)
COV_betaA.s1.1.gamma <- c(0.28, 0.6)
COV_betaB.s1.1.gamma <- c(0.73, 0.15)
frail.var.s1.1.gamma <- 0.27
frail.dist.s1.1.gamma <- "gamma"

pop.risks.s1.1.gamma <- est.pop.risks.frailty.func.vec(t.eval.in = t.eval.sim,
                                                     baselineA.in = baselineA.s1.1.gamma,
                                                     baselineB.in = baselineB.s1.1.gamma,
                                                     COV_betaA.in = COV_betaA.s1.1.gamma,
                                                     COV_betaB.in = COV_betaB.s1.1.gamma,
                                                     frail.dist.in = frail.dist.s1.1.gamma,
                                                     frail.var.in = frail.var.s1.1.gamma)

pop.risks.s1.1.gamma

###########
### MSM ###

shape12.s1.1.msm <- 1
scale12.s1.1.msm <- 118
shape13.s1.1.msm <- 1
scale13.s1.1.msm <- 123
shape24.s1.1.msm <- 1
scale24.s1.1.msm <- scale13.s1.1.msm*0.77
shape34.s1.1.msm <- 1
scale34.s1.1.msm <- scale12.s1.1.msm*0.77
beta12.cont.s1.1.msm <- 0.25
beta12.cont2.s1.1.msm <- 0.6
beta13.cont.s1.1.msm <- 0.73
beta13.cont2.s1.1.msm <- 0.15
beta24.cont.s1.1.msm <- 0.73
beta24.cont2.s1.1.msm <- 0.15
beta34.cont.s1.1.msm <- 0.25
beta34.cont2.s1.1.msm <- 0.6

pop.risks.s1.1.msm <- est.pop.risks.msm.func.vec(t.eval.in = t.eval.sim,                                   
                                               shape12.in = shape12.s1.1.msm, scale12.in = scale12.s1.1.msm,
                                               shape13.in = shape13.s1.1.msm, scale13.in = scale13.s1.1.msm,
                                               shape24.in = shape24.s1.1.msm, scale24.in = scale24.s1.1.msm,
                                               shape34.in = shape34.s1.1.msm, scale34.in = scale34.s1.1.msm,
                                               beta12.cont.in = beta12.cont.s1.1.msm, beta12.cont2.in = beta12.cont2.s1.1.msm, 
                                               beta13.cont.in = beta13.cont.s1.1.msm, beta13.cont2.in = beta13.cont2.s1.1.msm, 
                                               beta24.cont.in = beta24.cont.s1.1.msm, beta24.cont2.in = beta24.cont2.s1.1.msm,
                                               beta34.cont.in = beta34.cont.s1.1.msm, beta34.cont2.in = beta34.cont2.s1.1.msm)


######################
### No correlation ###

t.eval.sim <- 10
baselineA.s1.nocorr <- c(1,115)
baselineB.s1.nocorr <- c(1,120)
COV_betaA.s1.nocorr <- c(0.25, 0.6)
COV_betaB.s1.nocorr <- c(0.73, 0.15)
copula.param.s1.nocorr <- 0.000001
copula.s1.nocorr <- "clayton"
rotate.cop.s1.nocorr <- 0

pop.risks.s1.nocorr <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                     baselineA.in = baselineA.s1.nocorr,
                                                     baselineB.in = baselineB.s1.nocorr,
                                                     COV_betaA.in = COV_betaA.s1.nocorr,
                                                     COV_betaB.in = COV_betaB.s1.nocorr,
                                                     copula.param.in = copula.param.s1.nocorr,
                                                     copula.in = copula.s1.nocorr,
                                                     rotate.cop.in = rotate.cop.s1.nocorr)



###################################################
###################################################
### Calculate input parameters for scenario 1.2 ###
###################################################
###################################################

### This is the same marginal event rates as main scenario, but increasing the level of residual dependence


######################
### Clayton copula ###

t.eval.sim <- 10
baselineA.s1.2.clay <- c(1,115)
baselineB.s1.2.clay <- c(1,120)
COV_betaA.s1.2.clay <- c(0.25, 0.6)
COV_betaB.s1.2.clay <- c(0.73, 0.15)
copula.param.s1.2.clay <- 0.69
copula.s1.2.clay <- "clayton"
rotate.cop.s1.2.clay <- 0

pop.risks.s1.2.clay <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                     baselineA.in = baselineA.s1.2.clay,
                                                     baselineB.in = baselineB.s1.2.clay,
                                                     COV_betaA.in = COV_betaA.s1.2.clay,
                                                     COV_betaB.in = COV_betaB.s1.2.clay,
                                                     copula.param.in = copula.param.s1.2.clay,
                                                     copula.in = copula.s1.2.clay,
                                                     rotate.cop.in = rotate.cop.s1.2.clay)

pop.risks.s1.2.clay


#####################
### Gumbel copula ###

baselineA.s1.2.gumb <- c(1,115)
baselineB.s1.2.gumb <- c(1,120)
COV_betaA.s1.2.gumb <- c(0.25, 0.6)
COV_betaB.s1.2.gumb <- c(0.73, 0.15)
copula.param.s1.2.gumb <- 1.065
copula.s1.2.gumb <- "gumbel"
rotate.cop.s1.2.gumb <- 0

pop.risks.s1.2.gumb <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                     baselineA.in = baselineA.s1.2.gumb,
                                                     baselineB.in = baselineB.s1.2.gumb,
                                                     COV_betaA.in = COV_betaA.s1.2.gumb,
                                                     COV_betaB.in = COV_betaB.s1.2.gumb,
                                                     copula.param.in = copula.param.s1.2.gumb,
                                                     copula.in = copula.s1.2.gumb,
                                                     rotate.cop.in = rotate.cop.s1.2.gumb)

pop.risks.s1.2.gumb



####################
### Frank copula ###

baselineA.s1.2.frank <- c(1,115)
baselineB.s1.2.frank <- c(1,120)
COV_betaA.s1.2.frank <- c(0.25, 0.6)
COV_betaB.s1.2.frank <- c(0.73, 0.15)
copula.param.s1.2.frank <- 1.37
copula.s1.2.frank <- "frank"
rotate.cop.s1.2.frank <- 0

pop.risks.s1.2.frank <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                      baselineA.in = baselineA.s1.2.frank,
                                                      baselineB.in = baselineB.s1.2.frank,
                                                      COV_betaA.in = COV_betaA.s1.2.frank,
                                                      COV_betaB.in = COV_betaB.s1.2.frank,
                                                      copula.param.in = copula.param.s1.2.frank,
                                                      copula.in = copula.s1.2.frank,
                                                      rotate.cop.in = rotate.cop.s1.2.frank)

pop.risks.s1.2.frank


######################
### Normal frailty ###

baselineA.s1.2.norm <- c(1,144)
baselineB.s1.2.norm <- c(1,149)
COV_betaA.s1.2.norm <- c(0.25, 0.6)
COV_betaB.s1.2.norm <- c(0.73, 0.15)
frail.var.s1.2.norm <- 0.77
frail.dist.s1.2.norm <- "normal"

pop.risks.s1.2.norm <- est.pop.risks.frailty.func.vec(t.eval.in = t.eval.sim,
                                                      baselineA.in = baselineA.s1.2.norm,
                                                      baselineB.in = baselineB.s1.2.norm,
                                                      COV_betaA.in = COV_betaA.s1.2.norm,
                                                      COV_betaB.in = COV_betaB.s1.2.norm,
                                                      frail.dist.in = frail.dist.s1.2.norm,
                                                      frail.var.in = frail.var.s1.2.norm)

pop.risks.s1.2.norm


#####################
### Gamma frailty ###

baselineA.s1.2.gamma <- c(1,107)
baselineB.s1.2.gamma <- c(1,112)
COV_betaA.s1.2.gamma <- c(0.25, 0.6)
COV_betaB.s1.2.gamma <- c(0.73, 0.15)
frail.var.s1.2.gamma <- 0.69
frail.dist.s1.2.gamma <- "gamma"

pop.risks.s1.2.gamma <- est.pop.risks.frailty.func.vec(t.eval.in = t.eval.sim,
                                                       baselineA.in = baselineA.s1.2.gamma,
                                                       baselineB.in = baselineB.s1.2.gamma,
                                                       COV_betaA.in = COV_betaA.s1.2.gamma,
                                                       COV_betaB.in = COV_betaB.s1.2.gamma,
                                                       frail.dist.in = frail.dist.s1.2.gamma,
                                                       frail.var.in = frail.var.s1.2.gamma)
pop.risks.s1.2.gamma

###########
### MSM ###

shape12.s1.2.msm <- 1
scale12.s1.2.msm <- 122
shape13.s1.2.msm <- 1
scale13.s1.2.msm <- 126
shape24.s1.2.msm <- 1
scale24.s1.2.msm <- scale13.s1.2.msm*0.58
shape34.s1.2.msm <- 1
scale34.s1.2.msm <- scale12.s1.2.msm*0.58
beta12.cont.s1.2.msm <- 0.25
beta12.cont2.s1.2.msm <- 0.6
beta13.cont.s1.2.msm <- 0.73
beta13.cont2.s1.2.msm <- 0.15
beta24.cont.s1.2.msm <- 0.73
beta24.cont2.s1.2.msm <- 0.15
beta34.cont.s1.2.msm <- 0.25
beta34.cont2.s1.2.msm <- 0.6

pop.risks.s1.2.msm <- est.pop.risks.msm.func.vec(t.eval.in = t.eval.sim,
                                                 shape12.in = shape12.s1.2.msm, scale12.in = scale12.s1.2.msm,
                                                 shape13.in = shape13.s1.2.msm, scale13.in = scale13.s1.2.msm,
                                                 shape24.in = shape24.s1.2.msm, scale24.in = scale24.s1.2.msm,
                                                 shape34.in = shape34.s1.2.msm, scale34.in = scale34.s1.2.msm,
                                                 beta12.cont.in = beta12.cont.s1.2.msm, beta12.cont2.in = beta12.cont2.s1.2.msm,
                                                 beta13.cont.in = beta13.cont.s1.2.msm, beta13.cont2.in = beta13.cont2.s1.2.msm,
                                                 beta24.cont.in = beta24.cont.s1.2.msm, beta24.cont2.in = beta24.cont2.s1.2.msm,
                                                 beta34.cont.in = beta34.cont.s1.2.msm, beta34.cont2.in = beta34.cont2.s1.2.msm)


####################################################
####################################################
### Calculate input parameters for scenario s2.1 ###
####################################################
####################################################

### The mean marginal risk of each outcome will be 0.3, giving a mean naive joint risk of 0.09
### The difference from mean joint risk adjusted for X will be 1.25
### The difference from mean joint risk adjusted for X, and mean joint risk adjusted for X and residual correlation, will be 1.2.
### This is the same as scenario 1.1, with increased marginal risks.

######################
### Clayton copula ###

t.eval.sim <- 10
baselineA.s2.1.clay <- c(1,35)
baselineB.s2.1.clay <- c(1,35)
COV_betaA.s2.1.clay <- c(0.3, 0.85)
COV_betaB.s2.1.clay <- c(0.9, 0.22)
copula.param.s2.1.clay <- 0.55
copula.s2.1.clay <- "clayton"
rotate.cop.s2.1.clay <- 0

pop.risks.s2.1.clay <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                     baselineA.in = baselineA.s2.1.clay,
                                                     baselineB.in = baselineB.s2.1.clay,
                                                     COV_betaA.in = COV_betaA.s2.1.clay,
                                                     COV_betaB.in = COV_betaB.s2.1.clay,
                                                     copula.param.in = copula.param.s2.1.clay,
                                                     copula.in = copula.s2.1.clay,
                                                     rotate.cop.in = rotate.cop.s2.1.clay)

pop.risks.s2.1.clay


#####################
### Gumbel copula ###

baselineA.s2.1.gumb <- c(1,35)
baselineB.s2.1.gumb <- c(1,35)
COV_betaA.s2.1.gumb <- c(0.3, 0.85)
COV_betaB.s2.1.gumb <- c(0.9, 0.22)
copula.param.s2.1.gumb <- 1.16
copula.s2.1.gumb <- "gumbel"
rotate.cop.s2.1.gumb <- 0

pop.risks.s2.1.gumb <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                     baselineA.in = baselineA.s2.1.gumb,
                                                     baselineB.in = baselineB.s2.1.gumb,
                                                     COV_betaA.in = COV_betaA.s2.1.gumb,
                                                     COV_betaB.in = COV_betaB.s2.1.gumb,
                                                     copula.param.in = copula.param.s2.1.gumb,
                                                     copula.in = copula.s2.1.gumb,
                                                     rotate.cop.in = rotate.cop.s2.1.gumb)

pop.risks.s2.1.gumb



####################
### Frank copula ###

baselineA.s2.1.frank <- c(1,35)
baselineB.s2.1.frank <- c(1,35)
COV_betaA.s2.1.frank <- c(0.3, 0.85)
COV_betaB.s2.1.frank <- c(0.9, 0.22)
copula.param.s2.1.frank <- 1.52
copula.s2.1.frank <- "frank"
rotate.cop.s2.1.frank <- 0

pop.risks.s2.1.frank <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                      baselineA.in = baselineA.s2.1.frank,
                                                      baselineB.in = baselineB.s2.1.frank,
                                                      COV_betaA.in = COV_betaA.s2.1.frank,
                                                      COV_betaB.in = COV_betaB.s2.1.frank,
                                                      copula.param.in = copula.param.s2.1.frank,
                                                      copula.in = copula.s2.1.frank,
                                                      rotate.cop.in = rotate.cop.s2.1.frank)

pop.risks.s2.1.frank


######################
### Normal frailty ###

baselineA.s2.1.norm <- c(1,39.5)
baselineB.s2.1.norm <- c(1,39.5)
COV_betaA.s2.1.norm <- c(0.35, 0.85)
COV_betaB.s2.1.norm <- c(0.9, 0.3)
frail.var.s2.1.norm <- 0.75
frail.dist.s2.1.norm <- "normal"

pop.risks.s2.1.norm <- est.pop.risks.frailty.func.vec(t.eval.in = t.eval.sim,
                                                      baselineA.in = baselineA.s2.1.norm,
                                                      baselineB.in = baselineB.s2.1.norm,
                                                      COV_betaA.in = COV_betaA.s2.1.norm,
                                                      COV_betaB.in = COV_betaB.s2.1.norm,
                                                      frail.dist.in = frail.dist.s2.1.norm,
                                                      frail.var.in = frail.var.s2.1.norm)

pop.risks.s2.1.norm


#####################
### Gamma frailty ###

baselineA.s2.1.gamma <- c(1,31)
baselineB.s2.1.gamma <- c(1,31)
COV_betaA.s2.1.gamma <- c(0.35, 0.85)
COV_betaB.s2.1.gamma <- c(0.9, 0.3)
frail.var.s2.1.gamma <- 0.45
frail.dist.s2.1.gamma <- "gamma"

pop.risks.s2.1.gamma <- est.pop.risks.frailty.func.vec(t.eval.in = t.eval.sim,
                                                       baselineA.in = baselineA.s2.1.gamma,
                                                       baselineB.in = baselineB.s2.1.gamma,
                                                       COV_betaA.in = COV_betaA.s2.1.gamma,
                                                       COV_betaB.in = COV_betaB.s2.1.gamma,
                                                       frail.dist.in = frail.dist.s2.1.gamma,
                                                       frail.var.in = frail.var.s2.1.gamma)
pop.risks.s2.1.gamma

###########
### MSM ###
t.eval.sim <- 10
shape12.s2.1.msm <- 1
scale12.s2.1.msm <- 38
shape13.s2.1.msm <- 1
scale13.s2.1.msm <- 38
shape24.s2.1.msm <- 1
scale24.s2.1.msm <- scale13.s2.1.msm*0.61
shape34.s2.1.msm <- 1
scale34.s2.1.msm <- scale12.s2.1.msm*0.61
beta12.cont.s2.1.msm <- 0.22
beta12.cont2.s2.1.msm <- 0.80
beta13.cont.s2.1.msm <- 0.9
beta13.cont2.s2.1.msm <- 0.2
beta24.cont.s2.1.msm <- 0.9
beta24.cont2.s2.1.msm <- 0.2
beta34.cont.s2.1.msm <- 0.22
beta34.cont2.s2.1.msm <- 0.80

pop.risks.s2.1.msm <- est.pop.risks.msm.func.vec(t.eval.in = t.eval.sim,                                   
                                                 shape12.in = shape12.s2.1.msm, scale12.in = scale12.s2.1.msm,
                                                 shape13.in = shape13.s2.1.msm, scale13.in = scale13.s2.1.msm,
                                                 shape24.in = shape24.s2.1.msm, scale24.in = scale24.s2.1.msm,
                                                 shape34.in = shape34.s2.1.msm, scale34.in = scale34.s2.1.msm,
                                                 beta12.cont.in = beta12.cont.s2.1.msm, beta12.cont2.in = beta12.cont2.s2.1.msm, 
                                                 beta13.cont.in = beta13.cont.s2.1.msm, beta13.cont2.in = beta13.cont2.s2.1.msm, 
                                                 beta24.cont.in = beta24.cont.s2.1.msm, beta24.cont2.in = beta24.cont2.s2.1.msm,
                                                 beta34.cont.in = beta34.cont.s2.1.msm, beta34.cont2.in = beta34.cont2.s2.1.msm)
pop.risks.s2.1.msm


######################
### No correlation ###

t.eval.sim <- 10
baselineA.s2.nocorr <- c(1,35)
baselineB.s2.nocorr <- c(1,35)
COV_betaA.s2.nocorr <- c(0.3, 0.85)
COV_betaB.s2.nocorr <- c(0.9, 0.22)
copula.param.s2.nocorr <- 0.000001
copula.s2.nocorr <- "clayton"
rotate.cop.s2.nocorr <- 0

pop.risks.s2.nocorr <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                     baselineA.in = baselineA.s2.nocorr,
                                                     baselineB.in = baselineB.s2.nocorr,
                                                     COV_betaA.in = COV_betaA.s2.nocorr,
                                                     COV_betaB.in = COV_betaB.s2.nocorr,
                                                     copula.param.in = copula.param.s2.nocorr,
                                                     copula.in = copula.s2.nocorr,
                                                     rotate.cop.in = rotate.cop.s1.nocorr)


####################################################
####################################################
### Calculate input parameters for scenario s2.2 ###
####################################################
####################################################

### The mean marginal risk of each outcome will be 0.3, giving a mean naive joint risk of 0.09
### The difference from mean joint risk adjusted for X will be 1.25
### The difference from mean joint risk adjusted for X, and mean joint risk adjusted for X and residual correlation, will be 1.5.
### This is the same as scenario 1.2, with increased marginal risks.

######################
### Clayton copula ###

t.eval.sim <- 10
baselineA.s2.2.clay <- c(1,35)
baselineB.s2.2.clay <- c(1,35)
COV_betaA.s2.2.clay <- c(0.3, 0.85)
COV_betaB.s2.2.clay <- c(0.9, 0.22)
copula.param.s2.2.clay <- 2
copula.s2.2.clay <- "clayton"
rotate.cop.s2.2.clay <- 0

pop.risks.s2.2.clay <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                     baselineA.in = baselineA.s2.2.clay,
                                                     baselineB.in = baselineB.s2.2.clay,
                                                     COV_betaA.in = COV_betaA.s2.2.clay,
                                                     COV_betaB.in = COV_betaB.s2.2.clay,
                                                     copula.param.in = copula.param.s2.2.clay,
                                                     copula.in = copula.s2.2.clay,
                                                     rotate.cop.in = rotate.cop.s2.2.clay)

pop.risks.s2.2.clay


#####################
### Gumbel copula ###

baselineA.s2.2.gumb <- c(1,35)
baselineB.s2.2.gumb <- c(1,35)
COV_betaA.s2.2.gumb <- c(0.3, 0.85)
COV_betaB.s2.2.gumb <- c(0.9, 0.22)
copula.param.s2.2.gumb <- 1.54
copula.s2.2.gumb <- "gumbel"
rotate.cop.s2.2.gumb <- 0

pop.risks.s2.2.gumb <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                     baselineA.in = baselineA.s2.2.gumb,
                                                     baselineB.in = baselineB.s2.2.gumb,
                                                     COV_betaA.in = COV_betaA.s2.2.gumb,
                                                     COV_betaB.in = COV_betaB.s2.2.gumb,
                                                     copula.param.in = copula.param.s2.2.gumb,
                                                     copula.in = copula.s2.2.gumb,
                                                     rotate.cop.in = rotate.cop.s2.2.gumb)

pop.risks.s2.2.gumb



####################
### Frank copula ###

baselineA.s2.2.frank <- c(1,35)
baselineB.s2.2.frank <- c(1,35)
COV_betaA.s2.2.frank <- c(0.3, 0.85)
COV_betaB.s2.2.frank <- c(0.9, 0.22)
copula.param.s2.2.frank <- 4.2
copula.s2.2.frank <- "frank"
rotate.cop.s2.2.frank <- 0

pop.risks.s2.2.frank <- est.pop.risks.copula.func.vec(t.eval.in = t.eval.sim,
                                                      baselineA.in = baselineA.s2.2.frank,
                                                      baselineB.in = baselineB.s2.2.frank,
                                                      COV_betaA.in = COV_betaA.s2.2.frank,
                                                      COV_betaB.in = COV_betaB.s2.2.frank,
                                                      copula.param.in = copula.param.s2.2.frank,
                                                      copula.in = copula.s2.2.frank,
                                                      rotate.cop.in = rotate.cop.s2.2.frank)

pop.risks.s2.2.frank


######################
### Normal frailty ###

baselineA.s2.2.norm <- c(1,50.5)
baselineB.s2.2.norm <- c(1,51.5)
COV_betaA.s2.2.norm <- c(0.4, 1)
COV_betaB.s2.2.norm <- c(1.1, 0.35)
frail.var.s2.2.norm <- 1.38
frail.dist.s2.2.norm <- "normal"

pop.risks.s2.2.norm <- est.pop.risks.frailty.func.vec(t.eval.in = t.eval.sim,
                                                      baselineA.in = baselineA.s2.2.norm,
                                                      baselineB.in = baselineB.s2.2.norm,
                                                      COV_betaA.in = COV_betaA.s2.2.norm,
                                                      COV_betaB.in = COV_betaB.s2.2.norm,
                                                      frail.dist.in = frail.dist.s2.2.norm,
                                                      frail.var.in = frail.var.s2.2.norm)

pop.risks.s2.2.norm


#####################
### Gamma frailty ###

baselineA.s2.2.gamma <- c(1,25)
baselineB.s2.2.gamma <- c(1,25)
COV_betaA.s2.2.gamma <- c(0.4, 1)
COV_betaB.s2.2.gamma <- c(1.1, 0.35)
frail.var.s2.2.gamma <- 1.45
frail.dist.s2.2.gamma <- "gamma"

pop.risks.s2.2.gamma <- est.pop.risks.frailty.func.vec(t.eval.in = t.eval.sim,
                                                       baselineA.in = baselineA.s2.2.gamma,
                                                       baselineB.in = baselineB.s2.2.gamma,
                                                       COV_betaA.in = COV_betaA.s2.2.gamma,
                                                       COV_betaB.in = COV_betaB.s2.2.gamma,
                                                       frail.dist.in = frail.dist.s2.2.gamma,
                                                       frail.var.in = frail.var.s2.2.gamma)
pop.risks.s2.2.gamma


###########
### MSM ###
t.eval.sim <- 10
shape12.s2.2.msm <- 1
scale12.s2.2.msm <- 42
shape13.s2.2.msm <- 1
scale13.s2.2.msm <- 43
shape24.s2.2.msm <- 1
scale24.s2.2.msm <- scale13.s2.2.msm*0.33
shape34.s2.2.msm <- 1
scale34.s2.2.msm <- scale12.s2.2.msm*0.33
beta12.cont.s2.2.msm <- 0.22
beta12.cont2.s2.2.msm <- 0.7
beta13.cont.s2.2.msm <- 0.9
beta13.cont2.s2.2.msm <- 0.15
beta24.cont.s2.2.msm <- 0.9
beta24.cont2.s2.2.msm <- 0.15
beta34.cont.s2.2.msm <- 0.22
beta34.cont2.s2.2.msm <- 0.7

pop.risks.s2.2.msm <- est.pop.risks.msm.func.vec(t.eval.in = t.eval.sim,                                   
                                                 shape12.in = shape12.s2.2.msm, scale12.in = scale12.s2.2.msm,
                                                 shape13.in = shape13.s2.2.msm, scale13.in = scale13.s2.2.msm,
                                                 shape24.in = shape24.s2.2.msm, scale24.in = scale24.s2.2.msm,
                                                 shape34.in = shape34.s2.2.msm, scale34.in = scale34.s2.2.msm,
                                                 beta12.cont.in = beta12.cont.s2.2.msm, beta12.cont2.in = beta12.cont2.s2.2.msm, 
                                                 beta13.cont.in = beta13.cont.s2.2.msm, beta13.cont2.in = beta13.cont2.s2.2.msm, 
                                                 beta24.cont.in = beta24.cont.s2.2.msm, beta24.cont2.in = beta24.cont2.s2.2.msm,
                                                 beta34.cont.in = beta34.cont.s2.2.msm, beta34.cont2.in = beta34.cont2.s2.2.msm)
pop.risks.s2.2.msm


pop.risks.s1.nocorr

pop.risks.s1.1.clay
pop.risks.s1.1.gumb
pop.risks.s1.1.frank
pop.risks.s1.1.norm
pop.risks.s1.1.gamma
pop.risks.s1.1.msm

pop.risks.s1.2.clay
pop.risks.s1.2.gumb
pop.risks.s1.2.frank
pop.risks.s1.2.norm
pop.risks.s1.2.gamma
pop.risks.s1.2.msm

pop.risks.s2.nocorr

pop.risks.s2.1.clay
pop.risks.s2.1.gumb
pop.risks.s2.1.frank
pop.risks.s2.1.norm
pop.risks.s2.1.gamma
pop.risks.s2.1.msm

pop.risks.s2.2.clay
pop.risks.s2.2.gumb
pop.risks.s2.2.frank
pop.risks.s2.2.norm
pop.risks.s2.2.gamma
pop.risks.s2.2.msm


############################################################################################
############################################################################################
### Finally we want to calculate the range of predicted risks to assess calibration over ###
############################################################################################
############################################################################################

### Will calculate the risk at which only 1% of individuals will have a risk higher than
### And the risk at which only 1% of individuals will have a risk lower than

### Calc the 10th and 90th percentile for x1 and x2
p10 <- qnorm(0.1, 0, 1)
p90 <- qnorm(0.9, 0, 1)

### Then calculate the true risk at the point when x1 = x2 = p10, and x1 = x2 = p90

### Scenario s1.1
print("s1.1")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p10, x2.eval = p10, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.1.clay, rotate.cop = rotate.cop.s1.1.clay, 
                      baselineA = baselineA.s1.1.clay, COV_betaA = COV_betaA.s1.1.clay, 
                      baselineB = baselineB.s1.1.clay, COV_betaB = COV_betaB.s1.1.clay)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p90, x2.eval = p90, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.1.clay, rotate.cop = rotate.cop.s1.1.clay, 
                      baselineA = baselineA.s1.1.clay, COV_betaA = COV_betaA.s1.1.clay, 
                      baselineB = baselineB.s1.1.clay, COV_betaB = COV_betaB.s1.1.clay)


### Scenario s1.2
print("s1.2")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p10, x2.eval = p10, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.2.clay, rotate.cop = rotate.cop.s1.2.clay, 
                      baselineA = baselineA.s1.2.clay, COV_betaA = COV_betaA.s1.2.clay, 
                      baselineB = baselineB.s1.2.clay, COV_betaB = COV_betaB.s1.2.clay)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p90, x2.eval = p90, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.2.clay, rotate.cop = rotate.cop.s1.2.clay, 
                      baselineA = baselineA.s1.2.clay, COV_betaA = COV_betaA.s1.2.clay, 
                      baselineB = baselineB.s1.2.clay, COV_betaB = COV_betaB.s1.2.clay)

### Scenario s2.1
print("s2.1")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p10, x2.eval = p10, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.1.clay, rotate.cop = rotate.cop.s2.1.clay, 
                      baselineA = baselineA.s2.1.clay, COV_betaA = COV_betaA.s2.1.clay, 
                      baselineB = baselineB.s2.1.clay, COV_betaB = COV_betaB.s2.1.clay)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p90, x2.eval = p90, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.1.clay, rotate.cop = rotate.cop.s2.1.clay, 
                      baselineA = baselineA.s2.1.clay, COV_betaA = COV_betaA.s2.1.clay, 
                      baselineB = baselineB.s2.1.clay, COV_betaB = COV_betaB.s2.1.clay)

### Scenario s2.2
print("s2.2")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p10, x2.eval = p10, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.2.clay, rotate.cop = rotate.cop.s2.2.clay, 
                      baselineA = baselineA.s2.2.clay, COV_betaA = COV_betaA.s2.2.clay, 
                      baselineB = baselineB.s2.2.clay, COV_betaB = COV_betaB.s2.2.clay)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p90, x2.eval = p90, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.2.clay, rotate.cop = rotate.cop.s2.2.clay, 
                      baselineA = baselineA.s2.2.clay, COV_betaA = COV_betaA.s2.2.clay, 
                      baselineB = baselineB.s2.2.clay, COV_betaB = COV_betaB.s2.2.clay)


### Will calculate the risk at which only 5% of individuals will have a risk higher than
### And the risk at which only 5% of individuals will have a risk lower than

### Calc the 22.3th and 77.7th percentile for x1 and x2
p223 <- qnorm(0.223, 0, 1)
p777 <- qnorm(0.777, 0, 1)

### Then calculate the true risk at the point when x1 = x2 = p10, and x1 = x2 = p90

### Scenario s1.1
print("s1.1")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p223, x2.eval = p223, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.1.clay, rotate.cop = rotate.cop.s1.1.clay, 
                      baselineA = baselineA.s1.1.clay, COV_betaA = COV_betaA.s1.1.clay, 
                      baselineB = baselineB.s1.1.clay, COV_betaB = COV_betaB.s1.1.clay)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p777, x2.eval = p777, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.1.clay, rotate.cop = rotate.cop.s1.1.clay, 
                      baselineA = baselineA.s1.1.clay, COV_betaA = COV_betaA.s1.1.clay, 
                      baselineB = baselineB.s1.1.clay, COV_betaB = COV_betaB.s1.1.clay)

### Scenario s1.2
print("s1.2")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p223, x2.eval = p223, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.2.clay, rotate.cop = rotate.cop.s1.2.clay, 
                      baselineA = baselineA.s1.2.clay, COV_betaA = COV_betaA.s1.2.clay, 
                      baselineB = baselineB.s1.2.clay, COV_betaB = COV_betaB.s1.2.clay)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p777, x2.eval = p777, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.2.clay, rotate.cop = rotate.cop.s1.2.clay, 
                      baselineA = baselineA.s1.2.clay, COV_betaA = COV_betaA.s1.2.clay, 
                      baselineB = baselineB.s1.2.clay, COV_betaB = COV_betaB.s1.2.clay)

### Scenario s2.1
print("s2.1")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p223, x2.eval = p223, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.1.clay, rotate.cop = rotate.cop.s2.1.clay, 
                      baselineA = baselineA.s2.1.clay, COV_betaA = COV_betaA.s2.1.clay, 
                      baselineB = baselineB.s2.1.clay, COV_betaB = COV_betaB.s2.1.clay)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p777, x2.eval = p777, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.1.clay, rotate.cop = rotate.cop.s2.1.clay, 
                      baselineA = baselineA.s2.1.clay, COV_betaA = COV_betaA.s2.1.clay, 
                      baselineB = baselineB.s2.1.clay, COV_betaB = COV_betaB.s2.1.clay)

### Scenario s2.2
print("s2.2")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p223, x2.eval = p223, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.2.clay, rotate.cop = rotate.cop.s2.2.clay, 
                      baselineA = baselineA.s2.2.clay, COV_betaA = COV_betaA.s2.2.clay, 
                      baselineB = baselineB.s2.2.clay, COV_betaB = COV_betaB.s2.2.clay)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p777, x2.eval = p777, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.2.clay, rotate.cop = rotate.cop.s2.2.clay, 
                      baselineA = baselineA.s2.2.clay, COV_betaA = COV_betaA.s2.2.clay, 
                      baselineB = baselineB.s2.2.clay, COV_betaB = COV_betaB.s2.2.clay)



save.image("data/sim_calculate_input_parameters_2cont.RData")




### Scenario s1
print("s1")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p10, x2.eval = p10, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.nocorr, rotate.cop = rotate.cop.s1.nocorr, 
                      baselineA = baselineA.s1.nocorr, COV_betaA = COV_betaA.s1.nocorr, 
                      baselineB = baselineB.s1.nocorr, COV_betaB = COV_betaB.s1.nocorr)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p90, x2.eval = p90, 
                      copula = "clayton", 
                      copula.param = copula.param.s1.nocorr, rotate.cop = rotate.cop.s1.nocorr, 
                      baselineA = baselineA.s1.nocorr, COV_betaA = COV_betaA.s1.nocorr, 
                      baselineB = baselineB.s1.nocorr, COV_betaB = COV_betaB.s1.nocorr)


### Scenario s2.1
print("s2")
calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p10, x2.eval = p10, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.nocorr, rotate.cop = rotate.cop.s2.nocorr, 
                      baselineA = baselineA.s2.nocorr, COV_betaA = COV_betaA.s2.nocorr, 
                      baselineB = baselineB.s2.nocorr, COV_betaB = COV_betaB.s2.nocorr)

calc.true.risk.copula(t.eval = t.eval.sim, x1.eval = p90, x2.eval = p90, 
                      copula = "clayton", 
                      copula.param = copula.param.s2.nocorr, rotate.cop = rotate.cop.s2.nocorr, 
                      baselineA = baselineA.s2.nocorr, COV_betaA = COV_betaA.s2.nocorr, 
                      baselineB = baselineB.s2.nocorr, COV_betaB = COV_betaB.s2.nocorr)


# calc.true.risk.msm(
#   t.eval = t.eval.sim, x1.eval = p10, x2.eval = p10,
#   shape12.s1.1.msm, scale12.s1.1.msm, 
#   shape13.s1.1.msm, scale13.s1.1.msm, 
#   shape24.s1.1.msm, scale24.s1.1.msm, 
#   shape34.s1.1.msm, scale34.s1.1.msm,
#   beta12.cont.s1.1.msm, beta12.cont2.s1.1.msm, 
#   beta13.cont.s1.1.msm, beta13.cont2.s1.1.msm, 
#   beta24.cont.s1.1.msm, beta24.cont2.s1.1.msm, 
#   beta34.cont.s1.1.msm, beta34.cont2.s1.1.msm)
# 
# calc.true.risk.msm(
#   t.eval = t.eval.sim, x1.eval = p90, x2.eval = p90,
#   shape12.s1.1.msm, scale12.s1.1.msm, 
#   shape13.s1.1.msm, scale13.s1.1.msm, 
#   shape24.s1.1.msm, scale24.s1.1.msm, 
#   shape34.s1.1.msm, scale34.s1.1.msm,
#   beta12.cont.s1.1.msm, beta12.cont2.s1.1.msm, 
#   beta13.cont.s1.1.msm, beta13.cont2.s1.1.msm, 
#   beta24.cont.s1.1.msm, beta24.cont2.s1.1.msm, 
#   beta34.cont.s1.1.msm, beta34.cont2.s1.1.msm)
