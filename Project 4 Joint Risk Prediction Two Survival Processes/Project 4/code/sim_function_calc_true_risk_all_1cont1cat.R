#####################################################################################################
### This program will create functions to calculate the true risk for each  data generating mechanism
#####################################################################################################

### Clear workspace
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Load packages
# library(survival)
# library(gems)
# library(mstate)
# library(dplyr)
# library(tidyr)
# library(simsurv)
# library(coxme)
# library(frailtypack)
# library(frailtySurv)
# library(frailtyEM)
# library(CopulaCenR)
# library(copula)


##################
##################
### MSM ###
##################
##################

calc.true.risk.msm <- function(t.eval, x1.eval, x2.eval, #values of time, x1 and x2 to evaluate true risk at
                               shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                               shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                               shape24, scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                               shape34, scale34, #shape and scale for weibull baseline hazard for transition 3 -> 4
                               beta12.cont, beta12.cat, #covariate effects for transiion 12
                               beta13.cont, beta13.cat, #covariate effects for transiion 13
                               beta24.cont, beta24.cat, #covariate effects for transiion 24
                               beta34.cont, beta34.cat #covariate effects for transiion 34
                               ){
  
  
  # t.eval <- 1
  # x1.eval <- 0.09
  # x2.eval <- 0
  #   
  # shape12 <- 1
  # scale12 <- 1.5
  # shape13 <- 1
  # scale13 <- 1
  # shape24 <- 1
  # scale24 <- 0.8
  # shape34 <- 1
  # scale34 <- 1.3
  # 
  # beta12.cont <- 0.25
  # beta12.cat <- 0.5
  # 
  # beta13.cont <- 0.5
  # beta13.cat <- 0.5
  # 
  # beta24.cont <- 0.33
  # beta24.cat <- 0.25
  # 
  # beta34.cont <- 0.25
  # beta34.cat <- 0.25
  
  ## Hazard rate for each transition
  haz12 <- function(t.in, x1, x2){
    return(exp(beta12.cont*x1 + beta12.cat*x2)*(shape12/scale12)*((t.in)/scale12)^(shape12 - 1))
  }
  
  haz13 <- function(t.in, x1, x2){
    return(exp(beta13.cont*x1 + beta13.cat*x2)*(shape13/scale13)*((t.in)/scale13)^(shape13 - 1))
  }
  
  haz24 <- function(t.in, x1, x2){
    return(exp(beta24.cont*x1 + beta24.cat*x2)*(shape24/scale24)*((t.in)/scale24)^(shape24 - 1))
  }
  
  haz34 <- function(t.in, x1, x2){
    return(exp(beta34.cont*x1 + beta34.cat*x2)*(shape34/scale34)*((t.in)/scale34)^(shape12 - 1))
  }
  
  
  cumhaz12 <- function(u.in, t.in, x1, x2){
    return((exp(beta12.cont*x1 + beta12.cat*x2)*(t.in/scale12)^(shape12)) - (exp(beta12.cont*x1 + beta12.cat*x2)*(u.in/scale12)^(shape12)))
  }
  
  
  ## Cumulative hazard between u.in and t.in for each transition
  cumhaz13 <- function(u.in, t.in, x1, x2){
    return((exp(beta13.cont*x1 + beta13.cat*x2)*(t.in/scale13)^(shape13)) - (exp(beta13.cont*x1 + beta13.cat*x2)*(u.in/scale13)^(shape13)))
  }
  
  cumhaz24 <- function(u.in, t.in, x1, x2){
    return((exp(beta24.cont*x1 + beta24.cat*x2)*(t.in/scale24)^(shape24)) - (exp(beta24.cont*x1 + beta24.cat*x2)*(u.in/scale24)^(shape24)))
  }
  
  cumhaz34 <- function(u.in, t.in, x1, x2){
    return((exp(beta34.cont*x1 + beta34.cat*x2)*(t.in/scale34)^(shape34)) - (exp(beta34.cont*x1 + beta34.cat*x2)*(u.in/scale34)^(shape34)))
  }
  
  
  ## Probability of still being in initial state
  S.1 <- function(t.in, x1, x2){
    exp(-(cumhaz12(0, t.in, x1, x2) + cumhaz13(0, t.in, x1, x2)))
  }
  
  
  ## Probability of having transitioned from state 2 to 4 by time t.in, after entering state 2 at time u.in
  P.24 <- function(u.in, t.in, x1, x2){
    return(1-exp(cumhaz24(0, u.in, x1, x2) - cumhaz24(0, t.in, x1, x2)))
  }
  
  ## Probability of having transitioned from state 3 to 4 by time t.in, after entering state 2 at time u.in
  P.34 <- function(u.in, t.in, x1, x2){
    return(1-exp(cumhaz34(0, u.in, x1, x2) - cumhaz34(0, t.in, x1, x2)))
  }
  
  ## Probability of having not transitioned out of state 2 by time t.in, after entering state 2 at time u.in
  ## Given there is only one transition out of state 2, this is equal to 1 - P.24
  P.22 <- function(u.in, t.in, x1, x2){
    return(exp(cumhaz24(0, u.in, x1, x2) - cumhaz24(0, t.in, x1, x2)))
  }
  
  ## Probability of having not transitioned out of state 2 by time t.in, after entering state 2 at time u.in
  ## Given there is only one transition out of state 3, this is equal to 1 - P.34
  P.33 <- function(u.in, t.in, x1, x2){
    return(exp(cumhaz34(0, u.in, x1, x2) - cumhaz34(0, t.in, x1, x2)))
  }
  
  
  ### Now to calculate the true risks using these functions
  
  ###
  ### Lets start with the joint risk
  ###
  
  ## See equation (32) from Putter et al: tutorial in biostatistics
  ## Integrate this over 0 to some t value
  
  ## Calculate probability of having A then B before time t
  calc.P.A.B <- function(u.inn, t.inn, x1.inn, x2.inn){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r){
      return(haz12(t.in = r, x1 = x1.inn, x2 = x2.inn)*S.1(t.in = r, x1 = x1.inn, x2 = x2.inn)*
               P.24(u.in = r, t.in = t.inn, x1 = x1.inn, x2 = x2.inn))
    }
    
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.inn, upper = t.inn, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Calculate probability of having B then A before time t
  calc.P.B.A <- function(u.inn, t.inn, x1.inn, x2.inn){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r){
      return(haz13(t.in = r, x1 = x1.inn, x2 = x2.inn)*S.1(t.in = r, x1 = x1.inn, x2 = x2.inn)*
               P.34(u.in = r, t.in = t.inn, x1 = x1.inn, x2 = x2.inn))
    }
    
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.inn, upper = t.inn, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## The true joint risk is the sum of these
  risk.joint.true <- (calc.P.A.B(u.inn = 0, t.inn = t.eval, x1.inn = x1.eval, x2.inn = x2.eval) + 
                        calc.P.B.A(u.inn = 0, t.inn = t.eval, x1.inn = x1.eval, x2.inn = x2.eval))
  
  
  ###
  ### Now calculate the marignal risks
  ###
  
  ## The marginal risk of outcome A, is the sum of the probability of being in state 4, and state 2
  ## We have already calculated the probability of being in statee 4 (joint risk)
  
  ## Probability of being in state 2 (developing A but not B) is calculated using equation (32) from Putter et al
  calc.P.A.A <-function(u.inn, t.inn, x1.inn, x2.inn){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r){
      return(haz12(t.in = r, x1 = x1.inn, x2 = x2.inn)*S.1(t.in = r, x1 = x1.inn, x2 = x2.inn)*
               P.22(u.in = r, t.in = t.inn, x1 = x1.inn, x2 = x2.inn))
    }
    
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.inn, upper = t.inn, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  ## Probability of being in state 3 (developing B but not A) is calculated using equation (32) from Putter et al
  calc.P.B.B <-function(u.inn, t.inn, x1.inn, x2.inn){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r){
      return(haz13(t.in = r, x1 = x1.inn, x2 = x2.inn)*S.1(t.in = r, x1 = x1.inn, x2 = x2.inn)*
               P.33(u.in = r, t.in = t.inn, x1 = x1.inn, x2 = x2.inn))
    }
    
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.inn, upper = t.inn, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  ## Caclulate the marginal risks
  # Marginal risk of A, is probability of entering A and stay in A, or the probability of being in A and B
  risk.marg.true.A <- calc.P.A.A(u.inn = 0, t.inn = t.eval, x1.inn = x1.eval, x2.inn = x2.eval) + 
    calc.P.A.B(u.inn = 0, t.inn = t.eval, x1.inn = x1.eval, x2.inn = x2.eval) + 
    calc.P.B.A(u.inn = 0, t.inn = t.eval, x1.inn = x1.eval, x2.inn = x2.eval)
  
  # Marginal risk of B, is probability of entering B and stay in B, or the probability of being in A and B
  risk.marg.true.B <- calc.P.B.B(u.inn = 0, t.inn = t.eval, x1.inn = x1.eval, x2.inn = x2.eval) +
    calc.P.A.B(u.inn = 0, t.inn = t.eval, x1.inn = x1.eval, x2.inn = x2.eval) + 
    calc.P.B.A(u.inn = 0, t.inn = t.eval, x1.inn = x1.eval, x2.inn = x2.eval)
  
  ## Calculate joint survival
  surv.joint.true <- 1 - risk.marg.true.A - risk.marg.true.B + risk.joint.true
  
  ## Create output object
  output.obj <- c("surv.marg.true.A" = 1 - risk.marg.true.A,
                  "surv.marg.true.B" = 1 - risk.marg.true.B,
                  "surv.joint.true" = surv.joint.true,
                  "risk.marg.true.A" = risk.marg.true.A,
                  "risk.marg.true.B" = risk.marg.true.B,
                  "risk.joint.true" = risk.joint.true)
  
  return(output.obj)
}



##################
##################
### COPULA ###
##################
##################

calc.true.risk.copula <- function(t.eval, x1.eval, x2.eval, # time and baseline characteristics at which we want to assess risk
                                  copula, #type of copula, options are clayton, gumbel, frank, fgm, joe, amh
                                  copula.param, rotate.cop, #association parameter, for copula, and whether it sohuld be rotated or not
                                  baselineA, COV_betaA, #outcome1: baseline1 = shape and scale of weibull distribution, COV_beta1 = covariate effects
                                  baselineB, COV_betaB) #outcome2: baseline2 = shape and scale of weibull distribution, COV_beta2 = covariate effects
                                  {
  
  ## Start by calculating marginal risks
  surv.marg.true.A <- exp(-((t.eval/baselineA[2])^baselineA[1])*exp(x1.eval*COV_betaA[1] + x2.eval*COV_betaA[2]))
  surv.marg.true.B <- exp(-((t.eval/baselineB[2])^baselineB[1])*exp(x1.eval*COV_betaB[1] + x2.eval*COV_betaB[2]))
  
  risk.marg.true.A <- 1 - exp(-((t.eval/baselineA[2])^baselineA[1])*exp(x1.eval*COV_betaA[1] + x2.eval*COV_betaA[2]))
  risk.marg.true.B <- 1 - exp(-((t.eval/baselineB[2])^baselineB[1])*exp(x1.eval*COV_betaB[1] + x2.eval*COV_betaB[2]))
  
  ## Now create a copula object with the correct input parameter
  if (copula == "fgm"){cl.true <- fgmCopula(param = copula.param, dim = 2)}
  if (copula %in% c("clayton", "gumbel", "frank", "amh", "joe")){cl.true <- archmCopula(family = tolower(copula), param = copula.param, dim = 2)}
  
  ## Rotate if required
  if (rotate.cop == 90){cl.true <- rotCopula(cop1, flip = c(TRUE, FALSE))}
  if (rotate.cop == 180){cl.true <- rotCopula(cop1, flip = c(TRUE, TRUE))}
  
  ## Calculate joint survival probability
  surv.joint.true <- pCopula(c(surv.marg.true.A,surv.marg.true.B), cl.true)
  
  ## Calculate joint risk
  risk.joint.true <- 1 - surv.marg.true.A - surv.marg.true.B + surv.joint.true 
  
  ## Create output object
  output.obj <- c("surv.marg.true.A" = surv.marg.true.A,
                  "surv.marg.true.B" = surv.marg.true.B,
                  "surv.joint.true" = surv.joint.true,
                  "risk.marg.true.A" = risk.marg.true.A,
                  "risk.marg.true.B" = risk.marg.true.B,
                  "risk.joint.true" = risk.joint.true)
  
  return(output.obj)
}



##################
##################
### FRAILTY ###
##################
##################
 

### USE PCUBATURE FOR NORMAL, HCUBATURE FOR GAMMA
calc.true.risk.frailty <- function(t.eval, x1.eval, x2.eval, # time and baseline characteristics at which we want to assess risk
                                   frail.dist, frail.var, # frailty distribution, and variance of frailty term
                                   baselineA, COV_betaA, #1) shape and scale of weibull baseline hazard, 2) covariate effects
                                   baselineB, COV_betaB, #1) shape and scale of weibull baseline hazard, 2) covariate effects
                                   int.method) # Method of integration for integrating over the frailty term
  {
  
  #   t.eval <- 2
  #   x1.eval <- 0.1
  #   x2.eval <- 0
  #   shape.A.in <- 1
  #   scale.A.in <- 1.5
  #   shape.B.in <- 1
  #   scale.B.in <- 1
  #   beta1.A.in <- 0
  #   beta2.A.in <- 0
  #   beta1.B.in <- 0
  #   beta2.B.in <- 0
  #   frail.dist = "gamma"
  #   frail.var = 3
  
  
  ## Write a seperate function for frail.type = gamma or normal, as for normal the random effect is apllied on the exp() scale 
  if (frail.dist == "gamma"){## First write functions for marginal risk given the random effect (z)
    risk.marg.A.given.Z <- function(z){1-exp(-((t.eval/baselineA[2])^baselineA[1])*exp(x1.eval*COV_betaA[1] + x2.eval*COV_betaA[2])*z)}
    risk.marg.B.given.Z <- function(z){1-exp(-((t.eval/baselineB[2])^baselineB[1])*exp(x1.eval*COV_betaB[1] + x2.eval*COV_betaB[2])*z)}
    
    ## Create functions which we will integrate over the frailty term, to get risks/survival proabilities
    # Marg risk for A
    marg.risk.A.for.int <- function(z.in){
      risk.marg.A.given.Z(z.in)*dgamma(z.in, shape = 1/frail.var, scale = frail.var)
    } 
    
    # Marg risk for B
    marg.risk.B.for.int <- function(z.in){
      risk.marg.B.given.Z(z.in)*dgamma(z.in, shape = 1/frail.var, scale = frail.var)
    } 
    
    # Joint risk
    joint.risk.for.int <- function(z.in){
      risk.marg.A.given.Z(z.in)*risk.marg.B.given.Z(z.in)*dgamma(z.in, shape = 1/frail.var, scale = frail.var)
    } 
    
    # Marg surv for A (integrate over 1 - risk)
    marg.surv.A.for.int <- function(z.in){
      (1 - risk.marg.A.given.Z(z.in))*dgamma(z.in, shape = 1/frail.var, scale = frail.var)
    }
    
    # Marg surv for B (integrate over 1 - risk)
    marg.surv.B.for.int <- function(z.in){
      (1 - risk.marg.B.given.Z(z.in))*dgamma(z.in, shape = 1/frail.var, scale = frail.var)
    } 
    
    # Joint surv
    joint.surv.for.int <- function(z.in){
      (1 - risk.marg.A.given.Z(z.in))*(1 - risk.marg.B.given.Z(z.in))*dgamma(z.in, shape = 1/frail.var, scale = frail.var)
    } 
    
    ## Estimate the survival probabilities and risks by integrating these functions over the random effect
    
    ## Interate between 0 (lower limit) and 99.99th percentiles of the random effect distribution, otherwise we get some issues
    int.upper <- qgamma(0.99999, shape = 1/frail.var, scale = frail.var)
    
    risk.marg.true.A <- cubintegrate(f = marg.risk.A.for.int, lower = 0, upper = int.upper, method = int.method)$integral
    risk.marg.true.B <- cubintegrate(f = marg.risk.B.for.int, lower = 0, upper = int.upper, method = int.method)$integral
    risk.joint.true <- cubintegrate(f = joint.risk.for.int, lower = 0, upper = int.upper, method = int.method)$integral
    
    surv.marg.true.A <- cubintegrate(f = marg.surv.A.for.int, lower = 0, upper = int.upper, method = int.method)$integral
    surv.marg.true.B <- cubintegrate(f = marg.surv.B.for.int, lower = 0, upper = int.upper, method = int.method)$integral
    surv.joint.true <- cubintegrate(f = joint.surv.for.int, lower = 0, upper = int.upper, method = int.method)$integral                      
  }
  
  if (frail.dist == "normal"){## First write functions for marginal risk given the random effect (z)
    risk.marg.A.given.Z <- function(z){1-exp(-((t.eval/baselineA[2])^baselineA[1])*exp(x1.eval*COV_betaA[1] + x2.eval*COV_betaA[2])*exp(z))}
    risk.marg.B.given.Z <- function(z){1-exp(-((t.eval/baselineB[2])^baselineB[1])*exp(x1.eval*COV_betaB[1] + x2.eval*COV_betaB[2])*exp(z))}
    
    ## Create functions which we will integrate over the frailty term, to get risks/survival proabilities
    # Marg risk for A
    marg.risk.A.for.int <- function(z.in){
      risk.marg.A.given.Z(z.in)*dnorm(z.in, 0, frail.var)
    } 
    
    # Marg risk for B
    marg.risk.B.for.int <- function(z.in){
      risk.marg.B.given.Z(z.in)*dnorm(z.in, 0, frail.var)
    } 
    
    # Joint risk
    joint.risk.for.int <- function(z.in){
      risk.marg.A.given.Z(z.in)*risk.marg.B.given.Z(z.in)*dnorm(z.in, 0, frail.var)
    } 
    
    # Marg surv for A (integrate over 1 - risk)
    marg.surv.A.for.int <- function(z.in){
      (1 - risk.marg.A.given.Z(z.in))*dnorm(z.in, 0, frail.var)
    }
    
    # Marg surv for B (integrate over 1 - risk)
    marg.surv.B.for.int <- function(z.in){
      (1 - risk.marg.B.given.Z(z.in))*dnorm(z.in, 0, frail.var)
    } 
    
    # Joint surv
    joint.surv.for.int <- function(z.in){
      (1 - risk.marg.A.given.Z(z.in))*(1 - risk.marg.B.given.Z(z.in))*dnorm(z.in, 0, frail.var)
    } 
    
    ## Estimate the survival probabilities and risks by integrating these functions over the random effect
    
    ## Interate between 0 (lower limit) and 99.99th percentiles of the random effect distribution, otherwise we get some issues
    int.lower <- qnorm(0.00001, mean = 0, sd = frail.var)
    int.upper <- qnorm(0.99999, mean = 0, sd = frail.var)
    
    risk.marg.true.A <- cubintegrate(f = marg.risk.A.for.int, lower = int.lower, upper = int.upper, method = int.method)$integral
    risk.marg.true.B <- cubintegrate(f = marg.risk.B.for.int, lower = int.lower, upper = int.upper, method = int.method)$integral
    risk.joint.true <- cubintegrate(f = joint.risk.for.int, lower = int.lower, upper = int.upper, method = int.method)$integral
    
    surv.marg.true.A <- cubintegrate(f = marg.surv.A.for.int, lower = int.lower, upper = int.upper, method = int.method)$integral
    surv.marg.true.B <- cubintegrate(f = marg.surv.B.for.int, lower = int.lower, upper = int.upper, method = int.method)$integral
    surv.joint.true <- cubintegrate(f = joint.surv.for.int, lower = int.lower, upper = int.upper, method = int.method)$integral
  }
  
  ## Create output object
  output.obj <- c("surv.marg.true.A" = surv.marg.true.A,
                  "surv.marg.true.B" = surv.marg.true.B,
                  "surv.joint.true" = surv.joint.true,
                  "risk.marg.true.A" = risk.marg.true.A,
                  "risk.marg.true.B" = risk.marg.true.B,
                  "risk.joint.true" = risk.joint.true)
  
  return(output.obj)
}

##############
### Save image
##############
save.image("data/sim_function_calc_true_risk_all_1cont1cat.RData")



# ### Run some examples
# true.risk.msm <- calc.true.risk.msm(t.eval = 1, x1.eval = 0.3, x2.eval = 0,
#                                     shape12 = 1, scale12 = 1.5,
#                                     shape13 = 1, scale13 = 1,
#                                     shape24 = 1, scale24 = 0.8,
#                                     shape34 = 1, scale34 = 1.3,
#                                     beta12.cont = 0.25, beta12.cat = 0.5,
#                                     beta13.cont = 0.5, beta13.cat = 0.5,
#                                     beta24.cont = 0.33, beta24.cat = 0.25,
#                                     beta34.cont = 0.25, beta34.cat = 0.25)
# true.risk.msm
# 
# true.risk.copula <- calc.true.risk.copula(t.eval = 1, x1.eval = 0.1, x2.eval = 0,
#                                           copula = "clayton", copula.param = 3,
#                                           baselineA <- c(1,1.5), COV_betaA = c(0, 0), 
#                                           baselineB <- c(1,1), COV_betaB = c(0, 0))
# 
# true.risk.copula
# 
# true.risk.frailty.gamma <- calc.true.risk.frailty(
#   t.eval = 1, x1.eval = 0.1, x2.eval = 0,
#   frail.dist = "gamma", frail.var = 0.25,
#   baselineA <- c(1,1.5), COV_betaA = c(0, 0), 
#   baselineB <- c(1,1), COV_betaB = c(0, 0),
#   int.lower = 0.0001, int.upper = 20, int.method = "hcubature")
# 
# true.risk.frailty.normal <- calc.true.risk.frailty(
#   t.eval = 1, x1.eval = 0.1, x2.eval = 0,
#   frail.dist = "normal", frail.var = 0.25,
#   baselineA <- c(1,1.5), COV_betaA = c(0, 0), 
#   baselineB <- c(1,1), COV_betaB = c(0, 0),
#   int.lower = -25, int.upper = 25, int.method = "pcubature")
# 
# true.risk.frailty.gamma
# true.risk.frailty.normal