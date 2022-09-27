### This code is for the large sample validation, so the models are developed and calibrated on the same dataset

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.5/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions.R")
source("code/sim_functions_results.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")

### Define number of outcomes
K.sim <- 3

## Size of development datasets
N.devel.sim <- 500000


calc.sample.size <- function(coef.sim.multinomial, beta0.sim.multinomial, coef.sim.seqlog, beta0.sim.seqlog, N.devel.sim, K.sim, P.sim){
  
  ### Generate development datasets
  set.seed(1)
  
  dat.devel.list <- vector("list", 2)
  names(dat.devel.list) <- c("DGM.multinomial", "DGM.seqlog")
  dat.devel.list[[1]] <- generate.data.DGM.multinomial(K = K.sim, P = P.sim, coef = coef.sim.multinomial, beta0 = beta0.sim.multinomial, N = N.devel.sim)
  dat.devel.list[[2]] <- generate.data.DGM.seqlog(K = K.sim, P = P.sim, coef = coef.sim.seqlog, beta0 = beta0.sim.seqlog, N = N.devel.sim)
  
  ### Check prevalence of each outcome
  prop.table(table(dat.devel.list[[1]]$Y))
  prop.table(table(dat.devel.list[[2]]$Y))
  table(dat.devel.list[[1]]$Y)
  table(dat.devel.list[[2]]$Y)
  
  ### Need to fit pairwise models to get pairwise R2
  for (i in 1:2){
    dat.devel.list[[i]] <- dat.devel.list[[i]] %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                                                          Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                                          Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_))
  }
  
  fit.OvO.list <- vector("list", 2)
  for (i in 1:2){
    fit.OvO.list[[i]] <- vector("list", choose(K.sim, 2))
    
    fit.OvO.list[[i]][[1]] <- glm(Y1.2 ~ ., family = binomial(link = "logit"), data = subset(dat.devel.list[[i]], Y %in% c(1,2)) %>% 
                                    dplyr::select(-c(Y, Y.fac, Y1.3, Y2.3)))
    fit.OvO.list[[i]][[2]] <- glm(Y1.3 ~ ., family = binomial(link = "logit"), data = subset(dat.devel.list[[i]], Y %in% c(1,3)) %>% 
                                    dplyr::select(-c(Y, Y.fac, Y1.2, Y2.3)))
    fit.OvO.list[[i]][[3]] <- glm(Y2.3 ~ ., family = binomial(link = "logit"), data = subset(dat.devel.list[[i]], Y %in% c(2,3)) %>% 
                                    dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3)))
  }

  R2.list <- vector("list", 2)
  for (i in 1:2){
    R2.list[[i]] <- vector("list", choose(K.sim, 2))
    
    for (j in 1:choose(K.sim, 2)){
      R2.list[[i]][[j]] <- 1 - exp(-(fit.OvO.list[[i]][[j]]$null.deviance - fit.OvO.list[[i]][[j]]$deviance)/nrow(fit.OvO.list[[i]][[j]]$data))
    }
  }
  
  prop <- vector("list", 2)
  for (i in 1:2){
    prop[[i]] <- vector("list", choose(K.sim, 2))
    
    for (j in 1:choose(K.sim, 2)){
      prop[[i]][[j]] <- nrow(fit.OvO.list[[i]][[j]]$data)/N.devel.sim
    }
  }

  N.list <- vector("list", 2)
  for (i in 1:2){
    N.list[[i]] <- vector("list", choose(K.sim, 2))
    
    for (j in 1:choose(K.sim, 2)){
      N.list[[i]][[j]] <- (P.sim/((0.9-1)*log(1 - R2.list[[i]][[j]]/0.9)))/prop[[i]][[j]]
    }
  }
  
  N.1 <- c(max(unlist(N.list[[1]])), max(unlist(N.list[[2]])))
  
  ## Fit the multinomial logistic regression model
  fit.multinomial.list <- vector("list", 2)
  fit.multinomial.null.list <- vector("list", 2)
  
  fit.multinomial.list[[1]] <- vgam(Y.fac ~ x1 + x2 +x3 +x4 +x5, family = multinomial(ref = "1"), data = subset(dat.devel.list[[1]], select = -c(Y)))
  fit.multinomial.list[[2]] <- vgam(Y.fac ~ x1 + x2 +x3 +x4 +x5, family = multinomial(ref = "1"), data = subset(dat.devel.list[[2]], select = -c(Y)))
  
  fit.multinomial.null.list[[1]] <- vgam(Y.fac ~ 1, family = multinomial(ref = "1"), data = subset(dat.devel.list[[1]], select = -c(Y)))
  fit.multinomial.null.list[[2]] <- vgam(Y.fac ~ 1, family = multinomial(ref = "1"), data = subset(dat.devel.list[[2]], select = -c(Y)))
  
  
  R2.mult <- c(0, 0)
  R2.mult.max <- c(0, 0)
  for (i in 1:2){
    R2.mult[i] <- 1 - exp(-2*(logLik(fit.multinomial.list[[i]]) - logLik(fit.multinomial.null.list[[i]]))/N.devel.sim)
    R2.mult.max[i] <- 1 - exp(2*logLik(fit.multinomial.null.list[[i]])/N.devel.sim)
  }
  
  N.2 <- c(0, 0)
  for (i in 1:2){
    N.2[i] <- (P.sim*(K.sim -1))/(((R2.mult[i]/(R2.mult[i] + 0.05*R2.mult.max[i]))-1)*log(1 - R2.mult[i] - 0.05*R2.mult.max[i]))
  }
  
  N.3 <- c(0, 0)
  for (i in 1:2){
    probs <- prop.table(table(dat.devel.list[[i]]$Y))
    
    N_C3.1 <- qchisq(1-0.05/K.sim, 1)*probs[1]*(1-probs[1])/0.05^2 
    N_C3.2 <- qchisq(1-0.05/K.sim, 1)*probs[2]*(1-probs[2])/0.05^2 
    N_C3.3 <- qchisq(1-0.05/K.sim, 1)*probs[3]*(1-probs[3])/0.05^2
    
    N.3[i] <- max(N_C3.1, N_C3.2, N_C3.3)
  }
  
  return(list("DGM.mult" = c(N.1[1], N.2[1], N.3[1]), "DGM.seqlog" = c(N.1[2], N.2[2], N.3[2])))
}


##################
### Scenario 1 ###
##################

### Set seed.coef.sim
seed.coef.sim <- 1
set.seed(seed.coef.sim)

## And generate the simulation covariate effects
coef.sim <- rbind(c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)))

## Multinomial covariate effects
coef.sim.multinomial <- coef.sim
beta0.sim.multinomial <- c(-1, -1)*0.35

## Sequential logistic covariate effects
coef.sim.seqlog <- coef.sim
beta0.sim.seqlog <- c(0.5, 0)*2.25

samp.size.s1 <- calc.sample.size(coef.sim.multinomial = coef.sim.multinomial, 
                                 beta0.sim.multinomial = beta0.sim.multinomial, 
                                 coef.sim.seqlog = coef.sim.seqlog, 
                                 beta0.sim.seqlog = beta0.sim.seqlog, 
                                 N.devel.sim = N.devel.sim, K.sim = K.sim, P.sim = P.sim)
print("scen1")

##################
### Scenario 2 ###
##################

### Set seed.coef.sim
seed.coef.sim <- 1
set.seed(seed.coef.sim)

## And generate the simulation covariate effects
coef.sim <- rbind(c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))

## Multinomial covariate effects
coef.sim.multinomial <- coef.sim
beta0.sim.multinomial <- c(1, 1)*0.1

## Sequential logistic covariate effects
coef.sim.seqlog <- coef.sim
beta0.sim.seqlog <- c(0.5, 0)*2.25

samp.size.s2 <- calc.sample.size(coef.sim.multinomial = coef.sim.multinomial, 
                                 beta0.sim.multinomial = beta0.sim.multinomial, 
                                 coef.sim.seqlog = coef.sim.seqlog, 
                                 beta0.sim.seqlog = beta0.sim.seqlog, 
                                 N.devel.sim = N.devel.sim, K.sim = K.sim, P.sim = P.sim)
print("scen2")

##################
### Scenario 3 ###
##################

### Set seed.coef.sim
seed.coef.sim <- 1
set.seed(seed.coef.sim)

## And generate the simulation covariate effects
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)))

## Multinomial covariate effects
coef.sim.multinomial <- coef.sim
beta0.sim.multinomial <- c(-1, -1)*0.35

## Sequential logistic covariate effects
coef.sim.seqlog <- coef.sim
beta0.sim.seqlog <- c(0.5, 0)*2.25

samp.size.s3 <- calc.sample.size(coef.sim.multinomial = coef.sim.multinomial, 
                                 beta0.sim.multinomial = beta0.sim.multinomial, 
                                 coef.sim.seqlog = coef.sim.seqlog, 
                                 beta0.sim.seqlog = beta0.sim.seqlog, 
                                 N.devel.sim = N.devel.sim, K.sim = K.sim, P.sim = P.sim)
print("scen3")

##################
### Scenario 4 ###
##################

### Set seed.coef.sim
seed.coef.sim <- 1
set.seed(seed.coef.sim)

## And generate the simulation covariate effects
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))

## Multinomial covariate effects
coef.sim.multinomial <- coef.sim
beta0.sim.multinomial <- c(-1, -1)*0.35

## Sequential logistic covariate effects
coef.sim.seqlog <- coef.sim
beta0.sim.seqlog <- c(0.5, 0)*2.25

samp.size.s4 <- calc.sample.size(coef.sim.multinomial = coef.sim.multinomial, 
                                 beta0.sim.multinomial = beta0.sim.multinomial, 
                                 coef.sim.seqlog = coef.sim.seqlog, 
                                 beta0.sim.seqlog = beta0.sim.seqlog, 
                                 N.devel.sim = N.devel.sim, K.sim = K.sim, P.sim = P.sim)
print("scen4")

samp.size.s1
samp.size.s2
samp.size.s3
samp.size.s4

### Remove excess
rm(list=setdiff(ls(), list("samp.size.s1", "samp.size.s2", "samp.size.s3", "samp.size.s4")))

### Save image
save.image("data/sim_calc_min_sample_size.RData")
