### This code is for the large sample validation, so the models are developed and calibrated on the same dataset

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.4/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")

### Define number of outcomes
K.sim <- 3

## Size of development datasets
N.devel.sim <- 500000


calc.sample.size <- function(coef.sim.multinomial, beta0.sim.multinomial, coef.sim.seqlog, beta0.sim.seqlog, N.devel.sim, K.sim, P.sim){
  
  ### Set seed
  set.seed(1)
  
  ### Generate development datasets
  dat.devel.list <- vector("list", 2)
  names(dat.devel.list) <- c("DGM.multinomial", "DGM.seqlog")
  dat.devel.list[[1]] <- generate.data.DGM.mult(K = K.sim, P = P.sim, coef = coef.sim.multinomial, beta0 = beta0.sim.multinomial, N = N.devel.sim)
  dat.devel.list[[2]] <- generate.data.DGM.seqlog(K = K.sim, P = P.sim, coef = coef.sim.seqlog, beta0 = beta0.sim.seqlog, N = N.devel.sim)
  
  ### Check prevalence of each outcome
  prop.table(table(dat.devel.list[[1]]$Y))
  prop.table(table(dat.devel.list[[2]]$Y))
  
  ### Create outcome for binary model
  for (i in 1:2){
    dat.devel.list[[i]] <- mutate(dat.devel.list[[i]], Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                                          Y.fac != 1 ~ 0))
    dat.devel.list[[i]]$Y1.bin <- as.factor(dat.devel.list[[i]]$Y1.bin)
  }

  ### Fit BLR model and calc R2, and store in objects
  fit.list <- vector("list", 2)
  fit.null.list <- vector("list", 2)
  R2.list <- c(0, 0)
  R2.max.list <- c(0, 0)
  
  for (i in 1:2){
    fit.list[[i]] <- glm(Y1.bin ~ ., family = binomial(link = 'logit'), data = subset(dat.devel.list[[i]], select = -c(Y, Y.fac)))
    fit.null.list[[i]] <- glm(Y1.bin ~ 1, family = binomial(link = 'logit'), data = subset(dat.devel.list[[i]], select = -c(Y, Y.fac)))
    R2.list[i] <- 1 - exp(-(fit.list[[i]]$null.deviance - fit.list[[i]]$deviance)/N.devel.sim)
    R2.max.list[i] <- as.numeric(1 - exp(2*logLik(fit.null.list[[i]])/N.devel.sim))
  }
  
  ## Calculate N.1 using critirion (i)
  N.1 <- P.sim/((0.9-1)*log(1 - R2.list/0.9))
  
  ## Calcualte N.2 using criterion (ii)
  N.2 <- c(0, 0)
  for (i in 1:2){
    N.2[i] <- (P.sim*(K.sim -1))/(((R2.list[i]/(R2.list[i] + 0.05*R2.max.list[i]))-1)*log(1 - R2.list[i] - 0.05*R2.max.list[i]))
  }
  
  ## Calcualte N.3 using criterion (iii)
  N.3 <- c(0, 0)
  for (i in 1:2){
    probs <- prop.table(table(dat.devel.list[[i]]$Y))[1]
    
    N.3[i] <- probs*(1-probs)*(1.96/0.05)^2 
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


samp.size.s1
samp.size.s2
samp.size.s3
samp.size.s4

### Remove excess
rm(list=setdiff(ls(), list("samp.size.s1", "samp.size.s2", "samp.size.s3", "samp.size.s4")))

### Save image
save.image("data/sim_calc_min_sample_size.RData")