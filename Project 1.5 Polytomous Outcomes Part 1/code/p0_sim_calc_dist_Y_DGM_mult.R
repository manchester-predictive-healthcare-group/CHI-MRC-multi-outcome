######## What is the min and max of P(Y=1) for these different scenarios
### Run a little version of the simulation

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.5/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions.R")


### Scenario 1, K = 5
N.sim <- 500000
K.sim <- 5
P.sim <- 5
beta0.sim <- c(-1, -1, -1, -1)*0.35

#set.seed(101)
coef.sim <- rbind(c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)))

p1.probs <- generate.outcome.probs.DGM.mult(K = K.sim, #number of outcome categories
                                            P = P.sim, #number of predictors
                                            coef = coef.sim, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                            beta0 = beta0.sim, #vector of intercepts
                                            N = N.sim)
#hist(p1.probs[[1]]$p1)
colMeans(p1.probs[[1]])
quantile(p1.probs[[1]]$p1, probs = c(0.025, 0.975))


### Scenario 2, K = 5
N.sim <- 500000
K.sim <- 5
P.sim <- 5
beta0.sim <- c(1, 1, 1, 1)*0.1

#set.seed(101)
coef.sim <- rbind(c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))

p1.probs <- generate.outcome.probs.DGM.mult(K = K.sim, #number of outcome categories
                                            P = P.sim, #number of predictors
                                            coef = coef.sim, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                            beta0 = beta0.sim, #vector of intercepts
                                            N = N.sim)
#hist(p1.probs[[1]]$p1)
colMeans(p1.probs[[1]])
quantile(p1.probs$p1[[1]], probs = c(0.025, 0.975))



### Scenario 3, K = 5
N.sim <- 500000
K.sim <- 5
P.sim <- 5
beta0.sim <- c(-1, -1, -1, -1)*0.35

#set.seed(101)
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)))

p1.probs <- generate.outcome.probs.DGM.mult(K = K.sim, #number of outcome categories
                                            P = P.sim, #number of predictors
                                            coef = coef.sim, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                            beta0 = beta0.sim, #vector of intercepts
                                            N = N.sim)
#hist(p1.probs[[1]]$p1)
colMeans(p1.probs[[1]])
quantile(p1.probs[[1]]$p1, probs = c(0.025, 0.975))


### Scenario 4, K = 5
N.sim <- 500000
K.sim <- 5
P.sim <- 5
beta0.sim <- c(-1, -1, -1, -1)*0.35

#set.seed(101)
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))

p1.probs <- generate.outcome.probs.DGM.mult(K = K.sim, #number of outcome categories
                                            P = P.sim, #number of predictors
                                            coef = coef.sim, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                            beta0 = beta0.sim, #vector of intercepts
                                            N = N.sim)
#hist(p1.probs[[1]]$p1)
colMeans(p1.probs[[1]])
quantile(p1.probs[[1]]$p1, probs = c(0.025, 0.975))



### Scenario 1, K = 3
N.sim <- 500000
K.sim <- 3
P.sim <- 5
beta0.sim <- c(-1, -1)*0.35

#set.seed(101)
coef.sim <- rbind(c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)),
                  c(runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1), runif(1,-1,1)))

p1.probs <- generate.outcome.probs.DGM.mult(K = K.sim, #number of outcome categories
                                            P = P.sim, #number of predictors
                                            coef = coef.sim, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                            beta0 = beta0.sim, #vector of intercepts
                                            N = N.sim)
#hist(p1.probs[[1]]$p1)
colMeans(p1.probs[[1]])
quantile(p1.probs[[1]]$p1, probs = c(0.025, 0.975))


### Scenario 2, K = 3
N.sim <- 500000
K.sim <- 3
P.sim <- 5
beta0.sim <- c(1, 1)*0.1

#set.seed(101)
coef.sim <- rbind(c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))

p1.probs <- generate.outcome.probs.DGM.mult(K = K.sim, #number of outcome categories
                                            P = P.sim, #number of predictors
                                            coef = coef.sim, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                            beta0 = beta0.sim, #vector of intercepts
                                            N = N.sim)
#hist(p1.probs[[1]]$p1)
colMeans(p1.probs[[1]])
quantile(p1.probs$p1[[1]], probs = c(0.025, 0.975))



### Scenario 3, K = 3
N.sim <- 500000
K.sim <- 3
P.sim <- 5
beta0.sim <- c(-1, -1)*0.35

#set.seed(101)
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), runif(1,0,1), runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)))

p1.probs <- generate.outcome.probs.DGM.mult(K = K.sim, #number of outcome categories
                                            P = P.sim, #number of predictors
                                            coef = coef.sim, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                            beta0 = beta0.sim, #vector of intercepts
                                            N = N.sim)
#hist(p1.probs[[1]]$p1)
colMeans(p1.probs[[1]])
quantile(p1.probs[[1]]$p1, probs = c(0.025, 0.975))


### Scenario 4, K = 5
N.sim <- 500000
K.sim <- 3
P.sim <- 5
beta0.sim <- c(-1, -1)*0.35

#set.seed(101)
coef.sim <- rbind(c(-runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1), -runif(1,0,1)),
                  c(runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1), runif(1,0,1)))

p1.probs <- generate.outcome.probs.DGM.mult(K = K.sim, #number of outcome categories
                                            P = P.sim, #number of predictors
                                            coef = coef.sim, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                            beta0 = beta0.sim, #vector of intercepts
                                            N = N.sim)
#hist(p1.probs[[1]]$p1)
colMeans(p1.probs[[1]])
quantile(p1.probs[[1]]$p1, probs = c(0.025, 0.975))






