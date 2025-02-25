### Set working directory
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_8.2/")

### Source functions
R.func.sources <- list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load functions
source("R/sim_functions.R")

### Assign number for which repeat of this code we are using
args <- commandArgs(trailingOnly = TRUE)
set <- as.numeric(args[1])
print(set)
set.seed(set)

### Assign number of iterations for each simulation
n.iter <- 250

# ### Assign number of datasets we generate to estimate R2
# n.estR2 <- 500

### Validation dataset sample size
n.valid <- 1000000

### Assign number of bootstrap samples we will use to calculate S_boot
n.S.boot <- 200

### Assign dp
dp <- 3

### Assign number of simulations we are going to run
n.sim <- 50

### Create output objects
output.S.VH <- vector("list", n.sim)
output.S.boot <- vector("list", n.sim)
output.S.pop <- vector("list", n.sim)
output.R2.CS.app <- vector("list", n.sim)
output.R2.NAGEL.app <- vector("list", n.sim)
output.LR <- vector("list", n.sim)
output.L.full <- vector("list", n.sim)
output.L.null <- vector("list", n.sim)
output.D <- vector("list", n.sim)
output.C <- vector("list", n.sim)
input.data <- vector("list", n.sim)

### Run n.sim simulations
for (sim in 1:n.sim){
  
  ### Print progress
  print(paste("sim = ", sim, Sys.time()))  
  
  ### Generate number of predictors
  P1 <- sample(1:30, 1)
  
  ### Define number of extra variables in the linear predictor that we will not adjust for
  P2 <- sample(1:30, 1)
  
  ### Total number of predictors for DGM
  P.total <- P1 + P2
  
  ### Generate covariate effects
  coef <- runif(P.total, -0.5, 0.5)
  
  ### Generate covariance matrix
  covar_matrix <- diag(P.total)
  
  ### Assign beta (based around what outcome proportion it would equal to if all covariate effects were 0)
  temp.prop <- runif(1, 0.05, 0.95)
  beta0 <- -log((1-temp.prop)/temp.prop) 
  # This will be a beta0 between approximately (-3, 3)
  # Assigning beta0 uniformly between -3 and 3 would result in a lot of the probability distribution being
  # between 0.05 and 0.1, and 0.9 - 0.95.
  
  ## Pick random calculate required sample size
  nreq <- 5001
  while (nreq > 5000){
    nreq <- round(100 + rweibull(1, 0.5, 2500))
  }
  
  ### Create validation dataset
  dat.valid <- BLR.DGM(coef = coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                       beta0 = beta0, #vector of intercepts
                       covar_matrix = covar_matrix, #covariance matrix for patient covariates
                       P = P1,
                       N = n.valid)
  
  ### Fit model in validation dataset to estimate the C-statistic that could be obtained, fitting model in population
  BLR.model.lrm.valid <- lrm(Y.fac ~ ., data = dplyr::select(dat.valid, - c(Y)), x = TRUE, y = TRUE, maxit = 1000)
  BLR.model.glm.valid <- glm(Y.fac ~ . , data = dplyr::select(dat.valid, - c(Y)), family = binomial(link = "logit")) 
  C.pop <- as.numeric(BLR.model.lrm.valid$stats["C"])
  R2.CS.pop <- 1 - exp(-(BLR.model.glm.valid$null.deviance - BLR.model.glm.valid$deviance)/nrow(dat.valid))
  
  ### Save input information
  input.data[[sim]] <- vector("list", 9)
  names(input.data[[sim]]) <- c("P", "P.total", "coef", "beta0", "covar_matrix", "nreq", "prop", "C.pop", "R2.CS.pop")
  input.data[[sim]][["P"]] <- P1
  input.data[[sim]][["P.total"]] <- P.total
  input.data[[sim]][["coef"]] <- coef
  input.data[[sim]][["beta0"]] <- beta0
  input.data[[sim]][["covar_matrix"]] <- covar_matrix
  input.data[[sim]][["nreq"]] <- nreq
  input.data[[sim]][["prop"]] <- sum(dat.valid$Y == 1)/nrow(dat.valid)
  input.data[[sim]][["C.pop"]] <- C.pop
  input.data[[sim]][["R2.CS.pop"]] <- R2.CS.pop
  
  ### Create objects to store output
  S.VH.vec <- rep(NA, n.iter)
  S.boot.vec <- rep(NA, n.iter)
  S.pop.vec <- rep(NA, n.iter)
  R2.CS.app.vec <- rep(NA, n.iter)
  R2.NAGEL.app.vec <- rep(NA, n.iter)
  LR.vec <- rep(NA, n.iter)
  L.full.vec <- rep(NA, n.iter)
  L.null.vec <- rep(NA, n.iter)
  D.vec <- rep(NA, n.iter)
  C.vec <- rep(NA, n.iter)
  
  ### Print progress
  print(paste("sim = ", sim,
              "P = ", P1,
              "temp.prop = ", temp.prop,
              "nreq = ", nreq, 
              "time = ", Sys.time())) 
              
  ### Run simulation
  for (i in 1:n.iter){
    
    ### Print progress
    print(paste("iter = ", i, Sys.time()))  
    
    ### Create development dataset
    dat.devel <- BLR.DGM(coef = coef, 
                         beta0 = beta0, 
                         covar_matrix = covar_matrix,
                         P = P1,
                         N = nreq)
    
    ### Calculate heuristic shrinkage factor, bootstrapped shrinkage factor, and actual shrinkage factor
    S.out <- calc.S.VH.S.boot.S.pop(dat.devel = dat.devel, 
                             dat.valid = dat.valid,
                             n.S.boot = n.S.boot)
    
    ### Assign the vectors
    S.VH.vec[i] <- S.out["S.VH"]
    S.boot.vec[i] <- S.out["S.boot"]
    S.pop.vec[i] <- S.out["S.pop"]
    R2.CS.app.vec[i] <- S.out["R2.CS.app"]
    R2.NAGEL.app.vec[i] <- S.out["R2.NAGEL.app"]
    LR.vec[i] <- S.out["LR"]
    L.full.vec[i] <- S.out["L.full"]
    L.null.vec[i] <- S.out["L.null"]
    D.vec[i] <- S.out["D"]
    C.vec[i] <- S.out["C"]
    
  }
  
  ### Assign mean of S.VH and S.pop to output object
  output.S.VH[[sim]] <- S.VH.vec
  output.S.boot[[sim]] <- S.boot.vec
  output.S.pop[[sim]] <- S.pop.vec
  output.R2.CS.app[[sim]] <- R2.CS.app.vec
  output.R2.NAGEL.app[[sim]] <- R2.NAGEL.app.vec
  output.LR[[sim]] <- LR.vec
  output.L.full[[sim]] <- L.full.vec
  output.L.null[[sim]] <- L.null.vec
  output.D[[sim]] <- D.vec
  output.C[[sim]] <- C.vec
  
  rm(dat.valid, BLR.model.lrm.valid, BLR.model.glm.valid)
  
  ### Save output
  save.image(paste("data/run_sim_nocorr", set, ".RData", sep = ""))
  
}

### Save output
save.image(paste("data/run_sim_nocorr", set, ".RData", sep = ""))