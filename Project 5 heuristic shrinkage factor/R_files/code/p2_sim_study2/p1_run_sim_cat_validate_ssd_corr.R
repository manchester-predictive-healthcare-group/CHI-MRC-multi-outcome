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

### Assign number of iterations for each simulation
n.iter <- 250

### Validation dataset sample size
n.valid <- 1000000

### Assign number of bootstrap samples we will use to calculate S_boot
n.S.boot <- 200

### Assign dp
dp <- 3

### Assign number of simulations we are going to run
n.sim <- 25

### Create output objects
output.pavlou.S.VH <- vector("list", n.sim)
output.pavlou.S.boot <- vector("list", n.sim)
output.pavlou.S.pop <- vector("list", n.sim)
output.pavlou.R2.CS.app <- vector("list", n.sim)
output.pavlou.R2.NAGEL.app <- vector("list", n.sim)
output.pavlou.LR <- vector("list", n.sim)
output.pavlou.L.full <- vector("list", n.sim)
output.pavlou.L.null <- vector("list", n.sim)
output.pavlou.D <- vector("list", n.sim)
output.pavlou.C <- vector("list", n.sim)

output.riley.S.VH <- vector("list", n.sim)
output.riley.S.boot <- vector("list", n.sim)
output.riley.S.pop <- vector("list", n.sim)
output.riley.R2.CS.app <- vector("list", n.sim)
output.riley.R2.NAGEL.app <- vector("list", n.sim)
output.riley.LR <- vector("list", n.sim)
output.riley.L.full <- vector("list", n.sim)
output.riley.L.null <- vector("list", n.sim)
output.riley.D <- vector("list", n.sim)
output.riley.C <- vector("list", n.sim)

input.data <- vector("list", n.sim)

### Run n.sim simulations
for (sim in 1:n.sim){
  
  ### Print progress
  print(paste("set = ", set, "sim = ", sim, Sys.time()))  
  
  ### Set seed
  ### Having to do this manually, as the samplesizedev package resets the seed to set.seed(1) everytime it is ran
  set.seed(n.sim*(set-1) + sim)
  print(paste("the seed is ", n.sim*(set-1) + sim, Sys.time()))  
  
  ### Generate number of predictors
  P1 <- sample(1:30, 1)
  
  ### Define number of extra variables in the linear predictor that we will not adjust for
  P2 <- sample(1:30, 1)
  
  ### Total number of predictors for DGM
  P.total <- P1 + P2
  
  ### Generate covariate effects
  coef <- runif(P.total, -0.5, 0.5)
  
  ### Generate covariance matrix
  Q <- matrix(rep(FALSE, P.total*P.total), P.total)
  Lambda <- rep(1, P.total)
  covar_matrix <- rQ(Q, Lambda, f=rnorm)
  
  ### Generate marginal probabilities of each predictors
  marg_probs <- runif(length(coef), 0.05, 0.95)
  
  ### Assign beta (based around what outcome proportion it would equal to if all covariate effects were 0)
  temp.prop <- runif(1, 0.05, 0.95)
  beta0 <- -log((1-temp.prop)/temp.prop) 
  # This will be a beta0 between approximately (-3, 3)
  # Assigning beta0 uniformly between -3 and 3 would result in a lot of the probability distribution being
  # between 0.05 and 0.1, and 0.9 - 0.95.
  
  ### Create validation dataset
  dat.valid <- BLR.DGM.cat(coef = coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                           beta0 = beta0, #vector of intercepts
                           marg_probs = marg_probs, #marginal probabilities of each predictor
                           covar_matrix = covar_matrix, #covariance matrix for patient data
                           P = P1,
                           N = n.valid)
  
  ### Fit model in validation dataset to estimate the C-statistic that could be obtained, fitting model in population
  BLR.model.lrm.valid <- lrm(Y.fac ~ ., data = dplyr::select(dat.valid, - c(Y)), x = TRUE, y = TRUE, maxit = 1000)
  BLR.model.glm.valid <- glm(Y.fac ~ . , data = dplyr::select(dat.valid, - c(Y)), family = binomial(link = "logit")) 
  C.pop <- as.numeric(BLR.model.lrm.valid$stats["C"])
  
  ### Get the outcome proportion
  prop <- sum(dat.valid$Y)/nrow(dat.valid)
  
  ### Get the sample sizes
  # pavlou
  nreq.pavlou <- try(samplesizedev::samplesizedev(S = 0.9, phi = prop, c = C.pop, p = P1, nsim = 1000)$sim)
  
  # riley
  R2.CS.pop <- 1 - exp(-(BLR.model.glm.valid$null.deviance - BLR.model.glm.valid$deviance)/nrow(dat.valid))
  nreq.riley <- as.numeric(pmsampsize(type = "b", 
                                    csrsquared = R2.CS.pop, 
                                    parameters = P1, 
                                    shrinkage = 0.9, 
                                    prevalence = prop)$results_table["Criteria 1", "Samp_size"])

  ### Save input information
  input.data[[sim]] <- vector("list", 10)
  names(input.data[[sim]]) <- c("P", "P.total", "coef", "beta0", "covar_matrix", "nreq.pavlou", "nreq.riley", "prop", "C.pop", "R2.CS.pop")
  input.data[[sim]][["P"]] <- P1
  input.data[[sim]][["P.total"]] <- P.total
  input.data[[sim]][["coef"]] <- coef
  input.data[[sim]][["beta0"]] <- beta0
  input.data[[sim]][["covar_matrix"]] <- covar_matrix
  input.data[[sim]][["nreq.pavlou"]] <- nreq.pavlou
  input.data[[sim]][["nreq.riley"]] <- nreq.riley
  input.data[[sim]][["prop"]] <- sum(dat.valid$Y == 1)/nrow(dat.valid)
  input.data[[sim]][["C.pop"]] <- C.pop
  input.data[[sim]][["R2.CS.pop"]] <- R2.CS.pop
  
  ### Only run simulation if Riley sample size is < 10000, and nreq.pavlou could be estimated
  if (nreq.riley < 10000){
    
    ### Create objects to store output
    S.VH.vec.pavlou <- rep(NA, n.iter)
    S.boot.vec.pavlou <- rep(NA, n.iter)
    S.pop.vec.pavlou <- rep(NA, n.iter)
    R2.CS.app.vec.pavlou <- rep(NA, n.iter)
    R2.NAGEL.app.vec.pavlou <- rep(NA, n.iter)
    LR.vec.pavlou <- rep(NA, n.iter)
    L.full.vec.pavlou <- rep(NA, n.iter)
    L.null.vec.pavlou <- rep(NA, n.iter)
    D.vec.pavlou <- rep(NA, n.iter)
    C.vec.pavlou <- rep(NA, n.iter)
    
    S.VH.vec.riley <- rep(NA, n.iter)
    S.boot.vec.riley <- rep(NA, n.iter)
    S.pop.vec.riley <- rep(NA, n.iter)
    R2.CS.app.vec.riley <- rep(NA, n.iter)
    R2.NAGEL.app.vec.riley <- rep(NA, n.iter)
    LR.vec.riley <- rep(NA, n.iter)
    L.full.vec.riley <- rep(NA, n.iter)
    L.null.vec.riley <- rep(NA, n.iter)
    D.vec.riley <- rep(NA, n.iter)
    C.vec.riley <- rep(NA, n.iter)
    
    ### Print progress
    print(paste("sim = ", sim,
                "P = ", P1,
                "temp.prop = ", temp.prop,
                "nreq pavlou = ", nreq.pavlou, 
                "nreq riley = ", nreq.riley, 
                "time = ", Sys.time())) 
    
    ### Run simulation
    for (i in 1:n.iter){
      
      ### Print progress
      print(paste("iter = ", i, Sys.time()))  
      
      ###
      ### Run simulation for the sample sample calculate from the pavlou formula
      ###
      if (class(nreq.pavlou) != "try-error"){
        ### Create development dataset
        dat.devel.pavlou <- BLR.DGM.cat(coef = coef, 
                                    beta0 = beta0, 
                                    covar_matrix = covar_matrix,
                                    marg_probs = marg_probs,
                                    P = P1,
                                    N = nreq.pavlou)
        
        ### Calculate heuristic shrinkage factor, bootstrapped shrinkage factor, and actual shrinkage factor
        S.out.pavlou <- calc.S.pop(dat.devel = dat.devel.pavlou, 
                                               dat.valid = dat.valid,
                                               n.S.boot = n.S.boot)
        
        ### Assign the vectors
        S.pop.vec.pavlou[i] <- S.out.pavlou["S.pop"]
        R2.CS.app.vec.pavlou[i] <- S.out.pavlou["R2.CS.app"]
        R2.NAGEL.app.vec.pavlou[i] <- S.out.pavlou["R2.NAGEL.app"]
        LR.vec.pavlou[i] <- S.out.pavlou["LR"]
        L.full.vec.pavlou[i] <- S.out.pavlou["L.full"]
        L.null.vec.pavlou[i] <- S.out.pavlou["L.null"]
        D.vec.pavlou[i] <- S.out.pavlou["D"]
        C.vec.pavlou[i] <- S.out.pavlou["C"]
        
      } else {
        S.out.pavlou <- NA
      }
      
      ###
      ### Run simulation for the sample sample calculate from the pavlou formula
      ###
      
      ### Create development dataset
      dat.devel.riley <- BLR.DGM.cat(coef = coef, 
                                 beta0 = beta0, 
                                 covar_matrix = covar_matrix,
                                 marg_probs = marg_probs,
                                 P = P1,
                                 N = nreq.riley)
      
      ### Calculate heuristic shrinkage factor, bootstrapped shrinkage factor, and actual shrinkage factor
      S.out.riley <- calc.S.pop(dat.devel = dat.devel.riley, 
                                            dat.valid = dat.valid,
                                            n.S.boot = n.S.boot)
      
      ### Assign the vectors
      S.pop.vec.riley[i] <- S.out.riley["S.pop"]
      R2.CS.app.vec.riley[i] <- S.out.riley["R2.CS.app"]
      R2.NAGEL.app.vec.riley[i] <- S.out.riley["R2.NAGEL.app"]
      LR.vec.riley[i] <- S.out.riley["LR"]
      L.full.vec.riley[i] <- S.out.riley["L.full"]
      L.null.vec.riley[i] <- S.out.riley["L.null"]
      D.vec.riley[i] <- S.out.riley["D"]
      C.vec.riley[i] <- S.out.riley["C"]
      
    }
    
    ### Assign vectors of results to output object
    if (class(nreq.pavlou) != "try-error"){
      output.pavlou.S.pop[[sim]] <- S.pop.vec.pavlou
      output.pavlou.R2.CS.app[[sim]] <- R2.CS.app.vec.pavlou
      output.pavlou.R2.NAGEL.app[[sim]] <- R2.NAGEL.app.vec.pavlou
      output.pavlou.LR[[sim]] <- LR.vec.pavlou
      output.pavlou.L.full[[sim]] <- L.full.vec.pavlou
      output.pavlou.L.null[[sim]] <- L.null.vec.pavlou
      output.pavlou.D[[sim]] <- D.vec.pavlou
      output.pavlou.C[[sim]] <- C.vec.pavlou
    }
    
    output.riley.S.pop[[sim]] <- S.pop.vec.riley
    output.riley.R2.CS.app[[sim]] <- R2.CS.app.vec.riley
    output.riley.R2.NAGEL.app[[sim]] <- R2.NAGEL.app.vec.riley
    output.riley.LR[[sim]] <- LR.vec.riley
    output.riley.L.full[[sim]] <- L.full.vec.riley
    output.riley.L.null[[sim]] <- L.null.vec.riley
    output.riley.D[[sim]] <- D.vec.riley
    output.riley.C[[sim]] <- C.vec.riley
    
    rm(dat.valid, BLR.model.lrm.valid, BLR.model.glm.valid)
    
    ### Save output
    save.image(paste("data/run_sim_cat_validate_ssd_corr", set, ".RData", sep = ""))
    
  } else {
    
    rm(dat.valid, BLR.model.lrm.valid, BLR.model.glm.valid)
    
    ### Save output
    save.image(paste("data/run_sim_cat_validate_ssd_corr", set, ".RData", sep = ""))
    
  }
  
}

### Save output
save.image(paste("data/run_sim_cat_validate_ssd_corr", set, ".RData", sep = ""))