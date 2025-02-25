########################################################################
### 1.1A) Function to generate the data using the binary logistic DGM ###
########################################################################

# ### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
# BLR.DGM.OLD <- function(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
#                     beta0, #vector of intercepts
#                     covar_matrix, #covariance matrix for patient data
#                     N){
#   
#   ### Define number of predictors
#   P <- length(coef)
# 
#   ### Geneate covariate from a multivariate normal
#   dat.covar <- mvrnorm(N, mu = rep(0, P), Sigma = covar_matrix)
#   
#   ## Add variable names
#   colnames(dat.covar) <- paste("x", 1:P, sep = "")
#   
#   ### Create linear predictor
#   probs <- exp(beta0 + as.matrix(dat.covar)%*%coef)/(1 + exp(beta0 + as.matrix(dat.covar)%*%coef))
#   
#   ### Create outcome
#   Y <- rbinom(N, 1, probs)
#   
#   ### Add to an output dataset
#   dat.out <- data.frame("Y" = Y, 
#                         #"probs" = probs, 
#                         dat.covar)
#   dat.out$Y.fac <- factor(dat.out$Y)
#   
#   return(dat.out)
# }


#########################################################################
### 1.1B) Function to generate the data using the binary logistic DGM ###
### but we introduce extra variables which cause variation in the outcome 
### which we are not acconuting for in our models (with the OLD DGM, there)
### is basically a 1-1 relationship as number variables increases so does
### R2, etc 
#########################################################################

### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
BLR.DGM <- function(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                    beta0, #intercept
                    covar_matrix, #covariance matrix for patient data
                    P, #number of variables we will adjust for in model (variables we have access to)
                    N){
  
  ### Define number of predictors
  P.total <- length(coef)
 
  ### Geneate covariate from a multivariate normal
  dat.covar <- mvrnorm(N, mu = rep(0, P.total), Sigma = covar_matrix)
  
  ## Add variable names
  colnames(dat.covar) <- paste("x", 1:P.total, sep = "")
  
  ### Create linear predictor
  probs <- exp(beta0 + as.matrix(dat.covar)%*%coef)/(1 + exp(beta0 + as.matrix(dat.covar)%*%coef))
  
  ### Create outcome
  Y <- rbinom(N, 1, probs)
  
  ### Only retain first P columns of covar (i.e. variables we have access to/will adjust for)
  dat.covar <- dat.covar[, 1:P]
  
  ### Add to an output dataset
  dat.out <- data.frame("Y" = Y, 
                        #"probs" = probs, 
                        dat.covar)
  dat.out$Y.fac <- factor(dat.out$Y)
  
  return(dat.out)
}



########################################################################
### 1.2) Function to generate the data using the binary logistic DGM ###
########################################################################

# ### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
# BLR.DGM.cat.old <- function(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
#                     beta0, #vector of intercepts
#                     N){
#   
#   #   N <- 10
#   #   coef <- c(0.75, 0.75, 0.75, 0.75, 0.75)
#   #   beta0 <- 0
#   ### Define number of predictors
#   P <- length(coef)
#   
#   ### Create an input dataset of patient covariates
#   ## Create empty dataframe
#   dat.covar <- data.frame(matrix(, nrow = N, ncol = P))
#   ## Add variable names
#   colnames(dat.covar) <- paste("x", 1:P, sep = "")
#   ## Add covariates
#   for (i in 1:P){
#     dat.covar[, i] <- rbinom(N, 1, 0.5)
#   }
#   
#   ### Create linear predictor
#   probs <- exp(beta0 + as.matrix(dat.covar)%*%coef)/(1 + exp(beta0 + as.matrix(dat.covar)%*%coef))
#   
#   ### Create outcome
#   Y <- rbinom(N, 1, probs)
#   
#   ### Add to an output dataset
#   dat.out <- data.frame("Y" = Y, 
#                         #"probs" = probs, 
#                         dat.covar)
#   dat.out$Y.fac <- factor(dat.out$Y)
#   
#   return(dat.out)
# }


###
### DGL DGM Categorical no correlation
###

### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
BLR.DGM.cat <- function(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                        beta0, #intercept
                        marg_probs, #marginal probabilities of each predictor
                        covar_matrix, #covariance matrix for patient data
                        P, #number of variables we will adjust for in model (variables we have access to)
                        N){
  
  ### Define number of predictors
  P.total <- length(coef)
  
  ### Create an input dataset of patient covariates
  dat.covar.matrix <- bindata::rmvbin(N, margprob = marg_probs, sigma = covar_matrix)
  
  ### Add variable names
  colnames(dat.covar.matrix) <- paste("x", 1:P.total, sep = "")
  
  ### Turn into factors and create data frame
  dat.covar <- data.frame(dat.covar.matrix)
  dat.covar[] <- lapply(dat.covar, as.factor)
  
  ### Create linear predictor
  probs <- exp(beta0 + dat.covar.matrix%*%coef)/(1 + exp(beta0 + dat.covar.matrix%*%coef))
  
  ### Create outcome
  Y <- rbinom(N, 1, probs)
  
  ### Only retain first P columns of covar (i.e. variables we have access to/will adjust for)
  dat.covar <- dat.covar[, 1:P]
  
  ### Add to an output dataset
  dat.out <- data.frame("Y" = Y, 
                        #"probs" = probs, 
                        dat.covar)
  dat.out$Y.fac <- factor(dat.out$Y)
  
  return(dat.out)
}


# ###
# ### DGL DGM Categorical no correlation
# ###
# ### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
# BLR.DGM.cat.old2 <- function(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
#                     beta0, #intercept
#                     P, #number of variables we will adjust for in model (variables we have access to)
#                     N){
# 
#   ### Define number of predictors
#   P.total <- length(coef)
#   
#   ### Create an input dataset of patient covariates
#   ## Create empty dataframe
#   dat.covar <- data.frame(matrix(, nrow = N, ncol = P.total))
#   ## Add variable names
#   colnames(dat.covar) <- paste("x", 1:P.total, sep = "")
#   ## Add covariates
#   for (i in 1:P.total){
#     dat.covar[, i] <- as.factor(rbinom(N, 1, 0.5))
#   }
# 
#   ### Create linear predictor
#   probs <- exp(beta0 + apply(dat.covar, 2, as.numeric)%*%coef)/(1 + exp(beta0 + apply(dat.covar, 2, as.numeric)%*%coef))
#   
#   ### Create outcome
#   Y <- rbinom(N, 1, probs)
#   
#   ### Only retain first P columns of covar (i.e. variables we have access to/will adjust for)
#   dat.covar <- dat.covar[, 1:P]
#   
#   ### Add to an output dataset
#   dat.out <- data.frame("Y" = Y, 
#                         #"probs" = probs, 
#                         dat.covar)
#   dat.out$Y.fac <- factor(dat.out$Y)
#   
#   return(dat.out)
# }

####################################################################
### 1.3 Function to generate symetrical positive definite matrix ###
####################################################################
rQ <- function(Q, Lambda, f=rnorm) {
  normalize <- function(x) {
    v <- zapsmall(c(1, sqrt(sum(x * x))))[2]
    if (v == 0) v <- 1
    x / v
  }
  Q <- Q | t(Q)                    # Force symmetry by applying all constraints
  d <- nrow(Q) 
  if (missing(Lambda)) Lambda <- rep(1, d)
  R <- matrix(f(d^2), d, d)        # An array of column vectors
  for (i in seq_len(d)) {
    j <- which(Q[seq_len(i-1), i]) # Indices of the preceding orthogonal vectors
    R[, i] <- normalize(residuals(.lm.fit(R[, j, drop=FALSE], R[, i])))
  }
  R <- R %*% diag(Lambda)
  crossprod(R)
}



###
### 2.1) Fit a model and calculate heuristic shrinkage factor and actual shrinkage factor
###
calc.S.VH.S.pop <- function(dat.devel, dat.valid){
  
#   dat.devel <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
#                                    beta0, #vector of intercepts
#                                    N = 1000)
#   dat.valid <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
#                        beta0, #vector of intercepts
#                        N = 20000)
  
  ### Fit BLR in development dataset
  #BLR.model <- glm(Y.fac ~ . , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
  BLR.model.lrm <- lrm(Y.fac ~ ., data = dplyr::select(dat.devel, - c(Y)), x = TRUE, y = TRUE)
  #BLR.model.null <- glm(Y.fac ~ 1 , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
  
  if (BLR.model.lrm$fail == TRUE){
    return(c("S.VH" = NA, 
             "S.pop" = NA, 
             "R2.CS.app" = NA, 
             "R2.NAGEL.app" = NA,
             "D" = NA,
             "C" = NA,
             "LR" = NA,
             "L.full" = NA,
             "L.null" = NA,
             "P" = NA))
    
  }  else if (BLR.model.lrm$fail == FALSE){
    ### Calculate likelihood ratio statistic
    
    #LR <- BLR.model$null.deviance - BLR.model$deviance
    LR <- BLR.model.lrm$deviance[1] - BLR.model.lrm$deviance[2]
    
    ## Calculate log-likelihood
    loglik.full <- -0.5*BLR.model.lrm$deviance[2]
    loglik.null <- -0.5*BLR.model.lrm$deviance[1]
    
    
    ### Define P as number of predictor parameters estimated
    P <- length(BLR.model.lrm$coefficients)
    
    
    ###
    ### Calculate heuristic shrinkage factor
    ###
    
    S.VH <- as.numeric(1 - P/LR)
    
    
    ###
    ### Calculate shrinkage factor in population
    ###
    
    ### Generate predicted risks in validaiton dataset
    lp.valid <- predict(BLR.model.lrm, newdata = dat.valid, type = "lp")
    
    ### Add to validation dataset
    dat.valid$lp <- lp.valid
    
    ### Fit recalibration model to calculate shrinkage factor
    model.valid <- glm(Y.fac ~ lp, data = dat.valid, family = binomial(link = "logit")) 
    
    ### Return the slope as the shrinkage factor
    S.pop <- model.valid$coefficients["lp"]
    
    ###
    ### Calculate R2.CS.app, R2.CS.adj, R2.NAGEL.app, C, Dxy,
    ###
    R2.CS.app <- as.numeric(1 - exp(-(LR)/nrow(dat.devel)))
    R2.NAGEL.app <- R2.CS.app/(1 - exp(2*(-0.5*BLR.model.lrm$deviance[1])/nrow(dat.devel)))
    
    D <- as.numeric(BLR.model.lrm$stats["Dxy"])
    C <- as.numeric(BLR.model.lrm$stats["C"])
    
    
    ### Return heuristic and actual shrinkage factors
    return(c("S.VH" = S.VH, 
             #"S.boot" = S.boot,
             "S.pop" = as.numeric(S.pop), 
             "R2.CS.app" = R2.CS.app, 
             "R2.NAGEL.app" = R2.NAGEL.app,
             "D" = D,
             "C" = C,
             "LR" = LR,
             "L.full" = loglik.full,
             "L.null" = loglik.null,
             "P" = P))
  }
  
}
  

###
### 2.2) Fit a model and calculate heuristic shrinkage factor, shrinkage factor from bootstrapping,
### and actual shrinkage factor
###
calc.S.VH.S.boot.S.pop <- function(dat.devel, dat.valid, n.S.boot){

  #   dat.devel <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
  #                                    beta0, #vector of intercepts
  #                                    N = 1000)
  #   dat.valid <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
  #                        beta0, #vector of intercepts
  #                        N = 20000)
  
  ### Fit BLR in development dataset
  #BLR.model <- glm(Y.fac ~ . , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
  #BLR.model.null <- glm(Y.fac ~ 1 , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 

  BLR.model.lrm <- lrm(Y.fac ~ ., data = dplyr::select(dat.devel, - c(Y)), x = TRUE, y = TRUE, maxit = 1000)
  print("BLR.model.done")
  
  if (BLR.model.lrm$fail == TRUE){
    return(c("S.VH" = NA, 
             "S.boot" = NA,
             "S.pop" = NA, 
             "R2.CS.app" = NA, 
             "R2.NAGEL.app" = NA,
             "D" = NA,
             "C" = NA,
             "LR" = NA,
             "L.full" = NA,
             "L.null" = NA,
             "P" = NA))
    
  } else if (BLR.model.lrm$fail == FALSE){
    ### Calculate likelihood ratio statistic
    LR <- BLR.model.lrm$deviance[1] - BLR.model.lrm$deviance[2]
    
    ## Calculate log-likelihood
    loglik.full <- -0.5*BLR.model.lrm$deviance[2]
    loglik.null <- -0.5*BLR.model.lrm$deviance[1]
    
    ### Define P as number of predictor parameters estimated
    P <- length(BLR.model.lrm$coefficients)
    
    
    ###
    ### Calculate heuristic shrinkage factor
    ###
    S.VH <- as.numeric(1 - P/LR)
    
    
    ###
    ### Calculate S using bootstrapping
    ###
    
    ### First define a function to calculate S in a bootstrap sample
    function_calc_S_boot <- function(data, i){
      
      ### Take the bootstrap sample
      boot.dat <- data[i, ]
      
      ### Create model in bootstrapped dataset
      boot.model <- glm(Y.fac ~ . , data = dplyr::select(boot.dat, - c(Y)), family = binomial(link = "logit")) 
      
      ### Generate predictions using this model using the new dataset
      boot.lp <- predict(boot.model, newdata = dat.devel)
      
      ### Create a temporary dataset with both these things
      boot.dat.devel.temp <- data.frame(dat.devel, "boot.lp" = boot.lp)
      
      ### Calculate calibration slope
      boot.calib.model <- glm(Y.fac ~ boot.lp, data = boot.dat.devel.temp, family = binomial(link = "logit")) 
      
      ### Save slope
      return(boot.calib.model$coefficients["boot.lp"])
      
    }
    
    ### Run the bootstrapping
    ### Need to use try, specifically because of the scenarios with categorical predictors 
    ### When running the bootstrapping, for scenarios with small sample size
    ### and high number of predictors, we get bootstrap samples where the entire dataset has the same value
    ### for one of the predictors, which causes an error.
    boot.out <- try(boot(dat.devel, function_calc_S_boot, R = n.S.boot), silent = TRUE)
    
    ### Assign S.boot
    if (class(boot.out) == "try-error"){
      S.boot <- NA
    } else {
      S.boot <- mean(boot.out$t)
    }
    
    ###
    ### Calculate shrinkage factor in population
    ###
    
    ### Generate predicted risks in validaiton dataset
    lp.valid <- predict(BLR.model.lrm, newdata = dat.valid, type = "lp")
    
    ### Add to validation dataset
    dat.valid$lp <- lp.valid
    
    ### Fit recalibration model to calculate shrinkage factor
    model.valid <- glm(Y.fac ~ lp, data = dat.valid, family = binomial(link = "logit")) 
    
    ### Return the slope as the shrinkage factor
    S.pop <- model.valid$coefficients["lp"]
    
    ###
    ### Calculate R2.CS.app, R2.CS.adj, R2.NAGEL.app, C, Dxy,
    ###
    R2.CS.app <- as.numeric(1 - exp(-(LR)/nrow(dat.devel)))
    R2.NAGEL.app <- R2.CS.app/(1 - exp(2*(-0.5*BLR.model.lrm$deviance[1])/nrow(dat.devel)))
    
    D <- as.numeric(BLR.model.lrm$stats["Dxy"])
    C <- as.numeric(BLR.model.lrm$stats["C"])
    
    
    ### Return heuristic and actual shrinkage factors
    return(c("S.VH" = S.VH, 
             "S.boot" = S.boot,
             "S.pop" = as.numeric(S.pop), 
             "R2.CS.app" = R2.CS.app, 
             "R2.NAGEL.app" = R2.NAGEL.app,
             "D" = D,
             "C" = C,
             "LR" = LR,
             "L.full" = loglik.full,
             "L.null" = loglik.null,
             "P" = P))
  }
}


###
### 2.2) Fit a model and calculate required shrinkage shrinkage in validation cohort
###
calc.S.pop <- function(dat.devel, dat.valid, n.S.boot){
  
  #   dat.devel <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
  #                                    beta0, #vector of intercepts
  #                                    N = 1000)
  #   dat.valid <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
  #                        beta0, #vector of intercepts
  #                        N = 20000)
  
  ### Fit BLR in development dataset
  #BLR.model <- glm(Y.fac ~ . , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
  #BLR.model.null <- glm(Y.fac ~ 1 , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
  
  BLR.model.lrm <- lrm(Y.fac ~ ., data = dplyr::select(dat.devel, - c(Y)), x = TRUE, y = TRUE, maxit = 1000)
  print("BLR.model.done")
  
  if (BLR.model.lrm$fail == TRUE){
    return(c("S.VH" = NA, 
             "S.boot" = NA,
             "S.pop" = NA, 
             "R2.CS.app" = NA, 
             "R2.NAGEL.app" = NA,
             "D" = NA,
             "C" = NA,
             "LR" = NA,
             "L.full" = NA,
             "L.null" = NA,
             "P" = NA))
    
  } else if (BLR.model.lrm$fail == FALSE){
    ### Calculate likelihood ratio statistic
    LR <- BLR.model.lrm$deviance[1] - BLR.model.lrm$deviance[2]
    
    ## Calculate log-likelihood
    loglik.full <- -0.5*BLR.model.lrm$deviance[2]
    loglik.null <- -0.5*BLR.model.lrm$deviance[1]
    
    ### Define P as number of predictor parameters estimated
    P <- length(BLR.model.lrm$coefficients)
    
    ###
    ### Calculate shrinkage factor in population
    ###
    
    ### Generate predicted risks in validaiton dataset
    lp.valid <- predict(BLR.model.lrm, newdata = dat.valid, type = "lp")
    
    ### Add to validation dataset
    dat.valid$lp <- lp.valid
    
    ### Fit recalibration model to calculate shrinkage factor
    model.valid <- glm(Y.fac ~ lp, data = dat.valid, family = binomial(link = "logit")) 
    
    ### Return the slope as the shrinkage factor
    S.pop <- model.valid$coefficients["lp"]
    
    
    ###
    ### Calculate R2.CS.app, R2.CS.adj, R2.NAGEL.app, C, Dxy,
    ###
    R2.CS.app <- as.numeric(1 - exp(-(LR)/nrow(dat.devel)))
    R2.NAGEL.app <- R2.CS.app/(1 - exp(2*(-0.5*BLR.model.lrm$deviance[1])/nrow(dat.devel)))
    
    D <- as.numeric(BLR.model.lrm$stats["Dxy"])
    C <- as.numeric(BLR.model.lrm$stats["C"])
    
    
    ### Return heuristic and actual shrinkage factors
    return(c("S.pop" = as.numeric(S.pop), 
             "R2.CS.app" = R2.CS.app, 
             "R2.NAGEL.app" = R2.NAGEL.app,
             "D" = D,
             "C" = C,
             "LR" = LR,
             "L.full" = loglik.full,
             "L.null" = loglik.null,
             "P" = P))
  }
}

###
### 3.1) Function to calculate required sample size for a given model by calculating what R2 is in a population level model
###
calc.req.sampsize <- function(coef, beta0, covar_matrix, P, S){
  
  ### Fit a BLR on a large dataset so we will not suffer from overfitting
  dat.devel.temp <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                       beta0, #vector of intercepts
                       covar_matrix, #covariance matrix for patient data
                       P, #number of variables we will adjust for in model (variables we have access to)
                       N = 1000000)
  
  BLR.model <- glm(Y.fac ~ . , data = dplyr::select(dat.devel.temp, - c(Y)), family = binomial(link = "logit")) 
  #BLR.model.null <- glm(Y.fac ~ 1 , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
  
  R2.CS <- 1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel.temp))

  nreq <- pmsampsize(type = "b", 
             rsquared = R2.CS, 
             parameters = length(BLR.model$coefficients), 
             shrinkage = S, 
             prevalence = (sum(dat.devel.temp$Y == 1)/nrow(dat.devel.temp)))$results_table["Criteria 1", "Samp_size"]
  
  return(as.numeric(nreq))
  
}


###
### 3.2) Function to calculate required sample size for a given model by calculating what R2 is in datasets of the sample size
### that will be used in the model
###
calc.req.sampsize.realistic.V1 <- function(coef, beta0, covar_matrix, n.estR2){
  
  ### Fit a BLR on a large dataset so we will not suffer from overfitting
  dat.devel.pop <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                           beta0, #vector of intercepts
                           covar_matrix, #covariance matrix for patient data
                           N = 1000000)
  
  
  BLR.model.pop <- glm(Y.fac ~ . , data = dplyr::select(dat.devel.pop, - c(Y)), family = binomial(link = "logit")) 
  #BLR.model.null <- glm(Y.fac ~ 1 , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
  
  R2.CS.pop <- 1 - exp(-(BLR.model.pop$null.deviance - BLR.model.pop$deviance)/nrow(dat.devel.pop))
  
  nreq.pop <- as.numeric(pmsampsize(type = "b", 
                                    rsquared = R2.CS.pop, 
                                    parameters = length(BLR.model.pop$coefficients), 
                                    shrinkage = 0.9, 
                                    prevalence = (sum(dat.devel.pop$Y == 1)/nrow(dat.devel.pop)))$results_table["Criteria 1", "Samp_size"])
  #   nreq.pop
  #   R2.CS.pop
  
  R2.CS.vec <- rep(NA, n.estR2)
  S_VH.vec <- rep(NA, n.estR2)
  R2.CS.vec.adj1 <- rep(NA, n.estR2)
  R2.CS.vec.adj2 <- rep(NA, n.estR2)
  
  ### Loop through and calculate R2.CS for model's developed in datasets of size nreq
  for (j in 1:n.estR2){
    
    print(paste("calc n",j, Sys.time(), sep = " "))
    
    ### Create dataset
    dat.devel <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                         beta0, #vector of intercepts
                         covar_matrix, #covariance matrix for patient data
                         N = nreq.pop)
    
    ### Fit BLR
    BLR.model <- glm(Y.fac ~ . , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
    
    ### Calculate likelihood ratio statistic
    LR <- BLR.model$null.deviance - BLR.model$deviance
    
    ### Define P as number of predictor parameters estimated
    P <- length(BLR.model$coefficients)
    
    ### Calculate heuristic shrinkage factor
    #S_VH.vec[j] <- as.numeric(1 - P/LR)
    
    S_VH.vec[j] <- as.numeric(1 - length(BLR.model$coefficients)/(BLR.model$null.deviance - BLR.model$deviance))
    #R2.CS.vec[j] <- 1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel))
    #R2.CS.vec.adj1[j] <- 0.9*(1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel)))
    R2.CS.vec.adj2[j] <- S_VH.vec[j]*(1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel)))
  }
  
  #   mean(S_VH.vec)
  #   mean(R2.CS.vec)
  #   mean(R2.CS.vec.adj1)
  #   mean(R2.CS.vec.adj2)
  
  ### Assign actual R2.CS
  R2.CS.adj <- mean(R2.CS.vec.adj2)
  
  ### Calculate required sample size
  nreq <- as.numeric(pmsampsize(type = "b", 
                                rsquared = R2.CS.adj, 
                                parameters = length(BLR.model.pop$coefficients), 
                                shrinkage = 0.9, 
                                prevalence = (sum(dat.devel.pop$Y == 1)/nrow(dat.devel.pop)))$results_table["Criteria 1", "Samp_size"])
  
  return(c("nreq" = nreq, "nreq.pop" = nreq.pop))
  
}


###
### 3.3) Function to calculate required sample size for a given model by calculating what R2 is in datasets of the sample size
### that will be used in the model. This adjusts R2.CS_app to R2.CS_adj using an actual value of S (validated in population),
### rather than S_VH
###
calc.req.sampsize.realistic.V2 <- function(coef, beta0, covar_matrix, n.estR2){
 
  ### Fit a BLR on a large dataset so we will not suffer from overfitting
  dat.devel.pop <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                           beta0, #vector of intercepts
                           covar_matrix, #covariance matrix for patient data
                           N = 1000000)
  
  
  BLR.model.pop <- glm(Y.fac ~ . , data = dplyr::select(dat.devel.pop, - c(Y)), family = binomial(link = "logit")) 
  #BLR.model.null <- glm(Y.fac ~ 1 , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
  
  R2.CS.pop <- 1 - exp(-(BLR.model.pop$null.deviance - BLR.model.pop$deviance)/nrow(dat.devel.pop))
  
  nreq.pop <- as.numeric(pmsampsize(type = "b", 
                                    rsquared = R2.CS.pop, 
                                    parameters = length(BLR.model.pop$coefficients), 
                                    shrinkage = 0.9, 
                                    prevalence = (sum(dat.devel.pop$Y == 1)/nrow(dat.devel.pop)))$results_table["Criteria 1", "Samp_size"])
  #   nreq.pop
  #   R2.CS.pop
  
  R2.CS.vec <- rep(NA, n.estR2)
  S.vec <- rep(NA, n.estR2)
  R2.CS.vec.adj1 <- rep(NA, n.estR2)
  R2.CS.vec.adj2 <- rep(NA, n.estR2)
  
  ### Loop through and calculate R2.CS for model's developed in datasets of size nreq
  for (j in 1:n.estR2){
    
    print(paste("calc n",j, Sys.time(), sep = " "))
    
    ### Create dataset
    dat.devel <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                         beta0, #vector of intercepts
                         covar_matrix, #covariance matrix for patient data
                         N = nreq.pop)
    
    ### Fit BLR
    BLR.model <- glm(Y.fac ~ . , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
    
    ### Generate predicted risks in validaiton dataset
    lp.valid <- predict(BLR.model, newdata = dat.devel.pop)
    
    ### Add to validation dataset
    dat.devel.pop$lp <- lp.valid
    
    ### Fit recalibration model to calculate shrinkage factor
    model.valid <- glm(Y.fac ~ lp, data = dat.devel.pop, family = binomial(link = "logit")) 
    
    ### Return the slope as the shrinkage factor
    S.vec[j] <- as.numeric(model.valid$coefficients["lp"])
    #R2.CS.vec[j] <- 1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel))
    #R2.CS.vec.adj1[j] <- 0.9*(1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel)))
    R2.CS.vec.adj2[j] <- S.vec[j]*(1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel)))
  }
  
  #   mean(S_VH.vec)
  #   mean(R2.CS.vec)
  #   mean(R2.CS.vec.adj1)
  #   mean(R2.CS.vec.adj2)
  
  ### Assign actual R2.CS
  R2.CS.adj <- mean(R2.CS.vec.adj2)
  
  ### Calculate required sample size
  nreq <- as.numeric(pmsampsize(type = "b", 
                                rsquared = R2.CS.adj, 
                                parameters = length(BLR.model.pop$coefficients), 
                                shrinkage = 0.9, 
                                prevalence = (sum(dat.devel.pop$Y == 1)/nrow(dat.devel.pop)))$results_table["Criteria 1", "Samp_size"])
  
  return(c("nreq" = nreq, "nreq.pop" = nreq.pop, "R2.CS.pop" = R2.CS.pop))
  
}
#test <- calc.req.sampsize.realistic.V2(coef = c(0.5, 0.5, 0.5, 0.5, 0.5), beta0 = 0, n.estR2 = 10)


###
### 3.3) Function to calculate required sample size for a given model by calculating what R2 is in datasets of the sample size
### that will be used in the model. This adjusts R2.CS_app to R2.CS_adj using a fixed value of 0.9.
###
calc.req.sampsize.realistic.V3 <- function(coef, beta0, covar_matrix, n.estR2){
  
  ### Fit a BLR on a large dataset so we will not suffer from overfitting
  dat.devel.pop <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                           beta0, #vector of intercepts
                           covar_matrix, #covariance matrix for patient data
                           N = 1000000)
  
  
  BLR.model.pop <- glm(Y.fac ~ . , data = dplyr::select(dat.devel.pop, - c(Y)), family = binomial(link = "logit")) 
  #BLR.model.null <- glm(Y.fac ~ 1 , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
  
  R2.CS.pop <- 1 - exp(-(BLR.model.pop$null.deviance - BLR.model.pop$deviance)/nrow(dat.devel.pop))
  
  nreq.pop <- as.numeric(pmsampsize(type = "b", 
                                    rsquared = R2.CS.pop, 
                                    parameters = length(BLR.model.pop$coefficients), 
                                    shrinkage = 0.9, 
                                    prevalence = (sum(dat.devel.pop$Y == 1)/nrow(dat.devel.pop)))$results_table["Criteria 1", "Samp_size"])
  #   nreq.pop
  #   R2.CS.pop
  
  R2.CS.vec <- rep(NA, n.estR2)
  S.vec <- rep(NA, n.estR2)
  R2.CS.vec.adj1 <- rep(NA, n.estR2)
  R2.CS.vec.adj2 <- rep(NA, n.estR2)
  
  ### Loop through and calculate R2.CS for model's developed in datasets of size nreq
  for (j in 1:n.estR2){
    
    print(paste("calc n",j, Sys.time(), sep = " "))
    
    ### Create dataset
    dat.devel <- BLR.DGM(coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                         beta0, #vector of intercepts
                         covar_matrix, #covariance matrix for patient data
                         N = nreq.pop)
    
    ### Fit BLR
    BLR.model <- glm(Y.fac ~ . , data = dplyr::select(dat.devel, - c(Y)), family = binomial(link = "logit")) 
    
    #R2.CS.vec[j] <- 1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel))
    R2.CS.vec.adj1[j] <- 0.9*(1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel)))
    #R2.CS.vec.adj2[j] <- S.vec[j]*(1 - exp(-(BLR.model$null.deviance - BLR.model$deviance)/nrow(dat.devel)))
  }
  
  #   mean(S_VH.vec)
  #   mean(R2.CS.vec)
  #   mean(R2.CS.vec.adj1)
  #   mean(R2.CS.vec.adj2)
  
  ### Assign actual R2.CS
  R2.CS.adj <- mean(R2.CS.vec.adj1)
  
  ### Calculate required sample size
  nreq <- as.numeric(pmsampsize(type = "b", 
                                rsquared = R2.CS.adj, 
                                parameters = length(BLR.model.pop$coefficients), 
                                shrinkage = 0.9, 
                                prevalence = (sum(dat.devel.pop$Y == 1)/nrow(dat.devel.pop)))$results_table["Criteria 1", "Samp_size"])
  
  return(c("nreq" = nreq, "nreq.pop" = nreq.pop, "R2.CS.pop" = R2.CS.pop))
  
}
