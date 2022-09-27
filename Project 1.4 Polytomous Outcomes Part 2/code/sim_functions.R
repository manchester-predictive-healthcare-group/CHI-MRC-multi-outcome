####################################################################
### 1.1) Function to generate the data using the mulitnomial DGM ###
####################################################################

### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
generate.data.DGM.mult <- function(K, #number of outcome categories
                                   P, #number of predictors
                                   coef, #matrix of coefficients
                                   beta0, #vector of intercepts
                                   N){
  
  #   N <- 100
  #   K <- 4
  #   P <- 5
  #   coef <- rbind(c(0.75, 0.75, 0.75, 0.75, 0.75),
  #                 c(0.5, 0.5, 0.5, 0.5, 0.5),
  #                 c(0.25, 0.25, 0.25, 0.25, 0.25))
  #   beta0 <- c(0, 0, 0)
  
  stopifnot(length(beta0)+1 == K)
  stopifnot(nrow(coef)+1 == K)
  stopifnot(ncol(coef) == P)
  
  ### Create an input dataset of patient covariates
  ## Create empty dataframe
  dat.covar <- data.frame(matrix(, nrow = N, ncol = P))
  ## Add variable names
  colnames(dat.covar) <- paste("x", 1:P, sep = "")
  ## Add covariates
  for (i in 1:P){
    dat.covar[, i] <- rnorm(N, 0, 1)
  }
  
  ### Now add the probabilities of having each outcome event
  dat.probs <- data.frame(matrix(, nrow = N, ncol = K))
  ## Add variable names
  colnames(dat.probs) <- paste("p", 1:K, sep = "")
  ## Calculate the linear predictor for each equation in multinomial regression (P(Y=K)/(P(Y=1)))
  LP <- vector("list", K-1)
  for (i in 1:length(LP)){
    LP[[i]] <- exp(beta0[i] + as.matrix(dat.covar)%*%coef[i, ])
  }
  
  ## Use these to calculate probabilities
  # First outcome category
  dat.probs[,1] <- 1/(1 + rowSums(do.call("cbind", LP)))
  # Outcome categories 2 to K
  for (i in 2:K){
    dat.probs[,i] <- (LP[[i-1]])/(1 + rowSums(do.call("cbind", LP)))
  }
  
  ### Use these outcome probabilities to generate the outcomes from a multinomial distribution
  ## Define a function to calculate multinomial outcomes for a given probability vector
  generate.multinom <- function(probs.in){rmultinom(n = 1, size = 1, prob = probs.in)}
  ## Apply this to each row of the dataset
  multinom.outcomes <- t(apply(dat.probs, 1, function(x) generate.multinom(x)))
  
  ## Now for each row, need to convert this into a number depending on which row is equal to 1, and add it to a vector
  y.vec <- apply(multinom.outcomes, 1, function(x) which(x > 0))
  
  ## Combine this with dat.covar, into a development dataset
  dat.covar$Y <- y.vec
  dat.covar$Y.fac <- as.factor(dat.covar$Y)
  
  return(dat.covar)
}



############################################################################
### 1.2) Function to generate the data using the sequential logistic DGM ###
############################################################################

### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
generate.data.DGM.seqlog <- function(K, #number of outcome categories
                                     P, #number of predictors
                                     coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                     beta0, #vector of intercepts
                                     N){
  
  #   N <- 100
  #   K <- 4
  #   P <- 5
  #   coef <- rbind(c(0.75, 0.75, 0.75, 0.75, 0.75),
  #                 c(0.5, 0.5, 0.5, 0.5, 0.5),
  #                 c(0.25, 0.25, 0.25, 0.25, 0.25))
  #   beta0 <- c(0, 0, 0)
  
  stopifnot(length(beta0)+1 == K)
  stopifnot(nrow(coef)+1 == K)
  stopifnot(ncol(coef) == P)
  
  ### Create an input dataset of patient covariates
  ## Create empty dataframe
  dat.covar <- data.frame(matrix(, nrow = N, ncol = P))
  ## Add variable names
  colnames(dat.covar) <- paste("x", 1:P, sep = "")
  ## Add covariates
  for (i in 1:P){
    dat.covar[, i] <- rnorm(N, 0, 1)
  }
  
  ### Now add the probabilities of having each outcome event
  dat.probs <- data.frame(matrix(, nrow = N, ncol = K))
  ## Add variable names
  colnames(dat.probs) <- paste("p", 1:K, sep = "")
  ## Calculate the linear predictor for each equation in multinomial regression (P(Y=K)/(P(Y=1)))
  LP <- vector("list", K-1)
  for (i in 1:length(LP)){
    LP[[i]] <- exp(beta0[i] + as.matrix(dat.covar)%*%coef[i, ])
  }
  
  ## Use these to calculate probabilities
  ## First, calculate the probability of an event, for each sequential model
  dat.probs.individual.models <- vector("list", K-1)
  for (i in 1:(K-1)){
    dat.probs.individual.models[[i]] <- LP[[i]]/(1+LP[[i]])
  }
  
  ## Now calculate the probability of having each event
  # First outcome category
  dat.probs[,1] <- LP[[1]]/(1+LP[[1]])
  # Outcome categories 2 to K-1
  # It's the probability of i-1, divided by the probability of i-1 from the last seq log model, multiplied by the prob of
  # NOT having i-1 from the last seq log model, multiplied by the prob of having i from the next seqlog model
  for (i in 2:(K-1)){
    dat.probs[,i] <- (dat.probs[ , (i-1)])*((1+LP[[(i-1)]])/LP[[(i-1)]])*(1 - LP[[(i-1)]]/(1 + LP[[(i-1)]]))*LP[[i]]/(1+LP[[i]])
  }
  # Outcome category K
  dat.probs[,K] <- (dat.probs[ , (K-1)])*((1+LP[[(K-1)]])/LP[[(K-1)]])*(1 - LP[[(K-1)]]/(1 + LP[[(K-1)]]))
  
  ### Use these outcome probabilities to generate the outcomes from a multinomial distribution
  ## Define a function to calculate multinomial outcomes for a given probability vector
  generate.multinom <- function(probs.in){rmultinom(n = 1, size = 1, prob = probs.in)}
  ## Apply this to each row of the dataset
  multinom.outcomes <- t(apply(dat.probs, 1, function(x) generate.multinom(x)))
  
  ## Now for each row, need to convert this into a number depending on which row is equal to 1, and add it to a vector
  y.vec <- apply(multinom.outcomes, 1, function(x) which(x > 0))
  
  ## Combine this with dat.covar, into a development dataset
  dat.covar$Y <- y.vec
  dat.covar$Y.fac <- as.factor(dat.covar$Y)
  
  return(dat.covar)
}



##############################################################################################
### 2.1A) Function to fit multinomial model, and asses calibration in a validation dataset ###
##############################################################################################
calibrate.model.multinomial <- function(dat.devel, dat.valid){
  
  ### Fit a multinomial logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
  
  ### Fit the multinomial logistic regression model
  fit.multinomial <- vgam(Y.fac ~ ., family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
  
  ### Generate risks for individuals in the validation cohort
  ## First get the linear predictors for each submodel (log(P(Y=K)/P(Y=1)))
  pred.lp <- predict(fit.multinomial, newdata = dat.valid)
  ## Convert these into predicted probabilities
  ## For this simulation we only need predicted probability of outcome category 1
  pred.Y1 <- 1/(1 + rowSums(exp(pred.lp)))
  
  ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
  
  ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1.  Dataset also needs predicted probability and  LP in it.
  ## Create dataset
  dat.valid <- dat.valid
  ## Add predicted probability
  dat.valid$pred.Y1 <- pred.Y1
  ## Put predicted probability on scale of the linear predictor for a binary logistic regression
  dat.valid$LP.Y1 <- log(dat.valid$pred.Y1/(1-dat.valid$pred.Y1))
  ## Add binary definition of Y1
  dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                    Y.fac != 1 ~ 0))
  dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
  dat.valid <- arrange(dat.valid, pred.Y1)
  
  ### Now to do the calibration
  ### Going to do a bunch of different calibrations for comparison
  ## Create dataset to store all the output
  dat.valid.out <- dat.valid
  
  ### 1) Standard calibration plot, on LP
  ## Regress Y1.bin on LP.Y1
  fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
  
  ## Generate "predicted observed values, based off the calibration model
  dat.valid.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.valid, type = "response")
  
  #   ### 3) Standard calibration plot, on LP, with fractional polynomials
  #   ## Regress Y1.bin on LP.Y1
  #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
  #   
  #   ## Generate "predicted observed values, based off the calibration model
  #   dat.valid.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = dat.valid, type = "response")
  
  ### 6) Calibration plots with loess smoother, on probabilities
  ## Regress Y1.bin on pred.Y1
  fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
  
  ## Add to dataset
  dat.valid.out$pred.obs.Y1.loess <- fit.calib.loess$fitted
  
  ### Now to calculate discrimination
  Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
  
  ### Return output
  #return(select(dat.valid, -c(Y, Y.fac, LP.Y1, Y1.bin)))
  return(list(select(dat.valid.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
         fit.calib.lin$coefficients["LP.Y1"], Cstat))
}


###############################################################################################
### 2.1B) Function to fit multinomial model, and asses calibration in a validation dataset, ###
### but the oberved risks are estimated for a fixed vector of predicted risks, rather than  ###
### for the risks of individuals in the validation cohort                                   ###
###############################################################################################
calibrate.model.multinomial.predfix <- function(dat.devel, dat.valid, pred.eval){
  
  ### Fit a multinomial logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
  
  ### Fit the multinomial logistic regression model
  fit.multinomial <- vgam(Y.fac ~ ., family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
  
  ### Generate risks for individuals in the validation cohort
  ## First get the linear predictors for each submodel (log(P(Y=K)/P(Y=1)))
  pred.lp <- predict(fit.multinomial, newdata = dat.valid)
  ## Convert these into predicted probabilities
  ## For this simulation we only need predicted probability of outcome category 1
  pred.Y1 <- 1/(1 + rowSums(exp(pred.lp)))
  
  ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
  
  ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1.  Dataset also needs predicted probability and  LP in it.
  ## Create dataset
  dat.valid <- dat.valid
  ## Add predicted probability
  dat.valid$pred.Y1 <- pred.Y1
  ## Put predicted probability on scale of the linear predictor for a binary logistic regression
  dat.valid$LP.Y1 <- log(dat.valid$pred.Y1/(1-dat.valid$pred.Y1))
  ## Add binary definition of Y1
  dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                    Y.fac != 1 ~ 0))
  dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
  ## Arrange by predicted risk
  dat.valid <- arrange(dat.valid, pred.Y1)
  
  ### Now to do the calibration
  ### Going to do a bunch of different calibrations for comparison
  ## Create dataset to store all the output
  dat.valid.out <- dat.valid
  
  ## Also want to create a dataset, which contains "pred.eval" and "LP.ped.eval". So that I can pass it into the newdata arguments, 
  ## "LP.pred.eval" must be called "LP.Y1".
  dat.pred.eval.out <- data.frame("pred.Y1" = pred.eval, "LP.Y1" = log(pred.eval/(1-pred.eval)))
  
  ### 1) Standard calibration plot, on LP
  ## Regress Y1.bin on LP.Y1
  fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
  
  ## Generate "predicted observed values, based off the calibration model
  dat.pred.eval.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.pred.eval.out, type = "response")
  
  #   ### 3) Standard calibration plot, on LP, with fractional polynomials
  #   ## Regress Y1.bin on LP.Y1
  #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
  #   
  #   ## Generate "predicted observed values, based off the calibration model
  #   dat.pred.eval.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = dat.pred.eval.out, type = "response")
  
  
  ### 6) Calibration plots with loess smoother, on probabilities
  ## Regress Y1.bin on pred.Y1
  fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
  
  ## Generate "predicted observed values, based off the calibration model
  dat.pred.eval.out$pred.obs.Y1.loess <- predict(fit.calib.loess, newdata = dat.pred.eval.out$pred.Y1)
  
  
  ### Now to calculate discrimination
  Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
  
  ### Return output
  return(list(select(dat.pred.eval.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
              fit.calib.lin$coefficients["LP.Y1"], 
              Cstat))
}




#############################################################################################################################
### 2.2A) Function to fit multinomial model, with restricted cubic splines and asses calibration in a validation dataset, ###
### but the oberved risks are estimated for a fixed vector of predicted risks, rather than  ###
### for the risks of individuals in the validation cohort                                   ###
###############################################################################################


calibrate.model.multinomial.rcs <- function(dat.devel, dat.valid, n.knot.in){
  
  ### Fit a multinomial logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
  
  ### Fit the multinomial logistic regression model
  if (P.sim == 1){
    Knots.1 <- rcspline.eval(dat.devel$x1, nk = n.knot.in, knots.only = TRUE)
    fit.multinomial <- vgam(Y.fac ~ rcs(x1, Knots.1), family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
  } else if (P.sim == 3){
    Knots.1 <- rcspline.eval(dat.devel$x1, nk = n.knot.in, knots.only = TRUE)
    Knots.2 <- rcspline.eval(dat.devel$x2, nk = n.knot.in, knots.only = TRUE)
    Knots.3 <- rcspline.eval(dat.devel$x3, nk = n.knot.in, knots.only = TRUE)
    fit.multinomial <- vgam(Y.fac ~ rcs(x1, Knots.1) + rcs(x2, Knots.2) + rcs(x3, Knots.3), 
                            family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
  } else if (P.sim == 5){
    Knots.1 <- rcspline.eval(dat.devel$x1, nk = n.knot.in, knots.only = TRUE)
    Knots.2 <- rcspline.eval(dat.devel$x2, nk = n.knot.in, knots.only = TRUE)
    Knots.3 <- rcspline.eval(dat.devel$x3, nk = n.knot.in, knots.only = TRUE)
    Knots.4 <- rcspline.eval(dat.devel$x4, nk = n.knot.in, knots.only = TRUE)
    Knots.5 <- rcspline.eval(dat.devel$x5, nk = n.knot.in, knots.only = TRUE)
    fit.multinomial <- vgam(Y.fac ~ rcs(x1, Knots.1) + rcs(x2, Knots.2) + rcs(x3, Knots.3) + rcs(x4, Knots.4) + rcs(x5, Knots.5), 
                            family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
  } else if (P.sim == 10){
    Knots.1 <- rcspline.eval(dat.devel$x1, nk = n.knot.in, knots.only = TRUE)
    Knots.2 <- rcspline.eval(dat.devel$x2, nk = n.knot.in, knots.only = TRUE)
    Knots.3 <- rcspline.eval(dat.devel$x3, nk = n.knot.in, knots.only = TRUE)
    Knots.4 <- rcspline.eval(dat.devel$x4, nk = n.knot.in, knots.only = TRUE)
    Knots.5 <- rcspline.eval(dat.devel$x5, nk = n.knot.in, knots.only = TRUE)
    Knots.6 <- rcspline.eval(dat.devel$x6, nk = n.knot.in, knots.only = TRUE)
    Knots.7 <- rcspline.eval(dat.devel$x7, nk = n.knot.in, knots.only = TRUE)
    Knots.8 <- rcspline.eval(dat.devel$x8, nk = n.knot.in, knots.only = TRUE)
    Knots.9 <- rcspline.eval(dat.devel$x9, nk = n.knot.in, knots.only = TRUE)
    Knots.10 <- rcspline.eval(dat.devel$x10, nk = n.knot.in, knots.only = TRUE)
    fit.multinomial <- vgam(Y.fac ~ rcs(x1, Knots.1) + rcs(x2, Knots.2) + rcs(x3, Knots.3) + rcs(x4, Knots.4) + rcs(x5, Knots.5) +
                              rcs(x6, Knots.6) + rcs(x7, Knots.7) + rcs(x8, Knots.8) + rcs(x9, Knots.9) + rcs(x10, Knots.10), 
                            family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
  }
  
  ### Generate risks for individuals in the validation cohort
  ## First get the linear predictors for each submodel (log(P(Y=K)/P(Y=1)))
  pred.lp <- predict(fit.multinomial, newdata = dat.valid)
  
  ## Convert these into predicted probabilities
  ## For this simulation we only need predicted probability of outcome category 1
  pred.Y1 <- 1/(1 + rowSums(exp(pred.lp)))
  
  ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
  
  ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1.  Dataset also needs predicted probability and  LP in it.
  ## Create dataset
  dat.valid <- dat.valid
  ## Add predicted probability
  dat.valid$pred.Y1 <- pred.Y1
  ## Put predicted probability on scale of the linear predictor for a binary logistic regression
  dat.valid$LP.Y1 <- log(dat.valid$pred.Y1/(1-dat.valid$pred.Y1))
  ## Add binary definition of Y1
  dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                    Y.fac != 1 ~ 0))
  dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
  ## Arrange by predicted risk
  dat.valid <- arrange(dat.valid, pred.Y1)
  
  ### Now to do the calibration
  ### Going to do a bunch of different calibrations for comparison
  ## Create dataset to store all the output
  dat.valid.out <- dat.valid
  
  ### 1) Standard calibration plot, on LP
  ## Regress Y1.bin on LP.Y1
  fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
  
  ## Generate "predicted observed values, based off the calibration model
  dat.valid.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.valid, type = "response")
  
  #   ### 3) Standard calibration plot, on LP, with fractional polynomials
  #   ## Regress Y1.bin on LP.Y1
  #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
  #   
  #   ## Generate "predicted observed values, based off the calibration model
  #   dat.valid.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = dat.valid, type = "response")
  
  ### 6) Calibration plots with loess smoother, on probabilities
  ## Regress Y1.bin on pred.Y1
  fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
  
  ## Add to dataset
  dat.valid.out$pred.obs.Y1.loess <- fit.calib.loess$fitted
  
  ### Now to calculate discrimination
  Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
  
  ### Return output
  #return(select(dat.valid, -c(Y, Y.fac, LP.Y1, Y1.bin)))
  return(list(select(dat.valid.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
         fit.calib.lin$coefficients["LP.Y1"], Cstat))
}


###
### OLD VERSION WHERE KNOTS ARE NOT MAINTAINED BTWEEN DEVELOPMENT AND VALIDATION DATASETS
###
# calibrate.model.multinomial.rcs <- function(dat.devel, dat.valid, n.knot.in){
#   
#   ### Fit a multinomial logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
#   
#   ### Fit the multinomial logistic regression model
#   if (P.sim == 1){
#     fit.multinomial <- vgam(Y.fac ~ rcs(x1, n.knot.in), family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
#   } else if (P.sim == 3){
#     fit.multinomial <- vgam(Y.fac ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in), 
#                             family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
#   } else if (P.sim == 5){
#     fit.multinomial <- vgam(Y.fac ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in) + rcs(x4, n.knot.in) + rcs(x5, n.knot.in), 
#                             family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
#   } else if (P.sim == 10){
#     fit.multinomial <- vgam(Y.fac ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in) + rcs(x4, n.knot.in) + rcs(x5, n.knot.in) +
#                               rcs(x6, n.knot.in) + rcs(x7, n.knot.in) + rcs(x8, n.knot.in) + rcs(x9, n.knot.in) + rcs(x10, n.knot.in), 
#                             family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
#   }
#   
#   ### Generate risks for individuals in the validation cohort
#   ## First get the linear predictors for each submodel (log(P(Y=K)/P(Y=1)))
#   pred.lp <- predict(fit.multinomial, newdata = dat.valid)
#   
#   ## Convert these into predicted probabilities
#   ## For this simulation we only need predicted probability of outcome category 1
#   pred.Y1 <- 1/(1 + rowSums(exp(pred.lp)))
#   
#   ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
#   
#   ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1.  Dataset also needs predicted probability and  LP in it.
#   ## Create dataset
#   dat.valid <- dat.valid
#   ## Add predicted probability
#   dat.valid$pred.Y1 <- pred.Y1
#   ## Put predicted probability on scale of the linear predictor for a binary logistic regression
#   dat.valid$LP.Y1 <- log(dat.valid$pred.Y1/(1-dat.valid$pred.Y1))
#   ## Add binary definition of Y1
#   dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
#                                                     Y.fac != 1 ~ 0))
#   dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
#   ## Arrange by predicted risk
#   dat.valid <- arrange(dat.valid, pred.Y1)
#   
#   ### Now to do the calibration
#   ### Going to do a bunch of different calibrations for comparison
#   ## Create dataset to store all the output
#   dat.valid.out <- dat.valid
#   
#   ### 1) Standard calibration plot, on LP
#   ## Regress Y1.bin on LP.Y1
#   fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
#   
#   ## Generate "predicted observed values, based off the calibration model
#   dat.valid.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.valid, type = "response")
#   
#   #   ### 3) Standard calibration plot, on LP, with fractional polynomials
#   #   ## Regress Y1.bin on LP.Y1
#   #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
#   #   
#   #   ## Generate "predicted observed values, based off the calibration model
#   #   dat.valid.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = dat.valid, type = "response")
#   
#   ### 6) Calibration plots with loess smoother, on probabilities
#   ## Regress Y1.bin on pred.Y1
#   fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
#   
#   ## Add to dataset
#   dat.valid.out$pred.obs.Y1.loess <- fit.calib.loess$fitted
#   
#   ### Now to calculate discrimination
#   Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
#   
#   ### Return output
#   #return(select(dat.valid, -c(Y, Y.fac, LP.Y1, Y1.bin)))
#   return(list(select(dat.valid.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
#               fit.calib.lin$coefficients["LP.Y1"], Cstat))
# }



#############################################################################################################################
### 2.2B) Function to fit multinomial model, with restricted cubic splines and asses calibration in a validation dataset, ###
### but the oberved risks are estimated for a fixed vector of predicted risks, rather than  ###
### for the risks of individuals in the validation cohort                                   ###
###############################################################################################
###
### OLD VERSION WHERE KNOTS ARE NOT MAINTAINED BTWEEN DEVELOPMENT AND VALIDATION DATASETS
###
# calibrate.model.multinomial.predfix.rcs <- function(dat.devel, dat.valid, pred.eval, n.knot.in){
#   
#   ### Fit a multinomial logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
#   
#   ### Fit the multinomial logistic regression model
#   if (P.sim == 1){
#     fit.multinomial <- vgam(Y.fac ~ rcs(x1, n.knot.in), family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
#   } else if (P.sim == 3){
#     fit.multinomial <- vgam(Y.fac ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in), 
#                             family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
#   } else if (P.sim == 5){
#     fit.multinomial <- vgam(Y.fac ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in) + rcs(x4, n.knot.in) + rcs(x5, n.knot.in), 
#                             family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
#   } else if (P.sim == 10){
#     fit.multinomial <- vgam(Y.fac ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in) + rcs(x4, n.knot.in) + rcs(x5, n.knot.in) +
#                               rcs(x6, n.knot.in) + rcs(x7, n.knot.in) + rcs(x8, n.knot.in) + rcs(x9, n.knot.in) + rcs(x10, n.knot.in), 
#                             family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))
#   }
#   
#   ### Generate risks for individuals in the validation cohort
#   ## First get the linear predictors for each submodel (log(P(Y=K)/P(Y=1)))
#   pred.lp <- predict(fit.multinomial, newdata = dat.valid)
#   ## Convert these into predicted probabilities
#   ## For this simulation we only need predicted probability of outcome category 1
#   pred.Y1 <- 1/(1 + rowSums(exp(pred.lp)))
#   
#   ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
#   
#   ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1.  Dataset also needs predicted probability and  LP in it.
#   ## Create dataset
#   dat.valid <- dat.valid
#   ## Add predicted probability
#   dat.valid$pred.Y1 <- pred.Y1
#   ## Put predicted probability on scale of the linear predictor for a binary logistic regression
#   dat.valid$LP.Y1 <- log(dat.valid$pred.Y1/(1-dat.valid$pred.Y1))
#   ## Add binary definition of Y1
#   dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
#                                                     Y.fac != 1 ~ 0))
#   dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
#   ## Arrange by predicted risk
#   dat.valid <- arrange(dat.valid, pred.Y1)
#   
#   ### Now to do the calibration
#   ### Going to do a bunch of different calibrations for comparison
#   ## Create dataset to store all the output
#   dat.valid.out <- dat.valid
#   
#   ## Also want to create a dataset, which contains "pred.eval" and "LP.ped.eval". So that I can pass it into the newdata arguments, 
#   ## "LP.pred.eval" must be called "LP.Y1".
#   dat.pred.eval.out <- data.frame("pred.Y1" = pred.eval, "LP.Y1" = log(pred.eval/(1-pred.eval)))
#   
#   ### 1) Standard calibration plot, on LP
#   ## Regress Y1.bin on LP.Y1
#   fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
#   
#   ## Generate "predicted observed values, based off the calibration model
#   dat.pred.eval.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.pred.eval.out, type = "response")
#   
#   #   ### 3) Standard calibration plot, on LP, with fractional polynomials
#   #   ## Regress Y1.bin on LP.Y1
#   #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
#   #   
#   #   ## Generate "predicted observed values, based off the calibration model
#   #   dat.pred.eval.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = dat.pred.eval.out, type = "response")
#   
#   
#   ### 6) Calibration plots with loess smoother, on probabilities
#   ## Regress Y1.bin on pred.Y1
#   fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
#   
#   ## Generate "predicted observed values, based off the calibration model
#   dat.pred.eval.out$pred.obs.Y1.loess <- predict(fit.calib.loess, newdata = dat.pred.eval.out$pred.Y1)
#   
#   
#   ### Now to calculate discrimination
#   Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
#   
#   ### Return output
#   return(list(select(dat.pred.eval.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
#               fit.calib.lin$coefficients["LP.Y1"], 
#               Cstat))
# }



##################################################################################################
### 3.1A) Function to fit binary logistic model, and asses calibration in a validation dataset ###
##################################################################################################
calibrate.model.binary <- function(dat.devel, dat.valid){
  
  ### Start by transforming development dataset to have binary outcome definition for Y1
  dat.devel <- mutate(dat.devel, Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                    Y.fac != 1 ~ 0))
  dat.devel$Y1.bin <- as.factor(dat.devel$Y1.bin)
  
  ### Fit a binary logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
  
  ### Fit the binary logistic regression model
  fit.binary <- glm(Y1.bin ~ ., family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
  
  ### Generate risks for individuals in the validation cohort
  pred.Y1 <- predict(fit.binary, newdata = dat.valid, type = "response")
  ### Generate LP for individuals in the validation cohort
  LP.Y1 <- predict(fit.binary, newdata = dat.valid, type = "link")
  
  
  ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
  
  ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1. Dataset also needs predicted probability and LP in it.
  ## Add predicted probability
  dat.valid$pred.Y1 <- pred.Y1
  ## Add LP
  dat.valid$LP.Y1 <- LP.Y1
  ## Add binary definition of Y1
  dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                    Y.fac != 1 ~ 0))
  dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
  dat.valid <- arrange(dat.valid, pred.Y1)
  
  
  ### Now to do the calibration
  ### Going to do a bunch of different calibrations for comparison
  ## Create dataset to store all the output
  dat.valid.out <- dat.valid
  
  ### 1) Standard calibration plot, on LP
  ## Regress Y1.bin on LP.Y1
  fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
  
  ## Generate "predicted observed values, based off the calibration model
  dat.valid.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.valid, type = "response")
  
  #   ### 3) Standard calibration plot, on LP, with fractional polynomials
  #   ## Regress Y1.bin on LP.Y1
  #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
  #   
  #   ## Generate "predicted observed values, based off the calibration model
  #   dat.valid.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = pred.eval, type = "response")
  
  
  ### 6) Calibration plots with loess smoother, on probabilities
  ## Regress Y1.bin on pred.Y1
  fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
  
  ## Generate "predicted observed values, based off the calibration model
  dat.valid.out$pred.obs.Y1.loess <- predict(fit.calib.loess, newdata = dat.valid)
  
  
  ### Now to calculate discrimination
  Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
  
  ### Return output
  #return(select(dat.valid, -c(Y, Y.fac, LP.Y1, Y1.bin)))
  return(list(select(dat.valid.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
         fit.calib.lin$coefficients["LP.Y1"], Cstat))
}


###################################################################################################
### 3.1B) Function to fit binary logistic model, and asses calibration in a validation dataset, ###
### but the oberved risks are estimated for a fixed vector of predicted risks, rather than  ###
### for the risks of individuals in the validation cohort                                   ###
###############################################################################################
calibrate.model.binary.predfix <- function(dat.devel, dat.valid, pred.eval){
  
  ### Start by transforming development dataset to have binary outcome definition for Y1
  dat.devel <- mutate(dat.devel, Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                    Y.fac != 1 ~ 0))
  dat.devel$Y1.bin <- as.factor(dat.devel$Y1.bin)
  
  ### Fit a binary logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
  
  ### Fit the binary logistic regression model
  fit.binary <- glm(Y1.bin ~ ., family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
  
  ### Generate risks for individuals in the validation cohort
  pred.Y1 <- predict(fit.binary, newdata = dat.valid, type = "response")
  ### Generate LP for individuals in the validation cohort
  LP.Y1 <- predict(fit.binary, newdata = dat.valid, type = "link")
  
  
  ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
  
  ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1. Dataset also needs predicted probability and LP in it.
  ## Add predicted probability
  dat.valid$pred.Y1 <- pred.Y1
  ## Add LP
  dat.valid$LP.Y1 <- LP.Y1
  ## Add binary definition of Y1
  dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                    Y.fac != 1 ~ 0))
  dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
  dat.valid <- arrange(dat.valid, pred.Y1)
  
  ### Now to do the calibration
  ### Going to do a bunch of different calibrations for comparison
  ## Create dataset to store all the output
  dat.valid.out <- dat.valid
  
  ## Also want to create a dataset, which contains "pred.eval" and "LP.ped.eval". So that I can pass it into the newdata arguments, 
  ## "LP.pred.eval" must be called "LP.Y1".
  dat.pred.eval.out <- data.frame("pred.Y1" = pred.eval, "LP.Y1" = log(pred.eval/(1-pred.eval)))
  
  ### 1) Standard calibration plot, on LP
  ## Regress Y1.bin on LP.Y1
  fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
  
  ## Generate "predicted observed values, based off the calibration model
  dat.pred.eval.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.pred.eval.out, type = "response")
  
  #   ### 3) Standard calibration plot, on LP, with fractional polynomials
  #   ## Regress Y1.bin on LP.Y1
  #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
  #   
  #   ## Generate "predicted observed values, based off the calibration model
  #   dat.pred.eval.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = dat.pred.eval.out, type = "response")
  
  
  ### 6) Calibration plots with loess smoother, on probabilities
  ## Regress Y1.bin on pred.Y1
  fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
  
  ## Generate "predicted observed values, based off the calibration model
  dat.pred.eval.out$pred.obs.Y1.loess <- predict(fit.calib.loess, newdata = dat.pred.eval.out$pred.Y1)
  
  
  ### Now to calculate discrimination
  Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
  
  ### Return output
  return(list(select(dat.pred.eval.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
              fit.calib.lin$coefficients["LP.Y1"], 
              Cstat))
}


######################################################################################################
### 3.2A) Function to fit binary logistic  with rcs, and asses calibration in a validation dataset ###
######################################################################################################
calibrate.model.binary.rcs <- function(dat.devel, dat.valid, n.knot.in){

  ### Start by transforming development dataset to have binary outcome definition for Y1
  dat.devel <- mutate(dat.devel, Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                    Y.fac != 1 ~ 0))
  dat.devel$Y1.bin <- as.factor(dat.devel$Y1.bin)
  
  ### Fit a binary logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
  ### Fit the binary logistic regression model
  if (P.sim == 1){
    Knots.1 <- rcspline.eval(dat.devel$x1, nk = n.knot.in, knots.only = TRUE)
    fit.binary <- glm(Y1.bin ~ rcs(x1, Knots.1), 
                      family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
  } else if (P.sim == 3){
    Knots.1 <- rcspline.eval(dat.devel$x1, nk = n.knot.in, knots.only = TRUE)
    Knots.2 <- rcspline.eval(dat.devel$x2, nk = n.knot.in, knots.only = TRUE)
    Knots.3 <- rcspline.eval(dat.devel$x3, nk = n.knot.in, knots.only = TRUE)
    fit.binary <- glm(Y1.bin ~ rcs(x1, Knots.1) + rcs(x2, Knots.2) + rcs(x3, Knots.3), 
                      family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
  } else if (P.sim == 5){
    Knots.1 <- rcspline.eval(dat.devel$x1, nk = n.knot.in, knots.only = TRUE)
    Knots.2 <- rcspline.eval(dat.devel$x2, nk = n.knot.in, knots.only = TRUE)
    Knots.3 <- rcspline.eval(dat.devel$x3, nk = n.knot.in, knots.only = TRUE)
    Knots.4 <- rcspline.eval(dat.devel$x4, nk = n.knot.in, knots.only = TRUE)
    Knots.5 <- rcspline.eval(dat.devel$x5, nk = n.knot.in, knots.only = TRUE)
    fit.binary <- glm(Y1.bin ~ rcs(x1, Knots.1) + rcs(x2, Knots.2) + rcs(x3, Knots.3) + rcs(x4, Knots.4) + rcs(x5, Knots.5), 
                      family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
  } else if (P.sim == 10){
    Knots.1 <- rcspline.eval(dat.devel$x1, nk = n.knot.in, knots.only = TRUE)
    Knots.2 <- rcspline.eval(dat.devel$x2, nk = n.knot.in, knots.only = TRUE)
    Knots.3 <- rcspline.eval(dat.devel$x3, nk = n.knot.in, knots.only = TRUE)
    Knots.4 <- rcspline.eval(dat.devel$x4, nk = n.knot.in, knots.only = TRUE)
    Knots.5 <- rcspline.eval(dat.devel$x5, nk = n.knot.in, knots.only = TRUE)
    Knots.6 <- rcspline.eval(dat.devel$x6, nk = n.knot.in, knots.only = TRUE)
    Knots.7 <- rcspline.eval(dat.devel$x7, nk = n.knot.in, knots.only = TRUE)
    Knots.8 <- rcspline.eval(dat.devel$x8, nk = n.knot.in, knots.only = TRUE)
    Knots.9 <- rcspline.eval(dat.devel$x9, nk = n.knot.in, knots.only = TRUE)
    Knots.10 <- rcspline.eval(dat.devel$x10, nk = n.knot.in, knots.only = TRUE)
    fit.binary <- glm(Y1.bin ~ rcs(x1, Knots.1) + rcs(x2, Knots.2) + rcs(x3, Knots.3) + rcs(x4, Knots.4) + rcs(x5, Knots.5) + 
                        rcs(x6, Knots.6) + rcs(x7, Knots.7) + rcs(x8, Knots.8) + rcs(x9, Knots.9) + fp(x10, Knots.10), 
                      family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
  }
  
  ### Generate risks for individuals in the validation cohort
  pred.Y1 <- predict(fit.binary, newdata = dat.valid, type = "response")
  
  ### Generate LP for individuals in the validation cohort
  LP.Y1 <- predict(fit.binary, newdata = dat.valid, type = "link")
  
  
  ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
  
  ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1. Dataset also needs predicted probability and LP in it.
  ## Add predicted probability
  dat.valid$pred.Y1 <- pred.Y1
  ## Add LP
  dat.valid$LP.Y1 <- LP.Y1
  ## Add binary definition of Y1
  dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
                                                    Y.fac != 1 ~ 0))
  dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
  dat.valid <- arrange(dat.valid, pred.Y1)
  
  
  ### Now to do the calibration
  ### Going to do a bunch of different calibrations for comparison
  ## Create dataset to store all the output
  dat.valid.out <- dat.valid
  
  ### 1) Standard calibration plot, on LP
  ## Regress Y1.bin on LP.Y1
  fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
  
  ## Generate "predicted observed values, based off the calibration model
  dat.valid.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.valid, type = "response")
  
  #   ### 3) Standard calibration plot, on LP, with fractional polynomials
  #   ## Regress Y1.bin on LP.Y1
  #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
  #   
  #   ## Generate "predicted observed values, based off the calibration model
  #   dat.valid.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = pred.eval, type = "response")
  
  
  ### 6) Calibration plots with loess smoother, on probabilities
  ## Regress Y1.bin on pred.Y1
  fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
  
  ## Generate "predicted observed values, based off the calibration model
  dat.valid.out$pred.obs.Y1.loess <- predict(fit.calib.loess, newdata = dat.valid)
  
  
  ### Now to calculate discrimination
  Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
  
  ### Return output
  #return(select(dat.valid, -c(Y, Y.fac, LP.Y1, Y1.bin)))
  return(list(select(dat.valid.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
         fit.calib.lin$coefficients["LP.Y1"], Cstat))
}



###
### OLD VERSION WHERE KNOTS ARE NOT MAINTAINED BTWEEN DEVELOPMENT AND VALIDATION DATASETS
###
# calibrate.model.binary.rcs <- function(dat.devel, dat.valid, n.knot.in){
#   
#   ### Start by transforming development dataset to have binary outcome definition for Y1
#   dat.devel <- mutate(dat.devel, Y1.bin = case_when(Y.fac == 1 ~ 1,
#                                                     Y.fac != 1 ~ 0))
#   dat.devel$Y1.bin <- as.factor(dat.devel$Y1.bin)
#   
#   ### Fit a binary logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
#   ### Fit the binary logistic regression model
#   if (P.sim == 1){
#     fit.binary <- glm(Y1.bin ~ rcs(x1, n.knot.in), family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
#   } else if (P.sim == 3){
#     fit.binary <- glm(Y1.bin ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in), family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
#   } else if (P.sim == 5){
#     fit.binary <- glm(Y1.bin ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in) + rcs(x4, n.knot.in) + rcs(x5, n.knot.in), family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
#   } else if (P.sim == 10){
#     fit.binary <- glm(Y1.bin ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in) + rcs(x4, n.knot.in) + rcs(x5, n.knot.in) + 
#                         rcs(x6, n.knot.in) + rcs(x7, n.knot.in) + rcs(x8, n.knot.in) + rcs(x9, n.knot.in) + fp(x10, n.knot.in), family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
#   }
#   
#   ### Generate risks for individuals in the validation cohort
#   pred.Y1 <- predict(fit.binary, newdata = dat.valid, type = "response")
#   ### Generate LP for individuals in the validation cohort
#   LP.Y1 <- predict(fit.binary, newdata = dat.valid, type = "link")
#   
#   
#   ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
#   
#   ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1. Dataset also needs predicted probability and LP in it.
#   ## Add predicted probability
#   dat.valid$pred.Y1 <- pred.Y1
#   ## Add LP
#   dat.valid$LP.Y1 <- LP.Y1
#   ## Add binary definition of Y1
#   dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
#                                                     Y.fac != 1 ~ 0))
#   dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
#   dat.valid <- arrange(dat.valid, pred.Y1)
#   
#   
#   ### Now to do the calibration
#   ### Going to do a bunch of different calibrations for comparison
#   ## Create dataset to store all the output
#   dat.valid.out <- dat.valid
#   
#   ### 1) Standard calibration plot, on LP
#   ## Regress Y1.bin on LP.Y1
#   fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
#   
#   ## Generate "predicted observed values, based off the calibration model
#   dat.valid.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.valid, type = "response")
#   
#   #   ### 3) Standard calibration plot, on LP, with fractional polynomials
#   #   ## Regress Y1.bin on LP.Y1
#   #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
#   #   
#   #   ## Generate "predicted observed values, based off the calibration model
#   #   dat.valid.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = pred.eval, type = "response")
#   
#   
#   ### 6) Calibration plots with loess smoother, on probabilities
#   ## Regress Y1.bin on pred.Y1
#   fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
#   
#   ## Generate "predicted observed values, based off the calibration model
#   dat.valid.out$pred.obs.Y1.loess <- predict(fit.calib.loess, newdata = dat.valid)
#   
#   
#   ### Now to calculate discrimination
#   Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
#   
#   ### Return output
#   #return(select(dat.valid, -c(Y, Y.fac, LP.Y1, Y1.bin)))
#   return(list(select(dat.valid.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
#               fit.calib.lin$coefficients["LP.Y1"], Cstat))
# }


###################################################################################################
### 3.2B) Function to fit binary logistic model with rcs, and asses calibration in a validation dataset, ###
### but the oberved risks are estimated for a fixed vector of predicted risks, rather than  ###
### for the risks of individuals in the validation cohort                                   ###
###############################################################################################
###
### OLD VERSION WHERE KNOTS ARE NOT MAINTAINED BTWEEN DEVELOPMENT AND VALIDATION DATASETS
###
# calibrate.model.binary.predfix.rcs <- function(dat.devel, dat.valid, pred.eval, n.knot.in){
#   
#   ### Start by transforming development dataset to have binary outcome definition for Y1
#   dat.devel <- mutate(dat.devel, Y1.bin = case_when(Y.fac == 1 ~ 1,
#                                                     Y.fac != 1 ~ 0))
#   dat.devel$Y1.bin <- as.factor(dat.devel$Y1.bin)
#   
#   ### Fit a binary logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset, and predicted observed risks
#   ### Fit the binary logistic regression model
#   if (P.sim == 1){
#     fit.binary <- glm(Y1.bin ~ rcs(x1, n.knot.in), family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
#   } else if (P.sim == 3){
#     fit.binary <- glm(Y1.bin ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in), family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
#   } else if (P.sim == 5){
#     fit.binary <- glm(Y1.bin ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in) + rcs(x4, n.knot.in) + rcs(x5, n.knot.in), family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
#   } else if (P.sim == 10){
#     fit.binary <- glm(Y1.bin ~ rcs(x1, n.knot.in) + rcs(x2, n.knot.in) + rcs(x3, n.knot.in) + rcs(x4, n.knot.in) + rcs(x5, n.knot.in) + 
#                         rcs(x6, n.knot.in) + rcs(x7, n.knot.in) + rcs(x8, n.knot.in) + rcs(x9, n.knot.in) + fp(x10, n.knot.in), family = binomial(link = 'logit'), data = subset(dat.devel, select = -c(Y, Y.fac)))
#   }
#   
#   ### Generate risks for individuals in the validation cohort
#   pred.Y1 <- predict(fit.binary, newdata = dat.valid, type = "response")
#   ### Generate LP for individuals in the validation cohort
#   LP.Y1 <- predict(fit.binary, newdata = dat.valid, type = "link")
#   
#   
#   ### Now we want to produce a calibraiton plot for Y1, in the binary logistic sense
#   
#   ### To do this we need to create a new datset, with outcome Y1, which = 1 if Y = 1, but = 0 if Y != 1. Dataset also needs predicted probability and LP in it.
#   ## Add predicted probability
#   dat.valid$pred.Y1 <- pred.Y1
#   ## Add LP
#   dat.valid$LP.Y1 <- LP.Y1
#   ## Add binary definition of Y1
#   dat.valid <- mutate(dat.valid, Y1.bin = case_when(Y.fac == 1 ~ 1,
#                                                     Y.fac != 1 ~ 0))
#   dat.valid$Y1.bin.fac <- as.factor(dat.valid$Y1.bin)
#   dat.valid <- arrange(dat.valid, pred.Y1)
#   
#   ### Now to do the calibration
#   ### Going to do a bunch of different calibrations for comparison
#   ## Create dataset to store all the output
#   dat.valid.out <- dat.valid
#   
#   ## Also want to create a dataset, which contains "pred.eval" and "LP.ped.eval". So that I can pass it into the newdata arguments, 
#   ## "LP.pred.eval" must be called "LP.Y1".
#   dat.pred.eval.out <- data.frame("pred.Y1" = pred.eval, "LP.Y1" = log(pred.eval/(1-pred.eval)))
#   
#   ### 1) Standard calibration plot, on LP
#   ## Regress Y1.bin on LP.Y1
#   fit.calib.lin <- glm(Y1.bin ~ LP.Y1, family = binomial(link = 'logit'), data = dat.valid)
#   
#   ## Generate "predicted observed values, based off the calibration model
#   dat.pred.eval.out$pred.obs.Y1.lin <- predict(fit.calib.lin, newdata = dat.pred.eval.out, type = "response")
#   
#   #   ### 3) Standard calibration plot, on LP, with fractional polynomials
#   #   ## Regress Y1.bin on LP.Y1
#   #   fit.calib.rcs <- glm(Y1.bin ~ rcs(LP.Y1, 4), family = binomial(link = 'logit'), data = dat.valid)
#   #   
#   #   ## Generate "predicted observed values, based off the calibration model
#   #   dat.pred.eval.out$pred.obs.Y1.rcs <- predict(fit.calib.rcs, newdata = dat.pred.eval.out, type = "response")
#   
#   
#   ### 6) Calibration plots with loess smoother, on probabilities
#   ## Regress Y1.bin on pred.Y1
#   fit.calib.loess <- loess(Y1.bin ~ pred.Y1, data = dat.valid, method = 'loess')
#   
#   ## Generate "predicted observed values, based off the calibration model
#   dat.pred.eval.out$pred.obs.Y1.loess <- predict(fit.calib.loess, newdata = dat.pred.eval.out$pred.Y1)
#   
#   
#   ### Now to calculate discrimination
#   Cstat <- Cstat(dat.valid.out$pred.Y1, dat.valid.out$Y1.bin)
#   
#   ### Return output
#   return(list(select(dat.pred.eval.out, pred.Y1, pred.obs.Y1.lin, pred.obs.Y1.loess), 
#               fit.calib.lin$coefficients["LP.Y1"], 
#               Cstat))
# }


#####
### 4.1) Output predicted probabilities for a given set of input parameters
#####
### Write a function to output predicted probabilities
generate.outcome.probs.DGM.mult <- function(K, #number of outcome categories
                                            P, #number of predictors
                                            coef, #matrix of coefficients
                                            beta0, #vector of intercepts
                                            N){
  
  #   N <- 100
  #   K <- 4
  #   P <- 5
  #   coef <- c(0.5, -0.5, 0.5)
  #   beta0 <- c(0, 0, 0)
  
  #   N<- 500000
  #   K <- 5
  #   P <- 5
  #   beta0 <- c(0, 0, 0, 0)
  #   
  #   coef <- rbind(c(-0.5, -0.5, -0.5, -0.5, -0.5),
  #                     c(0.5, 0.5, 0.5, 0.5, 0.5),
  #                     c(0.5, 0.5, 0.5, 0.5, 0.5),
  #                     c(0.5, 0.5, 0.5, 0.5, 0.5))
  
  stopifnot(length(beta0)+1 == K)
  stopifnot(nrow(coef)+1 == K)
  stopifnot(ncol(coef) == P)
  
  ### Create an input dataset of patient covariates
  ## Create empty dataframe
  dat.covar <- data.frame(matrix(, nrow = N, ncol = P))
  ## Add variable names
  colnames(dat.covar) <- paste("x", 1:P, sep = "")
  ## Add covariates
  for (i in 1:P){
    dat.covar[, i] <- rnorm(N, 0, 1)
  }
  
  ### Now add the probabilities of having each outcome event
  dat.probs <- data.frame(matrix(, nrow = N, ncol = K))
  ## Add variable names
  colnames(dat.probs) <- paste("p", 1:K, sep = "")
  ## Calculate the linear predictor for each equation in multinomial regression (P(Y=K)/(P(Y=1)))
  LP <- vector("list", K-1)
  for (i in 1:length(LP)){
    LP[[i]] <- exp(beta0[i] + as.matrix(dat.covar)%*%coef[i, ])
  }
  
  ## Use these to calculate probabilities
  # First outcome category
  dat.probs[,1] <- 1/(1 + rowSums(do.call("cbind", LP)))
  # Outcome categories 2 to K
  for (i in 2:K){
    dat.probs[,i] <- (LP[[i-1]])/(1 + rowSums(do.call("cbind", LP)))
  }
  
  return(list(dat.probs, LP))
}



#####
### 4.2) Generate output predicted probabilities for a given set of input parameters
#####
### Write a function to output predicted probabilities
generate.outcome.probs.DGM.seqlog <- function(K, #number of outcome categories
                                              P, #number of predictors
                                              coef, #matrix of coefficients
                                              beta0, #vector of intercepts
                                              N){
  
  #   N <- 100
  #   K <- 4
  #   P <- 5
  #   coef <- c(0.5, -0.5, 0.5)
  #   beta0 <- c(0, 0, 0)
  
  #   N<- 500000
  #   K <- 5
  #   P <- 5
  #   beta0 <- c(0, 0, 0, 0)
  #   
  #   coef <- rbind(c(-0.5, -0.5, -0.5, -0.5, -0.5),
  #                     c(0.5, 0.5, 0.5, 0.5, 0.5),
  #                     c(0.5, 0.5, 0.5, 0.5, 0.5),
  #                     c(0.5, 0.5, 0.5, 0.5, 0.5))
  
  stopifnot(length(beta0)+1 == K)
  stopifnot(nrow(coef)+1 == K)
  stopifnot(ncol(coef) == P)
  
  ### Create an input dataset of patient covariates
  ## Create empty dataframe
  dat.covar <- data.frame(matrix(, nrow = N, ncol = P))
  ## Add variable names
  colnames(dat.covar) <- paste("x", 1:P, sep = "")
  ## Add covariates
  for (i in 1:P){
    dat.covar[, i] <- rnorm(N, 0, 1)
  }
  
  ### Now add the probabilities of having each outcome event
  dat.probs <- data.frame(matrix(, nrow = N, ncol = K))
  ## Add variable names
  colnames(dat.probs) <- paste("p", 1:K, sep = "")
  ## Calculate the linear predictor for each equation in multinomial regression (P(Y=K)/(P(Y=1)))
  LP <- vector("list", K-1)
  for (i in 1:length(LP)){
    LP[[i]] <- exp(beta0[i] + as.matrix(dat.covar)%*%coef[i, ])
  }
  
  ## Use these to calculate probabilities
  ## First, calculate the probability of an event, for each sequential model
  dat.probs.individual.models <- vector("list", K-1)
  for (i in 1:(K-1)){
    dat.probs.individual.models[[i]] <- LP[[i]]/(1+LP[[i]])
  }
  
  ## Now calculate the probability of having each event
  # First outcome category
  dat.probs[,1] <- LP[[1]]/(1+LP[[1]])
  # Outcome categories 2 to K-1
  # It's the probability of i-1, divided by the probability of i-1 from the last seq log model, multiplied by the prob of
  # NOT having i-1 from the last seq log model, multiplied by the prob of having i from the next seqlog model
  for (i in 2:(K-1)){
    dat.probs[,i] <- (dat.probs[ , (i-1)])*((1+LP[[(i-1)]])/LP[[(i-1)]])*(1 - LP[[(i-1)]]/(1 + LP[[(i-1)]]))*LP[[i]]/(1+LP[[i]])
  }
  # Outcome category K
  dat.probs[,K] <- (dat.probs[ , (K-1)])*((1+LP[[(K-1)]])/LP[[(K-1)]])*(1 - LP[[(K-1)]]/(1 + LP[[(K-1)]]))
  
  return(list(dat.probs, LP))
}


#####
### 4.3) Calculate quantiles for given input data
#####
calc.pred.quantiles.DGM.mult <- function(K, P, coef, beta0){
  p1.probs <- generate.outcome.probs.DGM.mult(K, #number of outcome categories
                                              P, #number of predictors
                                              coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                              beta0, #vector of intercepts
                                              N = 500000)
  
  return(quantile(p1.probs[[1]]$p1, probs = c(0.025, 0.975)))
}

calc.pred.quantiles.DGM.seqlog <- function(K, P, coef, beta0){
  p1.probs <- generate.outcome.probs.DGM.seqlog(K, #number of outcome categories
                                                P, #number of predictors
                                                coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
                                                beta0, #vector of intercepts
                                                N = 500000)
  
  return(quantile(p1.probs[[1]]$p1, probs = c(0.025, 0.975)))
}

