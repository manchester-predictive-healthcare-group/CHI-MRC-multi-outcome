####################################################################
### 1.1) Function to generate the data using the mulitnomial DGM ###
####################################################################

### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
generate.data.DGM.multinomial <- function(K, #number of outcome categories
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
### 1.2A) Function to generate the data using the sequential logistic DGM ###
############################################################################

### COMMENTED OUT, EITHER 1.2A) or 1.2B) MUST BE COMMENTED OUT, THEY BOTH HAVE SAME FUNCTION NAME

# ### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
# generate.data.DGM.seqlog <- function(K, #number of outcome categories
#                                    P, #number of predictors
#                                    coef, #vector of coefficients, assuming for each outcome, the effect of each predictor is the same
#                                    beta0, #vector of intercepts
#                                    N){
#   
#   #   N <- 100
#   #   K <- 4
#   #   P <- 5
#   #   coef <- rbind(c(0.75, 0.75, 0.75, 0.75, 0.75),
#   #                 c(0.5, 0.5, 0.5, 0.5, 0.5),
#   #                 c(0.25, 0.25, 0.25, 0.25, 0.25))
#   #   beta0 <- c(0, 0, 0)
#   
#   stopifnot(length(beta0)+1 == K)
#   stopifnot(nrow(coef)+1 == K)
#   stopifnot(ncol(coef) == P)
#   
#   ### Create an input dataset of patient covariates
#   ## Create empty dataframe
#   dat.covar <- data.frame(matrix(, nrow = N, ncol = P))
#   ## Add variable names
#   colnames(dat.covar) <- paste("x", 1:P, sep = "")
#   ## Add covariates
#   for (i in 1:P){
#     dat.covar[, i] <- rnorm(N, 0, 1)
#   }
#   
#   ### Now add the probabilities of having each outcome event
#   dat.probs <- data.frame(matrix(, nrow = N, ncol = K))
#   ## Add variable names
#   colnames(dat.probs) <- paste("p", 1:K, sep = "")
#   ## Calculate the linear predictor for each equation in multinomial regression (P(Y=K)/(P(Y=1)))
#   LP <- vector("list", K-1)
#   for (i in 1:length(LP)){
#     LP[[i]] <- exp(beta0[i] + as.matrix(dat.covar)%*%coef[i, ])
#   }
#   
#   ## Use these to calculate probabilities
#   ## First, calculate the probability of an event, for each sequential model
#   dat.probs.individual.models <- vector("list", K-1)
#   for (i in 1:(K-1)){
#     dat.probs.individual.models[[i]] <- LP[[i]]/(1+LP[[i]])
#   }
#   
#   ## Now calculate the probability of having each event
#   # First outcome category
#   dat.probs[,1] <- LP[[1]]/(1+LP[[1]])
#   # Outcome categories 2 to K-1
#   # It's the probability of i-1, divided by the probability of i-1 from the last seq log model, multiplied by the prob of
#   # NOT having i-1 from the last seq log model, multiplied by the prob of having i from the next seqlog model
#   for (i in 2:(K-1)){
#     dat.probs[,i] <- (dat.probs[ , (i-1)])*((1+LP[[(i-1)]])/LP[[(i-1)]])*(1 - LP[[(i-1)]]/(1 + LP[[(i-1)]]))*LP[[i]]/(1+LP[[i]])
#   }
#   # Outcome category K
#   dat.probs[,K] <- (dat.probs[ , (K-1)])*((1+LP[[(K-1)]])/LP[[(K-1)]])*(1 - LP[[(K-1)]]/(1 + LP[[(K-1)]]))
#   
#   ### Use these outcome probabilities to generate the outcomes from a multinomial distribution
#   ## Define a function to calculate multinomial outcomes for a given probability vector
#   generate.multinom <- function(probs.in){rmultinom(n = 1, size = 1, prob = probs.in)}
#   ## Apply this to each row of the dataset
#   multinom.outcomes <- t(apply(dat.probs, 1, function(x) generate.multinom(x)))
#   
#   ## Now for each row, need to convert this into a number depending on which row is equal to 1, and add it to a vector
#   y.vec <- apply(multinom.outcomes, 1, function(x) which(x > 0))
#   
#   ## Combine this with dat.covar, into a development dataset
#   dat.covar$Y <- y.vec
#   dat.covar$Y.fac <- as.factor(dat.covar$Y)
#   
#   return(dat.covar)
# }


###################################################################################################
### 1.2B) Function to generate the data using the sequential logistic DGM, reversed coefficients ###
###################################################################################################

### This generate data function used sequential logistic again, but rather than the coefficients giving the odds of having event i 
### (i.e. odds of event 1 for category 1), it gives the odds of not having event i (basically reversed). This is so the coefficients are 
### inline with the continuation ratio model from the VGAM package.

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
    dat.probs.individual.models[[i]] <- 1 - LP[[i]]/(1+LP[[i]])
  }
  
  ## Now calculate the probability of having each event
  # First outcome category
  dat.probs[,1] <- 1 - LP[[1]]/(1+LP[[1]])
  # Outcome categories 2 to K-1
  # It's the probability of i-1, divided by the probability of i-1 from the last seq log model, multiplied by the prob of
  # NOT having i-1 from the last seq log model, multiplied by the prob of having i from the next seqlog model
  for (i in 2:(K-1)){
    dat.probs[,i] <- (dat.probs[ , (i-1)])*(1/(1 - LP[[i-1]]/(1+LP[[i-1]])))*(LP[[(i-1)]]/(1 + LP[[(i-1)]]))*(1 - LP[[i]]/(1+LP[[i]]))
  }
  # Outcome category K
  dat.probs[,K] <- (dat.probs[ , (K-1)])*(1/(1 - LP[[K-1]]/(1+LP[[K-1]])))*(LP[[(K-1)]]/(1 + LP[[(K-1)]]))
  
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


###########################################################################################
### 1.3) Function to generate the data using OnevsAll models and standard normalisation ###
###########################################################################################

### Not sure if we actually need the prevalence command, this is determined through coef and beta0, which I will have to calculate in a seperate program?
generate.data.DGM.OvA <- function(K, #number of outcome categories
                                   P, #number of predictors
                                   coef, #matrix of coefficients
                                   beta0, #vector of intercepts
                                   N){
  
#     N <- 100
#     K <- 4
#     P <- 5
#     coef <- rbind(c(0.75, 0.75, 0.75, 0.75, 0.75),
#                   c(0.5, 0.5, 0.5, 0.5, 0.5),
#                   c(0.25, 0.25, 0.25, 0.25, 0.25),
#                   c(0.1, 0.1, 0.1, 0.1, 0.1))
#     beta0 <- c(0, 0, 0, 0)
  
  stopifnot(length(beta0) == K)
  stopifnot(nrow(coef) == K)
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
  ## Calculate the linear predictor for each One vs All equation (P(Y=K)/(P(Y!=K)))
  LP <- vector("list", K)
  for (i in 1:length(LP)){
    LP[[i]] <- exp(beta0[i] + as.matrix(dat.covar)%*%coef[i, ])
  }
  
  ## Use these to calculate probabilities
  # Outcome categories 1 to K
  for (i in 1:K){
    dat.probs[,i] <- LP[[i]]/(1+LP[[i]])
  }
  
  ## Normalise
  dat.probs <- dat.probs/rowSums(dat.probs)
  
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


#######################################################################################
### 1.4D) Function to generate the data using OnevsOne models and pairwise coupling ###
#######################################################################################

#####
### First need to define another function, which will do the pairwise coupling
### This function will do it for a vector of probabilities
#####
pairwise.couple <- function(probs.in, K){
  
  # Create array
  probs <- array(probs.in, dim = c(1, K*(K-1)/2))
  
  # Use code from following Github page, referenced in manuscript: https://gist.github.com/cheuerde/7c649749892c8623eee2#file-pairwise_coupling-r
  Q <- matrix(0,K,K)
  Q[lower.tri(Q)] <- 1 - probs
  Qt <- t(Q)
  Q[upper.tri(Q)] <- 1 - Qt[upper.tri(Qt)]
  diag(Q) <- rowSums(Q)
  Q <- Q / (K-1)
  
  # initial vector
  p <- rbeta(K,1,1)
  p <- p/sum(p)
  
  # updating the prob vector until equilibrium is reached
  for(i in 1:1000) p <- Q%*%p
  
  # return probs
  return(t(p))
}

#####
### Now write a function to generate the data using One vs One models and pairwise coupling
#####
generate.data.DGM.OvO.PC <- function(K, #number of outcome categories
                                  P, #number of predictors
                                  coef, #matrix of coefficients
                                  beta0, #vector of intercepts
                                  N){
  
#       N <- 100
#       K <- 4
#       P <- 5
#       coef <- rbind(c(0.75, 0.75, 0.75, 0.75, 0.75),
#                     c(0.5, 0.5, 0.5, 0.5, 0.5),
#                     c(0.25, 0.25, 0.25, 0.25, 0.25),
#                     c(0.1, 0.1, 0.1, 0.1, 0.1),
#                     c(0.7, 0.7, 0.7, 0.7, 0.7),
#                     c(0.4, 0.4, 0.4, 0.4, 0.4))
#       beta0 <- c(0, 0, 0, 0, 0, 0)
  
  stopifnot(length(beta0) == K*(K-1)/2)
  stopifnot(nrow(coef) == K*(K-1)/2)
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
  dat.probs <- data.frame(matrix(, nrow = N, ncol = (K*(K-1)/2)))
  ## Add variable names
  colnames(dat.probs) <- paste("p", 1:(K*(K-1)/2), sep = "")
  ## Calculate the linear predictor for each One vs One equation (P(Y=K1)/(P(Y=K2)))
  LP <- vector("list", K*(K-1)/2)
  for (i in 1:length(LP)){
    LP[[i]] <- exp(beta0[i] + as.matrix(dat.covar)%*%coef[i, ])
  }
  
  ## Use these to calculate probabilities for each equation
  # Outcome categories 1 to K
  for (i in 1:(K*(K-1)/2)){
    dat.probs[,i] <- LP[[i]]/(1+LP[[i]])
  }
  
  ## Use pairwise coupling to calculate risks of each outcome category
  ## First approach using couple function from kernlab package, which gives negative probabilities in some cases...
  #dat.probs <- couple(dat.probs)
  
  ## Second approach uses user-defined function
  dat.probs <- t(apply(dat.probs, 1, function(x) pairwise.couple(x, K)))
  
  ## Add variable names
  colnames(dat.probs) <- paste("p", 1:K, sep = "")
  
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


##############################################
### 2.1) Function to fit multinomial model ###
##############################################

### Fit a multinomial logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset
fit.model.multinomial <- function(dat.devel, dat.valid){
  
  ### Fit the multinomial logistic regression model
  fit.multinomial <- vgam(Y.fac ~ ., family = multinomial(ref = "1"), data = subset(dat.devel, select = -c(Y)))

  ### Generate risks for individuals in the validation cohort
  ## Extract the linear predictors, and the predicted risks
  pred.lp <- predict(fit.multinomial, newdata = dat.valid)
  pred.p <- predict(fit.multinomial, newdata = dat.valid, type = "response")
  
  ### Assign appropriate column names
  colnames(pred.lp) <- paste("lp", 1:ncol(pred.lp), sep = "")
  colnames(pred.p) <- paste("p", 1:ncol(pred.p), sep = "")
  
  ### Create output validation dataset
  valid.out <- cbind(select(dat.valid, Y, Y.fac), pred.lp, pred.p)
  
  ### Return output
  return(valid.out)
}




######################################################
### 2.2) Function to fit sequential logistic model ###
######################################################

### Fit a sequential logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset
fit.model.seqlog <- function(dat.devel, dat.valid, P){
  
  ### Create an ordinal variable version of the outcome
  dat.devel$Y.ord <- factor(dat.devel$Y.fac, order = TRUE)
  
  ### Fit the sequential logistic regression model
  fit.seqlog <- vgam(Y.ord ~ ., family = cratio(parallel = F), data = subset(dat.devel, select = -c(Y, Y.fac)))
  
  ### Generate risks for individuals in the validation cohort
  ## Extract the linear predictors, and the predicted risks
  pred.lp <- predict(fit.seqlog, newdata = dat.valid)
  pred.p <- predict(fit.seqlog, newdata = dat.valid, type = "response")
  
  ### Assign appropriate column names
  colnames(pred.lp) <- paste("lp", 1:ncol(pred.lp), sep = "")
  colnames(pred.p) <- paste("p", 1:ncol(pred.p), sep = "")
  
  ### Create output validation dataset
  valid.out <- cbind(dplyr::select(dat.valid, Y, Y.fac), pred.lp, pred.p)
  
  ### Return output
  return(valid.out)
}


##########################################################################
### 2.3) Function to fit One vs All models with standard normalisation ###
##########################################################################

### Fit one vs all models to the dataset and normalise the resulting risks
fit.model.OvA <- function(dat.devel, dat.valid){

  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.devel$Y))
  
  ### Create variables which can be used to do the One vs All models and add to datasets
  Y.dummies.devel <- model.matrix( ~ Y.fac - 1, dat.devel)
  Y.dummies.valid <- model.matrix( ~ Y.fac - 1, dat.valid)
  
  dat.devel <- cbind(dat.devel, Y.dummies.devel)
  dat.valid <- cbind(dat.valid, Y.dummies.valid)

  ### Fit a seperate one vs all model for each outcome category and store in a list
  fit.OvA.list <- vector("list", K)

  if (K == 3){
    fit.OvA.list[[1]] <- glm(Y.fac1 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac2, Y.fac3)))
    fit.OvA.list[[2]] <- glm(Y.fac2 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac1, Y.fac3)))
    fit.OvA.list[[3]] <- glm(Y.fac3 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac1, Y.fac2)))
  } else if (K == 4){
    fit.OvA.list[[1]] <- glm(Y.fac1 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac2, Y.fac3, Y.fac4)))
    fit.OvA.list[[2]] <- glm(Y.fac2 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac1, Y.fac3, Y.fac4)))
    fit.OvA.list[[3]] <- glm(Y.fac3 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac1, Y.fac2, Y.fac4)))
    fit.OvA.list[[4]] <- glm(Y.fac4 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac1, Y.fac2, Y.fac3)))
  } else if (K == 5){
    fit.OvA.list[[1]] <- glm(Y.fac1 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac2, Y.fac3, Y.fac4, Y.fac5)))
    fit.OvA.list[[2]] <- glm(Y.fac2 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac1, Y.fac3, Y.fac4, Y.fac5)))
    fit.OvA.list[[3]] <- glm(Y.fac3 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac1, Y.fac2, Y.fac4, Y.fac5)))
    fit.OvA.list[[4]] <- glm(Y.fac4 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac1, Y.fac2, Y.fac3, Y.fac5)))
    fit.OvA.list[[5]] <- glm(Y.fac5 ~ ., family = binomial(link = "logit"), data = dplyr::select(dat.devel, !c(Y, Y.fac, Y.fac1, Y.fac2, Y.fac3, Y.fac4)))
  }
  
  ### Use these models to create non-normalised predicted risks for each outcome category
  ## Create dataset to temporarily store these in
  risks.unnorm <- data.frame(matrix(NA, nrow = nrow(dat.valid), ncol = K))
  colnames(risks.unnorm) <- paste("p", 1:ncol(risks.unnorm), sep = "")
  
  ## Also create a dataset to store lp's in, which may be useful for checking weak calibration at a later point
  lp.unnorm <- data.frame(matrix(NA, nrow = nrow(dat.valid), ncol = K))
  colnames(lp.unnorm) <- paste("lp", 1:ncol(risks.unnorm), sep = "")
  
  ## Calculate the predicted risk for each outcome category using each model
  for (model in 1:K){
    lp.unnorm[, model] <- predict(fit.OvA.list[[model]], newdata = dat.valid)
    risks.unnorm[, model] <- predict(fit.OvA.list[[model]], newdata = dat.valid, type = "response")
  }
  
  ### Normalise the risks
  risks.norm <- risks.unnorm/rowSums(risks.unnorm)
  
  ### Create output dataset
  valid.out <- cbind(dplyr::select(dat.valid, Y, Y.fac), lp.unnorm, risks.norm)
 
  ### Return output
  return(valid.out)
}


#########################################################################
### 2.4) Function to fit One vs All models with softmax normalisation ###
#########################################################################

### Fit one vs all models to the dataset and normalise the resulting risks
fit.model.OvA.softmax <- function(dat.devel, dat.valid){
  
}



#####################################################################
### 2.5) Function to fit One vs One models with pairwise coupling ###
#####################################################################

#####
### Now creaate function to fit one vs all models to the dataset, produce predicted risks for validation dataset using pairwise coupling
#####
fit.model.OvO.PC <- function(dat.devel, dat.valid){

  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.devel$Y))
  
  ### Create variables for each pair of outcomes (I am actually not going to consider 5 outcome categories for this one)
  if (K == 3){
    dat.devel <- dat.devel %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                         Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                         Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_))
    
    dat.valid <- dat.valid %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                         Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                         Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_))
    
  } else if (K == 4){
    dat.devel <- dat.devel %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                         Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                         Y1.4 = case_when(Y == 1 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                         Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                         Y2.4 = case_when(Y == 2 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                         Y3.4 = case_when(Y == 3 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_))
    
    dat.valid <- dat.valid %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                         Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                         Y1.4 = case_when(Y == 1 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                         Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                         Y2.4 = case_when(Y == 2 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                         Y3.4 = case_when(Y == 3 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_))
  } else if (K == 5){
    dat.devel <- dat.devel %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                                      Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                      Y1.4 = case_when(Y == 1 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                      Y1.5 = case_when(Y == 1 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_),
                                      Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                      Y2.4 = case_when(Y == 2 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                      Y2.5 = case_when(Y == 2 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_),
                                      Y3.4 = case_when(Y == 3 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                      Y3.5 = case_when(Y == 3 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_),
                                      Y4.5 = case_when(Y == 4 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_))
    
    dat.valid <- dat.valid %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                                      Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                      Y1.4 = case_when(Y == 1 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                      Y1.5 = case_when(Y == 1 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_),
                                      Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                      Y2.4 = case_when(Y == 2 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                      Y2.5 = case_when(Y == 2 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_),
                                      Y3.4 = case_when(Y == 3 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                      Y3.5 = case_when(Y == 3 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_),
                                      Y4.5 = case_when(Y == 4 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_))
  }

  ### Fit the One vs One models and store in a list
  fit.OvO.list <- vector("list", choose(K, 2))

  if (K == 3){
    fit.OvO.list[[1]] <- glm(Y1.2 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(1,2)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.3, Y2.3)))
    fit.OvO.list[[2]] <- glm(Y1.3 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(1,3)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y2.3)))
    fit.OvO.list[[3]] <- glm(Y2.3 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(2,3)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3)))
  } else if (K == 4){
    fit.OvO.list[[1]] <- glm(Y1.2 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(1,2)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.3, Y1.4, Y2.3, Y2.4, Y3.4)))
    fit.OvO.list[[2]] <- glm(Y1.3 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(1,3)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.4, Y2.3, Y2.4, Y3.4)))
    fit.OvO.list[[3]] <- glm(Y1.4 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(1,4)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y2.3, Y2.4, Y3.4)))
    fit.OvO.list[[4]] <- glm(Y2.3 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(2,3)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y2.4, Y3.4)))
    fit.OvO.list[[5]] <- glm(Y2.4 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(2,4)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y2.3, Y3.4)))
    fit.OvO.list[[6]] <- glm(Y3.4 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(3,4)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y2.3, Y2.4)))
  } else if (K == 5){
    fit.OvO.list[[1]] <- glm(Y1.2 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(1,2)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.3, Y1.4, Y1.5, Y2.3, Y2.4, Y2.5, Y3.4, Y3.5, Y4.5)))
    fit.OvO.list[[2]] <- glm(Y1.3 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(1,3)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.4, Y1.5, Y2.3, Y2.4, Y2.5, Y3.4, Y3.5, Y4.5)))
    fit.OvO.list[[3]] <- glm(Y1.4 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(1,4)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.5, Y2.3, Y2.4, Y2.5, Y3.4, Y3.5, Y4.5)))
    fit.OvO.list[[4]] <- glm(Y1.5 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(1,5)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y2.3, Y2.4, Y2.5, Y3.4, Y3.5, Y4.5)))
    fit.OvO.list[[5]] <- glm(Y2.3 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(2,3)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y1.5, Y2.4, Y2.5, Y3.4, Y3.5, Y4.5)))
    fit.OvO.list[[6]] <- glm(Y2.4 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(2,4)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y1.5, Y2.3, Y2.5, Y3.4, Y3.5, Y4.5)))
    fit.OvO.list[[7]] <- glm(Y2.5 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(2,5)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y1.5, Y2.3, Y2.4, Y3.4, Y3.5, Y4.5)))
    fit.OvO.list[[8]] <- glm(Y3.4 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(3,4)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y1.5, Y2.3, Y2.4, Y2.5, Y3.5, Y4.5)))
    fit.OvO.list[[9]] <- glm(Y3.5 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(3,5)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y1.5, Y2.3, Y2.4, Y2.5, Y3.4, Y4.5)))
    fit.OvO.list[[10]] <- glm(Y4.5 ~ ., family = binomial(link = "logit"), data = subset(dat.devel, Y %in% c(4,5)) %>% 
                               dplyr::select(-c(Y, Y.fac, Y1.2, Y1.3, Y1.4, Y1.5, Y2.3, Y2.4, Y2.5, Y3.4, Y3.5)))
  } 

  ### Use these models to create pairwise risks of each category, also store linear predictors to enable weak calibration
  ## lp
  pairwise.lp <- lapply(fit.OvO.list, function(x){
    predict(x, newdata = dat.valid)
  })
  pairwise.lp <- do.call(cbind, pairwise.lp, 1)
  colnames(pairwise.lp) <- paste("lp", 1:ncol(pairwise.lp), sep = "")
  
  ## p
  pairwise.p <- lapply(fit.OvO.list, function(x){
    predict(x, newdata = dat.valid, type = "response")
  })
  pairwise.p <- do.call(cbind, pairwise.p, 1)
  
  ### Turn the pairwise risks into risks using pairwise coupling
  
  ## First approach using couple function from kernlab package, which gives negative probabilities in some cases...
  #pred.p <- couple(pairwise.p)
  
  ## Second approach uses user-defined function
  pred.p <- t(apply(pairwise.p, 1, function(x) pairwise.couple(x, K)))

  ## Assign colnames
  colnames(pred.p) <- paste("p", 1:ncol(pred.p), sep = "")
  
  ### Create output dataset
  valid.out <- cbind(dplyr::select(dat.valid, Y, Y.fac), pairwise.lp, pred.p)
  
  ### Return output
  return(valid.out)
}


###########################################################################################################
### 2.6) Function to fit sequential logistic model, using every possible ordering of outcome categories ###
###########################################################################################################

### Fit a sequential logistic regression model to the development dataset, and produce produce predicted risks for a validation dataset
fit.model.seqlog.all <- function(dat.devel, dat.valid, P){
  
#   dat.devel <- dat.devel.list[[1]]
#   dat.valid <- dat.devel.list[[1]]
#   P <- 5
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.devel$Y))
  
  ### Create list of permutations. Note we remove permuatations where the final two categories, 
  ### as these will result in the same predicted risks, and will reduce half the models we need to fit
  
  ## Start by generating all perms
  perm.list <- permn(c(1:K))
  
  ## Now identify which are duplicates with the final two elements swapped
  dups <- duplicated(lapply(perm.list, function(vec) vec[1:(K-2)]))
  
  ## Now subset the perm.list to only contain the non-duplicates
  perm.list <- perm.list[dups]
  
  ### Number of permuatations
  n.order <- length(perm.list)
  
  ### Want to repeat the following process, for every possible order of outcome categories, and store
  ### predicted risks, to then take an average over
  
  ### Create output object to store predicted risks
  predrisk.list <- vector("list", n.order)
  
  ### Cycle through and calculate predicted risk for every permuatation
  for (order in 1:n.order){
    
    ### Assign permuatation
    perm <- perm.list[[order]]
    print(paste(order, Sys.time()))
    print(perm)
    
    ### Create an ordinal variable version of the outcome
    dat.devel$Y.ord <- factor(dat.devel$Y.fac, order = TRUE, 
                              levels = c(1:K)[perm])
    
    ### Fit the sequential logistic regression model
    fit.seqlog <- vgam(Y.ord ~ ., family = cratio(parallel = F), data = subset(dat.devel, select = -c(Y, Y.fac)))
    
    ### Generate risks for individuals in the validation cohort
    ## Extract the linear predictors, and the predicted risks
    pred.p <- predict(fit.seqlog, newdata = dat.valid, type = "response")
    
    ### Put the predicted risks into a common order
    if (K == 3){
      pred.p <- pred.p[,c(which(perm == 1), which(perm == 2), which(perm == 3))]
    } else if (K == 5){
      pred.p <- pred.p[,c(which(perm == 1), which(perm == 2), which(perm == 3), which(perm == 4), which(perm == 5))]
    }
    
    #head(pred.p)
    
    ### Assign appropriate column names
    colnames(pred.p) <- paste("p", 1:ncol(pred.p), sep = "")
    
    ### Create output validation dataset
    predrisk.list[[order]] <- pred.p
    
  }
  
  ### Take average over all permuatations
  pred.p.mean <- apply(simplify2array(predrisk.list), 1:2, mean)
  
  ### Create output validation dataset
  valid.out <- cbind(dplyr::select(dat.valid, Y, Y.fac), pred.p.mean)
  
  ### Return output
  return(valid.out)
}

#################################################################################################################
### 3.1A) Function to assess overfitting via calibration slope and intercept using each model-specific methods ###
### Note these use the linear predictors calculated in the validation cohort, and can therefore only be applied
### to output from the appropriate model.

### Method = multinomial logistic regression
#################################################################################################################

calibrate.mod.specific.multinomial <- function(dat.valid.pred){
  
  ### dat.valid.pred must be output from fit.model.multinomial
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  
  ### Now fit model specific recalibration model
  if (K == 3){
    
    ## Add constraints
    i <- diag(2)
    i2 <- rbind(1, 0)
    i3 <- rbind(0, 1)
    clist <- list("(Intercept)" = i, "lp1" = i2, "lp2" = i3)
    clist
    
    calib.model <- vgam(Y.fac ~ lp1 + lp2, constraints = clist, 
                        data = dat.valid.pred, family = multinomial(refLevel = "1"))
    
    calib.model.offset <- vgam(dat.valid.pred$Y.fac ~ 1, offset = as.matrix(dat.valid.pred[, c("lp1", "lp2")]), 
                               family = multinomial(refLevel = "1"))
    
  } else if (K == 4){
    
    ## Add constraints
    i <- diag(3)
    i2 <- rbind(1, 0, 0)
    i3 <- rbind(0, 1, 0)
    i4 <- rbind(0, 0, 1)
    clist <- list("(Intercept)" = i, "lp1" = i2, "lp2" = i3, "lp3" = i4)
    clist
    
    calib.model <- vgam(Y.fac ~ lp1 + lp2 + lp3, constraints = clist, 
                        data = dat.valid.pred, family = multinomial(refLevel = "1"))
    
    calib.model.offset <- vgam(dat.valid.pred$Y.fac ~ 1, offset = as.matrix(dat.valid.pred[, c("lp1", "lp2", "lp3")]), 
                               family = multinomial(refLevel = "1"))
  } else if (K == 5){
    
    ## Add constraints
    i <- diag(4)
    i2 <- rbind(1, 0, 0, 0)
    i3 <- rbind(0, 1, 0, 0)
    i4 <- rbind(0, 0, 1, 0)
    i5 <- rbind(0, 0, 0, 1)
    clist <- list("(Intercept)" = i, "lp1" = i2, "lp2" = i3, "lp3" = i4, "lp4" = i5)
    clist
    
    calib.model <- vgam(Y.fac ~ lp1 + lp2 + lp3 + lp4, constraints = clist, 
                        data = dat.valid.pred, family = multinomial(refLevel = "1"))
    
    calib.model.offset <- vgam(dat.valid.pred$Y.fac ~ 1, offset = as.matrix(dat.valid.pred[, c("lp1", "lp2", "lp3", "lp4")]), 
                               family = multinomial(refLevel = "1"))
  }
  
  ### Extract relevant coefficients
  alpha <- calib.model.offset@coefficients
  names(alpha) <- paste("lp", 1:(K-1), sep = "")
  
  beta <- calib.model@coefficients[grep("lp", names(calib.model@coefficients))]
  names(beta) <- paste("lp", 1:(K-1), sep = "")
  
  ### Create output object and return
  output.obj <- list("alpha" = alpha, "beta" = beta)
  return(output.obj)
}



#################################################################################################################
### 3.1B) Function to assess overfitting via calibration slope and intercept using each model-specific methods ###
### Note these use the linear predictors calculated in the validation cohort, and can therefore only be applied
### to output from the appropriate model.

### Method = sequential logistic regression
#################################################################################################################

calibrate.mod.specific.seqlog <- function(dat.valid.pred){
  
  ### dat.valid.pred must be output from fit.model.seqlog
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  
  ### Create an ordinal variable version of the outcome
  dat.valid.pred$Y.ord <- factor(dat.valid.pred$Y.fac, order = TRUE)
  
  ### Now fit model specific recalibration model
  if (K == 3){
    
    ## Add constraints
    i <- diag(2)
    i2 <- rbind(1, 0)
    i3 <- rbind(0, 1)
    clist <- list("(Intercept)" = i, "lp1" = i2, "lp2" = i3)
    clist
    
    calib.model <- vgam(Y.ord ~ lp1 + lp2, constraints = clist, 
                        data = dat.valid.pred, family = cratio(parallel = F))
    
    calib.model.offset <- vgam(dat.valid.pred$Y.ord ~ 1, offset = as.matrix(dat.valid.pred[, c("lp1", "lp2")]), 
                               family = cratio(parallel = F))
    
  } else if (K == 4){
    
    ## Add constraints
    i <- diag(3)
    i2 <- rbind(1, 0, 0)
    i3 <- rbind(0, 1, 0)
    i4 <- rbind(0, 0, 1)
    clist <- list("(Intercept)" = i, "lp1" = i2, "lp2" = i3, "lp3" = i4)
    clist
    
    calib.model <- vgam(Y.ord ~ lp1 + lp2 + lp3, constraints = clist, 
                        data = dat.valid.pred, family = cratio(parallel = F))
    
    calib.model.offset <- vgam(dat.valid.pred$Y.ord ~ 1, offset = as.matrix(dat.valid.pred[, c("lp1", "lp2", "lp3")]), 
                               family = cratio(parallel = F))
  } else if (K == 5){
    
    ## Add constraints
    i <- diag(4)
    i2 <- rbind(1, 0, 0, 0)
    i3 <- rbind(0, 1, 0, 0)
    i4 <- rbind(0, 0, 1, 0)
    i5 <- rbind(0, 0, 0, 1)
    clist <- list("(Intercept)" = i, "lp1" = i2, "lp2" = i3, "lp3" = i4, "lp4" = i5)
    clist
    
    calib.model <- vgam(Y.ord ~ lp1 + lp2 + lp3 + lp4, constraints = clist, 
                        data = dat.valid.pred, family = cratio(parallel = F))
    
    calib.model.offset <- vgam(dat.valid.pred$Y.ord ~ 1, offset = as.matrix(dat.valid.pred[, c("lp1", "lp2", "lp3", "lp4")]), 
                               family = cratio(parallel = F))
  }
  
  ### Extract relevant coefficients
  # Note that for alpha, we take the negative of the intercept. This is because the c-ratio model predicts probability that Y>k for kth model,
  # This that means for k = 1, we are predicting probability that k!=1, and calibrate in tis direction too.
  # This is the opposite of calibration per outcome, where we will predict probability that Y = k (Prob Y = 1 for k = 1).
  # For k=1, the calibration should be the same for these two calibrations, so we take the negative of the intercept in the c-ratio model.
  # An alternative would be to change settings in cratio model to model probability !Y>k, (i.e. set reverse = TRUE)
  
  ## I'm actually not going to do this (use the -ve), because I have already run some programs without it like this. Instead I will add the negative
  ## sign when generating the table a tthe final step.
  #alpha <- -calib.model.offset@coefficients
  alpha <- calib.model.offset@coefficients
  names(alpha) <- paste("lp", 1:(K-1), sep = "")
  
  beta <- calib.model@coefficients[grep("lp", names(calib.model@coefficients))]
  names(beta) <- paste("lp", 1:(K-1), sep = "")
  
  ### Create output object and return
  output.obj <- list("alpha" = alpha, "beta" = beta)
  return(output.obj)
}


#################################################################################################################
### 3.1C) Function to assess overfitting via calibration slope and intercept using each model-specific methods ###
### Note these use the linear predictors calculated in the validation cohort, and can therefore only be applied
### to output from the appropriate model.

### Method = One vs All with normalisation
#################################################################################################################

calibrate.mod.specific.OvA <- function(dat.valid.pred){
  
  ### dat.valid.pred must be output from fit.model.OvA
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  
  ### Create variables which can be used to do the One vs All models and add to datasets
  Y.dummies <- model.matrix( ~ Y.fac - 1, dat.valid.pred)
  
  ### Combine with validation dataset
  dat.valid.pred <- cbind(dat.valid.pred, Y.dummies)
  
  ### Fit calibration model, extract intercept/slope per outcome
  if (K == 3){
    
    ### Fit recalibraiton model for slopes
    calib.model1 <- glm(Y.fac1 ~ lp1, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model2 <- glm(Y.fac2 ~ lp2, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model3 <- glm(Y.fac3 ~ lp3, family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Fit recalibration model for intercepts, where slope is fixed to 1
    calib.model1.offset <- glm(Y.fac1 ~ offset(lp1), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model2.offset <- glm(Y.fac2 ~ offset(lp2), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model3.offset <- glm(Y.fac3 ~ offset(lp3), family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Store output in a vector
    alpha.model <- c(calib.model1.offset$coefficients[1], 
                     calib.model2.offset$coefficients[1], 
                     calib.model3.offset$coefficients[1])
    names(alpha.model) <- paste("lp", 1:K, sep = "")
    beta.model <- c(calib.model1$coefficients[2], 
                    calib.model2$coefficients[2], 
                    calib.model3$coefficients[2])
    names(beta.model) <- paste("lp", 1:K, sep = "")
    
  } else if (K == 4){
    
    ### Fit recalibraiton model for slopes
    calib.model1 <- glm(Y.fac1 ~ lp1, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model2 <- glm(Y.fac2 ~ lp2, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model3 <- glm(Y.fac3 ~ lp3, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model4 <- glm(Y.fac4 ~ lp4, family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Fit recalibration model for intercepts, where slope is fixed to 1
    calib.model1.offset <- glm(Y.fac1 ~ offset(lp1), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model2.offset <- glm(Y.fac2 ~ offset(lp2), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model3.offset <- glm(Y.fac3 ~ offset(lp3), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model4.offset <- glm(Y.fac4 ~ offset(lp4), family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Store output in a vector
    ## Intercept
    alpha.model <- c(calib.model1.offset$coefficients[1], 
                     calib.model2.offset$coefficients[1], 
                     calib.model3.offset$coefficients[1], 
                     calib.model4.offset$coefficients[1])
    names(alpha.model) <- paste("lp", 1:K, sep = "")
    
    ## Slope
    beta.model <- c(calib.model1$coefficients[2], 
                    calib.model2$coefficients[2], 
                    calib.model3$coefficients[2], 
                    calib.model4$coefficients[2])
    names(beta.model) <- paste("lp", 1:K, sep = "")
  } else if (K == 5){
    
    ### Fit recalibraiton model for slopes
    calib.model1 <- glm(Y.fac1 ~ lp1, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model2 <- glm(Y.fac2 ~ lp2, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model3 <- glm(Y.fac3 ~ lp3, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model4 <- glm(Y.fac4 ~ lp4, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model5 <- glm(Y.fac5 ~ lp5, family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Fit recalibration model for intercepts, where slope is fixed to 1
    calib.model1.offset <- glm(Y.fac1 ~ offset(lp1), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model2.offset <- glm(Y.fac2 ~ offset(lp2), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model3.offset <- glm(Y.fac3 ~ offset(lp3), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model4.offset <- glm(Y.fac4 ~ offset(lp4), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.model5.offset <- glm(Y.fac5 ~ offset(lp5), family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Store output in a vector
    ## Intercept
    alpha.model <- c(calib.model1.offset$coefficients[1], 
                     calib.model2.offset$coefficients[1], 
                     calib.model3.offset$coefficients[1], 
                     calib.model4.offset$coefficients[1], 
                     calib.model5.offset$coefficients[1])
    names(alpha.model) <- paste("lp", 1:K, sep = "")
    
    ## Slope
    beta.model <- c(calib.model1$coefficients[2], 
                    calib.model2$coefficients[2], 
                    calib.model3$coefficients[2], 
                    calib.model4$coefficients[2], 
                    calib.model5$coefficients[2])
    names(beta.model) <- paste("lp", 1:K, sep = "")
  }
  
  ### Create output object and return
  output.obj <- list("alpha" = alpha.model, "beta" = beta.model)
  return(output.obj)
}


#################################################################################################################
### 3.1D) Function to assess overfitting via calibration slope and intercept using each model-specific methods ###
### Note these use the linear predictors calculated in the validation cohort, and can therefore only be applied
### to output from the appropriate model.

### Method = One vs One with pairwise coupling
#################################################################################################################

calibrate.mod.specific.OvO.PC <- function(dat.valid.pred){
  
  ### dat.valid.pred must be output from fit.model.OvO.PC
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  
  ### Create variables for each pair of outcomes (I am actually not going to consider 5 outcome categories for this one)
  if (K == 3){
    dat.valid.pred <- dat.valid.pred %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                                                Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                                Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_))
    
  } else if (K == 4){
    dat.valid.pred <- dat.valid.pred %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                                                Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                                Y1.4 = case_when(Y == 1 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                                Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                                Y2.4 = case_when(Y == 2 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                                Y3.4 = case_when(Y == 3 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_))
  } else if (K == 5){
    dat.valid.pred <- dat.valid.pred %>% mutate(Y1.2 = case_when(Y == 1 ~ 1, Y == 2 ~ 0, TRUE ~ NA_real_),
                                                Y1.3 = case_when(Y == 1 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                                Y1.4 = case_when(Y == 1 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                                Y1.5 = case_when(Y == 1 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_),
                                                Y2.3 = case_when(Y == 2 ~ 1, Y == 3 ~ 0, TRUE ~ NA_real_),
                                                Y2.4 = case_when(Y == 2 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                                Y2.5 = case_when(Y == 2 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_),
                                                Y3.4 = case_when(Y == 3 ~ 1, Y == 4 ~ 0, TRUE ~ NA_real_),
                                                Y3.5 = case_when(Y == 3 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_),
                                                Y4.5 = case_when(Y == 4 ~ 1, Y == 5 ~ 0, TRUE ~ NA_real_))
  }
  
  ### Fit the recalibration models and in a list
  ## Two lists, one for recalibraiton models, and one for recalibraiton model with slope as offset (to calculate intercept)
  calib.OvO.list <- vector("list", choose(K, 2))
  calib.OvO.list.offset <- vector("list", choose(K, 2))
  
  if (K == 3){
    ### Fit models
    calib.OvO.list[[1]] <- glm(Y1.2 ~ lp1, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,2)))
    calib.OvO.list[[2]] <- glm(Y1.3 ~ lp2, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,3)))
    calib.OvO.list[[3]] <- glm(Y2.3 ~ lp3, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,3)))
    
    calib.OvO.list.offset[[1]] <- glm(Y1.2 ~ offset(lp1), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,2)))
    calib.OvO.list.offset[[2]] <- glm(Y1.3 ~ offset(lp2), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,3)))
    calib.OvO.list.offset[[3]] <- glm(Y2.3 ~ offset(lp3), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,3)))
    
    ### Store output in a vector
    alpha.per.outcome <- c(calib.OvO.list.offset[[1]]$coefficients[1], 
                           calib.OvO.list.offset[[2]]$coefficients[1], 
                           calib.OvO.list.offset[[3]]$coefficients[1])
    names(alpha.per.outcome) <- paste("lp", 1:(K*(K-1)/2), sep = "")
    
    beta.per.outcome <- c(calib.OvO.list[[1]]$coefficients[2], 
                          calib.OvO.list[[2]]$coefficients[2], 
                          calib.OvO.list[[3]]$coefficients[2])
    names(beta.per.outcome) <- paste("lp", 1:(K*(K-1)/2), sep = "")
    
  } else if (K == 4){
    ### Fit models
    calib.OvO.list[[1]] <- glm(Y1.2 ~ lp1, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,2)))
    calib.OvO.list[[2]] <- glm(Y1.3 ~ lp2, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,3)))
    calib.OvO.list[[3]] <- glm(Y1.4 ~ lp3, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,4)))
    calib.OvO.list[[4]] <- glm(Y2.3 ~ lp4, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,3)))
    calib.OvO.list[[5]] <- glm(Y2.4 ~ lp5, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,4)))
    calib.OvO.list[[6]] <- glm(Y3.4 ~ lp6, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(3,4)))
    
    calib.OvO.list.offset[[1]] <- glm(Y1.2 ~ offset(lp1), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,2)))
    calib.OvO.list.offset[[2]] <- glm(Y1.3 ~ offset(lp2), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,3)))
    calib.OvO.list.offset[[3]] <- glm(Y1.4 ~ offset(lp3), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,4)))
    calib.OvO.list.offset[[4]] <- glm(Y2.3 ~ offset(lp4), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,3)))
    calib.OvO.list.offset[[5]] <- glm(Y2.4 ~ offset(lp5), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,4)))
    calib.OvO.list.offset[[6]] <- glm(Y3.4 ~ offset(lp6), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(3,4)))
    
    ### Store output in a vector
    alpha.per.outcome <- c(calib.OvO.list.offset[[1]]$coefficients[1], 
                           calib.OvO.list.offset[[2]]$coefficients[1], 
                           calib.OvO.list.offset[[3]]$coefficients[1],
                           calib.OvO.list.offset[[4]]$coefficients[1], 
                           calib.OvO.list.offset[[5]]$coefficients[1], 
                           calib.OvO.list.offset[[6]]$coefficients[1])
    names(alpha.per.outcome) <- paste("lp", 1:(K*(K-1)/2), sep = "")
    
    beta.per.outcome <- c(calib.OvO.list[[1]]$coefficients[2], 
                          calib.OvO.list[[2]]$coefficients[2], 
                          calib.OvO.list[[3]]$coefficients[2],
                          calib.OvO.list[[4]]$coefficients[2], 
                          calib.OvO.list[[5]]$coefficients[2], 
                          calib.OvO.list[[6]]$coefficients[2])
    names(beta.per.outcome) <- paste("lp", 1:(K*(K-1)/2), sep = "")
    
  } else if (K == 5){
    ### Fit models
    calib.OvO.list[[1]] <- glm(Y1.2 ~ lp1, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,2)))
    calib.OvO.list[[2]] <- glm(Y1.3 ~ lp2, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,3)))
    calib.OvO.list[[3]] <- glm(Y1.4 ~ lp3, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,4)))
    calib.OvO.list[[4]] <- glm(Y1.5 ~ lp4, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,5)))
    calib.OvO.list[[5]] <- glm(Y2.3 ~ lp5, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,3)))
    calib.OvO.list[[6]] <- glm(Y2.4 ~ lp6, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,4)))
    calib.OvO.list[[7]] <- glm(Y2.5 ~ lp7, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,5)))
    calib.OvO.list[[8]] <- glm(Y3.4 ~ lp8, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(3,4)))
    calib.OvO.list[[9]] <- glm(Y3.5 ~ lp9, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(3,5)))
    calib.OvO.list[[10]] <- glm(Y4.5 ~ lp10, family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(4,5)))
    
    calib.OvO.list.offset[[1]] <- glm(Y1.2 ~ offset(lp1), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,2)))
    calib.OvO.list.offset[[2]] <- glm(Y1.3 ~ offset(lp2), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,3)))
    calib.OvO.list.offset[[3]] <- glm(Y1.4 ~ offset(lp3), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,4)))
    calib.OvO.list.offset[[4]] <- glm(Y1.5 ~ offset(lp4), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(1,5)))
    calib.OvO.list.offset[[5]] <- glm(Y2.3 ~ offset(lp5), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,3)))
    calib.OvO.list.offset[[6]] <- glm(Y2.4 ~ offset(lp6), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,4)))
    calib.OvO.list.offset[[7]] <- glm(Y2.5 ~ offset(lp7), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(2,5)))
    calib.OvO.list.offset[[8]] <- glm(Y3.4 ~ offset(lp8), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(3,4)))
    calib.OvO.list.offset[[9]] <- glm(Y3.5 ~ offset(lp9), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(3,5)))
    calib.OvO.list.offset[[10]] <- glm(Y4.5 ~ offset(lp10), family = binomial(link = "logit"), data = subset(dat.valid.pred, Y %in% c(4,5)))
    
    ### Store output in a vector
    alpha.per.outcome <- c(calib.OvO.list.offset[[1]]$coefficients[1], 
                           calib.OvO.list.offset[[2]]$coefficients[1], 
                           calib.OvO.list.offset[[3]]$coefficients[1],
                           calib.OvO.list.offset[[4]]$coefficients[1], 
                           calib.OvO.list.offset[[5]]$coefficients[1], 
                           calib.OvO.list.offset[[6]]$coefficients[1], 
                           calib.OvO.list.offset[[7]]$coefficients[1],
                           calib.OvO.list.offset[[8]]$coefficients[1], 
                           calib.OvO.list.offset[[9]]$coefficients[1], 
                           calib.OvO.list.offset[[10]]$coefficients[1])
    names(alpha.per.outcome) <- paste("lp", 1:(K*(K-1)/2), sep = "")
    
    beta.per.outcome <- c(calib.OvO.list[[1]]$coefficients[2], 
                          calib.OvO.list[[2]]$coefficients[2], 
                          calib.OvO.list[[3]]$coefficients[2],
                          calib.OvO.list[[4]]$coefficients[2], 
                          calib.OvO.list[[5]]$coefficients[2], 
                          calib.OvO.list[[6]]$coefficients[2], 
                          calib.OvO.list[[7]]$coefficients[2],
                          calib.OvO.list[[8]]$coefficients[2], 
                          calib.OvO.list[[9]]$coefficients[2], 
                          calib.OvO.list[[10]]$coefficients[2])
    names(beta.per.outcome) <- paste("lp", 1:(K*(K-1)/2), sep = "")
    
  } 
  
  ### Create output object and return
  output.obj <- list("alpha" = alpha.per.outcome, "beta" = beta.per.outcome)
  return(output.obj)
}


######################################################################################
### 3.2) Function to produce calibration slope and intercept, per outcome category ###
######################################################################################

calibrate.weak.per.outcome <- function(dat.valid.pred){
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  
  ### Create variables which can be used to do the One vs All models and add to datasets
  Y.dummies <- model.matrix( ~ Y.fac - 1, dat.valid.pred)
  
  ### Combine with validation dataset
  dat.valid.pred <- cbind(dat.valid.pred, Y.dummies)
  
  ### Create logit(p_k) for each outcome category, fit calibration model, extract intercept/slope per outcome
  if (K == 3){
    ### Create logit(p_k)
    dat.valid.pred$logit_p1 <- log(dat.valid.pred$p1/(1 - dat.valid.pred$p1))
    dat.valid.pred$logit_p2 <- log(dat.valid.pred$p2/(1 - dat.valid.pred$p2))
    dat.valid.pred$logit_p3 <- log(dat.valid.pred$p3/(1 - dat.valid.pred$p3))
    
    ### Fit recalibraiton model for slopes
    calib.per.outcome1 <- glm(Y.fac1 ~ logit_p1, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome2 <- glm(Y.fac2 ~ logit_p2, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome3 <- glm(Y.fac3 ~ logit_p3, family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Fit recalibration model for intercepts, where slope is fixed to 1
    calib.per.outcome1.offset <- glm(Y.fac1 ~ offset(logit_p1), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome2.offset <- glm(Y.fac2 ~ offset(logit_p2), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome3.offset <- glm(Y.fac3 ~ offset(logit_p3), family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Store output in a vector
    alpha.per.outcome <- c(calib.per.outcome1.offset$coefficients[1], 
                           calib.per.outcome2.offset$coefficients[1], 
                           calib.per.outcome3.offset$coefficients[1])
    names(alpha.per.outcome) <- paste("p", 1:K, sep = "")
    
    beta.per.outcome <- c(calib.per.outcome1$coefficients[2], 
                          calib.per.outcome2$coefficients[2], 
                          calib.per.outcome3$coefficients[2])
    names(beta.per.outcome) <- paste("p", 1:K, sep = "")
    
  } else if (K == 4){
    ### Create logit(p_k)
    dat.valid.pred$logit_p1 <- log(dat.valid.pred$p1/(1 - dat.valid.pred$p1))
    dat.valid.pred$logit_p2 <- log(dat.valid.pred$p2/(1 - dat.valid.pred$p2))
    dat.valid.pred$logit_p3 <- log(dat.valid.pred$p3/(1 - dat.valid.pred$p3))
    dat.valid.pred$logit_p4 <- log(dat.valid.pred$p4/(1 - dat.valid.pred$p4))
    
    ### Fit recalibraiton model for slopes
    calib.per.outcome1 <- glm(Y.fac1 ~ logit_p1, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome2 <- glm(Y.fac2 ~ logit_p2, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome3 <- glm(Y.fac3 ~ logit_p3, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome4 <- glm(Y.fac4 ~ logit_p4, family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Fit recalibration model for intercepts, where slope is fixed to 1
    calib.per.outcome1.offset <- glm(Y.fac1 ~ offset(logit_p1), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome2.offset <- glm(Y.fac2 ~ offset(logit_p2), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome3.offset <- glm(Y.fac3 ~ offset(logit_p3), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome4.offset <- glm(Y.fac4 ~ offset(logit_p4), family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Store output in a vector
    ## Intercept
    alpha.per.outcome <- c(calib.per.outcome1.offset$coefficients[1], 
                           calib.per.outcome2.offset$coefficients[1], 
                           calib.per.outcome3.offset$coefficients[1], 
                           calib.per.outcome4.offset$coefficients[1])
    names(alpha.per.outcome) <- paste("p", 1:K, sep = "")
    
    ## Slope
    beta.per.outcome <- c(calib.per.outcome1$coefficients[2], 
                          calib.per.outcome2$coefficients[2], 
                          calib.per.outcome3$coefficients[2], 
                          calib.per.outcome4$coefficients[2])
    names(beta.per.outcome) <- paste("p", 1:K, sep = "")
  } else if (K == 5){
    ### Create logit(p_k)
    dat.valid.pred$logit_p1 <- log(dat.valid.pred$p1/(1 - dat.valid.pred$p1))
    dat.valid.pred$logit_p2 <- log(dat.valid.pred$p2/(1 - dat.valid.pred$p2))
    dat.valid.pred$logit_p3 <- log(dat.valid.pred$p3/(1 - dat.valid.pred$p3))
    dat.valid.pred$logit_p4 <- log(dat.valid.pred$p4/(1 - dat.valid.pred$p4))
    dat.valid.pred$logit_p5 <- log(dat.valid.pred$p5/(1 - dat.valid.pred$p5))
    
    ### Fit recalibraiton model for slopes
    calib.per.outcome1 <- glm(Y.fac1 ~ logit_p1, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome2 <- glm(Y.fac2 ~ logit_p2, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome3 <- glm(Y.fac3 ~ logit_p3, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome4 <- glm(Y.fac4 ~ logit_p4, family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome5 <- glm(Y.fac5 ~ logit_p5, family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Fit recalibration model for intercepts, where slope is fixed to 1
    calib.per.outcome1.offset <- glm(Y.fac1 ~ offset(logit_p1), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome2.offset <- glm(Y.fac2 ~ offset(logit_p2), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome3.offset <- glm(Y.fac3 ~ offset(logit_p3), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome4.offset <- glm(Y.fac4 ~ offset(logit_p4), family = binomial(link = "logit"), data = dat.valid.pred)
    calib.per.outcome5.offset <- glm(Y.fac4 ~ offset(logit_p5), family = binomial(link = "logit"), data = dat.valid.pred)
    
    ### Store output in a vector
    ## Intercept
    alpha.per.outcome <- c(calib.per.outcome1.offset$coefficients[1], 
                           calib.per.outcome2.offset$coefficients[1], 
                           calib.per.outcome3.offset$coefficients[1], 
                           calib.per.outcome4.offset$coefficients[1], 
                           calib.per.outcome5.offset$coefficients[1])
    names(alpha.per.outcome) <- paste("p", 1:K, sep = "")
    
    ## Slope
    beta.per.outcome <- c(calib.per.outcome1$coefficients[2], 
                          calib.per.outcome2$coefficients[2], 
                          calib.per.outcome3$coefficients[2], 
                          calib.per.outcome4$coefficients[2], 
                          calib.per.outcome5$coefficients[2])
    names(beta.per.outcome) <- paste("p", 1:K, sep = "")
  }
  
  ### Create output object and return
  output.obj <- list("alpha" = alpha.per.outcome, "beta" = beta.per.outcome)
  return(output.obj)
}



#####################################################################################################################################
### 3.3) Function to produce flexible calibration scatter plots using the polytomous recalibraiton framework of van Hoorde (2014) ###
#####################################################################################################################################

### Now want to produce function to evaluate caibration (get vector of predicted observed risks, to go with the predicted risks)
calibrate.model.flexible.polytomous <- function(dat.valid.pred){
  
  ### dat.valid.pred contains all the predicted risks for each individual in the validation cohort
  ### Want to follow the flexible calibration approach outlined by van Hoorde et al
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  
  ### Start by calculating Zj, for j = 2, 3, ..., K (code assumes K will not be bigger than 5)
  dat.valid.pred$Z2 <- log(dat.valid.pred$p2/dat.valid.pred$p1)
  dat.valid.pred$Z3 <- log(dat.valid.pred$p3/dat.valid.pred$p1)
  
  if (K >= 4){
    dat.valid.pred$Z4 <- log(dat.valid.pred$p4/dat.valid.pred$p1)
  }
  if (K >= 5){
    dat.valid.pred$Z5 <- log(dat.valid.pred$p5/dat.valid.pred$p1)
  }
  
  ### Now fit a flexible recalibration model using vector spline smoother
  if (K == 3){
    flex.recal.model <- vgam(Y.fac ~ sm.ps(Z2, df = 4) + sm.ps(Z3, df = 4), 
                             data = dat.valid.pred, family = multinomial(refLevel = "1"))
  } else if (K == 4){
    flex.recal.model <- vgam(Y.fac ~ sm.ps(Z2, df = 4) + sm.ps(Z3, df = 4) + sm.ps(Z4, df = 4), 
                             data = dat.valid.pred, family = multinomial(refLevel = "1"))
  } else if (K == 5){
    flex.recal.model <- vgam(Y.fac ~ sm.ps(Z2, df = 4) + sm.ps(Z3, df = 4) + sm.ps(Z4, df = 4) + sm.ps(Z5, df = 4), 
                             data = dat.valid.pred, family = multinomial(refLevel = "1"))
  }
  
  ### Now generate predicted risks for each individual either A), in the validation cohort, or B) fo rsome vector of predicted risks
  dat.valid.pred.obs <- predict(flex.recal.model, newdata = dat.valid.pred, type = "response")
  colnames(dat.valid.pred.obs) <- paste("p.obs", 1:ncol(dat.valid.pred.obs), sep = "")
  
  ### Add to the validation dataset
  dat.valid.pred <- cbind(dat.valid.pred, dat.valid.pred.obs)
  
  ### Remove unncesary variables
  if (K == 3){
    dat.valid.pred <- select(dat.valid.pred, c(Y, Y.fac, p1, p2, p3, p.obs1, p.obs2, p.obs3))
  } else if (K == 4){
    dat.valid.pred <- select(dat.valid.pred, c(Y, Y.fac, p1, p2, p3, p4, p.obs1, p.obs2, p.obs3, p.obs4))
  } else if (K == 5){
    dat.valid.pred <- select(dat.valid.pred, c(Y, Y.fac, p1, p2, p3, p4, p5, p.obs1, p.obs2, p.obs3, p.obs4, p.obs5))
  }
  
  ### Calculate the ECI from this data
  if (K == 3){
    # Calc event rates
    K1 <- sum(dat.valid.pred$Y == 1)/nrow(dat.valid.pred)
    K2 <- sum(dat.valid.pred$Y == 2)/nrow(dat.valid.pred)
    K3 <- sum(dat.valid.pred$Y == 3)/nrow(dat.valid.pred)
    # Calc numerator
    ECI.numer <- (sum((dat.valid.pred$p1 - dat.valid.pred$p.obs1)^2) + sum((dat.valid.pred$p2 - dat.valid.pred$p.obs2)^2) +
                    sum((dat.valid.pred$p3 - dat.valid.pred$p.obs3)^2))
    # Calc denominator
    ECI.denom <- (sum((dat.valid.pred$p1 - K1)^2) + sum((dat.valid.pred$p2 - K2)^2) +
                    sum((dat.valid.pred$p3 - K3)^2))
    # Calc ECI
    ECI <- ECI.numer/ECI.denom
    
  } else if (K == 4){
    
    # Calc event rates
    K1 <- sum(dat.valid.pred$Y == 1)/nrow(dat.valid.pred)
    K2 <- sum(dat.valid.pred$Y == 2)/nrow(dat.valid.pred)
    K3 <- sum(dat.valid.pred$Y == 3)/nrow(dat.valid.pred)
    K4 <- sum(dat.valid.pred$Y == 4)/nrow(dat.valid.pred)
    # Calc numerator
    ECI.numer <- (sum((dat.valid.pred$p1 - dat.valid.pred$p.obs1)^2) + sum((dat.valid.pred$p2 - dat.valid.pred$p.obs2)^2) +
                    sum((dat.valid.pred$p3 - dat.valid.pred$p.obs3)^2) + sum((dat.valid.pred$p4 - dat.valid.pred$p.obs4)^2))
    # Calc denominator
    ECI.denom <- (sum((dat.valid.pred$p1 - K1)^2) + sum((dat.valid.pred$p2 - K2)^2) +
                    sum((dat.valid.pred$p3 - K3)^2) + sum((dat.valid.pred$p4 - K4)^2))
    # Calc ECI
    ECI <- ECI.numer/ECI.denom
    
  } else if (K == 5){
    
    # Calc event rates
    K1 <- sum(dat.valid.pred$Y == 1)/nrow(dat.valid.pred)
    K2 <- sum(dat.valid.pred$Y == 2)/nrow(dat.valid.pred)
    K3 <- sum(dat.valid.pred$Y == 3)/nrow(dat.valid.pred)
    K4 <- sum(dat.valid.pred$Y == 4)/nrow(dat.valid.pred)
    K5 <- sum(dat.valid.pred$Y == 5)/nrow(dat.valid.pred)
    # Calc numerator
    ECI.numer <- (sum((dat.valid.pred$p1 - dat.valid.pred$p.obs1)^2) + sum((dat.valid.pred$p2 - dat.valid.pred$p.obs2)^2) +
                    sum((dat.valid.pred$p3 - dat.valid.pred$p.obs3)^2) + sum((dat.valid.pred$p4 - dat.valid.pred$p.obs4)^2) +
                    sum((dat.valid.pred$p5 - dat.valid.pred$p.obs5)^2))
    # Calc denominator
    ECI.denom <- (sum((dat.valid.pred$p1 - K1)^2) + sum((dat.valid.pred$p2 - K2)^2) +
                    sum((dat.valid.pred$p3 - K3)^2) + sum((dat.valid.pred$p4 - K4)^2) +
                    sum((dat.valid.pred$p5 - K5)^2))
    # Calc ECI
    ECI <- ECI.numer/ECI.denom
    
  }
  
  return(list("calib.data" = dat.valid.pred, "ECI" = ECI))
}





###########################################################################################
### 3.4) Function to produce flexible calibration curvess for each outcome in isolation ###
###########################################################################################

### Now want to produce function to evaluate caibration (get vector of predicted observed risks, to go with the predicted risks)
calibrate.model.flexible.per.outcome <- function(dat.valid.pred){
  
  ### dat.valid.pred contains all the predicted risks for each individual in the validation cohort
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  
  ### Create variables which can be used to do the One vs All models and add to datasets
  Y.dummies <- model.matrix( ~ Y.fac - 1, dat.valid.pred)
  
  ### Combine with validation dataset
  dat.valid.pred <- cbind(dat.valid.pred, Y.dummies)
  
  ### Create logit(p_k) for each outcome category, fit calibration model, extract intercept/slope per outcome
  if (K == 3){
 
    ### Fit recalibraiton model with loess smoother
    loess.outcome1 <- loess(Y.fac1 ~ p1, data = dat.valid.pred, method = 'loess')
    loess.outcome2 <- loess(Y.fac2 ~ p2, data = dat.valid.pred, method = 'loess')
    loess.outcome3 <- loess(Y.fac3 ~ p3, data = dat.valid.pred, method = 'loess')
    
    ## Add to dataset
    dat.valid.pred$p.obs1 <- loess.outcome1$fitted
    dat.valid.pred$p.obs2 <- loess.outcome2$fitted
    dat.valid.pred$p.obs3 <- loess.outcome3$fitted
    
    ## Remove unncesary variables
    dat.valid.pred <- select(dat.valid.pred, c(Y, Y.fac, p1, p2, p3, p.obs1, p.obs2, p.obs3))
    
  } else if (K == 4){
    
    ### Fit recalibraiton model with loess smoother
    loess.outcome1 <- loess(Y.fac1 ~ p1, data = dat.valid.pred, method = 'loess')
    loess.outcome2 <- loess(Y.fac2 ~ p2, data = dat.valid.pred, method = 'loess')
    loess.outcome3 <- loess(Y.fac3 ~ p3, data = dat.valid.pred, method = 'loess')
    loess.outcome4 <- loess(Y.fac4 ~ p4, data = dat.valid.pred, method = 'loess')
    
    ## Add to dataset
    dat.valid.pred$p.obs1 <- loess.outcome1$fitted
    dat.valid.pred$p.obs2 <- loess.outcome2$fitted
    dat.valid.pred$p.obs3 <- loess.outcome3$fitted
    dat.valid.pred$p.obs4 <- loess.outcome4$fitted
    
    ## Remove unncesary variables
    dat.valid.pred <- select(dat.valid.pred, c(Y, Y.fac, p1, p2, p3, p4, p.obs1, p.obs2, p.obs3, p.obs4))
    
  } else if (K == 5){
    
    ### Fit recalibraiton model with loess smoother
    loess.outcome1 <- loess(Y.fac1 ~ p1, data = dat.valid.pred, method = 'loess')
    loess.outcome2 <- loess(Y.fac2 ~ p2, data = dat.valid.pred, method = 'loess')
    loess.outcome3 <- loess(Y.fac3 ~ p3, data = dat.valid.pred, method = 'loess')
    loess.outcome4 <- loess(Y.fac4 ~ p4, data = dat.valid.pred, method = 'loess')
    loess.outcome5 <- loess(Y.fac5 ~ p5, data = dat.valid.pred, method = 'loess')
    
    ## Add to dataset
    dat.valid.pred$p.obs1 <- loess.outcome1$fitted
    dat.valid.pred$p.obs2 <- loess.outcome2$fitted
    dat.valid.pred$p.obs3 <- loess.outcome3$fitted
    dat.valid.pred$p.obs4 <- loess.outcome4$fitted
    dat.valid.pred$p.obs5 <- loess.outcome5$fitted
    
    ## Remove unncesary variables
    dat.valid.pred <- select(dat.valid.pred, c(Y, Y.fac, p1, p2, p3, p4, p5, p.obs1, p.obs2, p.obs3, p.obs4, p.obs5))
  }
  
  ## Return output object
  return(dat.valid.pred)
}


#####################################################################################
### 3.5) Function to produce observed vs expected ratios in each outcome category ###
#####################################################################################

### We calculate the in the validation cohort, the average predicted risk for each category, 
### and then the observd risk as the proportion of individuals who had that event
calibrate.model.expected.observed <- function(dat.valid.pred){
  
  ### dat.valid.pred contains all the predicted risks for each individual in the validation cohort
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  
  ### Create a vector as output object
  output.object <- rep(NA, K)
  names(output.object) <- paste("p", 1:K, sep = "")
  
  ### Run through outcome categories
  for (i in 1:K){
    ## Calc obs
    observed <- sum(dat.valid.pred$Y == i)/length(dat.valid.pred$Y)
    ## Calc expected
    expected <- mean(dat.valid.pred[,paste("p", i, sep = "")])
    ## Calc ratio
    output.object[i] <- observed/expected
  }
  
  ### Return output object
  return(output.object)
}


#########################################################################################################################
### 3.6) Function to calculate PDI. Taken from literature. Referenced in manuscript and more details on authors below ###
### I have added a few lines at the top on the function to make the dataset from these simulations suitable for use in this function.
#########################################################################################################################

#2020-11-20
###########################################################
#Anamaria Savu, PhD                                                                                         #
#Canadian VIGOUR Centre                                                                                #
#University of Alberta                                                                                        #
# savu@ualberta.ca                                                                                            #
#November 2020                                                                                                #
###########################################################
#PDI – polytomous discrimination index                                                        #
###########################################################

###############################################################
# data = input dataset of data.frame type
#             required
#              
# nbs = number of bootstrap samples to perform of numeric type
#           required
###############################################################
#  Usage notes
#  data must have one of its components named $outcome for storing 
# the outcome values. The values of the outcome must cover all integer 
# values between a lower and an upper bound. The remaining components 
# of the data frame # ($p1, $p2, …) must record the estimated probabilities
# in a specific order: $p1 probabilities for the lowest value of the outcome, 
# $p2 probabilities for the second lowest value of the outcome, and so on.
###############################################################


#Estimates of PDI and its components
pdiest<-function(dat.valid.pred){
  
  ### My additions
  ## Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  ## Create outcome variable
  dat.valid.pred$outcome <- dat.valid.pred$Y
  ## Retain only neccesary variables and assign to dataset called 'data'
  if (K == 3){
    data <- dplyr::select(dat.valid.pred, c(outcome, p1, p2, p3))
  } else if (K == 4){
    data <- dplyr::select(dat.valid.pred, c(outcome, p1, p2, p3, p4))
  } else if (K == 5){
    data <- dplyr::select(dat.valid.pred, c(outcome, p1, p2, p3, p4, p5))
  }
  ## Reduce to 100,000 observations max for computational reasons
  data <- data[1:min(nrow(data), 100000), ]
  
  ### ORIGINAL CODE
  y<-data$outcome
  ymin<-min(y)
  ymax<-max(y)
  noutcome<-ymax-ymin
  p<-prod(table(y))
  pdi<-c()
  
  for (i in 1:(noutcome+1)){
    
    predprob<-data[,(i+1)]  #READ predicted probabilities for level i
    t0<-table(predprob,y)   #CALCULATE frequencies of predicted probabilities for level i by outcome
    
    dim1<-dim(t0)[1]
    dim2<-dim(t0)[2]
    t<-cbind(t0[,i],t0[,-i]) #REORDER columns
    restrictt<- if (noutcome == 1){matrix(t[,2:(noutcome+1)],ncol=1)} else {t[,2:(noutcome+1)] } #REMOVE first column of t
    
    c<-apply(restrictt,2,cumsum) #CALCULATE cumulative frequencies of predicted probabilities for level i by outcome
    cnew<- if (noutcome == 1) {rbind(rep(0,noutcome),matrix(c[1:(dim(c)[1]-1),],ncol=))} else {rbind(rep(0,noutcome),c[1:(dim(c)[1]-1),])} #INTRODUCE a row of zeros at the begining of c
    
    mat<-c()                     #MATRIX of 0s and 1s of dimension 2^(noutcome) x noutcome
    for (j in 1:noutcome){
      mat0<-cbind(mat,0)
      mat1<-cbind(mat,1)
      mat<-rbind(mat0,mat1)}
    
    r<-0
    for (k in 1:dim(mat)[1]){
      dt<-t(apply(restrictt, 1, function(x) mat[k,]*x))
      dcnew<-t(apply(cnew, 1, function(x) (1-mat[k,])*x))
      dfinal<-if (noutcome == 1) {cbind(t[,1],t(dt+dcnew))} else {cbind(t[,1],dt+dcnew)} #TAKE all combinations of frequencies and cumulative frequencies
      r<-r+sum(apply(dfinal,1,prod))/(1+sum(mat[k,]))}                                   #MULTIPLYIES across rows
    
    r<-r/p     #PDI component for outcome i
    pdi<-rbind(pdi,r)
  }
  pdi<-rbind(mean(pdi),pdi)
  
  ## Only return the PDI (rather than category specific)
  return(as.numeric(pdi[1]))}

# #Estimates and bootstrap 95% confidence intervals for PDI and its components
# pdifunction<-function(data,nbs){
#   #PDI estimate
#   estimate<-pdiest(data)
#   #BOOTSTRAP
#   samplesize<-dim(data)[1]
#   for (i in 1:nbs)
#   {vec<-sample.int(samplesize,size=samplesize, replace=TRUE)
#    mydatabs<-data[vec,]
#    if (i<2) {pdibs<-pdiest(mydatabs)
#    } else {
#      pdibs<-cbind(pdibs,pdiest(mydatabs))
#    }
#   }
#   
#   stderr <- sqrt(apply(pdibs, 1, var))
#   lowerci<- pmax(0, estimate - 1.96*stderr)
#   upperci<- pmin(1, estimate + 1.96*stderr)
#   
#   estci<-cbind(estimate, lowerci, upperci)
#   estci
# }



###########################################################################
### 3.7) Function to calculate C-statistic for each outcome vs the rest ###
###########################################################################

### We calculate the in the validation cohort, the average predicted risk for each category, 
### and then the observd risk as the proportion of individuals who had that event
calc.Cstat.per.outcome <- function(dat.valid.pred){
  
  ### dat.valid.pred contains all the predicted risks for each individual in the validation cohort
  
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(dat.valid.pred$Y))
  
  ### Create variables which can be used to do the One vs All models and add to datasets
  Y.dummies <- model.matrix( ~ Y.fac - 1, dat.valid.pred)
  
  ### Combine with validation dataset
  dat.valid.pred <- cbind(dat.valid.pred, Y.dummies)
  
  ### Create a vector as output object
  output.object <- rep(NA, K)
  names(output.object) <- paste("p", 1:K, sep = "")
  
  Cstat <- Cstat(dat.valid.pred[ ,paste("p", 1, sep = "")], dat.valid.pred[ ,paste("Y.fac", 1, sep = "")])
  
  ### Run through outcome categories
  for (i in 1:K){
    ## Calc Cstatistic and assign
    output.object[i] <- Cstat(dat.valid.pred[ ,paste("p", i, sep = "")], dat.valid.pred[ ,paste("Y.fac", i, sep = "")])
  }
  
  ### Return output object
  return(output.object)
}


##############################################################
### 3.7A) Function to run simulation for large sample size ###
##############################################################

run.sim.large.sample.OLD <- function(dat.devel){
  
  pred.multinomial <- fit.model.multinomial(dat.devel, dat.devel)
  pred.seqlog <- fit.model.seqlog(dat.devel, dat.devel)
  pred.OvA <- fit.model.OvA(dat.devel, dat.devel)
  pred.OvO.PC <- fit.model.OvO.PC(dat.devel, dat.devel)
  
#   wc.ms.multinomial <- calibrate.mod.specific.multinomial(pred.multinomial)
#   wc.ms.seqlog <- calibrate.mod.specific.seqlog(pred.seqlog)
#   wc.ms.OvA <- calibrate.mod.specific.OvA(pred.OvA)
#   wc.ms.OvO.PC <- calibrate.mod.specific.OvO.PC(pred.OvO.PC)
  
  wc.po.multinomial <- calibrate.weak.per.outcome(pred.multinomial)
  wc.po.seqlog <- calibrate.weak.per.outcome(pred.seqlog)
  wc.po.OvA <- calibrate.weak.per.outcome(pred.OvA)
  wc.po.OvO.PC <- calibrate.weak.per.outcome(pred.OvO.PC)
  
  OvE.multinomial <- calibrate.model.expected.observed(pred.multinomial)
  OvE.seqlog <- calibrate.model.expected.observed(pred.seqlog)
  OvE.OvA <- calibrate.model.expected.observed(pred.OvA)
  OvE.OvO.PC <- calibrate.model.expected.observed(pred.OvO.PC)
  
  Cstat.multinomial <- calc.Cstat.per.outcome(pred.multinomial)
  Cstat.seqlog <- calc.Cstat.per.outcome(pred.seqlog)
  Cstat.OvA <- calc.Cstat.per.outcome(pred.OvA)
  Cstat.OvO.PC <- calc.Cstat.per.outcome(pred.OvO.PC)
  
  PDI.multinomial <- pdiest(pred.multinomial)
  PDI.seqlog <- pdiest(pred.seqlog)
  PDI.OvA <- pdiest(pred.OvA)
  PDI.OvO.PC <- pdiest(pred.OvO.PC)
  
  calib.flex.mlr.multinomial <- calibrate.model.flexible.polytomous(pred.multinomial)
  calib.flex.mlr.seqlog <- calibrate.model.flexible.polytomous(pred.seqlog)
  calib.flex.mlr.OvA <- calibrate.model.flexible.polytomous(pred.OvA)
  calib.flex.mlr.OvO.PC <- calibrate.model.flexible.polytomous(pred.OvO.PC)
  
  calib.flex.po.multinomial <- calibrate.model.flexible.per.outcome(pred.multinomial)
  calib.flex.po.seqlog <- calibrate.model.flexible.per.outcome(pred.seqlog)
  calib.flex.po.OvA <- calibrate.model.flexible.per.outcome(pred.OvA)
  calib.flex.po.OvO.PC <- calibrate.model.flexible.per.outcome(pred.OvO.PC)
  
  return(list(#"wc.ms.multinomial" = wc.ms.multinomial, "wc.ms.seqlog" = wc.ms.seqlog, "wc.ms.OvA" = wc.ms.OvA, "wc.ms.OvO.PC" = wc.ms.OvO.PC,
              "wc.po.multinomial" = wc.po.multinomial, "wc.po.seqlog" = wc.po.seqlog, "wc.po.OvA" = wc.po.OvA, "wc.po.OvO.PC" = wc.po.OvO.PC,
              "OvE.multinomial" = OvE.multinomial, "OvE.seqlog" = OvE.seqlog, "OvE.OvA" = OvE.OvA, "OvE.OvO.PC" = OvE.OvO.PC,
              "Cstat.multinomial" = Cstat.multinomial, "Cstat.seqlog" = Cstat.seqlog, "Cstat.OvA" = Cstat.OvA, "Cstat.OvO.PC" = Cstat.OvO.PC,
              "PDI.multinomial" = PDI.multinomial, "PDI.seqlog" = PDI.seqlog, "PDI.OvA" = PDI.OvA, "PDI.OvO.PC" = PDI.OvO.PC,
              "calib.flex.mlr.multinomial" = calib.flex.mlr.multinomial, "calib.flex.mlr.seqlog" = calib.flex.mlr.seqlog, 
              "calib.flex.mlr.OvA" = calib.flex.mlr.OvA, "calib.flex.mlr.OvO.PC" = calib.flex.mlr.OvO.PC,
              "calib.flex.po.multinomial" = calib.flex.po.multinomial, "calib.flex.po.seqlog" = calib.flex.po.seqlog, 
              "calib.flex.po.OvA" = calib.flex.po.OvA, "calib.flex.po.OvO.PC" = calib.flex.po.OvO.PC
              ))
}


#############################################################################################################################
### 3.7B) Function to run simulation for large sample size, including doing sequential logistic with every possible order ###
#############################################################################################################################

run.sim.large.sample <- function(dat.devel){
  
  print(paste("multinomial", Sys.time()))
  pred.multinomial <- fit.model.multinomial(dat.devel, dat.devel)
  print(paste("seqlog", Sys.time()))
  pred.seqlog <- fit.model.seqlog(dat.devel, dat.devel)
  print(paste("seqlog.all", Sys.time()))
  pred.seqlog.all <- fit.model.seqlog.all(dat.devel, dat.devel)
  print(paste("OvA", Sys.time()))
  pred.OvA <- fit.model.OvA(dat.devel, dat.devel)
  print(paste("OvO.PC", Sys.time()))
  pred.OvO.PC <- fit.model.OvO.PC(dat.devel, dat.devel)
  
#   print(paste("wc.ms", Sys.time()))
#   wc.ms.multinomial <- calibrate.mod.specific.multinomial(pred.multinomial)
#   wc.ms.seqlog <- calibrate.mod.specific.seqlog(pred.seqlog)
#   wc.ms.OvA <- calibrate.mod.specific.OvA(pred.OvA)
#   wc.ms.OvO.PC <- calibrate.mod.specific.OvO.PC(pred.OvO.PC)
  
  print(paste("wc.po", Sys.time()))
  wc.po.multinomial <- calibrate.weak.per.outcome(pred.multinomial)
  wc.po.seqlog <- calibrate.weak.per.outcome(pred.seqlog)
  wc.po.seqlog.all <- calibrate.weak.per.outcome(pred.seqlog.all)
  wc.po.OvA <- calibrate.weak.per.outcome(pred.OvA)
  wc.po.OvO.PC <- calibrate.weak.per.outcome(pred.OvO.PC)
  
  print(paste("OvE", Sys.time()))
  OvE.multinomial <- calibrate.model.expected.observed(pred.multinomial)
  OvE.seqlog <- calibrate.model.expected.observed(pred.seqlog)
  OvE.seqlog.all <- calibrate.model.expected.observed(pred.seqlog.all)
  OvE.OvA <- calibrate.model.expected.observed(pred.OvA)
  OvE.OvO.PC <- calibrate.model.expected.observed(pred.OvO.PC)
  
  print(paste("Cstat", Sys.time()))
  Cstat.multinomial <- calc.Cstat.per.outcome(pred.multinomial)
  Cstat.seqlog <- calc.Cstat.per.outcome(pred.seqlog)
  Cstat.seqlog.all <- calc.Cstat.per.outcome(pred.seqlog.all)
  Cstat.OvA <- calc.Cstat.per.outcome(pred.OvA)
  Cstat.OvO.PC <- calc.Cstat.per.outcome(pred.OvO.PC)
  
  print(paste("PDI", Sys.time()))
  PDI.multinomial <- pdiest(pred.multinomial)
  PDI.seqlog <- pdiest(pred.seqlog)
  PDI.seqlog.all <- pdiest(pred.seqlog.all)
  PDI.OvA <- pdiest(pred.OvA)
  PDI.OvO.PC <- pdiest(pred.OvO.PC)
  
  print(paste("calib.flex.mlr", Sys.time()))
  calib.flex.mlr.multinomial <- calibrate.model.flexible.polytomous(pred.multinomial)
  calib.flex.mlr.seqlog <- calibrate.model.flexible.polytomous(pred.seqlog)
  calib.flex.mlr.seqlog.all <- calibrate.model.flexible.polytomous(pred.seqlog.all)
  calib.flex.mlr.OvA <- calibrate.model.flexible.polytomous(pred.OvA)
  calib.flex.mlr.OvO.PC <- calibrate.model.flexible.polytomous(pred.OvO.PC)
  
  print(paste("calib.flex.po", Sys.time()))
  calib.flex.po.multinomial <- calibrate.model.flexible.per.outcome(pred.multinomial)
  calib.flex.po.seqlog <- calibrate.model.flexible.per.outcome(pred.seqlog)
  calib.flex.po.seqlog.all <- calibrate.model.flexible.per.outcome(pred.seqlog.all)
  calib.flex.po.OvA <- calibrate.model.flexible.per.outcome(pred.OvA)
  calib.flex.po.OvO.PC <- calibrate.model.flexible.per.outcome(pred.OvO.PC)
  
  return(list(#"wc.ms.multinomial" = wc.ms.multinomial, "wc.ms.seqlog" = wc.ms.seqlog, "wc.ms.OvA" = wc.ms.OvA, "wc.ms.OvO.PC" = wc.ms.OvO.PC,
              "wc.po.multinomial" = wc.po.multinomial, "wc.po.seqlog" = wc.po.seqlog, "wc.po.seqlog.all" = wc.po.seqlog.all, "wc.po.OvA" = wc.po.OvA, "wc.po.OvO.PC" = wc.po.OvO.PC,
              "OvE.multinomial" = OvE.multinomial, "OvE.seqlog" = OvE.seqlog, "OvE.seqlog.all" = OvE.seqlog.all, "OvE.OvA" = OvE.OvA, "OvE.OvO.PC" = OvE.OvO.PC,
              "Cstat.multinomial" = Cstat.multinomial, "Cstat.seqlog" = Cstat.seqlog, "Cstat.seqlog.all" = Cstat.seqlog.all, "Cstat.OvA" = Cstat.OvA, "Cstat.OvO.PC" = Cstat.OvO.PC,
              "PDI.multinomial" = PDI.multinomial, "PDI.seqlog" = PDI.seqlog, "PDI.seqlog.all" = PDI.seqlog.all, "PDI.OvA" = PDI.OvA, "PDI.OvO.PC" = PDI.OvO.PC,
              "calib.flex.mlr.multinomial" = calib.flex.mlr.multinomial, "calib.flex.mlr.seqlog" = calib.flex.mlr.seqlog,
              "calib.flex.mlr.seqlog.all" = calib.flex.mlr.seqlog.all, 
              "calib.flex.mlr.OvA" = calib.flex.mlr.OvA, "calib.flex.mlr.OvO.PC" = calib.flex.mlr.OvO.PC,
              "calib.flex.po.multinomial" = calib.flex.po.multinomial, "calib.flex.po.seqlog" = calib.flex.po.seqlog, 
              "calib.flex.po.seqlog.all" = calib.flex.po.seqlog.all, 
              "calib.flex.po.OvA" = calib.flex.po.OvA, "calib.flex.po.OvO.PC" = calib.flex.po.OvO.PC
  ))
}



#########################################################################
### 3.7C) Function to run simulation for large sample size, for K = 5 ###
#########################################################################

run.sim.large.sample.K5.OLD <- function(dat.devel){
  
  pred.multinomial <- fit.model.multinomial(dat.devel, dat.devel)
  pred.seqlog <- fit.model.seqlog(dat.devel, dat.devel)
  pred.OvA <- fit.model.OvA(dat.devel, dat.devel)
  pred.OvO.PC <- fit.model.OvO.PC(dat.devel, dat.devel)
#   
#   wc.ms.multinomial <- calibrate.mod.specific.multinomial(pred.multinomial)
#   wc.ms.seqlog <- calibrate.mod.specific.seqlog(pred.seqlog)
#   wc.ms.OvA <- calibrate.mod.specific.OvA(pred.OvA)
#   wc.ms.OvO.PC <- calibrate.mod.specific.OvO.PC(pred.OvO.PC)
#   
  wc.po.multinomial <- calibrate.weak.per.outcome(pred.multinomial)
  wc.po.seqlog <- calibrate.weak.per.outcome(pred.seqlog)
  wc.po.OvA <- calibrate.weak.per.outcome(pred.OvA)
  wc.po.OvO.PC <- calibrate.weak.per.outcome(pred.OvO.PC)
#   
  OvE.multinomial <- calibrate.model.expected.observed(pred.multinomial)
  OvE.seqlog <- calibrate.model.expected.observed(pred.seqlog)
  OvE.OvA <- calibrate.model.expected.observed(pred.OvA)
  OvE.OvO.PC <- calibrate.model.expected.observed(pred.OvO.PC)
#   
  Cstat.multinomial <- calc.Cstat.per.outcome(pred.multinomial)
  Cstat.seqlog <- calc.Cstat.per.outcome(pred.seqlog)
  Cstat.OvA <- calc.Cstat.per.outcome(pred.OvA)
  Cstat.OvO.PC <- calc.Cstat.per.outcome(pred.OvO.PC)
  
  PDI.multinomial <- pdiest(pred.multinomial)
  PDI.seqlog <- pdiest(pred.seqlog)
  PDI.OvA <- pdiest(pred.OvA)
  PDI.OvO.PC <- pdiest(pred.OvO.PC)
  
  calib.flex.mlr.multinomial <- calibrate.model.flexible.polytomous(pred.multinomial)
  calib.flex.mlr.seqlog <- calibrate.model.flexible.polytomous(pred.seqlog)
  calib.flex.mlr.OvA <- calibrate.model.flexible.polytomous(pred.OvA)
  calib.flex.mlr.OvO.PC <- calibrate.model.flexible.polytomous(pred.OvO.PC)
  
  calib.flex.po.multinomial <- calibrate.model.flexible.per.outcome(pred.multinomial)
  calib.flex.po.seqlog <- calibrate.model.flexible.per.outcome(pred.seqlog)
  calib.flex.po.OvA <- calibrate.model.flexible.per.outcome(pred.OvA)
  calib.flex.po.OvO.PC <- calibrate.model.flexible.per.outcome(pred.OvO.PC)
  
  return(list(#"wc.ms.multinomial" = wc.ms.multinomial, "wc.ms.seqlog" = wc.ms.seqlog, "wc.ms.OvA" = wc.ms.OvA, "wc.ms.OvO.PC" = wc.ms.OvO.PC,
              "wc.po.multinomial" = wc.po.multinomial, "wc.po.seqlog" = wc.po.seqlog, "wc.po.OvA" = wc.po.OvA, "wc.po.OvO.PC" = wc.po.OvO.PC,
              "OvE.multinomial" = OvE.multinomial, "OvE.seqlog" = OvE.seqlog, "OvE.OvA" = OvE.OvA, "OvE.OvO.PC" = OvE.OvO.PC,
              "Cstat.multinomial" = Cstat.multinomial, "Cstat.seqlog" = Cstat.seqlog, "Cstat.OvA" = Cstat.OvA, "Cstat.OvO.PC" = Cstat.OvO.PC,
              "PDI.multinomial" = PDI.multinomial, "PDI.seqlog" = PDI.seqlog, "PDI.OvA" = PDI.OvA, "PDI.OvO.PC" = PDI.OvO.PC,
              "calib.flex.mlr.multinomial" = calib.flex.mlr.multinomial, "calib.flex.mlr.seqlog" = calib.flex.mlr.seqlog, 
              "calib.flex.mlr.OvA" = calib.flex.mlr.OvA, "calib.flex.mlr.OvO.PC" = calib.flex.mlr.OvO.PC,
              "calib.flex.po.multinomial" = calib.flex.po.multinomial, "calib.flex.po.seqlog" = calib.flex.po.seqlog, 
              "calib.flex.po.OvA" = calib.flex.po.OvA, "calib.flex.po.OvO.PC" = calib.flex.po.OvO.PC
  ))
}


#############################################################################################################################
### 3.7D) Function to run simulation for large sample size, including doing sequential logistic with every possible order ###
### For K = 5
#############################################################################################################################

run.sim.large.sample.K5 <- function(dat.devel){
  
  print(paste("multinomial", Sys.time()))
  pred.multinomial <- fit.model.multinomial(dat.devel, dat.devel)
  print(paste("seqlog", Sys.time()))
  pred.seqlog <- fit.model.seqlog(dat.devel, dat.devel)
  print(paste("seqlog.all", Sys.time()))
  pred.seqlog.all <- fit.model.seqlog.all(dat.devel, dat.devel)
  print(paste("OvA", Sys.time()))
  pred.OvA <- fit.model.OvA(dat.devel, dat.devel)
  print(paste("OvO.PC", Sys.time()))
  pred.OvO.PC <- fit.model.OvO.PC(dat.devel, dat.devel)
  
#   print(paste("wc.ms", Sys.time()))
#   wc.ms.multinomial <- calibrate.mod.specific.multinomial(pred.multinomial)
#   wc.ms.seqlog <- calibrate.mod.specific.seqlog(pred.seqlog)
#   wc.ms.OvA <- calibrate.mod.specific.OvA(pred.OvA)
#   wc.ms.OvO.PC <- calibrate.mod.specific.OvO.PC(pred.OvO.PC)
  
  print(paste("wc.po", Sys.time()))
  wc.po.multinomial <- calibrate.weak.per.outcome(pred.multinomial)
  wc.po.seqlog <- calibrate.weak.per.outcome(pred.seqlog)
  wc.po.seqlog.all <- calibrate.weak.per.outcome(pred.seqlog.all)
  wc.po.OvA <- calibrate.weak.per.outcome(pred.OvA)
  wc.po.OvO.PC <- calibrate.weak.per.outcome(pred.OvO.PC)
  
  print(paste("OvE", Sys.time()))
  OvE.multinomial <- calibrate.model.expected.observed(pred.multinomial)
  OvE.seqlog <- calibrate.model.expected.observed(pred.seqlog)
  OvE.seqlog.all <- calibrate.model.expected.observed(pred.seqlog.all)
  OvE.OvA <- calibrate.model.expected.observed(pred.OvA)
  OvE.OvO.PC <- calibrate.model.expected.observed(pred.OvO.PC)
  
  print(paste("Cstat", Sys.time()))
  Cstat.multinomial <- calc.Cstat.per.outcome(pred.multinomial)
  Cstat.seqlog <- calc.Cstat.per.outcome(pred.seqlog)
  Cstat.seqlog.all <- calc.Cstat.per.outcome(pred.seqlog.all)
  Cstat.OvA <- calc.Cstat.per.outcome(pred.OvA)
  Cstat.OvO.PC <- calc.Cstat.per.outcome(pred.OvO.PC)
  
  print(paste("PDI", Sys.time()))
  PDI.multinomial <- pdiest(pred.multinomial)
  PDI.seqlog <- pdiest(pred.seqlog)
  PDI.seqlog.all <- pdiest(pred.seqlog.all)
  PDI.OvA <- pdiest(pred.OvA)
  PDI.OvO.PC <- pdiest(pred.OvO.PC)
  
  print(paste("calib.flex.mlr", Sys.time()))
  calib.flex.mlr.multinomial <- calibrate.model.flexible.polytomous(pred.multinomial)
  calib.flex.mlr.seqlog <- calibrate.model.flexible.polytomous(pred.seqlog)
  calib.flex.mlr.seqlog.all <- calibrate.model.flexible.polytomous(pred.seqlog.all)
  calib.flex.mlr.OvA <- calibrate.model.flexible.polytomous(pred.OvA)
  calib.flex.mlr.OvO.PC <- calibrate.model.flexible.polytomous(pred.OvO.PC)
  
  print(paste("calib.flex.po", Sys.time()))
  calib.flex.po.multinomial <- calibrate.model.flexible.per.outcome(pred.multinomial)
  calib.flex.po.seqlog <- calibrate.model.flexible.per.outcome(pred.seqlog)
  calib.flex.po.seqlog.all <- calibrate.model.flexible.per.outcome(pred.seqlog.all)
  calib.flex.po.OvA <- calibrate.model.flexible.per.outcome(pred.OvA)
  calib.flex.po.OvO.PC <- calibrate.model.flexible.per.outcome(pred.OvO.PC)
  
  return(list("wc.ms.multinomial" = wc.ms.multinomial, "wc.ms.seqlog" = wc.ms.seqlog, "wc.ms.OvA" = wc.ms.OvA, "wc.ms.OvO.PC" = wc.ms.OvO.PC,
              "wc.po.multinomial" = wc.po.multinomial, "wc.po.seqlog" = wc.po.seqlog, "wc.po.seqlog.all" = wc.po.seqlog.all, "wc.po.OvA" = wc.po.OvA, "wc.po.OvO.PC" = wc.po.OvO.PC,
              "OvE.multinomial" = OvE.multinomial, "OvE.seqlog" = OvE.seqlog, "OvE.seqlog.all" = OvE.seqlog.all, "OvE.OvA" = OvE.OvA, "OvE.OvO.PC" = OvE.OvO.PC,
              "Cstat.multinomial" = Cstat.multinomial, "Cstat.seqlog" = Cstat.seqlog, "Cstat.seqlog.all" = Cstat.seqlog.all, "Cstat.OvA" = Cstat.OvA, "Cstat.OvO.PC" = Cstat.OvO.PC,
              "PDI.multinomial" = PDI.multinomial, "PDI.seqlog" = PDI.seqlog, "PDI.seqlog.all" = PDI.seqlog.all, "PDI.OvA" = PDI.OvA, "PDI.OvO.PC" = PDI.OvO.PC,
              "calib.flex.mlr.multinomial" = calib.flex.mlr.multinomial, "calib.flex.mlr.seqlog" = calib.flex.mlr.seqlog,
              "calib.flex.mlr.seqlog.all" = calib.flex.mlr.seqlog.all, 
              "calib.flex.mlr.OvA" = calib.flex.mlr.OvA, "calib.flex.mlr.OvO.PC" = calib.flex.mlr.OvO.PC,
              "calib.flex.po.multinomial" = calib.flex.po.multinomial, "calib.flex.po.seqlog" = calib.flex.po.seqlog, 
              "calib.flex.po.seqlog.all" = calib.flex.po.seqlog.all, 
              "calib.flex.po.OvA" = calib.flex.po.OvA, "calib.flex.po.OvO.PC" = calib.flex.po.OvO.PC
  ))
}

#############################################################
### 3.8A) Function to run simulation for small sample size ###
#############################################################

run.sim.small.sample.OLD <- function(dat.devel, dat.valid){
  
  pred.multinomial <- fit.model.multinomial(dat.devel, dat.valid)
  pred.seqlog <- fit.model.seqlog(dat.devel, dat.valid)
  pred.OvA <- fit.model.OvA(dat.devel, dat.valid)
  pred.OvO.PC <- fit.model.OvO.PC(dat.devel, dat.valid)
  
  wc.ms.multinomial <- calibrate.mod.specific.multinomial(pred.multinomial)
  wc.ms.seqlog <- calibrate.mod.specific.seqlog(pred.seqlog)
  wc.ms.OvA <- calibrate.mod.specific.OvA(pred.OvA)
  wc.ms.OvO.PC <- calibrate.mod.specific.OvO.PC(pred.OvO.PC)
  
  wc.po.multinomial <- calibrate.weak.per.outcome(pred.multinomial)
  wc.po.seqlog <- calibrate.weak.per.outcome(pred.seqlog)
  wc.po.OvA <- calibrate.weak.per.outcome(pred.OvA)
  wc.po.OvO.PC <- calibrate.weak.per.outcome(pred.OvO.PC)
  
  OvE.multinomial <- calibrate.model.expected.observed(pred.multinomial)
  OvE.seqlog <- calibrate.model.expected.observed(pred.seqlog)
  OvE.OvA <- calibrate.model.expected.observed(pred.OvA)
  OvE.OvO.PC <- calibrate.model.expected.observed(pred.OvO.PC)
  
  Cstat.multinomial <- calc.Cstat.per.outcome(pred.multinomial)
  Cstat.seqlog <- calc.Cstat.per.outcome(pred.seqlog)
  Cstat.OvA <- calc.Cstat.per.outcome(pred.OvA)
  Cstat.OvO.PC <- calc.Cstat.per.outcome(pred.OvO.PC)
  
  PDI.multinomial <- pdiest(pred.multinomial)
  PDI.seqlog <- pdiest(pred.seqlog)
  PDI.OvA <- pdiest(pred.OvA)
  PDI.OvO.PC <- pdiest(pred.OvO.PC)
  
  ECI.multinomial <- calibrate.model.flexible.polytomous(pred.multinomial)$ECI
  ECI.seqlog <- calibrate.model.flexible.polytomous(pred.seqlog)$ECI
  ECI.OvA <- calibrate.model.flexible.polytomous(pred.OvA)$ECI
  ECI.OvO.PC <- calibrate.model.flexible.polytomous(pred.OvO.PC)$ECI
    
  return(list("wc.ms.multinomial" = wc.ms.multinomial, "wc.ms.seqlog" = wc.ms.seqlog, "wc.ms.OvA" = wc.ms.OvA, "wc.ms.OvO.PC" = wc.ms.OvO.PC,
              "wc.po.multinomial" = wc.po.multinomial, "wc.po.seqlog" = wc.po.seqlog, "wc.po.OvA" = wc.po.OvA, "wc.po.OvO.PC" = wc.po.OvO.PC,
              "OvE.multinomial" = OvE.multinomial, "OvE.seqlog" = OvE.seqlog, "OvE.OvA" = OvE.OvA, "OvE.OvO.PC" = OvE.OvO.PC,
              "Cstat.multinomial" = Cstat.multinomial, "Cstat.seqlog" = Cstat.seqlog, "Cstat.OvA" = Cstat.OvA, "Cstat.OvO.PC" = Cstat.OvO.PC,
              "PDI.multinomial" = PDI.multinomial, "PDI.seqlog" = PDI.seqlog, "PDI.OvA" = PDI.OvA, "PDI.OvO.PC" = PDI.OvO.PC,
              "ECI.multinomial" = ECI.multinomial, "ECI.seqlog" = ECI.seqlog, "ECI.OvA" = ECI.OvA, "ECI.OvO.PC" = ECI.OvO.PC
  ))
}


#############################################################
### 3.8B) Function to run simulation for small sample size ###
### Includes seqlog.all
#############################################################

run.sim.small.sample <- function(dat.devel, dat.valid){
  
  print(paste("multinomial", Sys.time()))
  pred.multinomial <- fit.model.multinomial(dat.devel, dat.valid)
  print(paste("seqlog", Sys.time()))
  pred.seqlog <- fit.model.seqlog(dat.devel, dat.valid)
  print(paste("seqlog.all", Sys.time()))
  pred.seqlog.all <- fit.model.seqlog.all(dat.devel, dat.valid)
  print(paste("OvA", Sys.time()))
  pred.OvA <- fit.model.OvA(dat.devel, dat.valid)
  print(paste("OvO.PC", Sys.time()))
  pred.OvO.PC <- fit.model.OvO.PC(dat.devel, dat.valid)
  
  print(paste("wc.ms", Sys.time()))
  wc.ms.multinomial <- calibrate.mod.specific.multinomial(pred.multinomial)
  wc.ms.seqlog <- calibrate.mod.specific.seqlog(pred.seqlog)
  wc.ms.OvA <- calibrate.mod.specific.OvA(pred.OvA)
  wc.ms.OvO.PC <- calibrate.mod.specific.OvO.PC(pred.OvO.PC)
  
  print(paste("wc.po", Sys.time()))
  wc.po.multinomial <- calibrate.weak.per.outcome(pred.multinomial)
  wc.po.seqlog <- calibrate.weak.per.outcome(pred.seqlog)
  wc.po.seqlog.all <- calibrate.weak.per.outcome(pred.seqlog.all)
  wc.po.OvA <- calibrate.weak.per.outcome(pred.OvA)
  wc.po.OvO.PC <- calibrate.weak.per.outcome(pred.OvO.PC)
  
  print(paste("OvE", Sys.time()))
  OvE.multinomial <- calibrate.model.expected.observed(pred.multinomial)
  OvE.seqlog <- calibrate.model.expected.observed(pred.seqlog)
  OvE.seqlog.all <- calibrate.model.expected.observed(pred.seqlog.all)
  OvE.OvA <- calibrate.model.expected.observed(pred.OvA)
  OvE.OvO.PC <- calibrate.model.expected.observed(pred.OvO.PC)
  
  print(paste("Cstat", Sys.time()))
  Cstat.multinomial <- calc.Cstat.per.outcome(pred.multinomial)
  Cstat.seqlog <- calc.Cstat.per.outcome(pred.seqlog)
  Cstat.seqlog.all <- calc.Cstat.per.outcome(pred.seqlog.all)
  Cstat.OvA <- calc.Cstat.per.outcome(pred.OvA)
  Cstat.OvO.PC <- calc.Cstat.per.outcome(pred.OvO.PC)
  
  print(paste("PDI", Sys.time()))
  PDI.multinomial <- pdiest(pred.multinomial)
  PDI.seqlog <- pdiest(pred.seqlog)
  PDI.seqlog.all <- pdiest(pred.seqlog.all)
  PDI.OvA <- pdiest(pred.OvA)
  PDI.OvO.PC <- pdiest(pred.OvO.PC)
  
  print(paste("ECI", Sys.time()))
  ECI.multinomial <- calibrate.model.flexible.polytomous(pred.multinomial)$ECI
  ECI.seqlog <- calibrate.model.flexible.polytomous(pred.seqlog)$ECI
  ECI.seqlog.all <- calibrate.model.flexible.polytomous(pred.seqlog.all)$ECI
  ECI.OvA <- calibrate.model.flexible.polytomous(pred.OvA)$ECI
  ECI.OvO.PC <- calibrate.model.flexible.polytomous(pred.OvO.PC)$ECI
  
  return(list("wc.ms.multinomial" = wc.ms.multinomial, "wc.ms.seqlog" = wc.ms.seqlog, "wc.ms.OvA" = wc.ms.OvA, "wc.ms.OvO.PC" = wc.ms.OvO.PC,
              "wc.po.multinomial" = wc.po.multinomial, "wc.po.seqlog" = wc.po.seqlog, "wc.po.seqlog.all" = wc.po.seqlog.all, "wc.po.OvA" = wc.po.OvA, "wc.po.OvO.PC" = wc.po.OvO.PC,
              "OvE.multinomial" = OvE.multinomial, "OvE.seqlog" = OvE.seqlog, "OvE.seqlog.all" = OvE.seqlog.all, "OvE.OvA" = OvE.OvA, "OvE.OvO.PC" = OvE.OvO.PC,
              "Cstat.multinomial" = Cstat.multinomial, "Cstat.seqlog" = Cstat.seqlog, "Cstat.seqlog.all" = Cstat.seqlog.all, "Cstat.OvA" = Cstat.OvA, "Cstat.OvO.PC" = Cstat.OvO.PC,
              "PDI.multinomial" = PDI.multinomial, "PDI.seqlog" = PDI.seqlog, "PDI.seqlog.all" = PDI.seqlog.all, "PDI.OvA" = PDI.OvA, "PDI.OvO.PC" = PDI.OvO.PC,
              "ECI.multinomial" = ECI.multinomial, "ECI.seqlog" = ECI.seqlog, "ECI.seqlog.all" = ECI.seqlog.all, "ECI.OvA" = ECI.OvA, "ECI.OvO.PC" = ECI.OvO.PC
  ))
}


#############################################################
### 3.8B) Function to run simulation for small sample size for K = 5 ###
### We don't calculate ECI, its too computationally intensive
#############################################################

run.sim.small.sample.K5 <- function(dat.devel, dat.valid){
  
  print(paste("multinomial", Sys.time()))
  pred.multinomial <- fit.model.multinomial(dat.devel, dat.valid)
  print(paste("seqlog", Sys.time()))
  pred.seqlog <- fit.model.seqlog(dat.devel, dat.valid)
  print(paste("seqlog.all", Sys.time()))
  pred.seqlog.all <- fit.model.seqlog.all(dat.devel, dat.valid)
  print(paste("OvA", Sys.time()))
  pred.OvA <- fit.model.OvA(dat.devel, dat.valid)
  print(paste("OvO.PC", Sys.time()))
  pred.OvO.PC <- fit.model.OvO.PC(dat.devel, dat.valid)
  
  print(paste("wc.ms", Sys.time()))
  wc.ms.multinomial <- calibrate.mod.specific.multinomial(pred.multinomial)
  wc.ms.seqlog <- calibrate.mod.specific.seqlog(pred.seqlog)
  wc.ms.OvA <- calibrate.mod.specific.OvA(pred.OvA)
  wc.ms.OvO.PC <- calibrate.mod.specific.OvO.PC(pred.OvO.PC)
  
  print(paste("wc.po", Sys.time()))
  wc.po.multinomial <- calibrate.weak.per.outcome(pred.multinomial)
  wc.po.seqlog <- calibrate.weak.per.outcome(pred.seqlog)
  wc.po.seqlog.all <- calibrate.weak.per.outcome(pred.seqlog.all)
  wc.po.OvA <- calibrate.weak.per.outcome(pred.OvA)
  wc.po.OvO.PC <- calibrate.weak.per.outcome(pred.OvO.PC)
  
  print(paste("OvE", Sys.time()))
  OvE.multinomial <- calibrate.model.expected.observed(pred.multinomial)
  OvE.seqlog <- calibrate.model.expected.observed(pred.seqlog)
  OvE.seqlog.all <- calibrate.model.expected.observed(pred.seqlog.all)
  OvE.OvA <- calibrate.model.expected.observed(pred.OvA)
  OvE.OvO.PC <- calibrate.model.expected.observed(pred.OvO.PC)
  
  print(paste("Cstat", Sys.time()))
  Cstat.multinomial <- calc.Cstat.per.outcome(pred.multinomial)
  Cstat.seqlog <- calc.Cstat.per.outcome(pred.seqlog)
  Cstat.seqlog.all <- calc.Cstat.per.outcome(pred.seqlog.all)
  Cstat.OvA <- calc.Cstat.per.outcome(pred.OvA)
  Cstat.OvO.PC <- calc.Cstat.per.outcome(pred.OvO.PC)
  
  print(paste("PDI", Sys.time()))
  PDI.multinomial <- pdiest(pred.multinomial)
  PDI.seqlog <- pdiest(pred.seqlog)
  PDI.seqlog.all <- pdiest(pred.seqlog.all)
  PDI.OvA <- pdiest(pred.OvA)
  PDI.OvO.PC <- pdiest(pred.OvO.PC)
  
#   print(paste("ECI", Sys.time()))
#   ECI.multinomial <- calibrate.model.flexible.polytomous(pred.multinomial)$ECI
#   ECI.seqlog <- calibrate.model.flexible.polytomous(pred.seqlog)$ECI
#   ECI.seqlog.all <- calibrate.model.flexible.polytomous(pred.seqlog.all)$ECI
#   ECI.OvA <- calibrate.model.flexible.polytomous(pred.OvA)$ECI
#   ECI.OvO.PC <- calibrate.model.flexible.polytomous(pred.OvO.PC)$ECI
  
  return(list("wc.ms.multinomial" = wc.ms.multinomial, "wc.ms.seqlog" = wc.ms.seqlog, "wc.ms.OvA" = wc.ms.OvA, "wc.ms.OvO.PC" = wc.ms.OvO.PC,
              "wc.po.multinomial" = wc.po.multinomial, "wc.po.seqlog" = wc.po.seqlog, "wc.po.seqlog.all" = wc.po.seqlog.all, "wc.po.OvA" = wc.po.OvA, "wc.po.OvO.PC" = wc.po.OvO.PC,
              "OvE.multinomial" = OvE.multinomial, "OvE.seqlog" = OvE.seqlog, "OvE.seqlog.all" = OvE.seqlog.all, "OvE.OvA" = OvE.OvA, "OvE.OvO.PC" = OvE.OvO.PC,
              "Cstat.multinomial" = Cstat.multinomial, "Cstat.seqlog" = Cstat.seqlog, "Cstat.seqlog.all" = Cstat.seqlog.all, "Cstat.OvA" = Cstat.OvA, "Cstat.OvO.PC" = Cstat.OvO.PC,
              "PDI.multinomial" = PDI.multinomial, "PDI.seqlog" = PDI.seqlog, "PDI.seqlog.all" = PDI.seqlog.all, "PDI.OvA" = PDI.OvA, "PDI.OvO.PC" = PDI.OvO.PC
              #"ECI.multinomial" = ECI.multinomial, "ECI.seqlog" = ECI.seqlog, "ECI.seqlog.all" = ECI.seqlog.all, "ECI.OvA" = ECI.OvA, "ECI.OvO.PC" = ECI.OvO.PC
  ))
}




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
    dat.probs.individual.models[[i]] <- 1 - LP[[i]]/(1+LP[[i]])
  }
  
  ## Now calculate the probability of having each event
  # First outcome category
  dat.probs[,1] <- 1 - LP[[1]]/(1+LP[[1]])
  # Outcome categories 2 to K-1
  # It's the probability of i-1, divided by the probability of i-1 from the last seq log model, multiplied by the prob of
  # NOT having i-1 from the last seq log model, multiplied by the prob of having i from the next seqlog model
  for (i in 2:(K-1)){
    dat.probs[,i] <- (dat.probs[ , (i-1)])*(1/(1 - LP[[i-1]]/(1+LP[[i-1]])))*(LP[[(i-1)]]/(1 + LP[[(i-1)]]))*(1 - LP[[i]]/(1+LP[[i]]))
  }
  # Outcome category K
  dat.probs[,K] <- (dat.probs[ , (K-1)])*(1/(1 - LP[[K-1]]/(1+LP[[K-1]])))*(LP[[(K-1)]]/(1 + LP[[(K-1)]]))
  
  return(list(dat.probs, LP))
}

