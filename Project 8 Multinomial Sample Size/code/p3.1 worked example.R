##### This R code will estimate R2_CS for each distinct logistic regression model, following the process of Riley et al (doi: 10.1002/sim.8806)
##### The C statistic values on which the calculation will be based, are extracted from van Calster et al (doi: 10.1186/1471-2288-10-96)

set.seed(101)

### Load relevant packages
library(pROC)

###################################################################################
###################################################################################
##### STEP'S 1 and 2: Identify values for Q, p_k, p_k_r, max(R2_CS), R2_CS_adj and 
##### R2_CS_adj_k_r
###################################################################################
###################################################################################

### First define the number of events in each category
EV1 <- 2557 # benign
EV2 <- 186 # Borderline
EV3 <- 176 # Stage 1
EV4 <- 467 # Stage 2 - 4
EV5 <- 120 # Metastatic

############
### Define Q
############
Q <- 17

########################
### Define p_k and p_k_r
########################

## Define total number of events
n.total <- EV1 + EV2 + EV3 + EV4 + EV5

## Calculate p_k
p.1 <- EV1/n.total
p.2 <- EV2/n.total
p.3 <- EV3/n.total
p.4 <- EV4/n.total
p.5 <- EV5/n.total

p1
p2
p3
p4
p5

## Calculate p_k_r
p.1.2 <- (EV1 + EV2)/n.total
p.1.3 <- (EV1 + EV3)/n.total
p.1.4 <- (EV1 + EV4)/n.total
p.1.5 <- (EV1 + EV5)/n.total

p.2.3 <- (EV2 + EV3)/n.total
p.2.4 <- (EV2 + EV4)/n.total
p.2.5 <- (EV2 + EV5)/n.total

p.3.4 <- (EV3 + EV4)/n.total
p.3.5 <- (EV3 + EV5)/n.total

p.4.5 <- (EV4 + EV5)/n.total

########################
### Calculate max(R2_CS)
########################
max_R2_CS <- 1 - (p.1^p.1*p.2^p.2*p.3^p.3*p.4^p.4*p.5^p.5)^2
max_R2_CS


########################
### Calculate R2_CS_adj
########################

### Calculate an estimte of R2_CS_app, based off R2_NAGEL = 0.15
R2_CS_adj <- 0.15*max_R2_CS
R2_CS_adj


###########################
### Calculate R2_CS_adj_k_r
###########################

### Define the C-statistic values for each model
### These were calculated in a temporal validation, so will give estimates of R2_CS
### that do not need to be adjusted for optimism (i.e. R2_CS_adj)

### The pairwise C-statistic for Benign vs Borderline malignent is 0.85
### The pairwise C-statistic for Benign vs Stage 1 0.92
### The pairwise C-statistic for Benign vs Stage 2 - 4 is 0.99
### The pairwise C-statistic for Benign vs Metastatic is 0.95
C.1.2 <- 0.85
C.1.3 <- 0.92
C.1.4 <- 0.99
C.1.5 <- 0.95

### The pairwise C-statistic for borderline vs Stage 1 is 0.75
### The pairwise C-statistic for borderline vs Stage 2 - 4 is 0.95
### The pairwise C-statistic for borderline vs Metastatic is 0.87
C.2.3 <- 0.75
C.2.4 <- 0.95
C.2.5 <- 0.87

### The pairwise C-statistic for stage 1 vs Stage 2 - 4 is 0.87
### The pairwise C-statistic for stage 1 vs Metastatic is 0.71
C.3.4 <- 0.87
C.3.5 <- 0.71

### The pairwise C-statistic for Stage 2 - 4 vs Metastatic is 0.82
C.4.5 <- 0.82

### Calculate pairwise outcome proportions (phi), of category k relative to category i
phi.1.2 <- EV2/(EV1 + EV2)
phi.1.3 <- EV3/(EV1 + EV3)
phi.1.4 <- EV4/(EV1 + EV4)
phi.1.5 <- EV5/(EV1 + EV5)

phi.2.3 <- EV3/(EV2 + EV3)
phi.2.4 <- EV4/(EV2 + EV4)
phi.2.5 <- EV5/(EV2 + EV5)

phi.3.4 <- EV4/(EV3 + EV4)
phi.3.5 <- EV5/(EV3 + EV5)

phi.4.5 <- EV5/(EV4 + EV5)

phi.1.2
phi.1.3
phi.1.4
phi.1.5

phi.2.3
phi.2.4
phi.2.5

phi.3.4
phi.3.5

phi.4.5


### Create a function to simulate R2 from a C-statistic and pairwise outcome proportion
simulate.R2 <- function(N, prop.in, C.in){

  ## Create an empty dataset
  output.dat <- data.frame(matrix(ncol = 2, nrow = N))
  colnames(output.dat) <- c("Y", "LP")
  
  ## Create the outcome variable
  Y.vec <- rbinom(N, 1, prop.in)
  
  ## Create the vector of mean values for the linear predictor data generation
  Y.vec.mu <- Y.vec*sqrt(2)*qnorm(C.in, 0, 1)
  
  ## Generate the linear predictor
  LP.vec <- rnorm(N,Y.vec.mu,1)

  ## Assign these into an output dataset
  output.dat$Y <- as.integer(Y.vec)
  output.dat$LP <- LP.vec
  
  ## Fit a logistic regression to this dataset
  model.out <- glm(Y ~ LP.vec, family = binomial(link = "logit"), data = output.dat)
  model.out
  
  ## Check the AUC is correct, matches input
  #C.stat.sim <- as.numeric(roc(Y ~ LP.vec, data = output.dat)$auc)
  
  ## Fit a null model also, to calculate likelihood ratio
  model.null <- glm(Y ~ 1, family = binomial(link = "logit"), data = output.dat)
  
  ## Calculate likelihood ratio statistics
  LR <- as.numeric(-2*(logLik(model.null) - logLik(model.out)))

  ## Calculate R2_CS_APP
  R2_CS_APP <- 1 - exp(-LR/N)
  R2_CS_APP

  ## Output object
  return(R2_CS_APP)
}

### Calculate R2_CS_adj.k.i according to the simulation approach of Riley.
R2_CS_adj.1.2 <- simulate.R2(1000000, phi.1.2, C.1.2)
R2_CS_adj.1.3 <- simulate.R2(1000000, phi.1.3, C.1.3)
R2_CS_adj.1.4 <- simulate.R2(1000000, phi.1.4, C.1.4)
R2_CS_adj.1.5 <- simulate.R2(1000000, phi.1.5, C.1.5)

R2_CS_adj.2.3 <- simulate.R2(1000000, phi.2.3, C.2.3)
R2_CS_adj.2.4 <- simulate.R2(1000000, phi.2.4, C.2.4)
R2_CS_adj.2.5 <- simulate.R2(1000000, phi.2.5, C.2.5)

R2_CS_adj.3.4 <- simulate.R2(1000000, phi.3.4, C.3.4)
R2_CS_adj.3.5 <- simulate.R2(1000000, phi.3.5, C.3.5)

R2_CS_adj.4.5 <- simulate.R2(1000000, phi.4.5, C.4.5)


R2_CS_adj.1.2
R2_CS_adj.1.3
R2_CS_adj.1.4
R2_CS_adj.1.5

R2_CS_adj.2.3
R2_CS_adj.2.4
R2_CS_adj.2.5

R2_CS_adj.3.4
R2_CS_adj.3.5

R2_CS_adj.4.5



###########################
###########################
##### STEP 3: Criterion (i)
###########################
###########################

## Let S be the level of shrinkage we are targeting
S <- 0.9

## Calculate m_k_r
m.1.2 <- Q/((S - 1)*log(1 - R2_CS_adj.1.2/S))
m.1.3 <- Q/((S - 1)*log(1 - R2_CS_adj.1.3/S))
m.1.4 <- Q/((S - 1)*log(1 - R2_CS_adj.1.4/S))
m.1.5 <- Q/((S - 1)*log(1 - R2_CS_adj.1.5/S))

m.2.3 <- Q/((S - 1)*log(1 - R2_CS_adj.2.3/S))
m.2.4 <- Q/((S - 1)*log(1 - R2_CS_adj.2.4/S))
m.2.5 <- Q/((S - 1)*log(1 - R2_CS_adj.2.5/S))

m.3.4 <- Q/((S - 1)*log(1 - R2_CS_adj.3.4/S))
m.3.5 <- Q/((S - 1)*log(1 - R2_CS_adj.3.5/S))

m.4.5 <- Q/((S - 1)*log(1 - R2_CS_adj.4.5/S))


### Calculate n_k_r for criterion (i) for each submodel
N_C1.1.2 <- m.1.2/p.1.2
N_C1.1.3 <- m.1.3/p.1.3
N_C1.1.4 <- m.1.4/p.1.4
N_C1.1.5 <- m.1.5/p.1.5

N_C1.2.3 <- m.2.3/p.2.3
N_C1.2.4 <- m.2.4/p.2.4
N_C1.2.5 <- m.2.5/p.2.5

N_C1.3.4 <- m.3.4/p.3.4
N_C1.3.5 <- m.3.5/p.3.5

N_C1.4.5 <- m.4.5/p.4.5

N_C1.1.2
N_C1.1.3
N_C1.1.4
N_C1.1.5

N_C1.2.3
N_C1.2.4
N_C1.2.5

N_C1.3.4
N_C1.3.5

N_C1.4.5


### Take the ceiling of the maximum of these as the sample size for criteiron (i)
N_C1 <- ceiling(max(N_C1.1.2, N_C1.1.3, N_C1.1.4, N_C1.1.5, N_C1.2.3, N_C1.2.4, N_C1.2.5, 
            N_C1.3.4, N_C1.3.5, N_C1.4.5))
N_C1


### Now calculate number of each event we expect to see in a datast of this size
N_C1*p.1
N_C1*p.2
N_C1*p.3
N_C1*p.4
N_C1*p.5


############################
############################
##### STEP 4: Criterion (ii)
############################
############################

N_C2 <- 4*Q/((R2_CS_adj/(R2_CS_adj + 0.05*max_R2_CS) - 1)*log(1 - R2_CS_adj - 0.05*max_R2_CS))
N_C2 <- ceiling(N_C2)
N_C2

### Now calculate number of each event we expect to see in a datast of this size
N_C2*p.1
N_C2*p.2
N_C2*p.3
N_C2*p.4
N_C2*p.5


#############################
#############################
##### STEP 5: Criterion (iii)
#############################
#############################

N_C3.1 <- qchisq(1-0.05/5, 1)*p.1*(1-p.1)/0.05^2 
N_C3.2 <- qchisq(1-0.05/5, 1)*p.2*(1-p.2)/0.05^2 
N_C3.3 <- qchisq(1-0.05/5, 1)*p.3*(1-p.3)/0.05^2 
N_C3.4 <- qchisq(1-0.05/5, 1)*p.4*(1-p.4)/0.05^2 
N_C3.5 <- qchisq(1-0.05/5, 1)*p.5*(1-p.5)/0.05^2

N_C3.1
N_C3.2
N_C3.3
N_C3.4
N_C3.5

N_C3 <- ceiling(max(N_C3.1, N_C3.2, N_C3.3, N_C3.4, N_C3.5))
N_C3

### Now calculate number of each event we expect to see in a datast of this size
N_C3*p.1
N_C3*p.2
N_C3*p.3
N_C3*p.4
N_C3*p.5


#####################################################################
#####################################################################
##### STEP 6: Take the maximum sample size across all three criteria
#####################################################################
#####################################################################

N_C1
N_C2
N_C3

N <- max(N_C1, N_C2, N_C3)
N
