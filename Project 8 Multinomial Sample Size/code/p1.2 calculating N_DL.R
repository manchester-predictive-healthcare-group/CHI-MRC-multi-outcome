### This program will calculate the required sample sizes N_DL

### This is the sample size applying Riley's criteria to distinct logistic models for each comparison that is made in the multinomial framework
### and taking the largest one.

#install.packages("VGAM")

library(VGAM)

### Set seed
set.seed(101)

### First define a function that will cerate a dataset according to various beta values
create.dataset <- function(npatt, beta02, beta12, beta22, beta32, beta42, beta52, 
                           beta03, beta13, beta23, beta33, beta43, beta53){
  x1 <- rnorm(npatt, 0, 1)
  x2 <- rnorm(npatt, 0, 1)
  x3 <- rnorm(npatt, 0, 1)
  x4 <- rnorm(npatt, 0, 1)
  x5 <- rnorm(npatt, 0, 1)
  p1 <- 1/(1 + exp(beta02 + x1*beta12 + x2*beta22 + x3*beta32 + x4*beta42 + x5*beta52) + 
             exp(beta03 + x1*beta13 + x2*beta23 + x3*beta33 + x4*beta43 + x5*beta53))
  p2 <- exp(beta02 + x1*beta12 + x2*beta22 + x3*beta32 + x4*beta42 + x5*beta52)/(1 + exp(beta02 + x1*beta12 + x2*beta22 + x3*beta32 + x4*beta42 + x5*beta52) + 
                                                                                   exp(beta03 + x1*beta13 + x2*beta23 + x3*beta33 + x4*beta43 + x5*beta53))
  p3 <- exp(beta03 + x1*beta13 + x2*beta23 + x3*beta33 + x4*beta43 + x5*beta53)/(1 + exp(beta02 + x1*beta12 + x2*beta22 + x3*beta32 + x4*beta42 + x5*beta52) + 
                                                                                   exp(beta03 + x1*beta13 + x2*beta23 + x3*beta33 + x4*beta43 + x5*beta53))
  multinom.probs <- data.frame("p1" = p1, "p2" = p2, "p3" = p3)
  
  ### Now to generate outcome data from a multinomial distribution for each patient
  multinom.outcomes <- matrix(nrow = 3, ncol = npatt)
  for (i in 1:npatt){
    multinom.outcomes[, i] <- rmultinom(n = 1, size = 1, prob = multinom.probs[i, ])
  }
  
  ### Now for each column, need to convert this into a number depending on which row is equal to 1, and add it to a vector
  y.vec <- rep(0, npatt)
  for (i in 1:npatt){
    y.vec[i] <- which(multinom.outcomes[,i] > 0)
  }
  
  
  ### Now combine the outcome data, and the predictor data into a dataset
  devel.data <- data.frame("x1" = x1, "x2" = x2, "x3" = x3, "x4" = x4, "x5" = x5, "y.num" = y.vec)
  
  ### Create a non-numeric vector for y, so it is easier to deal with
  devel.data$y.cat <- rep(0,npatt)
  devel.data$y.cat[devel.data$y.num == 1] <- "cat1"
  devel.data$y.cat[devel.data$y.num == 2] <- "cat2"
  devel.data$y.cat[devel.data$y.num == 3] <- "cat3"
  return(devel.data)
}


### Now generate a dataset for each scenario
data.valid.s1 <- create.dataset(500000, 0, .5, -0.25, -0.125, 0.25, 0.375,
                                0, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s2 <- create.dataset(500000, -0, .5, -0.25, -0.125, 0.25, 0.375,
                                -0.75, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s3 <- create.dataset(500000, -0.35, .5, -0.25, -0.125, 0.25, 0.375,
                                -0.85, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s4 <- create.dataset(500000, -0.4, .5, -0.25, -0.125, 0.25, 0.375,
                                -1.7, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s5 <- create.dataset(500000, -0.53, .5, -0.25, -0.125, 0.25, 0.375,
                                -2.4, 0.375, -0.5, -0.25, -0.375, 0.125)


data.valid.s7 <- create.dataset(500000, 0, 1, -0.5, -0.25, 0.5, 0.75,
                                0, 0.75, -1, -0.5, -0.75, 0.25)
data.valid.s8 <- create.dataset(500000, 0, 1, -0.5, -0.25, 0.5, 0.75,
                                -1, 0.75, -1, -0.5, -0.75, 0.25)
data.valid.s9 <- create.dataset(500000, -0.4, 1, -0.5, -0.25, 0.5, 0.75,
                                -1, 0.75, -1, -0.5, -0.75, 0.25)
data.valid.s10 <- create.dataset(500000, -0.5, 1, -0.5, -0.25, 0.5, 0.75,
                                -2, 0.75, -1, -0.5, -0.75, 0.25)
data.valid.s11 <- create.dataset(500000, -0.65, 1, -0.5, -0.25, 0.5, 0.75,
                                -2.85, 0.75, -1, -0.5, -0.75, 0.25)



### Write a function to calculate N_DL
calculate.N.DL <- function(data.in){
  
### Split into two datasets for the independent regressions
data.in.lp1 <- data.in[data.in$y.cat %in% c("cat1","cat2"), ]
data.in.lp2 <- data.in[data.in$y.cat %in% c("cat1","cat3"), ]

### Create a binary outcome variable in each of these datasets
data.in.lp1$Y.bin <- as.numeric(data.in.lp1$y.cat == "cat1")
data.in.lp2$Y.bin <- as.numeric(data.in.lp2$y.cat == "cat1")

### Fit a logistic regression to each dataset
logis.mod.1 <- glm(Y.bin ~ x1 + x2 + x3 + x4 + x5, family = binomial(link = "logit"), 
                   data = data.in.lp1)
logis.mod.2 <- glm(Y.bin ~ x1 + x2 + x3 + x4 + x5, family = binomial(link = "logit"), 
                   data = data.in.lp2)

### Also fit null models
logis.mod.null.1 <- glm(Y.bin ~ 1, family = binomial(link = "logit"), 
                   data = data.in.lp1)
logis.mod.null.2 <- glm(Y.bin ~ 1, family = binomial(link = "logit"), 
                   data = data.in.lp2)

### Then calculate the likelihood ratio for each model
logis.LR.1 <- -2*(logLik(logis.mod.null.1) - logLik(logis.mod.1))
logis.LR.2 <- -2*(logLik(logis.mod.null.2) - logLik(logis.mod.2))

### Calculaet R2_APP
logis.R2_APP.1 <- 1 - exp(-logis.LR.1/nrow(data.in.lp1))
logis.R2_APP.2 <- 1 - exp(-logis.LR.2/nrow(data.in.lp2))

### Calculate S_VH
logis.S_VH.1 <- 1 + (5/(nrow(data.in.lp1)*log(1-logis.R2_APP.1)))
logis.S_VH.2 <- 1 + (5/(nrow(data.in.lp2)*log(1-logis.R2_APP.2)))

### Calculate R2_ADJ
logis.R2_ADJ.1 <- logis.S_VH.1*logis.R2_APP.1
logis.R2_ADJ.2 <- logis.S_VH.2*logis.R2_APP.2

### Calculate required n for each model
logis.N.init.1 <- 5/((0.9-1)*log(1-(logis.R2_ADJ.1/0.9)))
logis.N.init.2 <- 5/((0.9-1)*log(1-(logis.R2_ADJ.2/0.9)))

### However only a proportion of total individuals in the model will be used to calculate this linear predictor, therefore need to dvidie it
### by the proportion of patients that actually have the outcome category in the ratio of the linear predictor
prop.lp1 <- nrow(data.in.lp1)/nrow(data.in)
prop.lp2 <- nrow(data.in.lp2)/nrow(data.in)

logis.N.DL.1 <- logis.N.init.1/prop.lp1
logis.N.DL.2 <- logis.N.init.2/prop.lp2

return(c(logis.N.DL.1, logis.N.DL.2))}

### Calculate N.DL for each scenario (will need to take the max of the output)
N.DL.1 <- calculate.N.DL(data.valid.s1) 
N.DL.2 <- calculate.N.DL(data.valid.s2) 
N.DL.3 <- calculate.N.DL(data.valid.s3) 
N.DL.4 <- calculate.N.DL(data.valid.s4) 
N.DL.5 <- calculate.N.DL(data.valid.s5) 

N.DL.7 <- calculate.N.DL(data.valid.s7) 
N.DL.8 <- calculate.N.DL(data.valid.s8) 
N.DL.9 <- calculate.N.DL(data.valid.s9) 
N.DL.10 <- calculate.N.DL(data.valid.s10) 
N.DL.11 <- calculate.N.DL(data.valid.s11) 


### These two scenarios come later as they were introduced as scenarios at a later
### time
data.valid.s6 <- create.dataset(500000, -2.9, .5, -0.25, -0.125, 0.25, 0.375,
                                -2.9, 0.375, -0.5, -0.25, -0.375, 0.125)
data.valid.s12 <- create.dataset(500000, -3.5, 1, -0.5, -0.25, 0.5, 0.75,
                                 -3.5, 0.75, -1, -0.5, -0.75, 0.25)

N.DL.6 <- calculate.N.DL(data.valid.s6) 
N.DL.12 <- calculate.N.DL(data.valid.s12) 


N.DL.1
N.DL.2 
N.DL.3
N.DL.4
N.DL.5
N.DL.6
N.DL.7
N.DL.8
N.DL.9 
N.DL.10 
N.DL.11 
N.DL.12

rm(list=setdiff(ls(),list("N.DL.1", "N.DL.2", "N.DL.3", 
                          "N.DL.4", "N.DL.5", "N.DL.6",
                          "N.DL.7", "N.DL.8", "N.DL.9", 
                          "N.DL.10", "N.DL.11", "N.DL.12")))

save.image("R_out/calculating N_DL.RData")
