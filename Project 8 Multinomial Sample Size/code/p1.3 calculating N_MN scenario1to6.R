### This program will calculate the required sample sizes N_MN

### These are the sample sizes obtained by applying Riley et al's criteria directly in the multinomial framework (N_MN)

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
data.valid.s6 <- create.dataset(500000, -2.9, .5, -0.25, -0.125, 0.25, 0.375,
                                -2.9, 0.375, -0.5, -0.25, -0.375, 0.125)


### Now fit a model in each of these and get the loglikelihood
model.s1 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s1)
model.s2 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s2)
model.s3 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s3)
model.s4 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s4)
model.s5 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s5)
model.s6 <- vgam(y.cat ~ x1 + x2 + x3 + x4 + x5, family = multinomial(ref = "cat1"), data = data.valid.s6)

## Fit null models
model.null.s1 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s1)
model.null.s2 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s2)
model.null.s3 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s3)
model.null.s4 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s4)
model.null.s5 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s5)
model.null.s6 <- vgam(y.cat ~ 1, family = multinomial(ref = "cat1"), data = data.valid.s6)


### Extract  loglikelihoods
lr.1 <- -2*(model.null.s1@criterion$loglikelihood - model.s1@criterion$loglikelihood)
lr.2 <- -2*(model.null.s2@criterion$loglikelihood - model.s2@criterion$loglikelihood)
lr.3 <- -2*(model.null.s3@criterion$loglikelihood - model.s3@criterion$loglikelihood)
lr.4 <- -2*(model.null.s4@criterion$loglikelihood - model.s4@criterion$loglikelihood)
lr.5 <- -2*(model.null.s5@criterion$loglikelihood - model.s5@criterion$loglikelihood)
lr.6 <- -2*(model.null.s6@criterion$loglikelihood - model.s6@criterion$loglikelihood)


### Calculaet R2_APP
R2_APP.1 <- 1 - exp(-lr.1/nrow(data.valid.s1))
R2_APP.2 <- 1 - exp(-lr.2/nrow(data.valid.s2))
R2_APP.3 <- 1 - exp(-lr.3/nrow(data.valid.s3))
R2_APP.4 <- 1 - exp(-lr.4/nrow(data.valid.s4))
R2_APP.5 <- 1 - exp(-lr.5/nrow(data.valid.s5))
R2_APP.6 <- 1 - exp(-lr.6/nrow(data.valid.s6))

### Calculate S_VH
S_VH.1 <- 1 + (10/(nrow(data.valid.s1)*log(1-R2_APP.1)))
S_VH.2 <- 1 + (10/(nrow(data.valid.s2)*log(1-R2_APP.2)))
S_VH.3 <- 1 + (10/(nrow(data.valid.s3)*log(1-R2_APP.3)))
S_VH.4 <- 1 + (10/(nrow(data.valid.s4)*log(1-R2_APP.4)))
S_VH.5 <- 1 + (10/(nrow(data.valid.s5)*log(1-R2_APP.5)))
S_VH.6 <- 1 + (10/(nrow(data.valid.s6)*log(1-R2_APP.6)))

### Calculate R2_ADJ
R2_ADJ.1 <- S_VH.1*R2_APP.1
R2_ADJ.2 <- S_VH.2*R2_APP.2
R2_ADJ.3 <- S_VH.3*R2_APP.3
R2_ADJ.4 <- S_VH.4*R2_APP.4
R2_ADJ.5 <- S_VH.5*R2_APP.5
R2_ADJ.6 <- S_VH.6*R2_APP.6

### Calculate required n
n.req.MN.1 <- 10/((0.9-1)*log(1-(R2_ADJ.1/0.9)))
n.req.MN.2 <- 10/((0.9-1)*log(1-(R2_ADJ.2/0.9)))
n.req.MN.3 <- 10/((0.9-1)*log(1-(R2_ADJ.3/0.9)))
n.req.MN.4 <- 10/((0.9-1)*log(1-(R2_ADJ.4/0.9)))
n.req.MN.5 <- 10/((0.9-1)*log(1-(R2_ADJ.5/0.9)))
n.req.MN.6 <- 10/((0.9-1)*log(1-(R2_ADJ.6/0.9)))

### Print n.req
n.req.MN.1
n.req.MN.2
n.req.MN.3
n.req.MN.4
n.req.MN.5
n.req.MN.6

rm(list=setdiff(ls(),list("n.req.MN.1", "n.req.MN.2", "n.req.MN.3", 
                          "n.req.MN.4", "n.req.MN.5", "n.req.MN.6")))

save.image("R_out/calculating N_MN scen 1to6.RData")


