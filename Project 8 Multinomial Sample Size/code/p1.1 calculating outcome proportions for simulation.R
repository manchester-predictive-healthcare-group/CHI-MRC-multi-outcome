### Load required packages
library("VGAM")

### This code will return the outcome proportions when using the coefficients defining the different scenarios
### in the main siulation

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

## Scenario 1
print("start scenario1")
test <- create.dataset(500000, 0, .5, -0.25, -0.125, 0.25, 0.375,
                       0, 0.375, -0.5, -0.25, -0.375, 0.125)

scen1.outcome.dist <- prop.table(table(test$y.cat))
scen1.outcome.dist

rm(test)
## Results in cat1 = 0.33, cat2 = 0.33, cat3 = 0.34

## Scenario 2
print("start scenario2")
test <- create.dataset(500000, -0, .5, -0.25, -0.125, 0.25, 0.375,
                       -0.75, 0.375, -0.5, -0.25, -0.375, 0.125)

scen2.outcome.dist <- prop.table(table(test$y.cat))
scen2.outcome.dist

rm(test)
## Results in cat1 = 0.4, cat2 = 0.4, cat3 = 0.2

## Scenario 3
print("start scenario3")
test <- create.dataset(500000, -0.35, .5, -0.25, -0.125, 0.25, 0.375,
                       -0.85, 0.375, -0.5, -0.25, -0.375, 0.125)

scen3.outcome.dist <- prop.table(table(test$y.cat))
scen3.outcome.dist


rm(test)
## Results in cat1 = 0.46, cat2 = 0.33, cat3 = 0.21

## Scenario 4
print("start scenario4")
test <- create.dataset(500000, -0.4, .5, -0.25, -0.125, 0.25, 0.375,
                       -1.7, 0.375, -0.5, -0.25, -0.375, 0.125)

scen4.outcome.dist <- prop.table(table(test$y.cat))
scen4.outcome.dist

rm(test)
## Results in cat1 = 0.53, cat2 = 0.36, cat3 = 0.11


## Scenario 5
print("start scenario5")
test <- create.dataset(500000, -0.53, .5, -0.25, -0.125, 0.25, 0.375,
                       -2.4, 0.375, -0.5, -0.25, -0.375, 0.125)

scen5.outcome.dist <- prop.table(table(test$y.cat))
scen5.outcome.dist

rm(test)
## Results in cat1 = 0.58, cat2 = 0.36, cat3 = 0.06


## Scenario 6
print("start scneario6")
test <- create.dataset(500000, -2.9, .5, -0.25, -0.125, 0.25, 0.375,
                       -2.9, 0.375, -0.5, -0.25, -0.375, 0.125)

scen6.outcome.dist <- prop.table(table(test$y.cat))
scen6.outcome.dist

rm(test)
## Results in cat1 = 0.88, cat2 = 0.06, cat3 = 0.06

## Scenario 7
print("start scenario7")
test <- create.dataset(500000, 0, 1, -0.5, -0.25, 0.5, 0.75,
                       0, 0.75, -1, -0.5, -0.75, 0.25)

scen7.outcome.dist <- prop.table(table(test$y.cat))
scen7.outcome.dist

rm(test)
## Results in cat1 = 0.33, cat2 = 0.33, cat3 = 0.34


## Scenario 8
print("start scenario8")
test <- create.dataset(500000, 0, 1, -0.5, -0.25, 0.5, 0.75,
                       -1, 0.75, -1, -0.5, -0.75, 0.25)

scen8.outcome.dist <- prop.table(table(test$y.cat))
scen8.outcome.dist

rm(test)
## Results in cat1 = 0.4, cat2 = 0.4, cat3 = 0.19



## Scenario 9
print("start scenario9")
test <- create.dataset(500000, -0.4, 1, -0.5, -0.25, 0.5, 0.75,
                       -1, 0.75, -1, -0.5, -0.75, 0.25)

scen9.outcome.dist <- prop.table(table(test$y.cat))
scen9.outcome.dist

rm(test)
## Results in cat1 = 0.45, cat2 = 0.33, cat3 = 0.21


## Scenario 10
print("start scenario10")
test <- create.dataset(500000, -0.5, 1, -0.5, -0.25, 0.5, 0.75,
                       -2, 0.75, -1, -0.5, -0.75, 0.25)

scen10.outcome.dist <- prop.table(table(test$y.cat))
scen10.outcome.dist

rm(test)
## Results in cat1 = 0.52, cat2 = 0.36, cat3 = 0.11


## Test 11
print("start test11")
test <- create.dataset(500000, -0.65, 1, -0.5, -0.25, 0.5, 0.75,
                       -2.85, 0.75, -1, -0.5, -0.75, 0.25)

scen11.outcome.dist <- prop.table(table(test$y.cat))
scen11.outcome.dist

rm(test)
## Results in cat1 = 0.58, cat2 = 0.36, cat3 = 0.06



## Scenario 12
print("start scneario12")
test <- create.dataset(500000, -3.5, 1, -0.5, -0.25, 0.5, 0.75,
                       -3.5, 0.75, -1, -0.5, -0.75, 0.25)

scen12.outcome.dist <- prop.table(table(test$y.cat))
scen12.outcome.dist

rm(test)
## Results in cat1 = 0.88, cat2 = 0.06, cat3 = 0.06

save.image("R_out/outcome_proportions.RData")
