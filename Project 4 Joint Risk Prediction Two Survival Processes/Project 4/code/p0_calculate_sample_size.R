#install.packages("pmsampsize")
library(pmsampsize)


### Estimate R2_ADJ assuming as R2_NAGEL of 0.15

## Let n be 1000
n <- 1000
## follow up time t is 10 years for each individual (admittedly slightly less due to censoring)
t <- 10*n
## Number of events is mean 10 year joint risk (0.0125), mulitplied by n 0.0125
E <- n*0.0125

## Calculate lnull and maxR2
lnull <- E*log(E/t) + E #eqn 13 from Riley et al
maxR2 <- 1 - exp(2*lnull/n) 
maxR2

## Note this is identical if we increase n, which makes sense
n <- 5000
t <- 10*n
E <- n*0.0125
lnull <- E*log(E/t) + E
maxR2 <- 1 - exp(2*lnull/n)
maxR2


### Assume R2_ADJ equates to R2_NAGEL of 0.15
R2_ADJ <- 0.15*maxR2

### Calculate sample size manually
s_vh <- 0.9
p <- 2
n.req <- p/((s_vh - 1)*log(1 - R2_ADJ/s_vh))
n.req


### Calculate using sample size package
pmsampsize(type = "s",
           rsquared = R2_ADJ,
           parameters = 2,
           shrinkage = 0.9,
           rate = 0.00125,
           timepoint = 10,
           meanfup = 9)


### They give same answer of 896





pmsampsize()


  
  
  
  
  
  
  
2/()