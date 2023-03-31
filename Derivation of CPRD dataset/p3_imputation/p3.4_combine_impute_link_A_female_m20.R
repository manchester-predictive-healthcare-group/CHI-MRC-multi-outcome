### Combine the imputed datasets into one object

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Load packages
library(mice)

### Load image
load("data_Aurum_65plus/data_intermediate/impute_link_A_female_m20.RData")

### Combine into one object
mice.comb <- ibind(mice.out[[1]], mice.out[[2]])
mice.comb <- ibind(mice.comb, mice.out[[3]])
mice.comb <- ibind(mice.comb, mice.out[[4]])
mice.comb <- ibind(mice.comb, mice.out[[5]])
mice.comb <- ibind(mice.comb, mice.out[[6]])
mice.comb <- ibind(mice.comb, mice.out[[7]])
mice.comb <- ibind(mice.comb, mice.out[[8]])
mice.comb <- ibind(mice.comb, mice.out[[9]])
mice.comb <- ibind(mice.comb, mice.out[[10]])

### Remove everything except what we need
rm(list=setdiff(ls(), list("mice.comb", "data.for.imp", "gender.var")))

### Save image
save.image("data_Aurum_65plus/data_intermediate/imputed_datesets_link_A_female_m20.RData")
print("image saved")

