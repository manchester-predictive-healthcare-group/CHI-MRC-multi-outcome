### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")
source("code/z_functions_ce.R")

### Choose t.eval
t.eval <- ceiling(10*365.25)

### Read in seed/iter/n.cohort
args <- commandArgs(trailingOnly = T)
n.cohort <- as.numeric(args[1])
print(paste("n.cohort = ", n.cohort, sep = ""))

### There are a different number of .RData files depending on the size of n.cohort
if (n.cohort == 5000){
  n.set <- 1000
} else if (n.cohort == 3000){
  n.set <- 2000
} else if (n.cohort == 100000){
  n.set <- 100000
}

### Load and combine predicted risks data into a list
probtrans.all.list <- vector("list", n.set)
false <- 0

for (set in 1:n.set){
  ## Load data
  print(set)
  if (file.exists(paste("data/ce_calc_probtrans_N", n.cohort, "_set", set, ".RData", sep = "")) == TRUE){
    load(paste("data/ce_calc_probtrans_N", n.cohort, "_set", set, ".RData", sep = ""))
    probtrans.all.list[[set]] <- probtrans.set
  } else if (file.exists(paste("data/ce_calc_probtrans_N", n.cohort, "_set", set, ".RData", sep = "")) == FALSE){
    false <- c(false, set)
  }
}


length(false)
print("FALSE")
false
write.table(false[-1], file = "code/ce_fail.txt", row.names = FALSE, col.names = FALSE)

### Error with scientific notation, load set = 100000 manually
load(paste("data/ce_calc_probtrans_N", 100000, "_set", 100000, ".RData", sep = ""))
probtrans.all.list[[100000]] <- probtrans.set

### Each list contains 10 datasets (one for each time point)
### Combine the 1000 datasets for each time point, and arrange them by person_id
probtrans.all.list2 <- vector("list", 10)
for (i in 1:10){
  print(paste(i, Sys.time()))
  probtrans.all.list2[[i]] <- do.call("rbind", lapply(probtrans.all.list, function(x) x[[i]]))
  probtrans.all.list2[[i]] <- arrange(probtrans.all.list2[[i]], person_id)
  colnames(probtrans.all.list2[[i]])[2:7] <- paste("pstate", 1:6, sep = "")
}
probtrans.all.list <- probtrans.all.list2
rm(probtrans.all.list2)


### Get a list of patids for which generating the predicted risks was succesful
patids.success <- probtrans.all.list[[1]]$person_id
length(patids.success)

### Load the data on which the models were fitted
load(paste("data/ce_fit_csh_msm_N", n.cohort, ".RData", sep = ""))

### Create a "raw" (not mstate format) dataset
data.raw <- subset(complete.data, person_id %in% patids.success)

### Arrange data.raw by person_id
data.raw <- arrange(data.raw, person_id)
### Arrange predicted risks by patid
probtrans.all.list[[10]] <- arrange(probtrans.all.list[[10]], person_id)

### Combine with predicted risks
str(data.raw)
str(probtrans.all.list[[10]])
data.raw <- cbind(data.raw, probtrans.all.list[[10]][2:13])

### Reduce the mstate datset to jut those individuals in this analysis (currently not all ran succesfully)
data.mstate <- subset(complete.data.prep.valid, person_id %in% patids.success)

### Reduce data.raw to just have variables of interst
colnames(data.raw)
data.raw <- data.raw[, c(1:13, 17, 28, 75:103)]
str(data.raw)

### Save image
rm(list = setdiff(ls(), list("data.raw", "data.mstate", "t.eval", "patids.success", "probtrans.all.list", "tmat", "n.cohort", "covs")))
save.image(paste("data/ce_combine_probtrans_N", n.cohort, ".RData", sep = ""))
print("IMAGE SAVED")



