###
### This programs calculates transition probabilities for individuals in the validation cohort, according to the model developed
### on a sample size of n.cohort
###

### Note that for each model developed (we developed multiple models on different sample sizes), this has to be run for a 
### range of values of 'set', which should be parallelised
### i.e. the program calculates transition probabilities for groups of set.size individuals from the validation cohort.
### 'set' should take the values 1:100000/set.size

### In this simulation, we used set.size of 100 for N = 5000, and a set.size of 1 for N = 100,000

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")

### Read in n.cohort
args <- commandArgs(trailingOnly = T)
n.cohort <- as.numeric(args[1])
set <- as.numeric(args[2])
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("set = ", set, sep = ""))

### Load in raw dataset which contain censoring dates (will be required, as currently people can have events on date of censoring,
### which causes an issue)
load(paste("data/ce_fit_csh_msm_N", n.cohort, ".RData", sep = ""))

### Pick number of individuals we want to calculate risks for in one go
set.size <- 50

### Create datasets for each time point we may want to validate (yearly intervals)
probtrans.set <- vector("list", 10)

### Turn each one into an empty data frame with appropriate column names
for (i in 1:10){
  ### Create data frame with correct number of columns (number of states x2)
  probtrans.set[[i]] <- data.frame(matrix(rep(1,13), ncol = 13))
  ### Name the columns
  colnames(probtrans.set[[i]]) <- c("person_id", paste("pt.",1:6, sep = ""), paste("se.", 1:6, sep = ""))
  ### Make it empty
  probtrans.set[[i]] <- probtrans.set[[i]][FALSE, ]
}

### Generate transition probs for the set/batch of individuals defined by the variable "set"
for (id in ((set-1)*set.size + 1):(set*set.size)){
  
  ### Paste progress
  print(paste("id =", id, Sys.time(), sep = " "))
  
  ### Define person_id
  person_id_temp <- person.ids.valid[id]
  
  ### Create temporary dataset to make predictors with
  ## Get location of individual
  person_id_loc <- which(complete.data.prep.valid$person_id == person_id_temp)
  ## Extract dataset with 27 rows of this individual
  temp.data <- complete.data.prep.valid[rep(person_id_loc[1], max(tmat[!is.na(tmat)])), c("person_id", covs)]
  ## Add trans
  temp.data$trans <- 1:max(tmat[!is.na(tmat)])
  ## Attribute transition matrix
  attr(temp.data, "trans") <- tmat
  ## Expand the covariates
  temp.data <- expand.covs(temp.data, covs, longnames = FALSE)
  ## Add strata
  temp.data$strata <- temp.data$trans
  print(Sys.time())
  
  ### Calculate transition specific hazards
  msm.fit <- msfit(object = msm.cox.fit, trans = tmat, 
                   newdata = temp.data)
  print(Sys.time())
  
  ### Calculate probtrans
  probtrans <- probtrans(msm.fit, predt = 0)
  print(Sys.time())
  
  ### Extract probtrans and se for each yearly interval and put into appropriate datasets
  for (i in 1:10){
    
    ### Extract prob trans
    pt.se <- probtrans[[1]] %>% 
      ## Extract row corresponding to crrecot time point
      slice(min(which(probtrans[[1]]$time > (365.25*i)))) %>% 
      ## Extract all the correct columns
      select(c(paste("pstate",1:6, sep = ""), paste("se", 1:6, sep = "")))
    
    ### Create the row to add by adding patid to pt.se
    pt.se <- data.frame("person_id" = person_id_temp, pt.se)
    
    ### Bind into output dataset
    probtrans.set[[i]] <- rbind(probtrans.set[[i]], pt.se)
  }
  
  gc()
  
}
 
### Subset the validation cohort
print("SAVING IMAGE")
print(Sys.time())
rm(list=setdiff(ls(), list("probtrans.set", "set", "n.cohort")))
save.image(paste("data/ce_calc_probtrans_N", n.cohort, "_set", set, ".RData", sep = ""))


