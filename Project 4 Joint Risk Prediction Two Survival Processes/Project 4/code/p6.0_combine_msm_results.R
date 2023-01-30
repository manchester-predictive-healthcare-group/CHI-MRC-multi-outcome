### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Source packages
source("Project 4/code/sim_function_load_packages.R")
### Load functions
source("Project 4/code/sim_function_clinical_example.R")


### Cycle through model.type and data.size
# for (model.type in c("A", "B", "C")){
#   for (data.size in c(200000, 400000)){
for (model.type in c("C")){
  for (data.size in c(200000)){
    
    ### Clear results from previous scenario to avoid contamination
    rm(list = setdiff(ls(), list("model.type", "data.size")))
    
    ### Load output for iter = 1
    load(paste("Project 4/data/clinical_example_msm_fit_model", model.type, "_n", data.size, "_iter", 1, ".RData", sep = ""))
    
    ### Assign output risks to a vector
    predrisk.msm <- msm.risk.joint.est
    
    ### Now look through and get output risks for each value of iter >2 and concatenate into a vector
    for (iter in 2:1000){
      ### Load dataset
      load(paste("Project 4/data/clinical_example_msm_fit_model", model.type, "_n", data.size, "_iter", iter, ".RData", sep = ""))
      
      ### Combine output
      predrisk.msm <- c(predrisk.msm, msm.risk.joint.est)
      
      print(iter)
    }
    
    ### Remove everything we don't need
    rm(list = setdiff(ls(), list("model.type", "data.size", "predrisk.msm")))
    
    save.image(paste("Project 4/data/clinical_example_msm_fit_model", model.type, "_n", data.size, "_combined.RData", sep = ""))
    
  }
}

