### This program will load the output from the frailty model, which is really large, and just save the predicted risks, so that these can be loaded
### quickly for other programs which calculate performance metric. The full output will remain saved for checking model convergence, etc.

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Source packages
source("Project 4/code/sim_function_load_packages.R")
### Load functions
source("Project 4/code/sim_function_clinical_example.R")

### Run through model types and sample sizes
# for (model.type in c("A", "B", "C")){
#   for (data.size in c(200000, 400000)){
for (model.type in c("C")){
  for (data.size in c(200000)){  
    
    print(model.type)
    print(data.size)
    
    ###
    ### Load normal model
    ###
    load(paste("Project 4/data/clinical_example_fit_frailty_normal_model", model.type, "_n", data.size, ".RData", sep = ""))
    
    ### Remove all the stan model information
    predrisk.frailty.normal[[1]] <- NA
    
    ### Remove all other excess objects
    rm(list=setdiff(ls(), list("predrisk.frailty.normal", "model.type", "data.size")))
    
    ### Save image
    save.image(paste("Project 4/data/clinical_example_fit_frailty_normal_model", model.type, "_n", data.size, "_risksonly.RData", sep = ""))
    print("normal dome")
    
    ###
    ### Load gamma model
    ###
    load(paste("Project 4/data/clinical_example_fit_frailty_gamma_model", model.type, "_n", data.size, ".RData", sep = ""))
    
    ### Remove all the stan model information
    predrisk.frailty.gamma[[1]] <- NA
    
    ### Remove all other excess objects
    rm(list=setdiff(ls(), list("predrisk.frailty.gamma", "model.type", "data.size")))
    
    ### Save image
    save.image(paste("Project 4/data/clinical_example_fit_frailty_gamma_model", model.type, "_n", data.size, "_risksonly.RData", sep = ""))
    print("gamma done")
  }
}

