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
    
    ###
    ### Load normal model
    ###
    
    load(paste("Project 4/data/clinical_example_fit_frailty_normal_model", model.type, "_n", data.size, ".RData", sep = ""))
    
    traceplot.betas <- traceplot(predrisk.frailty.normal[["rstan.model"]], pars = c("betas_A", "betas_B"), inc_warmup = TRUE)
    traceplot.other <- traceplot(predrisk.frailty.normal[["rstan.model"]], pars = c("shape_A", "intercept_A", "shape_B", "intercept_B", "frail_param"), inc_warmup = TRUE)
    
    ggsave(paste("Project 4/figures/clin.example.traceplot.f.normal.model", model.type, ".n", data.size, ".betas.png", sep = ""), traceplot.betas, type = "cairo-png")
    ggsave(paste("Project 4/figures/clin.example.traceplot.f.normal.model", model.type, ".n", data.size, ".other.png", sep = ""), traceplot.other, type = "cairo-png")
    
    
    ###
    ### Load gamma model
    ###
    
    load(paste("Project 4/data/clinical_example_fit_frailty_gamma_model", model.type, "_n", data.size, ".RData", sep = ""))
    
    traceplot.betas <- traceplot(predrisk.frailty.gamma[["rstan.model"]], pars = c("betas_A", "betas_B"), inc_warmup = TRUE)
    traceplot.other <- traceplot(predrisk.frailty.gamma[["rstan.model"]], pars = c("shape_A", "intercept_A", "shape_B", "intercept_B", "frail_param"), inc_warmup = TRUE)
    
    ggsave(paste("Project 4/figures/clin.example.traceplot.f.gamma.model", model.type, ".n", data.size, ".betas.png", sep = ""), traceplot.betas, type = "cairo-png")
    ggsave(paste("Project 4/figures/clin.example.traceplot.f.gamma.model", model.type, ".n", data.size, ".other.png", sep = ""), traceplot.other, type = "cairo-png")
    
  }
}
