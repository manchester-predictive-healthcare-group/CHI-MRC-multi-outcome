### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Source packages
source("Project 4/code/sim_function_load_packages.R")
### Load functions
source("Project 4/code/sim_function_clinical_example.R")
### Load data
variables.vec <- c("Age", "gender", "Smoking", "SBP", "Cholhdl_ratio", "IMD", "BMI", "Ethnicity6")
model.type <- "C"
source("Project 4/code/p4_clinical_example_load_data_200000.R")

##################
##################
### FIT MODELS ###
##################
##################

print("fit frailty gamma")
Sys.time()
predrisk.frailty.gamma <- calc.predrisk.frailty.parallel.5.init(data.devel = data.devel, 
                                                                 data.valid = data.valid, 
                                                                 t.eval = t.eval, 
                                                                 variables.vec = variables.vec, 
                                                                 frail.dist = "gamma", 
                                                                 baseline.dist = "weibull", 
                                                                 n.iter = 4000)
print("frailty gamma fitted")
Sys.time()
save.image(paste("Project 4/data/clinical_example_fit_frailty_gamma_model", model.type, "_n", data.size, ".RData", sep = ""))

