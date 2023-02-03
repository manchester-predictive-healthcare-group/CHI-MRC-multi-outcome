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

### Load iter number
args <- commandArgs(trailingOnly = T)
iter <- as.numeric(args[1])
print(iter)

Sys.time()
predrisk.msm <- calc.predrisk.msm(data.devel = data.devel, 
                                  data.valid = data.valid,
                                  t.eval = t.eval, 
                                  variables.vec = variables.vec,
                                  msm.iter = iter,
                                  msm.chunk = 100)
Sys.time()

msm.risk.joint.est <-predrisk.msm[["risk.joint.est"]]
msm.fit <- predrisk.msm[["msm.fit"]]

if (iter == 1){
rm(list = setdiff(ls(), list("msm.risk.joint.est", "msm.fit", "data.size", "iter", "t.eval", "variables.vec", "data.devel", "data.valid", "model.type")))
} else if (iter != 1){
rm(list = setdiff(ls(), list("msm.risk.joint.est", "data.size", "iter", "model.type")))
}

save.image(paste("Project 4/data/clinical_example_msm_fit_model", model.type, "_n", data.size, "_iter", iter, ".RData", sep = ""))
print("FINISHED")