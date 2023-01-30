### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Source packages
source("code/sim_function_load_packages.R")

### Source input parameters for calibration plots
source("code/sim_function_calibration_input_parameters_s2.R")

### Define size of development and validation cohorts
n.devel.fix <- 5000
n.valid.fix <- 1000

### Load results functions
load("data/sim_function_results.RData")


### Combine the results for each dgm
dgm.nocorr.out.comb <- create.combined.output.objects.nocorr()
print("dgm nocorr fail")
dgm.nocorr.out.comb[[4]]

save.image(paste("data/sim_results_", scen, "_n", n.devel.fix, "v", n.valid.fix, ".RData", sep = ""))

###########################
### Analyse the results ###
###########################

res.DGM.nocorr <- analyse.results(dgm.nocorr.out.comb)

### Remove combined data, and just save data we will use to build the ggplots
rm(dgm.nocorr.out.comb)
save.image(paste("data/sim_results_", scen, "_n", n.devel.fix, "v", n.valid.fix, ".RData", sep = ""))




