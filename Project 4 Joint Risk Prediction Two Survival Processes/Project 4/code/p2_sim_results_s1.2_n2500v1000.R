### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Source packages
source("code/sim_function_load_packages.R")

### Source input parameters for calibration plots
source("code/sim_function_calibration_input_parameters_s1.2.R")

### Define size of development and validation cohorts
n.devel.fix <- 2500
n.valid.fix <- 1000

### Load results functions
load("data/sim_function_results.RData")


### Combine the results for each dgm
dgm.clay.out.comb <- create.combined.output.objects("clay")
print("dgm clay fail")
dgm.clay.out.comb[[4]]

dgm.gumb.out.comb <- create.combined.output.objects("gumb")
print("dgm gumb fail")
dgm.gumb.out.comb[[4]]

dgm.frank.out.comb <- create.combined.output.objects("frank")
print("dgm frank fail")
dgm.frank.out.comb[[4]]

dgm.gamma.out.comb <- create.combined.output.objects("gamma")
print("dgm gamma fail")
dgm.gamma.out.comb[[4]]

dgm.normal.out.comb <- create.combined.output.objects("normal")
print("dgm normal fail")
dgm.normal.out.comb[[4]]

dgm.msm.out.comb <- create.combined.output.objects("msm")
print("dgm msm fail")
dgm.msm.out.comb[[4]]

save.image(paste("data/sim_results_", scen, "_n", n.devel.fix, "v", n.valid.fix, ".RData", sep = ""))

###########################
### Analyse the results ###
###########################

print("ANALYSE DGM CLAY")
res.DGM.clay <- analyse.results(dgm.clay.out.comb)
print("ANALYSE DGM GUMB")
res.DGM.gumb <- analyse.results(dgm.gumb.out.comb)
print("ANALYSE DGM FRANK")
res.DGM.frank <- analyse.results(dgm.frank.out.comb)
print("ANALYSE DGM MSM")
res.DGM.msm <- analyse.results(dgm.msm.out.comb)
print("ANALYSE DGM NORMAL")
res.DGM.normal <- analyse.results(dgm.normal.out.comb)
print("ANALYSE DGM GAMMA")
res.DGM.gamma <- analyse.results(dgm.gamma.out.comb)

### Remove combined data, and just save data we will use to build the ggplots
rm(dgm.clay.out.comb, dgm.gumb.out.comb, dgm.frank.out.comb, dgm.msm.out.comb, dgm.normal.out.comb, dgm.gamma.out.comb)
save.image(paste("data/sim_results_", scen, "_n", n.devel.fix, "v", n.valid.fix, ".RData", sep = ""))
