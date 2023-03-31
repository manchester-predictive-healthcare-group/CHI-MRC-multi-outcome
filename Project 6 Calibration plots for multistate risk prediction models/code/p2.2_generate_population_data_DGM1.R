###
### This program will generate a cohort of size 1,000,000 data according to DGM1
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Choose cohort size
n.cohort <- 200000

### Calculate the scale for each p, taken from Putter et al 2006, Estimation and prediction in a multistate model for breast cancer, Figure 2
csh.survival.probs.desired <- c("12" = 0.9, "13" = 0.8, "15" = 0.99, "24" = 0.55, "25" = 0.95, "34" = 0.7, "35" = 0.15, "45" = 0.05)

### Apply function to calculate the desired scales
scales.sim <- sapply(csh.survival.probs.desired, calc.scale)
scales.sim

### Defne the max follow up to be 1 day after 10 years
max.follow <- ceiling(10*365.25) + 1

### Generate the data, but paralellise the process to improve speed (note I could move to CSF at a later date if this process is highly time consuming,
### but at the moment I think other things require the CSF more (fitting the actual multistate models))
print("START DATA GEN")
cl <- makeCluster(5)
registerDoParallel(5)
start_time_fit_parallel <- Sys.time()
data_parallel_list<-(foreach(input=1:5, .combine=list, .multicombine=TRUE, 
                             .packages=c("gems", "dplyr")) %dopar%{
                       ## Set seed
                       set.seed(input)
                       
                       ## Generate baseline data
                       x.baseline <- data.frame("x1" = rnorm(n.cohort, 0, 1), "x2" = rnorm(n.cohort, 0, 1))
                       
                       ## Generate data
                       data.out <- gen.dat.DGM1(n = n.cohort, #number of patients to simulate
                                                max.follow = max.follow, #max follow up, informative censoring 1 day after 7 years
                                                shape12 = 1, scale12 = scales.sim["12"], #shape and scale for weibull baseline hazard for transition 1 -> 2
                                                shape13 = 1, scale13 = scales.sim["13"], #shape and scale for weibull baseline hazard for transition 1 -> 3
                                                shape15 = 1, scale15 = scales.sim["15"], #shape and scale for weibull baseline hazard for transition 1 -> 5
                                                shape24 = 1, scale24 = scales.sim["24"], #shape and scale for weibull baseline hazard for transition 2 -> 4
                                                shape25 = 1, scale25 = scales.sim["25"], #shape and scale for weibull baseline hazard for transition 2 -> 5
                                                shape34 = 1, scale34 = scales.sim["34"], #shape and scale for weibull baseline hazard for transition 3 -> 4
                                                shape35 = 1, scale35 = scales.sim["35"], #shape and scale for weibull baseline hazard for transition 3 -> 5
                                                shape45 = 1, scale45 = scales.sim["45"], #shape and scale for weibull baseline hazard for transition 3 -> 5
                                                beta12.x1 = 0.5, beta12.x2 = -0.5, #covariate effects for transiion 12
                                                beta13.x1 = 0.5, beta13.x2 = -0.5, #covariate effects for transiion 13
                                                beta15.x1 = 0.5, beta15.x2 = -0.5, #covariate effects for transiion 15
                                                beta24.x1 = 0.5, beta24.x2 = -0.5, #covariate effects for transiion 24
                                                beta25.x1 = 0.5, beta25.x2 = -0.5, #covariate effects for transiion 25
                                                beta34.x1 = 0.5, beta34.x2 = -0.5, #covariate effects for transiion 34
                                                beta35.x1 = 0.5, beta35.x2 = -0.5, #covariate effects for transiion 35
                                                beta45.x1 = 0.5, beta45.x2 = -0.5, #covariate effects for transiion 45
                                                x.in = x.baseline, #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
                                                numsteps = max.follow)
                       
                       ## Return data
                       data.out
                     })
end_time_fit_parallel <- Sys.time()
diff_fit_parallel <- start_time_fit_parallel - end_time_fit_parallel
stopCluster(cl)
print("FINISH DATA GEN")
diff_fit_parallel



###
### Calculate true transition probabilities for each individual
###

### Write a function to it for a single individual
print(paste("CALC TRUE RISKS ", Sys.time(), sep = ""))
calc.true.tp.ind <- function(row){
  return(calc.true.transition.probs.DGM1(0, ceiling(7*365.25), row[1], row[2],
                                         shape12 = 1, scale12 = scales.sim["12"], #shape and scale for weibull baseline hazard for transition 1 -> 2
                                         shape13 = 1, scale13 = scales.sim["13"], #shape and scale for weibull baseline hazard for transition 1 -> 3
                                         shape15 = 1, scale15 = scales.sim["15"], #shape and scale for weibull baseline hazard for transition 1 -> 5
                                         shape24 = 1, scale24 = scales.sim["24"], #shape and scale for weibull baseline hazard for transition 2 -> 4
                                         shape25 = 1, scale25 = scales.sim["25"], #shape and scale for weibull baseline hazard for transition 2 -> 5
                                         shape34 = 1, scale34 = scales.sim["34"], #shape and scale for weibull baseline hazard for transition 3 -> 4
                                         shape35 = 1, scale35 = scales.sim["35"], #shape and scale for weibull baseline hazard for transition 3 -> 5
                                         shape45 = 1, scale45 = scales.sim["45"], #shape and scale for weibull baseline hazard for transition 3 -> 5
                                         beta12.x1 = 0.5, beta12.x2 = -0.5, #covariate effects for transiion 12
                                         beta13.x1 = 0.5, beta13.x2 = -0.5, #covariate effects for transiion 13
                                         beta15.x1 = 0.5, beta15.x2 = -0.5, #covariate effects for transiion 15
                                         beta24.x1 = 0.5, beta24.x2 = -0.5, #covariate effects for transiion 24
                                         beta25.x1 = 0.5, beta25.x2 = -0.5, #covariate effects for transiion 25
                                         beta34.x1 = 0.5, beta34.x2 = -0.5, #covariate effects for transiion 34
                                         beta35.x1 = 0.5, beta35.x2 = -0.5, #covariate effects for transiion 35
                                         beta45.x1 = 0.5, beta45.x2 = -0.5 #covariate effects for transiion 45
  ))
}

print("START CALC TRUE RISKS")
cl <- makeCluster(5)
registerDoParallel(5)
start_time_ctr_parallel <- Sys.time()
p.true_parallel_list<-(foreach(input=data_parallel_list, .combine=list, .multicombine=TRUE, 
                               .packages=c("gems", "dplyr")) %dopar%{
                                 
                                 ### Calc true risk for each individual in this dataset
                                 p.true <- apply(select(input[["cohort"]], x1, x2), 1, calc.true.tp.ind)
                                 p.true <- data.frame(t(p.true))
                                 colnames(p.true) <- paste("p.true", 1:5, sep = "")
                                 print(paste("FINISH TRUE RISKS ", Sys.time(), sep = ""))
                                 
                                 ## Return data
                                 p.true
                               })
end_time_ctr_parallel <- Sys.time()
diff_ctr_parallel <- start_time_ctr_parallel - end_time_ctr_parallel
stopCluster(cl)
print("FINISH CALC TRUE RISKS")
diff_ctr_parallel


### Combine data into one data frame
temp.data.cohort <- do.call(rbind, lapply(data_parallel_list, function(x) {x[[1]]}))
temp.data.p.true <- do.call(rbind, p.true_parallel_list)

### Assign patient ID's and remove rownames
temp.data.cohort$patid <- 1:nrow(temp.data.cohort)
rownames(temp.data.cohort) <- NULL
rownames(temp.data.p.true) <- NULL

### Combine
temp.data <- cbind(temp.data.cohort, temp.data.p.true)

### Put into an oject that can be uesd by subsequent functions that convert into mstate form
data.pop <- list("cohort" = temp.data, "max.follow" = max.follow)

### Save workspace
rm(list=setdiff(ls(), list("data.pop", "scales.sim")))
save.image("data/generate_population_data_DGM1.RData")