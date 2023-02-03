### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Source packages
source("code/sim_function_load_packages.R")

### Load data
load("data_Aurum_65plus/data_intermediate/create_cohort_cdm_link.RData")

### We are interseted in incidence of CHD_MI and Diabetes type 2
### Therefore want to remove individuals who have got ahistory of either of these things
dat.pre.imp.cdm.link.A <- dat.pre.imp.cdm.link.A[(dat.pre.imp.cdm.link.A$CHD_MI_hist == 0) & 
                                                   (dat.pre.imp.cdm.link.A$Diab_t2_hist == 0) , ]


### Create male and female cohorts
data.anal.female <- filter(dat.pre.imp.cdm.link.A, gender == 1)
data.anal.male <- filter(dat.pre.imp.cdm.link.A, gender == 0)
data.anal <- dat.pre.imp.cdm.link.A
rm(dat.pre.imp.cdm.link.A, dat.pre.imp.cdm.link.B)



#####################
### ENTIRE COHORT ###
#####################


### Also need to create new variable called CHD_MI_Diab (time until both have occured)
data.anal$CHD_MI_Diab_ev_t <- pmax(data.anal$CHD_MI_ev_t, data.anal$Diab_t2_ev_t)
data.anal$CHD_MI_Diab_ev_c <- as.factor(pmin(as.numeric(data.anal$CHD_MI_ev_c == 1), as.numeric(data.anal$Diab_t2_ev_c == 1)))


### Fit survival models
surv.obj.CHD_MI <- survfit(Surv(CHD_MI_ev_t, CHD_MI_ev_c) ~ 1, data = data.anal)
surv.obj.Diab_t2 <- survfit(Surv(Diab_t2_ev_t, Diab_t2_ev_c) ~ 1, data = data.anal)
surv.obj.CHD_MI_Diab <- survfit(Surv(CHD_MI_Diab_ev_t, CHD_MI_Diab_ev_c) ~ 1, data = data.anal)

### Calculate survival probabilities and risks
surv.10y.CHD_MI <- exp(-surv.obj.CHD_MI$cumhaz[min(which(surv.obj.CHD_MI$time > 3652.5))])
surv.10y.Diab_t2 <- exp(-surv.obj.Diab_t2$cumhaz[min(which(surv.obj.Diab_t2$time > 3652.5))])
surv.10y.CHD_MI_Diab <- exp(-surv.obj.CHD_MI_Diab$cumhaz[min(which(surv.obj.CHD_MI_Diab$time > 3652.5))])

risk.10y.CHD_MI <- 1 - surv.10y.CHD_MI
risk.10y.Diab_t2 <- 1 - surv.10y.Diab_t2
risk.10y.CHD_MI_Diab <- 1 - surv.10y.CHD_MI_Diab

print("risk CHD_MI")
risk.10y.CHD_MI
print("risk Diab_t2")
risk.10y.Diab_t2
print("risk CHD_MI Diab joint")
risk.10y.CHD_MI_Diab
print("risk CHD_MI Diab naive")
risk.10y.CHD_MI*risk.10y.Diab_t2

Sys.time()
save.image("Project 4/data/sim_calculate_overall_risk_CPRD_CHD_MI_Diab_cohortA.RData")

##############
### FEMALE ###
##############

### Also need to create new variable called CHD_MI_Diab (time until both have occured)
data.anal.female$CHD_MI_Diab_ev_t <- pmax(data.anal.female$CHD_MI_ev_t, data.anal.female$Diab_t2_ev_t)
data.anal.female$CHD_MI_Diab_ev_c <- as.factor(pmin(as.numeric(data.anal.female$CHD_MI_ev_c == 1), as.numeric(data.anal.female$Diab_t2_ev_c == 1)))


### Fit survival models
surv.obj.CHD_MI <- survfit(Surv(CHD_MI_ev_t, CHD_MI_ev_c) ~ 1, data = data.anal.female)
surv.obj.Diab_t2 <- survfit(Surv(Diab_t2_ev_t, Diab_t2_ev_c) ~ 1, data = data.anal.female)
surv.obj.CHD_MI_Diab <- survfit(Surv(CHD_MI_Diab_ev_t, CHD_MI_Diab_ev_c) ~ 1, data = data.anal.female)

### Calculate survival probabilities and risks
surv.10y.CHD_MI.female <- exp(-surv.obj.CHD_MI$cumhaz[min(which(surv.obj.CHD_MI$time > 3652.5))])
surv.10y.Diab_t2.female <- exp(-surv.obj.Diab_t2$cumhaz[min(which(surv.obj.Diab_t2$time > 3652.5))])
surv.10y.CHD_MI_Diab.female <- exp(-surv.obj.CHD_MI_Diab$cumhaz[min(which(surv.obj.CHD_MI_Diab$time > 3652.5))])

risk.10y.CHD_MI.female <- 1 - surv.10y.CHD_MI.female
risk.10y.Diab_t2.female <- 1 - surv.10y.Diab_t2.female
risk.10y.CHD_MI_Diab.female <- 1 - surv.10y.CHD_MI_Diab.female

print("risk CHD_MI")
risk.10y.CHD_MI.female
print("risk Diab_t2")
risk.10y.Diab_t2.female
print("risk CHD_MI Diab joint")
risk.10y.CHD_MI_Diab.female
print("risk CHD_MI Diab naive")
risk.10y.CHD_MI.female*risk.10y.Diab_t2.female

Sys.time()
save.image("Project 4/data/sim_calculate_overall_risk_CPRD_CHD_MI_Diab_cohortA.RData")

############
### MALE ###
############


### Also need to create new variable called CHD_MI_Diab (time until both have occured)
data.anal.male$CHD_MI_Diab_ev_t <- pmax(data.anal.male$CHD_MI_ev_t, data.anal.male$Diab_t2_ev_t)
data.anal.male$CHD_MI_Diab_ev_c <- as.factor(pmin(as.numeric(data.anal.male$CHD_MI_ev_c == 1), as.numeric(data.anal.male$Diab_t2_ev_c == 1)))


### Fit survival models
surv.obj.CHD_MI <- survfit(Surv(CHD_MI_ev_t, CHD_MI_ev_c) ~ 1, data = data.anal.male)
surv.obj.Diab_t2 <- survfit(Surv(Diab_t2_ev_t, Diab_t2_ev_c) ~ 1, data = data.anal.male)
surv.obj.CHD_MI_Diab <- survfit(Surv(CHD_MI_Diab_ev_t, CHD_MI_Diab_ev_c) ~ 1, data = data.anal.male)

### Calculate survival probabilities and risks
surv.10y.CHD_MI.male <- exp(-surv.obj.CHD_MI$cumhaz[min(which(surv.obj.CHD_MI$time > 3652.5))])
surv.10y.Diab_t2.male <- exp(-surv.obj.Diab_t2$cumhaz[min(which(surv.obj.Diab_t2$time > 3652.5))])
surv.10y.CHD_MI_Diab.male <- exp(-surv.obj.CHD_MI_Diab$cumhaz[min(which(surv.obj.CHD_MI_Diab$time > 3652.5))])

risk.10y.CHD_MI.male <- 1 - surv.10y.CHD_MI.male
risk.10y.Diab_t2.male <- 1 - surv.10y.Diab_t2.male
risk.10y.CHD_MI_Diab.male <- 1 - surv.10y.CHD_MI_Diab.male

print("risk CHD_MI")
risk.10y.CHD_MI.male
print("risk Diab_t2")
risk.10y.Diab_t2.male
print("risk CHD_MI Diab joint")
risk.10y.CHD_MI_Diab.male
print("risk CHD_MI Diab naive")
risk.10y.CHD_MI.male*risk.10y.Diab_t2.male

Sys.time()
save.image("Project 4/data/sim_calculate_overall_risk_CPRD_CHD_MI_Diab_cohortA.RData")
