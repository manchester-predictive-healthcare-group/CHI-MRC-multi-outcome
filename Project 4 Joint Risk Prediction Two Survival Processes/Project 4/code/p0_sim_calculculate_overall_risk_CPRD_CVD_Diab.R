### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Source packages
source("code/sim_function_load_packages.R")

### Load data
load("data_Aurum_65plus/data_intermediate/create_cohort_cdm_link.RData")

### We are interseted in incidence of CHD_MI + Stroke_TIA = CVD and Diabetes type 2
### Therefore want to remove individuals who have got ahistory of either of these things
dat.pre.imp.cdm.link.A <- dat.pre.imp.cdm.link.A[(dat.pre.imp.cdm.link.A$CHD_MI_hist == 0) & 
                                                   (dat.pre.imp.cdm.link.A$Diab_t2_hist == 0) &
                                                   (dat.pre.imp.cdm.link.A$Stroke_TIA_hist == 0) , ]


### Create male and female cohorts
data.anal.female <- filter(dat.pre.imp.cdm.link.A, gender == 1)
data.anal.male <- filter(dat.pre.imp.cdm.link.A, gender == 0)
data.anal <- dat.pre.imp.cdm.link.A
rm(dat.pre.imp.cdm.link.A, dat.pre.imp.cdm.link.B)



#####################
### ENTIRE COHORT ###
#####################

### Also need to create new variable called CVD  (time until either have occured)
data.anal$CVD_ev_t <- pmin(data.anal$CHD_MI_ev_t, data.anal$Stroke_TIA_ev_t)
data.anal$CVD_ev_c <- as.factor(pmax(as.numeric(data.anal$CHD_MI_ev_c == 1), as.numeric(data.anal$Stroke_TIA_ev_c == 1)))

### Also need to create new variable called CVD_Diab  (time until both have occured)
data.anal$CVD_Diab_ev_t <- pmax(data.anal$CVD_ev_t, data.anal$Diab_t2_ev_t)
data.anal$CVD_Diab_ev_c <- as.factor(pmin(as.numeric(data.anal$CVD_ev_c == 1), as.numeric(data.anal$Diab_t2_ev_c == 1)))


### Fit survival models
surv.obj.CVD <- survfit(Surv(CVD_ev_t, CVD_ev_c) ~ 1, data = data.anal)
surv.obj.CHD_MI <- survfit(Surv(CHD_MI_ev_t, CHD_MI_ev_c) ~ 1, data = data.anal)
surv.obj.Stroke_TIA <- survfit(Surv(Stroke_TIA_ev_t, Stroke_TIA_ev_c) ~ 1, data = data.anal)
surv.obj.Diab_t2 <- survfit(Surv(Diab_t2_ev_t, Diab_t2_ev_c) ~ 1, data = data.anal)
surv.obj.CVD_Diab <- survfit(Surv(CVD_Diab_ev_t, CVD_Diab_ev_c) ~ 1, data = data.anal)


### Calculate survival probabilities and risks
surv.10y.CVD <- exp(-surv.obj.CVD$cumhaz[min(which(surv.obj.CVD$time > 3652.5))])
surv.10y.CHD_MI <- exp(-surv.obj.CHD_MI$cumhaz[min(which(surv.obj.CHD_MI$time > 3652.5))])
surv.10y.Stroke_TIA <- exp(-surv.obj.Stroke_TIA$cumhaz[min(which(surv.obj.Stroke_TIA$time > 3652.5))])
surv.10y.Diab_t2 <- exp(-surv.obj.Diab_t2$cumhaz[min(which(surv.obj.Diab_t2$time > 3652.5))])
surv.10y.CVD_Diab <- exp(-surv.obj.CVD_Diab$cumhaz[min(which(surv.obj.CVD_Diab$time > 3652.5))])

risk.10y.CVD <- 1 - surv.10y.CVD
risk.10y.CHD_MI <- 1 - surv.10y.CHD_MI
risk.10y.Stroke_TIA <- 1 - surv.10y.Stroke_TIA
risk.10y.Diab_t2 <- 1 - surv.10y.Diab_t2
risk.10y.CVD_Diab <- 1 - surv.10y.CVD_Diab

print("risk CHD MI")
risk.10y.CHD_MI
print("risk Stroke TIA")
risk.10y.Stroke_TIA
print("risk CVD")
risk.10y.CVD
print("risk Diab")
risk.10y.Diab_t2
print("risk CVD Diab joint")
risk.10y.CVD_Diab
print("risk CVD Diab naive")
risk.10y.CVD*risk.10y.Diab_t2

Sys.time()
save.image("Project 4/data/sim_calculate_overall_risk_CPRD_CVD_Diab_cohortA.RData")

##############
### FEMALE ###
##############

### Also need to create new variable called CVD  (time until either have occured)
data.anal.female$CVD_ev_t <- pmin(data.anal.female$CHD_MI_ev_t, data.anal.female$Stroke_TIA_ev_t)
data.anal.female$CVD_ev_c <- as.factor(pmax(as.numeric(data.anal.female$CHD_MI_ev_c == 1), as.numeric(data.anal.female$Stroke_TIA_ev_c == 1)))

### Also need to create new variable called CVD_Diab  (time until both have occured)
data.anal.female$CVD_Diab_ev_t <- pmax(data.anal.female$CVD_ev_t, data.anal.female$Diab_t2_ev_t)
data.anal.female$CVD_Diab_ev_c <- as.factor(pmin(as.numeric(data.anal.female$CVD_ev_c == 1), as.numeric(data.anal.female$Diab_t2_ev_c == 1)))


### Fit survival models
surv.obj.CVD <- survfit(Surv(CVD_ev_t, CVD_ev_c) ~ 1, data = data.anal.female)
surv.obj.CHD_MI <- survfit(Surv(CHD_MI_ev_t, CHD_MI_ev_c) ~ 1, data = data.anal.female)
surv.obj.Stroke_TIA <- survfit(Surv(Stroke_TIA_ev_t, Stroke_TIA_ev_c) ~ 1, data = data.anal.female)
surv.obj.Diab_t2 <- survfit(Surv(Diab_t2_ev_t, Diab_t2_ev_c) ~ 1, data = data.anal.female)
surv.obj.CVD_Diab <- survfit(Surv(CVD_Diab_ev_t, CVD_Diab_ev_c) ~ 1, data = data.anal.female)


### Calculate survival probabilities and risks
surv.10y.CVD.female <- exp(-surv.obj.CVD$cumhaz[min(which(surv.obj.CVD$time > 3652.5))])
surv.10y.CHD_MI.female <- exp(-surv.obj.CHD_MI$cumhaz[min(which(surv.obj.CHD_MI$time > 3652.5))])
surv.10y.Stroke_TIA.female <- exp(-surv.obj.Stroke_TIA$cumhaz[min(which(surv.obj.Stroke_TIA$time > 3652.5))])
surv.10y.Diab_t2.female <- exp(-surv.obj.Diab_t2$cumhaz[min(which(surv.obj.Diab_t2$time > 3652.5))])
surv.10y.CVD_Diab.female <- exp(-surv.obj.CVD_Diab$cumhaz[min(which(surv.obj.CVD_Diab$time > 3652.5))])

risk.10y.CVD.female <- 1 - surv.10y.CVD.female
risk.10y.CHD_MI.female <- 1 - surv.10y.CHD_MI.female
risk.10y.Stroke_TIA.female <- 1 - surv.10y.Stroke_TIA.female
risk.10y.Diab_t2.female <- 1 - surv.10y.Diab_t2.female
risk.10y.CVD_Diab.female <- 1 - surv.10y.CVD_Diab.female

print("risk CHD MI")
risk.10y.CHD_MI.female
print("risk Stroke TIA")
risk.10y.Stroke_TIA.female
print("risk CVD")
risk.10y.CVD.female
print("risk Diab")
risk.10y.Diab_t2.female
print("risk CVD Diab joint")
risk.10y.CVD_Diab.female
print("risk CVD Diab naive")
risk.10y.CVD.female*risk.10y.Diab_t2.female

Sys.time()
save.image("Project 4/data/sim_calculate_overall_risk_CPRD_CVD_Diab_cohortA.RData")

############
### MALE ###
############

### Also need to create new variable called CVD  (time until either have occured)
data.anal.male$CVD_ev_t <- pmin(data.anal.male$CHD_MI_ev_t, data.anal.male$Stroke_TIA_ev_t)
data.anal.male$CVD_ev_c <- as.factor(pmax(as.numeric(data.anal.male$CHD_MI_ev_c == 1), as.numeric(data.anal.male$Stroke_TIA_ev_c == 1)))

### Also need to create new variable called CVD_Diab  (time until both have occured)
data.anal.male$CVD_Diab_ev_t <- pmax(data.anal.male$CVD_ev_t, data.anal.male$Diab_t2_ev_t)
data.anal.male$CVD_Diab_ev_c <- as.factor(pmin(as.numeric(data.anal.male$CVD_ev_c == 1), as.numeric(data.anal.male$Diab_t2_ev_c == 1)))


### Fit survival models
surv.obj.CVD <- survfit(Surv(CVD_ev_t, CVD_ev_c) ~ 1, data = data.anal.male)
surv.obj.CHD_MI <- survfit(Surv(CHD_MI_ev_t, CHD_MI_ev_c) ~ 1, data = data.anal.male)
surv.obj.Stroke_TIA <- survfit(Surv(Stroke_TIA_ev_t, Stroke_TIA_ev_c) ~ 1, data = data.anal.male)
surv.obj.Diab_t2 <- survfit(Surv(Diab_t2_ev_t, Diab_t2_ev_c) ~ 1, data = data.anal.male)
surv.obj.CVD_Diab <- survfit(Surv(CVD_Diab_ev_t, CVD_Diab_ev_c) ~ 1, data = data.anal.male)


### Calculate survival probabilities and risks
surv.10y.CVD.male <- exp(-surv.obj.CVD$cumhaz[min(which(surv.obj.CVD$time > 3652.5))])
surv.10y.CHD_MI.male <- exp(-surv.obj.CHD_MI$cumhaz[min(which(surv.obj.CHD_MI$time > 3652.5))])
surv.10y.Stroke_TIA.male <- exp(-surv.obj.Stroke_TIA$cumhaz[min(which(surv.obj.Stroke_TIA$time > 3652.5))])
surv.10y.Diab_t2.male <- exp(-surv.obj.Diab_t2$cumhaz[min(which(surv.obj.Diab_t2$time > 3652.5))])
surv.10y.CVD_Diab.male <- exp(-surv.obj.CVD_Diab$cumhaz[min(which(surv.obj.CVD_Diab$time > 3652.5))])

risk.10y.CVD.male <- 1 - surv.10y.CVD.male
risk.10y.CHD_MI.male <- 1 - surv.10y.CHD_MI.male
risk.10y.Stroke_TIA.male <- 1 - surv.10y.Stroke_TIA.male
risk.10y.Diab_t2.male <- 1 - surv.10y.Diab_t2.male
risk.10y.CVD_Diab.male <- 1 - surv.10y.CVD_Diab.male

print("risk CHD MI")
risk.10y.CHD_MI.male
print("risk Stroke TIA")
risk.10y.Stroke_TIA.male
print("risk CVD")
risk.10y.CVD.male
print("risk Diab")
risk.10y.Diab_t2.male
print("risk CVD Diab joint")
risk.10y.CVD_Diab.male
print("risk CVD Diab naive")
risk.10y.CVD.male*risk.10y.Diab_t2.male

Sys.time()
save.image("Project 4/data/sim_calculate_overall_risk_CPRD_CVD_Diab_cohortA.RData")
