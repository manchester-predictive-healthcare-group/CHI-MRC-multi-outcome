### Run the imputation for female cohort, index date of tpye A, and no linkage required for the cohort (primary care data only)

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Load packages
library(dplyr)
library(mice)
library(doParallel)
library(foreach)

### Load the dataset of interest
load("data_Aurum_65plus/data_intermediate/create_cohort_cdm_link.RData")

### Set a gender variable, so we can differentiate between code for men and women, without needing to alter code further down
### female = 1, male = 0
gender.var <- 1

### Assign the dataset we want to impute to a common name so I can re-use functions all programs
data.for.imp <- filter(dat.pre.imp.cdm.link.A, gender == gender.var)
str(data.for.imp)


################################################
### Preliminary stuff: create imputation objects 
################################################

### Run a empty imputation to extract input objects and alter them
imp0<-mice(data.for.imp,maxit=0,print=FALSE)
str(imp0)

### Check logged event warnings
## Variables of issue are gender (constant) and study_dtindex_r (colinear). Neither will be used in imputation, so this is fine.
imp0$loggedEvents

### Extract methods matrix
meth<-imp0$meth
length(meth)
## Using predictive mean matching for continuous variables, and polytomous regression for categorical
meth

### Extract predictor matrix
pred<-imp0$pred
colnames(pred)
## Set variables we don't want to use imputation procedure
## For outcomes, we are just using the presence of the outcome at baseline in the imputation procedure
## For death, we use the censoring indicator, and the NelsenAalen estimator of death at the time of death/censored
pred[,"person_id"] <- 0
pred[,"care_site_id"] <- 0
pred[,"study_dtindex"] <- 0
pred[,"study_dtindex_r"] <- 0
pred[,"gender"] <- 0
pred[,"CKD_hist_t"] <- 0
pred[,"CKD_ev_c"] <- 0
pred[,"CKD_ev_t"] <- 0
pred[,"Diab_t2_hist_t"] <- 0
pred[,"Diab_t2_ev_c"] <- 0
pred[,"Diab_t2_ev_t"] <- 0
pred[,"AF_hist_t"] <- 0
pred[,"AF_ev_c"] <- 0
pred[,"AF_ev_t"] <- 0
pred[,"HF_hist_t"] <- 0
pred[,"HF_ev_c"] <- 0
pred[,"HF_ev_t"] <- 0
pred[,"CHD_MI_hist_t"] <- 0
pred[,"CHD_MI_ev_c"] <- 0
pred[,"CHD_MI_ev_t"] <- 0
pred[,"Stroke_TIA_hist_t"] <- 0
pred[,"Stroke_TIA_ev_c"] <- 0
pred[,"Stroke_TIA_ev_t"] <- 0
pred[,"Death_t"] <- 0

### The prevalence of eating disorders is so low for men and women, we arenot going to use it
pred[,"Eating_disorders"] <- 0
pred[,"Hepatic_failure"] <- 0

### For womena we don't want to use prostate either
if (gender.var == 1){pred[,"Prostate"] <- 0}
pred


### Do another dry run with the new meth and predictor objects to check it's working
imp00 <- mice(data.for.imp, maxit = 0, print = FALSE, meth = meth, pred = pred)

### Check the visiting sequence (only really neccesary when imputing interaction terms using exact values, based on imputed values, etc)
imp00$vis


#######################
### Run the imputation
#######################
time.in <- Sys.time()
print("start imp")

cl <- makeCluster(10)
registerDoParallel(10)
mice.out<-(foreach(seed.num=c(10001,10002,10003,10004,10005,10006,10007,10008,10009,10010), .combine=list, .multicombine=TRUE, 
                .packages=c("dplyr","mice","tidyr"))
        %dopar%{mice(data.for.imp,m=2,meth=meth,pred=pred,seed=seed.num,maxit=20,print=FALSE)
        })
stopCluster(cl)

time.out <- Sys.time()

time.diff <- time.out - time.in

time.in
time.out
time.diff

rm(imp0, imp00)

### Save image
save.image("data_Aurum_65plus/data_intermediate/impute_link_A_female_m20.RData")
print("image saved")


