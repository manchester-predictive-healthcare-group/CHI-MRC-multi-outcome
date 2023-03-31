### This program will create cohorts under the data common model, with variables based only on both primary care an dlinked HES data
### Male/female cohorts, and ones with study_dtindex_A and study_dtindex_B

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome")

### Load packages
library(dplyr)

### First read in a subset of the raw data to work with
dat.pre.imp.A <- read.table("data_Aurum_65plus/data_intermediate/sas_cohort_var_A.csv", 
                          sep="," , header=TRUE)
         
str(dat.pre.imp.A, list.len = ncol(dat.pre.imp.A))
                 
dat.pre.imp.B <- read.table("data_Aurum_65plus/data_intermediate/sas_cohort_var_B.csv", 
                          sep="," , header=TRUE)
                          
str(dat.pre.imp.B, list.len = ncol(dat.pre.imp.B))


### Start by filtering these cohorts to patients with eligible linkage, and actual linkage to HES, LSOA and ONS
nrow(dat.pre.imp.A)
dat.pre.imp.A <- filter(dat.pre.imp.A, linkage_elig == 1, hes_e == 1, death_e == 1, lsoa_e == 1)
nrow(dat.pre.imp.A)

nrow(dat.pre.imp.B)
dat.pre.imp.B <- filter(dat.pre.imp.B, linkage_elig == 1, hes_e == 1, death_e == 1, lsoa_e == 1)
nrow(dat.pre.imp.B)


### Also remove patients without a gender = 0 or 1
dat.pre.imp.A <- filter(dat.pre.imp.A, gender != 2)
dat.pre.imp.B <- filter(dat.pre.imp.B, gender != 2)

### Write a function to convert the data into the common data model
### Note although we extracted t1 diabetes to potentially be an outcome, we are going to use it as a predictor 
### Its name (diab_t1_primhist) has same format as other comorbidities

### Turn this into a function
create_cdm_link <- function(dat.in, cohort.in){
  ### Step 1) Calculate outcomes under common data model names
  
  ### Note that we take the max of primhist_t and heshist_t, as if both events occured, we want whichever one happened first,
  ### which is the bigger number. This also works with the fact that if either didn't occur, then primhist_t = -1, which will
  ### always be the smaller entity
  
  ### On the contrary, time until first event, primev_t and hesev_t, we want to minimum (whichever happens first)
  dat.out <- mutate(dat.in,
                        ## Outcomes
                        CKD_hist = as.factor(pmax(CKD_primhist, CKD_heshist)),
                        CKD_hist_t = pmax(CKD_primhist_t, CKD_heshist_t),
                        CKD_ev_c = as.factor(pmax(CKD_primev_c, CKD_hesev_c)),
                        CKD_ev_t = pmin(CKD_primev_t, CKD_hesev_t),
                        
                        Diab_t2_hist = as.factor(pmax(Diab_t2_primhist, Diab_t2_heshist)),
                        Diab_t2_hist_t = pmax(Diab_t2_primhist_t, Diab_t2_heshist_t),
                        Diab_t2_ev_c = as.factor(pmax(Diab_t2_primev_c, Diab_t2_hesev_c)),
                        Diab_t2_ev_t = pmin(Diab_t2_primev_t, Diab_t2_hesev_t),
                        
                        AF_hist = as.factor(pmax(AF_primhist, AF_heshist)),
                        AF_hist_t = pmax(AF_primhist_t, AF_heshist_t),
                        AF_ev_c = as.factor(pmax(AF_primev_c, AF_hesev_c)),
                        AF_ev_t = pmin(AF_primev_t, AF_hesev_t),
                        
                        HF_hist = as.factor(pmax(HF_primhist, HF_heshist)),
                        HF_hist_t = pmax(HF_primhist_t, HF_heshist_t),
                        HF_ev_c = as.factor(pmax(HF_primev_c, HF_hesev_c)),
                        HF_ev_t = pmin(HF_primev_t, HF_hesev_t),
                        
                        CHD_MI_hist = as.factor(pmax(CHD_MI_primhist, CHD_MI_heshist)),
                        CHD_MI_hist_t = pmax(CHD_MI_primhist_t, CHD_MI_heshist_t),
                        CHD_MI_ev_c = as.factor(pmax(CHD_MI_primev_c, CHD_MI_hesev_c)),
                        CHD_MI_ev_t = pmin(CHD_MI_primev_t, CHD_MI_hesev_t),
                        
                        Stroke_TIA_hist = as.factor(pmax(Stroke_TIA_primhist, Stroke_TIA_heshist)),
                        Stroke_TIA_hist_t = pmax(Stroke_TIA_primhist_t, Stroke_TIA_heshist_t),
                        Stroke_TIA_ev_c = as.factor(pmax(Stroke_TIA_primev_c, Stroke_TIA_hesev_c)),
                        Stroke_TIA_ev_t = pmin(Stroke_TIA_primev_t, Stroke_TIA_hesev_t),
                        
                        Diab_t2_hist = as.factor(pmax(Diab_t2_primhist, Diab_t2_heshist)),
                        Diab_t2_hist_t = pmax(Diab_t2_primhist_t, Diab_t2_heshist_t),
                        Diab_t2_ev_c = as.factor(pmax(Diab_t2_primev_c, Diab_t2_hesev_c)),
                        Diab_t2_ev_t = pmin(Diab_t2_primev_t, Diab_t2_hesev_t),
                        
                        Death_c = as.factor(Death_c),
                        
                        ## Predictors
                        gender = as.factor(gender),
                        Alcohol_misuse = as.factor(Alcohol_misuse_primhist), 
                        Eating_disorders = as.factor(Eating_disorders_primhist), 
                        Asthma = as.factor(Asthma_primhist), 
                        Anxiety_disorders = as.factor(Anxiety_disorders_primhist), 
                        Depression = as.factor(Depression_primhist), 
                        Visual_impairment = as.factor(Visual_impairment_primhist), 
                        Bronchiectasis = as.factor(Bronchiectasis_primhist), 
                        Hepatic_failure = as.factor(Hepatic_failure_primhist), 
                        Viral_hepatitis = as.factor(Viral_hepatitis_primhist), 
                        Sinusitis = as.factor(Sinusitis_primhist), 
                        COPD = as.factor(COPD_primhist), 
                        Dementia = as.factor(Dementia_primhist),
                        Diverticular = as.factor(Diverticular_primhist), 
                        Epilepsy = as.factor(Epilepsy_primhist), 
                        Hearing_loss = as.factor(Hearing_loss_primhist), 
                        Hypertension = as.factor(Hypertension_primhist), 
                        IBS = as.factor(IBS_primhist), 
                        Intellectual_dis = as.factor(Intellectual_dis_primhist), 
                        MS = as.factor(MS_primhist), 
                        Parkinsons = as.factor(Parkinsons_primhist), 
                        Perip_vascular = as.factor(Perip_vascular_primhist), 
                        Psoriasis = as.factor(Psoriasis_primhist), 
                        Substance_misuse = as.factor(Substance_misuse_primhist), 
                        RA = as.factor(RA_primhist),
                        Schizophrenia = as.factor(Schizophrenia_primhist), 
                        Bipolar = as.factor(Bipolar_primhist), 
                        Thyroid = as.factor(Thyroid_primhist), 
                        Peptic_ulcer = as.factor(Peptic_ulcer_primhist), 
                        IBD = as.factor(IBD_primhist), 
                        Prostate = as.factor(Prostate_primhist),
                        Diab_t1 = as.factor(Diab_t1_primhist),
                        Smoking = as.factor(Smoking),
                        Smoking_anyhist = as.factor(Smoking_anyhist),
                        IMD = as.factor(IMD),
                        Ethnicity6 = as.factor(Ethnicity6),
                        Ethnicity16 = as.factor(Ethnicity16) 
                        
  )
  
  ### Step 2) Rename all study_dtindex to a name not dependent on cohort
  if (cohort.in == "A"){
  dat.out <- rename(dat.out, study_dtindex = study_dtindex_A, study_dtindex_r = study_dtindex_A_r)
  }
  
  if (cohort.in == "B"){
  dat.out <- rename(dat.out, study_dtindex = study_dtindex_B, study_dtindex_r = study_dtindex_B_r)
  }
  
  
  ### Step 3) Retain only variables of interest to us
  dat.out <- select(dat.out,
                        ## Not for imputation
                        person_id,
                        care_site_id,
                        study_dtindex,
                        study_dtindex_r,
                        ## Demographic and test
                        Age,
                        gender,
                        BMI, 
                        Cholhdl_ratio, 
                        Ethnicity6, 
                        #Ethnicity16, 
                        SBP, 
                        Smoking, 
                        #Smoking_anyhist, 
                        IMD,
                        ## Medical history
                        Alcohol_misuse, 
                        Eating_disorders, 
                        Asthma, 
                        Anxiety_disorders, 
                        Depression, 
                        Visual_impairment, 
                        Bronchiectasis, 
                        Hepatic_failure, 
                        Viral_hepatitis, 
                        Sinusitis, 
                        COPD, 
                        Dementia,
                        Diverticular, 
                        Epilepsy, 
                        Hearing_loss, 
                        Hypertension, 
                        IBS, 
                        Intellectual_dis, 
                        MS, 
                        Parkinsons, 
                        Perip_vascular, 
                        Psoriasis, 
                        Substance_misuse, 
                        RA,
                        Schizophrenia, 
                        Bipolar, 
                        Thyroid, 
                        Peptic_ulcer, 
                        IBD, 
                        Prostate,
                        Diab_t1,
                        ## Outcomes
                        ## Death
                        Death_t, Death_c, Death_NelsonAalen, Death_NelsonAalen_link,
                        ## CKD 
                        CKD_hist, CKD_hist_t, CKD_ev_c, CKD_ev_t, 
                        ## Diabetes
                        Diab_t2_hist, Diab_t2_hist_t, Diab_t2_ev_c, Diab_t2_ev_t,
                        ## AF
                        AF_hist, AF_hist_t, AF_ev_c, AF_ev_t, 
                        ## HF
                        HF_hist, HF_hist_t, HF_ev_c, HF_ev_t, 
                        ## CHD/MI
                        CHD_MI_hist, CHD_MI_hist_t, CHD_MI_ev_c, CHD_MI_ev_t,
                        ## Stroke/TIA
                        Stroke_TIA_hist, Stroke_TIA_hist_t, Stroke_TIA_ev_c, Stroke_TIA_ev_t)
  
  ### Step 4) Retain only the correct NelsonAalen estimator of death
  ### Given these patients do require linkage, we drop the NelsonAalenestimate of death which was calculated in all patients
  dat.out <- select(dat.out, -Death_NelsonAalen)
  
  print("Returning data")
  return(dat.out)
}


### Apply the function to the cohort
dat.pre.imp.cdm.link.A <- create_cdm_link(dat.in = dat.pre.imp.A, cohort.in = "A")
dat.pre.imp.cdm.link.B <- create_cdm_link(dat.in = dat.pre.imp.B, cohort.in = "B")

print("functions done")
str(dat.pre.imp.cdm.link.A)
str(dat.pre.imp.cdm.link.B)

head(dat.pre.imp.A, n = 5)
head(dat.pre.imp.cdm.link.A, n = 5)

### Remove unnecesary stuff
rm(list=setdiff(ls(),list("dat.pre.imp.cdm.link.A","dat.pre.imp.cdm.link.B")))


### Save image
save.image("data_Aurum_65plus/data_intermediate/create_cohort_cdm_link.RData")


colnames(dat.pre.imp.A)
testdat <- dat.pre.imp.A[,c("person_id", "linkage_elig", "lsoa_e", "hes_e", "death_e")]

n <- 2000000
n - sum(testdat$linkage_elig[1:n])
sum(testdat$linkage_elig[1:n])

n <- 2500000
n - sum(testdat$linkage_elig[1:n])
sum(testdat$linkage_elig[1:n])

n <- 2600000
n - sum(testdat$linkage_elig[1:n])
sum(testdat$linkage_elig[1:n])

n <- 2700000
n - sum(testdat$linkage_elig[1:n])
sum(testdat$linkage_elig[1:n])

sum((testdat$linkage_elig == 0)[1:2700000])
sum((testdat$linkage_elig == 0)[2700001:nrow(testdat)])/(nrow(testdat)-2700001)

sum((testdat$linkage_elig == 0)[1:2700000])/(sum((testdat$linkage_elig == 0)[1:2700000]) + sum((testdat$linkage_elig == 0)[2700001:nrow(testdat)]))

head(testdat, n = 1000)

testdat[100000:101000, ]
n <- 2800000
n - sum(testdat$linkage_elig[1:n])
sum(testdat$linkage_elig[1:n])

n <- 2900000
n - sum(testdat$linkage_elig[1:n])
sum(testdat$linkage_elig[1:n])

n <- 3000000
n - sum(testdat$linkage_elig[1:n])
sum(testdat$linkage_elig[1:n])

