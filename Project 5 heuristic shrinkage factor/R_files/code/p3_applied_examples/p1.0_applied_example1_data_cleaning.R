# ##############################################################################

# Author of code: Glen P. Martin.

# This script takes the mimic-iii cohort that was extract from 
# the database (see sql files), applies some data cleaning steps, 
# defines the outcomes to be modelled in the analysis,
# and then imputes missing data using multiple imputation.

# ##############################################################################

### setwd
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_8.2")

library(tidyr)
library(dplyr)
library(readr)

MIMIC_Cohort <- read_csv("data/mimiccohort_withMV.csv",
                         col_names = TRUE,
                         na = c("", "NA", "null", "NULL"))

####----------------------------------------------------------------------------
## Data clean
####----------------------------------------------------------------------------

#Clean age - define a categorical version of age (not normally recommended but 
# used here as MIMIC masked ages over 89 years)
MIMIC_Cohort <- MIMIC_Cohort %>%
  mutate(age_grouped = factor(ifelse(age < 30, "<30",
                                     ifelse(age < 40, "<40",
                                            ifelse(age < 50, "<50",
                                                   ifelse(age < 60, "<60",
                                                          ifelse(age < 70, "<70",
                                                                 ifelse(age < 80, "<80", 
                                                                        ">80")))))),
                              levels = c("<30", "<40", "<50", "<60", "<70", "<80", ">80")),
         .after = age)

#Define factor levels
MIMIC_Cohort <- MIMIC_Cohort %>% 
  mutate(gender = case_when(gender == "M" ~ 1,
                            gender == "F" ~ 0,
                            is.na(gender) ~ NA),
         
         admission_type = case_when(admission_type == "EMERGENCY" ~ 1,
                                    is.na(admission_type) ~ NA,
                                    TRUE ~ 0))


#Clean glucose values
MIMIC_Cohort$glucose_min[which(MIMIC_Cohort$glucose_mean > 140000)] <- NA
MIMIC_Cohort$glucose_max[which(MIMIC_Cohort$glucose_mean > 140000)] <- NA
MIMIC_Cohort$glucose_mean[which(MIMIC_Cohort$glucose_mean > 140000)] <- NA


####----------------------------------------------------------------------------
## Define the outcomes
####----------------------------------------------------------------------------

MIMIC_Cohort <- MIMIC_Cohort %>%
  group_by(subject_id, hadm_id, icustay_id) %>%
  mutate(
    #define AKI using criteria in national guidelines
    AKI = case_when(is.na(creatinine_max_day2_3) | 
                      is.na(creatinine_min) ~ NA,
                    creatinine_max_day2_3 >= (1.5 * creatinine_min) |
                        abs(creatinine_max_day2_3 - creatinine_min) >= 0.3 ~ 1, 
                    TRUE ~ 0),
    
    #MV outcome = more than 48 hours on MV after time zero (so subtract the time 
    #  on MV before index date from the cumulative time on MV across the admission):
    MV = case_when(is.na(cumulative_mv_duration) |
                     is.na(time_on_mv_within_first24hr) ~ NA,
                   cumulative_mv_duration - 
                     time_on_mv_within_first24hr >= (48*60) ~ 1,
                   TRUE ~ 0),
    
    Long_ICU_LOS = case_when(is.na(icu_los_days) ~ NA,
                             #subtract 1 as long LOS is >48 hrs after end of 24hr
                             icu_los_days - 1 >= 2 ~ 1 
                             , TRUE ~ 0),
    
    Mort = case_when(is.na(hospital_expire_flag) ~ NA,
                     hospital_expire_flag == 1 ~ 1,
                     TRUE ~ 0)) %>%
  ungroup() 


####----------------------------------------------------------------------------
## Apply exclusion criteria
####----------------------------------------------------------------------------

Analysis_Cohort <- MIMIC_Cohort %>%
  filter(age >= 18) %>% #patients aged over 18 years (i.e. exclude neonatal)
  filter(icu_los_days > 1) %>% #exclude any ICU admission with a LOS < 24 hours
  filter(icustay_seq == 1) %>% #only include pt first ICU within each hospitalisation
  filter(!is.na(baseline_egfr)) %>%
  filter(baseline_egfr >= 60) %>% #exclude anyone with baseline eGFR <60 
  filter(mv_before_icu != 1) %>% #exclude those who were on MV before admission
  filter(!is.na(AKI)) %>% #exclude cases with missing outcome
  filter(!is.na(MV)) %>%
  filter(!is.na(Long_ICU_LOS)) %>%
  filter(!is.na(Mort))
  

####----------------------------------------------------------------------------
## Subset the variables that we wish to use for modelling 
####----------------------------------------------------------------------------
Analysis_Cohort <- Analysis_Cohort %>%
  select(subject_id, hadm_id, icustay_id,
         age,
         age_grouped,
         gender,
         admission_type,
         ethnicity_grouped,
         ends_with("_mean"),
         AKI,
         MV,
         Long_ICU_LOS,
         Mort) %>%
  select(-meanbp_mean) #remove this since we have both systolic and diastolic BP
  

# #####---------------------------------------------------------------------------
# ### Handle Missing data
# #####---------------------------------------------------------------------------

# Consider complete case analysis on the mimic iii cohort, for simplicity of
#   modelling and bootstrap:
Analysis_Cohort <- Analysis_Cohort %>% 
  filter(complete.cases(.))

saveRDS(Analysis_Cohort, "data/ce_df_mimic.rds")

