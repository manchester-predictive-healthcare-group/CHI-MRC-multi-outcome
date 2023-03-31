### Run the imputation for female cohort, index date of tpye A, and no linkage required for the cohort (primary care data only)

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome")

### Get rid of scientific notation
options(scipen = 999)

### Load packages
library(table1)
library(knitr)
library(dplyr)

### Load the data under the common data model
load("data_Aurum_65plus/data_intermediate/create_cohort_cdm_nolink.RData")
load("data_Aurum_65plus/data_intermediate/create_cohort_cdm_link.RData")

### Datasets are dat.pre.imp.cdm.nolink.A and dat.pre.imp.cdm.nolink.B
#str(dat.pre.imp.cdm.nolink.A)

### Create a function to only calculate mean and sd of continuous variables
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(round(x, 2)), digits = 3, digits.pct = 2), c("", 
       "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}


# ### Apply labels to variables
# label(dat.pre.imp.cdm.nolink.A$Smoking) <- "Smoking status"
# label(dat.pre.imp.cdm.nolink.A$Cholhdl_ratio) <- "Cholesterol/HDL ratio"
# label(dat.pre.imp.cdm.nolink.A$Ethnicity6) <- "Ethnicity"
# ## Medical history
# label(dat.pre.imp.cdm.nolink.A$Alcohol_misuse) <- "Alcohol misuse"
# label(dat.pre.imp.cdm.nolink.A$Eating_disorders) <- "Eating disorders"
# label(dat.pre.imp.cdm.nolink.A$Anxiety_disorders) <- "Anxiety"
# label(dat.pre.imp.cdm.nolink.A$Visual_impairment) <- "Visual impairment"
# label(dat.pre.imp.cdm.nolink.A$Hepatic_failure) <- "Hepatic failure"
# label(dat.pre.imp.cdm.nolink.A$Viral_hepatitis) <- "Viral hepatitis"
# label(dat.pre.imp.cdm.nolink.A$Diverticular) <- "Diverticular disease"
# label(dat.pre.imp.cdm.nolink.A$Hearing_loss) <- "Hearing loss"
# label(dat.pre.imp.cdm.nolink.A$Intellectual_dis) <- "Learning disabilities"
# label(dat.pre.imp.cdm.nolink.A$Perip_vascular) <- "Peripheral vascular disease"
# label(dat.pre.imp.cdm.nolink.A$Substance_misuse) <- "Psychoactive substance misuse"
# label(dat.pre.imp.cdm.nolink.A$Bipolar) <- "Bipolar disorders"
# label(dat.pre.imp.cdm.nolink.A$Thyroid) <- "Thyroid disorders"
# label(dat.pre.imp.cdm.nolink.A$Peptic_ulcer) <- "Peptic ulcer disease"
# label(dat.pre.imp.cdm.nolink.A$Prostate) <- "Prostate disorders"
# label(dat.pre.imp.cdm.nolink.A$Diab_t1) <- "Diabetes tpye 1"
# ## History of outcomes
# label(dat.pre.imp.cdm.nolink.A$CKD_hist) <- "H/O: CKD"
# label(dat.pre.imp.cdm.nolink.A$Diab_t2_hist) <- "H/O: Diabetes type 2"
# label(dat.pre.imp.cdm.nolink.A$AF_hist) <- "H/O: Atrial fibrillation" 
# label(dat.pre.imp.cdm.nolink.A$HF_hist) <- "H/O: Heart Failure"
# label(dat.pre.imp.cdm.nolink.A$CHD_MI_hist) <- "H/O: CHD or MI"
# label(dat.pre.imp.cdm.nolink.A$Stroke_TIA_hist) <- "H/O: Stroke or TIA"
# 
# 
# ### Define levels of categorical variables
# levels(dat.pre.imp.cdm.nolink.A$Smoking) <- c("Never", "Ex", "Current")
# levels(dat.pre.imp.cdm.nolink.A$Cholhdl_ratio) <- "Cholesterol/HDL ratio"
# levels(dat.pre.imp.cdm.nolink.A$Ethnicity6) <- c("White", "Mixed race", "South asian", "Black", "Chinese and other", "Other")
# levels(dat.pre.imp.cdm.nolink.A$IMD) <- c("1 (most deprived)", "2", "3", "4", "5 (least deprived)")
# ## Medical history
# levels(dat.pre.imp.cdm.nolink.A$Alcohol_misuse) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Eating_disorders) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Asthma) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Anxiety_disorders) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Depression) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Visual_impairment) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Bronchiectasis) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Hepatic_failure) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Viral_hepatitis) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Sinusitis) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$COPD) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Dementia) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Diverticular) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Epilepsy) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Hearing_loss) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Hypertension) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$IBS) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Intellectual_dis) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$MS) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Parkinsons) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Perip_vascular) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Psoriasis) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Substance_misuse) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$RA) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Schizophrenia) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Bipolar) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Thyroid) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Peptic_ulcer) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$IBD) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Prostate) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Diab_t1) <- c("Absent", "Present")
# ## History of outcomes
# levels(dat.pre.imp.cdm.nolink.A$CKD_hist) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Diab_t2_hist) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$AF_hist) <- c("Absent", "Present") 
# levels(dat.pre.imp.cdm.nolink.A$HF_hist) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$CHD_MI_hist) <- c("Absent", "Present")
# levels(dat.pre.imp.cdm.nolink.A$Stroke_TIA_hist) <- c("Absent", "Present")

### Want to do continuous and categorical variables separately
create.table1.hist <- function(data.in){
  
  ### Apply labels to variables
  label(data.in$Smoking) <- "Smoking status"
  label(data.in$Cholhdl_ratio) <- "Cholesterol/HDL ratio"
  label(data.in$Ethnicity6) <- "Ethnicity"
  ## Medical history
  label(data.in$Alcohol_misuse) <- "Alcohol misuse"
  label(data.in$Eating_disorders) <- "Eating disorders"
  label(data.in$Anxiety_disorders) <- "Anxiety"
  label(data.in$Visual_impairment) <- "Visual impairment"
  label(data.in$Hepatic_failure) <- "Hepatic failure"
  label(data.in$Viral_hepatitis) <- "Viral hepatitis"
  label(data.in$Diverticular) <- "Diverticular disease"
  label(data.in$Hearing_loss) <- "Hearing loss"
  label(data.in$Intellectual_dis) <- "Learning disabilities"
  label(data.in$Perip_vascular) <- "Peripheral vascular disease"
  label(data.in$Substance_misuse) <- "Psychoactive substance misuse"
  label(data.in$Bipolar) <- "Bipolar disorders"
  label(data.in$Thyroid) <- "Thyroid disorders"
  label(data.in$Peptic_ulcer) <- "Peptic ulcer disease"
  label(data.in$Prostate) <- "Prostate disorders"
  label(data.in$Diab_t1) <- "Diabetes tpye 1"
  ## History of outcomes
  label(data.in$CKD_hist) <- "H/O: CKD"
  label(data.in$Diab_t2_hist) <- "H/O: Diabetes type 2"
  label(data.in$AF_hist) <- "H/O: Atrial fibrillation" 
  label(data.in$HF_hist) <- "H/O: Heart Failure"
  label(data.in$CHD_MI_hist) <- "H/O: CHD or MI"
  label(data.in$Stroke_TIA_hist) <- "H/O: Stroke or TIA"
  
  
  ### Define levels of categorical variables
  levels(data.in$Smoking) <- c("Never", "Ex", "Current")
  levels(data.in$Ethnicity6) <- c("White", "Mixed race", "South asian", "Black", "Chinese and other", "Other")
  levels(data.in$IMD) <- c("1 (most deprived)", "2", "3", "4", "5 (least deprived)")
  ## Medical history
  levels(data.in$Alcohol_misuse) <- c("Absent", "Present")
  levels(data.in$Eating_disorders) <- c("Absent", "Present")
  levels(data.in$Asthma) <- c("Absent", "Present")
  levels(data.in$Anxiety_disorders) <- c("Absent", "Present")
  levels(data.in$Depression) <- c("Absent", "Present")
  levels(data.in$Visual_impairment) <- c("Absent", "Present")
  levels(data.in$Bronchiectasis) <- c("Absent", "Present")
  levels(data.in$Hepatic_failure) <- c("Absent", "Present")
  levels(data.in$Viral_hepatitis) <- c("Absent", "Present")
  levels(data.in$Sinusitis) <- c("Absent", "Present")
  levels(data.in$COPD) <- c("Absent", "Present")
  levels(data.in$Dementia) <- c("Absent", "Present")
  levels(data.in$Diverticular) <- c("Absent", "Present")
  levels(data.in$Epilepsy) <- c("Absent", "Present")
  levels(data.in$Hearing_loss) <- c("Absent", "Present")
  levels(data.in$Hypertension) <- c("Absent", "Present")
  levels(data.in$IBS) <- c("Absent", "Present")
  levels(data.in$Intellectual_dis) <- c("Absent", "Present")
  levels(data.in$MS) <- c("Absent", "Present")
  levels(data.in$Parkinsons) <- c("Absent", "Present")
  levels(data.in$Perip_vascular) <- c("Absent", "Present")
  levels(data.in$Psoriasis) <- c("Absent", "Present")
  levels(data.in$Substance_misuse) <- c("Absent", "Present")
  levels(data.in$RA) <- c("Absent", "Present")
  levels(data.in$Schizophrenia) <- c("Absent", "Present")
  levels(data.in$Bipolar) <- c("Absent", "Present")
  levels(data.in$Thyroid) <- c("Absent", "Present")
  levels(data.in$Peptic_ulcer) <- c("Absent", "Present")
  levels(data.in$IBD) <- c("Absent", "Present")
  levels(data.in$Prostate) <- c("Absent", "Present")
  levels(data.in$Diab_t1) <- c("Absent", "Present")
  ## History of outcomes
  levels(data.in$CKD_hist) <- c("Absent", "Present")
  levels(data.in$Diab_t2_hist) <- c("Absent", "Present")
  levels(data.in$AF_hist) <- c("Absent", "Present") 
  levels(data.in$HF_hist) <- c("Absent", "Present")
  levels(data.in$CHD_MI_hist) <- c("Absent", "Present")
  levels(data.in$Stroke_TIA_hist) <- c("Absent", "Present")
  
  ### Create the table 1
  table1.out.hist <- table1(~ Age +
                         BMI +
                         Cholhdl_ratio +
                         Ethnicity6 +
                         #Ethnicity16 +
                         SBP +
                         Smoking +
                         #Smoking_anyhist
                         IMD +
                         ## Medical history
                         Alcohol_misuse +
                         Eating_disorders +
                         Asthma +
                         Anxiety_disorders +
                         Depression +
                         Visual_impairment +
                         Bronchiectasis +
                         Hepatic_failure +
                         Viral_hepatitis +
                         Sinusitis +
                         COPD +
                         Dementia +
                         Diverticular +
                         Epilepsy +
                         Hearing_loss +
                         Hypertension +
                         IBS +
                         Intellectual_dis +
                         MS +
                         Parkinsons +
                         Perip_vascular +
                         Psoriasis +
                         Substance_misuse +
                         RA +
                         Schizophrenia +
                         Bipolar +
                         Thyroid +
                         Peptic_ulcer +
                         IBD +
                         Prostate +
                         Diab_t1 +
                         ## History of outcomes
                         CKD_hist +
                         Diab_t2_hist +
                         AF_hist + 
                         HF_hist +
                         CHD_MI_hist +
                         Stroke_TIA_hist| gender, 
                          data = data.in, 
                          render.missing = NULL,
                          render.categorical = "FREQ (PCTnoNA%)",
                          render.cont = my.render.cont)
  return(table1.out.hist)
}


### Calculate the tables for predictor variables
table1.nolink.A.hist <- create.table1.hist(dat.pre.imp.cdm.nolink.A)
table1.nolink.B.hist <- create.table1.hist(dat.pre.imp.cdm.nolink.B)  

table1.link.A.hist <- create.table1.hist(dat.pre.imp.cdm.link.A)
table1.link.B.hist <- create.table1.hist(dat.pre.imp.cdm.link.B) 


### Need to create a separate table for outcomes
### Want to calculate total number of events
### Total follow up time
### An incidence rate per 1000 person years
### For each outcome, need to do this  only within the subset of individuals who don't
### have the outcome at baseline.

### Try and turn this into a function
calc.rate <- function(data.in, varname, var.in.time, var.in.ind, var.in.hist, gender.in){
  time.at.risk <- sum(data.in[ , var.in.time][data.in[ , var.in.hist] == 0 & data.in[ , "gender"] == gender.in])
  number.events <- sum((data.in[ , var.in.ind] == 1)[data.in[ , var.in.hist] == 0 & data.in[ , "gender"] == gender.in])
  rate <- number.events*365.25*1000/time.at.risk
  n.at.risk <- sum(data.in[ , var.in.hist] == 0 & data.in[ , "gender"] == gender.in)
  
  output <- c(n.at.risk, time.at.risk, number.events, round(rate,2))
  names(output) <- paste(varname, c(".n.at.risk", ".time.at.risk", ".number.events", ".rate"), sep = "")
  return(output)
}


### Want to do death separately, as this will be calculated for all individuals
calc.rate.death <- function(data.in, gender.in){
  time.at.risk <- sum(data.in[ , "Death_t"][data.in[ , "gender"] == gender.in])
  number.events <- sum((data.in[ , "Death_c"] == 1)[data.in[ , "gender"] == gender.in])
  rate <- number.events*365.25*1000/time.at.risk
  n.at.risk <- sum(data.in[ , "gender"] == gender.in)
  
  output <- c(n.at.risk, time.at.risk, number.events, round(rate,2))
  names(output) <- paste("Death", c(".n.at.risk", ".time.at.risk", ".number.events", ".rate"), sep = "")
  return(output)
}

### Now use these functions to create a table for the output variables
create.table1.events <- function(data.in){

#   data.in <- dat.pre.imp.cdm.link.A
#   data.in <- dat.pre.imp.cdm.nolink.A
# 
#   dat.nolink.1000 <- dat.pre.imp.cdm.nolink.A[1:1000, ]
#   dat.link.1000 <- dat.pre.imp.cdm.link.A[1:1000, ]
  
CKD.rate.male <- calc.rate(data.in, "CKD", "CKD_ev_t", "CKD_ev_c", "CKD_hist", 0)
CKD.rate.female <- calc.rate(data.in, "CKD", "CKD_ev_t", "CKD_ev_c", "CKD_hist", 1)

Diab_t2.rate.male <- calc.rate(data.in, "Diab_t2", "Diab_t2_ev_t", "Diab_t2_ev_c", "Diab_t2_hist", 0)
Diab_t2.rate.female <- calc.rate(data.in, "Diab_t2", "Diab_t2_ev_t", "Diab_t2_ev_c", "Diab_t2_hist", 1)

AF.rate.male <- calc.rate(data.in, "AF", "AF_ev_t", "AF_ev_c", "AF_hist", 0)
AF.rate.female <- calc.rate(data.in, "AF", "AF_ev_t", "AF_ev_c", "AF_hist", 1)

HF.rate.male <- calc.rate(data.in, "HF", "HF_ev_t", "HF_ev_c", "HF_hist", 0)
HF.rate.female <- calc.rate(data.in, "HF", "HF_ev_t", "HF_ev_c", "HF_hist", 1)

CHD_MI.rate.male <- calc.rate(data.in, "CHD_MI", "CHD_MI_ev_t", "CHD_MI_ev_c", "CHD_MI_hist", 0)
CHD_MI.rate.female <- calc.rate(data.in, "CHD_MI", "CHD_MI_ev_t", "CHD_MI_ev_c", "CHD_MI_hist", 1)

Stroke_TIA.rate.male <- calc.rate(data.in, "Stroke_TIA", "Stroke_TIA_ev_t", "Stroke_TIA_ev_c", "Stroke_TIA_hist", 0)
Stroke_TIA.rate.female <- calc.rate(data.in, "Stroke_TIA", "Stroke_TIA_ev_t", "Stroke_TIA_ev_c", "Stroke_TIA_hist", 1)

death.rate.male <- calc.rate.death(data.in, 0)
death.rate.female <- calc.rate.death(data.in, 1)

table1.out.event <- cbind(c(CKD.rate.male,
                            Diab_t2.rate.male,
                            AF.rate.male,
                            HF.rate.male,
                            CHD_MI.rate.male,
                            Stroke_TIA.rate.male,
                            death.rate.male),
                          c(CKD.rate.female,
                            Diab_t2.rate.female,
                            AF.rate.female,
                            HF.rate.female,
                            CHD_MI.rate.female,
                            Stroke_TIA.rate.female,
                            death.rate.female))

colnames(table1.out.event) <- c("male", "female")

return(table1.out.event)}

### Calculate the tables for the outcome events
table1.nolink.A.events <- create.table1.events(dat.pre.imp.cdm.nolink.A)
table1.nolink.B.events <- create.table1.events(dat.pre.imp.cdm.nolink.B)  

table1.link.A.events <- create.table1.events(dat.pre.imp.cdm.link.A)
table1.link.B.events <- create.table1.events(dat.pre.imp.cdm.link.B) 



### Now to do the missing data tables
create.miss <- function(data.in){
  miss.out.male <- data.in[data.in[ , "gender"] == 0, ] %>%
    summarise_all(funs(sum(is.na(.)) / length(.)))
  miss.out.male2 <- miss.out.male[miss.out.male != 0]
  names(miss.out.male2) <- names(miss.out.male)[miss.out.male != 0]
  
  miss.out.female <- data.in[data.in[ , "gender"] == 1, ] %>%
    summarise_all(funs(sum(is.na(.)) / length(.)))
  miss.out.female2 <- miss.out.female[miss.out.female != 0]
  names(miss.out.female2) <- names(miss.out.female)[miss.out.female != 0]
  
  miss.out.final <- rbind(miss.out.male2, miss.out.female2)
  rownames(miss.out.final) <- c("male", "female")
  return(miss.out.final)
}



table.miss.nolink.A <- create.miss(dat.pre.imp.cdm.nolink.A)
table.miss.nolink.B <- create.miss(dat.pre.imp.cdm.nolink.B)  

table.miss.link.A <- create.miss(dat.pre.imp.cdm.link.A)
table.miss.link.B <- create.miss(dat.pre.imp.cdm.link.B) 

table.miss.nolink.A
table.miss.nolink.B
table.miss.link.A
table.miss.link.B


## Remove unncesary stuff
rm(dat.pre.imp.cdm.nolink.A, dat.pre.imp.cdm.nolink.B, dat.pre.imp.cdm.link.A, dat.pre.imp.cdm.link.B)

### Load the data under the common data model
save.image("data_Aurum_65plus/data_intermediate/create_table1.RData")
print("image saved")



# test.1000 <- dat.pre.imp.cdm.nolink.A[1:10000, ]
# ## Outcomes
#                         ## Death
#                         Death_t, Death_c, Death_NelsonAalen, Death_NelsonAalen_link,
#                         ## CKD 
#                         CKD_ev_c, CKD_ev_t, 
#                         ## Diabetes
#                         Diab_t2_ev_c, Diab_t2_ev_t,
#                         ## AF
#                         AF_ev_c, AF_ev_t, 
#                         ## HF
#                         HF_ev_c, HF_ev_t, 
#                         ## CHD/MI
#                         CHD_MI_ev_c, CHD_MI_ev_t,
#                         ## Stroke/TIA
#                        Stroke_TIA_ev_c, Stroke_TIA_ev_t
#                         
#                         
# ## Demographic and test
#                         Age,
#                         gender,
#                         BMI, 
#                         Cholhdl_ratio, 
#                         Ethnicity6, 
#                         #Ethnicity16, 
#                         SBP, 
#                         Smoking, 
#                         #Smoking_anyhist, 
#                         IMD,
#                         ## Medical history
#                         Alcohol_misuse, 
#                         Eating_disorders, 
#                         Asthma, 
#                         Anxiety_disorders, 
#                         Depression, 
#                         Visual_impairment, 
#                         Bronchiectasis, 
#                         Hepatic_failure, 
#                         Viral_hepatitis, 
#                         Sinusitis, 
#                         COPD, 
#                         Dementia,
#                         Diverticular, 
#                         Epilepsy, 
#                         Hearing_loss, 
#                         Hypertension, 
#                         IBS, 
#                         Intellectual_dis, 
#                         MS, 
#                         Parkinsons, 
#                         Perip_vascular, 
#                         Psoriasis, 
#                         Substance_misuse, 
#                         RA,
#                         Schizophrenia, 
#                         Bipolar, 
#                         Thyroid, 
#                         Peptic_ulcer, 
#                         IBD, 
#                         Prostate,
#                         Diab_t1,
#                         ## Outcomes
#                         ## Death
#                         Death_t, Death_c, Death_NelsonAalen, Death_NelsonAalen_link,
#                         ## CKD 
#                         CKD_hist, CKD_hist_t, CKD_ev_c, CKD_ev_t, 
#                         ## Diabetes
#                         Diab_t2_hist, Diab_t2_hist_t, Diab_t2_ev_c, Diab_t2_ev_t,
#                         ## AF
#                         AF_hist, AF_hist_t, AF_ev_c, AF_ev_t, 
#                         ## HF
#                         HF_hist, HF_hist_t, HF_ev_c, HF_ev_t, 
#                         ## CHD/MI
#                         CHD_MI_hist, CHD_MI_hist_t, CHD_MI_ev_c, CHD_MI_ev_t,
#                         ## Stroke/TIA
#                         Stroke_TIA_hist, Stroke_TIA_hist_t, Stroke_TIA_ev_c, Stroke_TIA_ev_t
# 
# rm()


