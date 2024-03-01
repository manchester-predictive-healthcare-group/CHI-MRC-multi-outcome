### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Load packages
source("Project_6/code/sim_function_load_packages.R")

### Read in n.cohort
args <- commandArgs(trailingOnly = T)
n.cohort <- as.numeric(args[1])
print(paste("n.cohort = ", n.cohort, sep = ""))

### Set seed
set.seed(101)

### Load female and male data
load("data_Aurum_65plus/data_intermediate/imputed_datesets_link_A_female_m20.RData")
complete.data.female <- complete(mice.comb)

load("data_Aurum_65plus/data_intermediate/imputed_datesets_link_A_male_m20.RData")
complete.data.male <- complete(mice.comb)

### Combine into single dataset
complete.data <- rbind(complete.data.female, complete.data.male)

### Load in raw dataset which contain censoring dates (will be required, as currently people can have events on date of censoring,
### which causes an issue)
load("data_Aurum_65plus/data_intermediate/sas_cohorts_raw.RData")

### Reduce dataset to just person_id and date of censoring (date of censoring includes date of death)
dat.pre.imp.A.dtcens <- select(dat.pre.imp.A, person_id, dtcens_combdeath_r)

### Merge with complete.data
complete.data <- left_join(complete.data, dat.pre.imp.A.dtcens, by = c("person_id"))

### Save raw dataset, as going to have to do lots of further data steps to get dataset in correct format
complete.data.raw <- complete.data
complete.data <- complete.data.raw
### Remove old datasets clogging up wotkspace
rm(list = setdiff(ls(), list("complete.data", "complete.data.raw", "n.cohort")))


### The issue is that we want to model events that happen on the day of censoring
### If we push the date of censoring to be one day after an events that happen on day of censoring,
### we will pick up on these events and they will be included in the incidence rates. Patients will then be censored
### after 1 day of being in the subsequent state. Individuals that died on the same day as an event, will still move
### into the death state on the day of censoring.

# ### Remove all individuals who have already had a CVD or T2D event
# complete.data <- subset(complete.data,
#                           CHD_MI_hist == 0 & Stroke_TIA_hist == 0 & HF_hist == 0 & Diab_t2_hist == 0)

### Create validation and development datasets by sampling at random from this dataset
# data.size <- 200000
# sample.rows <- sample(1:nrow(complete.data), data.size, replace = FALSE)
# complete.data <- complete.data[sample.rows, ]

### Also need to create new variable called CVD  (time until either have occured)
complete.data$CVD_ev_t <- pmin(complete.data$CHD_MI_ev_t, complete.data$Stroke_TIA_ev_t, complete.data$HF_ev_t)
complete.data$CVD_ev_c <- as.factor(pmax(as.numeric(complete.data$CHD_MI_ev_c == 1), as.numeric(complete.data$Stroke_TIA_ev_c == 1), as.numeric(complete.data$HF_ev_c == 1)))
complete.data$CVD_hist <- as.factor(pmax(as.numeric(complete.data$CHD_MI_hist == 1), 
                                         as.numeric(complete.data$Stroke_TIA_hist == 1), 
                                         as.numeric(complete.data$HF_hist == 1)))

### Currently there are CKD events which happen on day of censoring (approx 1%)
### Wheres there are no CVD or Diab_t2 events happening on the day of censoring (or death)...
### Have double checked, this is because we only look for medical codes prior to date of death (and censoring).
### We get CKD events that are on or beyond date of death, because we didn't apply the same restriction to test data,
### from which some CKD events are identified.

### I am therefore going to censoring all events which happen on the day of censoring (or death) for consistency with the other variables

### Start by creating a censorng variable which is time from study_dtindex
complete.data <- complete.data %>% mutate(dtcens_var = dtcens_combdeath_r - study_dtindex_r)

### Check number of events that happen on same day as censoring or death
sum(complete.data$CKD_ev_t >= complete.data$dtcens_var & complete.data$CKD_ev_c == 1)/sum(complete.data$CKD_ev_c == 1)
sum(complete.data$CKD_ev_t >= complete.data$Death_t & complete.data$CKD_ev_c == 1 & complete.data$Death_c == 1)/sum(complete.data$CKD_ev_c == 1)

### If an event happens on the same day as dtcens_var, set to censored (not this will also cover individuals who have CKD event on same day as death)
complete.data$CKD_ev_c[complete.data$CKD_ev_t >= complete.data$dtcens_var & complete.data$CKD_ev_c == 1] <- 0


### Create each of the outcome variables for each state
complete.data$o2 <- complete.data$CVD_ev_t
complete.data$c2 <- as.numeric(complete.data$CVD_ev_c) - 1

complete.data$o3 <- complete.data$Diab_t2_ev_t
complete.data$c3 <- as.numeric(complete.data$Diab_t2_ev_c) - 1

complete.data$o4 <- complete.data$CKD_ev_t
complete.data$c4 <- as.numeric(complete.data$CKD_ev_c) - 1


## State 5 is CVD + Diab_t2
complete.data$o5 <- pmax(complete.data$o2, complete.data$o3)
complete.data$c5 <- pmin(as.numeric(complete.data$c2 == 1), as.numeric(complete.data$c3 == 1))
#complete.data$c5 <- as.factor(pmin(as.numeric(complete.data$c2 == 1), as.numeric(complete.data$c3 == 1)))


## State 6 is CVD + CKD
complete.data$o6 <- pmax(complete.data$o2, complete.data$o4)
complete.data$c6 <- pmin(as.numeric(complete.data$c2 == 1), as.numeric(complete.data$c4 == 1))
#complete.data$c6 <- as.factor(pmin(as.numeric(complete.data$c2 == 1), as.numeric(complete.data$c4 == 1)))


## State 7 is Diab_t2 + CKD
complete.data$o7 <- pmax(complete.data$o3, complete.data$o4)
complete.data$c7 <- pmin(as.numeric(complete.data$c3 == 1), as.numeric(complete.data$c4 == 1))
#complete.data$c7 <- as.factor(pmin(as.numeric(complete.data$c3 == 1), as.numeric(complete.data$c4 == 1)))


## State 8 is CVD + Diab_t2 + CKD
complete.data$o8 <- pmax(complete.data$o2, complete.data$o3, complete.data$o4)
complete.data$c8 <- pmin(as.numeric(complete.data$c2 == 1), as.numeric(complete.data$c3 == 1), as.numeric(complete.data$c4 == 1))
#complete.data$c8 <- as.factor(pmin(as.numeric(complete.data$c2 == 1), as.numeric(complete.data$c3 == 1), as.numeric(complete.data$c4 == 1)))

## Create death
complete.data$o9 <- complete.data$Death_t
complete.data$c9 <- as.numeric(complete.data$Death_c) - 1
#complete.data$c9 <- complete.data$Death_c

# complete.data <- complete.data[, c(1:5, 43:ncol(complete.data))]
# complete.data.raw <- complete.data
str(complete.data)

### If:
### events 2 + 3 happen same time as 5
### events 2 + 4 happen same time as 6
### events 3 + 4 happen same time as 7
### then set to non-event
complete.data$c2[(complete.data$o2 == complete.data$o5 & 
                    complete.data$c2 == 1 & complete.data$c5 == 1)] <- 0
complete.data$c3[(complete.data$o3 == complete.data$o5 & 
                    complete.data$c3 == 1 & complete.data$c5 == 1)] <- 0

complete.data$c2[(complete.data$o2 == complete.data$o6 & 
                    complete.data$c2 == 1 & complete.data$c6 == 1)] <- 0
complete.data$c4[(complete.data$o4 == complete.data$o6 & 
                    complete.data$c4 == 1 & complete.data$c6 == 1)] <- 0

complete.data$c3[(complete.data$o3 == complete.data$o7 & 
                    complete.data$c3 == 1 & complete.data$c7 == 1)] <- 0
complete.data$c4[(complete.data$o4 == complete.data$o7 & 
                    complete.data$c4 == 1  & complete.data$c7 == 1)] <- 0


### If:
### events 2 + 3 + 4 happen same time as 8
### then set to non-event
# complete.data$c2[(complete.data$o2 == complete.data$o8 &
#                complete.data$c2 == 1 & complete.data$c8 == 1)] <- 0
# complete.data$c3[(complete.data$o3 == complete.data$o8 &
#                     complete.data$c3 == 1 & complete.data$c8 == 1)] <- 0
# complete.data$c4[(complete.data$o4 == complete.data$o8 &
#                     complete.data$c4 == 1 & complete.data$c8 == 1)] <- 0
complete.data$c5[(complete.data$o5 == complete.data$o8 &
                    complete.data$c5 == 1 & complete.data$c8 == 1)] <- 0
complete.data$c6[(complete.data$o6 == complete.data$o8 &
                    complete.data$c6 == 1 & complete.data$c8 == 1)] <- 0
complete.data$c7[(complete.data$o7 == complete.data$o8 &
                    complete.data$c7 == 1 & complete.data$c8 == 1)] <- 0


### Now reduce states 5, 6, 7 and 8 into a single multimorbidity state
complete.data$o5 <- pmin(complete.data$o5, complete.data$o6, complete.data$o7, complete.data$o8)
complete.data$c5 <- pmax(complete.data$c5, complete.data$c6, complete.data$c7, complete.data$c8)

### Turn state 6 into hte death staet
complete.data$o6 <- complete.data$o9
complete.data$c6 <- complete.data$c9

### First create date of censoring
### Create transition matrix
tmat <- transMat(x = list(c(2, 3, 4, 5, 6), #1
                          c(5, 6), #2
                          c(5, 6), #3
                          c(5, 6), #4
                          c(6), #5
                          c()), #6
                 names = 1:6)
tmat

### Define covariates we want to use 
covs <- c("Age", "gender", "Cholhdl_ratio", "BMI", "SBP", "Smoking", "IMD", "Hypertension", "Depression", "Alcohol_misuse")

### Create msprep object

### Select patients we are going to extract (develop model on 100000, validate on 100000)
ids <- sample(1:nrow(complete.data), 1000000, replace = FALSE)
ids.devel <- ids[1:n.cohort]
ids.valid <- ids[100001:200000]
person.ids.devel <- complete.data$person_id[ids.devel]
person.ids.valid <- complete.data$person_id[ids.valid]

complete.data.prep.devel <- msprep(data = subset(complete.data, person_id %in% person.ids.devel), 
                             trans = tmat,
                             time = c(NA, paste("o", 2:6, sep = "")),
                             status = c(NA, paste("c", 2:6, sep = "")),
                             keep = c(covs, "person_id"))

complete.data.prep.valid <- msprep(data = subset(complete.data, person_id %in% person.ids.valid), 
                                   trans = tmat,
                                   time = c(NA, paste("o", 2:6, sep = "")),
                                   status = c(NA, paste("c", 2:6, sep = "")),
                                   keep = c(covs, "person_id"))
paste("msprep", Sys.time())

### Check how many events there are for each transition
events(complete.data.prep.devel)
Sys.time()

### Expand covariates and create Tstart, Tstop variables, etc
complete.data.prep.devel <- expand.covs(complete.data.prep.devel, covs, longnames = FALSE)
complete.data.prep.valid <- expand.covs(complete.data.prep.valid, covs, longnames = FALSE)
str(complete.data.prep.devel)
paste("expand.covs", Sys.time())

### Fit the Cox model for each transition
### These are the cause specific hazards (rather than competing risks/transition hazards of the msm)
### Note we don't want to model the transitions 4, 5, 6, 7, 11, 15, 19
### These are the transition that are double jumps (i.e. not in our model, and only happen due to peculiarities in the data)
adj.msm <- paste(paste("Age.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("Cholhdl_ratio.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("BMI.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("SBP.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("gender.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("Hypertension.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("Depression.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("Alcohol_misuse.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("Smoking1.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("Smoking2.", (1:12)[-c(4)], sep = "", collapse = "+"),  
                 paste("IMD1.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("IMD2.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("IMD3.", (1:12)[-c(4)], sep = "", collapse = "+"), 
                 paste("IMD4.", (1:12)[-c(4)], sep = "", collapse = "+"),  
                 sep = "+")
formula.msm <- paste("Surv(Tstart, Tstop, status) ~ ", adj.msm, " + strata(trans)", sep = "")
msm.cox.fit <- coxph(as.formula(formula.msm), data = complete.data.prep.devel)
paste("cox.fit", Sys.time())

### Add a column called strata (required when making predictors for new patients
complete.data.prep.devel$strata <- complete.data.prep.devel$trans
complete.data.prep.valid$strata <- complete.data.prep.valid$trans

### Save the work space with the cause specific hazards, want to parallelise the process of generating transition hazards
### and transition probabilities as it takes ages
rm(list=setdiff(ls(), list("complete.data", "complete.data.prep.devel", "complete.data.prep.valid", 
                           "msm.cox.fit", "covs", "ids.devel", "ids.valid", "ids", 
                           "person.ids.devel", "person.ids.valid", "tmat", "n.cohort")))
save.image(paste("Project_6/data/ce/ce_fit_csh_msm_N", n.cohort, ".RData", sep = ""))
print("IMAGE SAVED")