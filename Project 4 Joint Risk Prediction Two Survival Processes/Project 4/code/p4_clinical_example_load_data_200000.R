####################################################
####################################################
####################################################
### START OF FITTING MODELS AND GENERATING RISKS ###
####################################################
####################################################
####################################################

### Load female and male data
load("data_Aurum_65plus/data_intermediate/imputed_datesets_link_A_female_m20.RData")
complete.data.1.female <- complete(mice.comb)

load("data_Aurum_65plus/data_intermediate/imputed_datesets_link_A_male_m20.RData")
complete.data.1.male <- complete(mice.comb)

### Combine into single dataset
complete.data.1 <- rbind(complete.data.1.female, complete.data.1.male)
rm(complete.data.1.female, complete.data.1.male, dat.pre.imp.cdm.link.A, dat.pre.imp.cdm.link.B, data.for.imp, mice.comb)

### Set seed
set.seed(10101)

### Remove all individuals who have already had a CVD or T2D event
complete.data.1 <- subset(complete.data.1,
                          CHD_MI_hist == 0 & Stroke_TIA_hist == 0 & HF_hist == 0 & Diab_t2_hist == 0)

### Create validation and development datasets by sampling at random from this dataset
data.size <- 200000
devel.prop <- 0.5
sample.rows <- sample(1:nrow(complete.data.1), data.size, replace = FALSE)
data.devel <- complete.data.1[sample.rows[1:(data.size*devel.prop)], ]
data.valid <- complete.data.1[sample.rows[(data.size*devel.prop + 1):data.size], ]
rm(complete.data.1)

### Choose time point to evaluate calibration
t.eval <- 3652.50

### Transform datasets ready for analysis
data.devel <- transform.data.cop(data.devel, variables.vec = variables.vec)
data.valid <- transform.data.cop(data.valid, variables.vec = variables.vec)