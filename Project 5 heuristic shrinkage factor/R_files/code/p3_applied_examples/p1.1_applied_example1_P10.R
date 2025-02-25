### setwd
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_8.2")

### Load libraries
library(rms)
library(boot)

### Set seed
set.seed(555)

### Rename
df <- readRDS("data/ce_df_mimic.rds")

###
colnames(df)

###
str(df$AKI)

### Build a model to predict AKI

### Set n.iter
n.iter <- 1000

### set n.vec
n.vec <- c(100, 250, 500, 1000)

### Set number of iterations for bootstrapped
n.S.boot <- 500

### Set predictors to adjust for in models 
predictors <- c("age", "gender", "bicarbonate_mean", "creatinine_mean", "chloride_mean", "hemoglobin_mean", "platelet_mean", "potassium_mean", "ptt_mean", "inr_mean")
#,"pt_mean" , "bun_mean", "wbc_mean", "heartrate_mean", "sysbp_mean", "diasbp_mean", "resprate_mean", "tempc_mean", "spo2_mean", "glucose_mean"

### Create output
data.output <- vector("list", length(n.vec))
names(data.output) <- n.vec
for (i in 1:length(n.vec)){
  data.output[[i]] <- data.frame("S.VH" = rep(NA, n.iter), "S.boot" = rep(NA, n.iter), "R2.CS.app" = rep(NA, n.iter), "C.app" = rep(NA, n.iter))
}

###
### Define a function to calculate S using bootstrapping
###
function_calc_S_boot <- function(data, i){
  
  ### Take the bootstrap sample
  boot.dat <- data[i, ]
  
  ### Create model in bootstrapped dataset
  boot.model <- glm(as.formula(paste("AKI ~ ", paste(predictors, collapse = "+"), sep = "")), data = boot.dat, family = binomial(link = "logit")) 
  
  ### Generate predictions using this model using the new dataset
  boot.lp <- predict(boot.model, newdata = data.temp)
  
  ### Create a temporary dataset with both these things
  boot.dat.devel.temp <- data.frame(data.temp, "boot.lp" = boot.lp)
  
  ### Calculate calibration slope
  boot.calib.model <- glm(AKI ~ boot.lp, data = boot.dat.devel.temp, family = binomial(link = "logit")) 
  
  ### Save slope
  return(boot.calib.model$coefficients["boot.lp"])
  
}

### Now run the simulation
### 1) Sample n individuals
### 2) Fit a model in this dataset
### 3) Extract data, record S.VH and S.boot, C.app and R2.CS.app

### For each sample size
for (i in 1:length(n.vec)){
  
  print(paste("samp size =", n.vec[i], Sys.time(), sep = " "))
  
  ### Assign sample sizze
  n.samp <- n.vec[i]
  
  ### Run the sim
  for (iter in 1:n.iter){
   
    print(paste("iter =", iter, Sys.time(), sep = " "))
    
    ### Sample df at random
    data.temp <- df[sample(1:nrow(df), n.samp, replace = FALSE), ]
    
    ### Build model
    model.temp <- lrm(as.formula(paste("AKI ~ ", paste(predictors, collapse = "+"), sep = "")), data = data.temp)
    
    ### Record required statistics for prediction of S.pop
    LR <- model.temp$deviance[1] - model.temp$deviance[2]
    C.app <- model.temp$stats["C"]
    R2.CS.app <- as.numeric(1 - exp(-(LR)/nrow(data.temp)))
    P <- length(model.temp$coefficients)
    S.VH <- as.numeric(1 - P/LR)
    
    ###
    ### Calculate S.boot
    
    ### Run the bootstrapping
    boot.out <- boot(data.temp, function_calc_S_boot, R = n.S.boot)
    
    ### Assign S.boot
    S.boot <- mean(boot.out$t)
    
    ###
    ### Save the output
    data.output[[i]][iter, ] <- c(as.numeric(S.VH), as.numeric(S.boot),  as.numeric(R2.CS.app), as.numeric(C.app))
    
  }
  
  save.image("data/applied_example1_P10.RData")
}


### Calculate means
data.output.means <- lapply(data.output, colMeans, na.rm = TRUE)

### Finally, get C-statistic of model in entire dataset
biggest.model <- lrm(as.formula(paste("AKI ~ ", paste(predictors, collapse = "+"), sep = "")), data = df)
biggest.C.app <- biggest.model$stats["C"]
biggest.C.app

### Get mean C-statistics across the bootstrap iterations
lapply(1:7, function(x){mean(data.output[[x]][,"C.app"])})

rm(list=setdiff(ls(), list("data.output", "data.output.means", "biggest.C.app")))
save.image("data/applied_example1_P10.RData")
print(paste("FINISHED", Sys.time()))

