###
### This program runs the assessment of mean calibration in the small sample simulation
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
scen <- args[1]
n.cohort <- as.numeric(args[2])
n.pctls <- as.numeric(args[3])

print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))

### Load simulated data
load(paste("data/apply_cens_", scen, ".RData", sep = ""))
source("code/z_functions.R")

### Choose number of simulations
n.sim <- 1000

### Set seed
set.seed(1001)

###
### Calibration using binary logistic regression at time point t (with and without IPCW)
###
print(paste("START CALIB BLR ", Sys.time(), sep = ""))

### Create list to store calibration of each estimates
### (include a list to store true calibration in each simulation iteration)
calib.blr.int.list <- vector("list", 4)
calib.blr.slopes.list <- vector("list", 4)
calib.blr.diff.list <- vector("list", 4)
calib.mlr.int.list <- vector("list", 4)
calib.mlr.slopes.list <- vector("list", 4)
calib.mlr.diff.list <- vector("list", 4)
calib.aj.int.list <- vector("list", 4)
calib.true.list <- vector("list", 4)

### Calculate the calibration intecepts and slopes
for (i in 1:3){
  
  ### Create dataframe to store intercepts and slopes for BLR/MLR (No weights, IPCW MPSEC, IPCW PSPEC, IPCW DGM SPEC)
  calib.blr.int.list[[i]] <- vector("list", 4)
  names(calib.blr.int.list[[i]]) <- c("NO.WEIGHT", "IPCW.MSPEC", "IPCW.PSPEC", "IPCW.DGMSPEC")
  
  calib.blr.slopes.list[[i]] <- vector("list", 4)
  names(calib.blr.slopes.list[[i]]) <- c("NO.WEIGHT", "IPCW.MSPEC", "IPCW.PSPEC", "IPCW.DGMSPEC")
  
  calib.blr.diff.list[[i]] <- vector("list", 4)
  names(calib.blr.diff.list[[i]]) <- c("NO.WEIGHT", "IPCW.MSPEC", "IPCW.PSPEC", "IPCW.DGMSPEC")
  
  calib.mlr.int.list[[i]] <- vector("list", 4)
  names(calib.mlr.int.list[[i]]) <- c("NO.WEIGHT", "IPCW.MSPEC", "IPCW.PSPEC", "IPCW.DGMSPEC")
  
  calib.mlr.slopes.list[[i]] <- vector("list", 4)
  names(calib.mlr.slopes.list[[i]]) <- c("NO.WEIGHT", "IPCW.MSPEC", "IPCW.PSPEC", "IPCW.DGMSPEC")
  
  calib.mlr.diff.list[[i]] <- vector("list", 4)
  names(calib.mlr.diff.list[[i]]) <- c("NO.WEIGHT", "IPCW.MSPEC", "IPCW.PSPEC", "IPCW.DGMSPEC")
  
  for (j in 1:4){
    calib.blr.int.list[[i]][[j]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.blr.int.list[[i]][[j]]) <- paste("state", 1:5, sep = "")
    calib.blr.slopes.list[[i]][[j]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.blr.slopes.list[[i]][[j]]) <- paste("state", 1:5, sep = "")
    calib.blr.diff.list[[i]][[j]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.blr.diff.list[[i]][[j]]) <- paste("state", 1:5, sep = "")
    calib.mlr.int.list[[i]][[j]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.mlr.int.list[[i]][[j]]) <- paste("state", 1:5, sep = "")
    calib.mlr.slopes.list[[i]][[j]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.mlr.slopes.list[[i]][[j]]) <- paste("state", 1:5, sep = "")
    calib.mlr.diff.list[[i]][[j]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.mlr.diff.list[[i]][[j]]) <- paste("state", 1:5, sep = "")
  }
  
  ### Create data frames for Aalen-Johansen estimator and true calibration
  calib.aj.int.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
  colnames(calib.aj.int.list[[i]]) <- paste("state", 1:5, sep = "")
  
  calib.true.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
  colnames(calib.true.list[[i]]) <- paste("state", 1:5, sep = "")
}

###
### Run the simulation loop
###
for (sim in 1:n.sim){
  
  print(paste("cohort size = ", n.cohort, ", sim =", sim, Sys.time(), sep = " "))
  
  ### Choose patids at random
  patids.vec <- sample((1:1000000)[-390375], n.cohort, replace = FALSE)
  # We exclude individual 390375, as this person was removed from M1C2 dataset as they caused an error

  ### Extract 100,000 individuals, which is what we're doing the large sample analysis with
  data.mstate.reduc <- data.mstate.obj[["data.mstate"]][data.mstate.obj[["data.mstate"]]$patid %in% patids.vec, ]
  data.raw.reduc <- data.mstate.obj[["data.raw"]][data.mstate.obj[["data.raw"]]$patid %in% patids.vec, ]
  tmat <- data.mstate.obj[["tmat"]]
  
  ###
  ### Extract the true risks into their own object
  ###
  p.true <- data.raw.reduc %>% select(paste("p.true", 1:5, sep = ""))
  
  ###
  ### Create miss-calibrated risks (although they still need to sum to 1...)
  ### Need to get creative for this
  ###
  p.est.perf <- p.true
  p.est.over <- exp((log(p.true/(1-p.true)) + 0.5))/(1 + exp((log(p.true/(1-p.true)) + 0.5)))
  p.est.under <- exp((log(p.true/(1-p.true)) - 0.5))/(1 + exp((log(p.true/(1-p.true)) - 0.5)))
  
  ### ALSO WANT TO HAVE SEPERATE DEVIATIONS WHEN CONSIDERING CHANGES TO MEAN CALIBRATION AND MODERATE CALIBRATION
  ### Put these into a list
  p.est <- list(p.est.perf, p.est.over, p.est.under)
  
  ###
  ### Calculate AJ within percentiles of individuals grouped by predicted risk
  ###
  
  ### Define max.state
  max.state <- 5
  
  ### Create a list to store the sorted data, the Aalen-Johansen estimator for each group, 
  data.sort.pctls <- vector("list", max.state)
  obs.aj.pctls <- vector("list", max.state)
  
  ### Loop through the states of interest
  for (state in 1:max.state){
    print(state)
    ###
    ### Create n.pctls datasets, grouped by predicted risk
    ###
    
    ### Arrange data by risk of category of interest (note predicted risks always share same order with true risks, we can can sort by these)
    p.est.data.frame <- arrange(data.raw.reduc, !!sym(paste("p.true", state, sep = "")))
    
    ### Create a list of length n.pctls
    data.sort.pctls[[state]] <- vector("list", n.pctls)
    group.size <- nrow(p.est.data.frame)/n.pctls
    for (i in 1:n.pctls){
      data.sort.pctls[[state]][[i]] <- slice(p.est.data.frame, (((i-1)*group.size)+1):(i*group.size))
    }
    
    ###
    ### Calculate the Aalen-Johansen estimator within each group
    ###
    
    ### Create a list to store them
    obs.aj.pctls[[state]] <- vector("list", n.pctls)
    
    ### Calculate obs.aj for each group
    for (i in 1:n.pctls){
      print(paste("Calc obs.aj = ", state, "pctl = ", i, Sys.time()))
      obs.aj.pctls[[state]][[i]] <- calc.calib.aj(data.mstate = subset(data.mstate.reduc, patid %in% data.sort.pctls[[state]][[i]]$patid), 
                                                  tmat = tmat, t.eval = t.eval)$obs.aj
    }
  }
  
  ### Create an object with mean of the Aalen-Johansen estimator when ordered by risk of each category
  obs.aj.pctls.mean.df <- lapply(obs.aj.pctls, function(x) {colMeans(do.call("rbind", x))}) 
  
  ### Put into a vector
  obs.aj.pctls.mean <- c(obs.aj.pctls.mean.df[[1]][1],
                         obs.aj.pctls.mean.df[[2]][2],
                         obs.aj.pctls.mean.df[[3]][3],
                         obs.aj.pctls.mean.df[[4]][4],
                         obs.aj.pctls.mean.df[[5]][5])
  ###
  ### Calculate calibration for each type of predicted risk
  ###
  for (i in 1:3){
    
    ###
    ### Calculate the intercepts and slopes using mlr
    ###
    temp.calib.blr <- calc.calib.blr(data.mstate.reduc,
                                     data.raw.reduc,
                                     t.eval = t.eval,
                                     p.est = p.est[[i]])
    
    temp.calib.blr.ipcw <- calc.calib.blr.ipcw(data.mstate.reduc,
                                               data.raw.reduc,
                                               t.eval = t.eval,
                                               p.est = p.est[[i]])
    
    ### Assign to output
    calib.blr.int.list[[i]][["NO.WEIGHT"]][sim, ] <- temp.calib.blr$int
    calib.blr.int.list[[i]][["IPCW.MSPEC"]][sim, ] <- temp.calib.blr.ipcw$int.mspec
    calib.blr.int.list[[i]][["IPCW.PSPEC"]][sim, ] <- temp.calib.blr.ipcw$int.pspec
    calib.blr.int.list[[i]][["IPCW.DGMSPEC"]][sim, ] <- temp.calib.blr.ipcw$int.DGMspec
    
    calib.blr.slopes.list[[i]][["NO.WEIGHT"]][sim, ] <- temp.calib.blr$slopes
    calib.blr.slopes.list[[i]][["IPCW.MSPEC"]][sim, ] <- temp.calib.blr.ipcw$slopes.mspec
    calib.blr.slopes.list[[i]][["IPCW.PSPEC"]][sim, ] <- temp.calib.blr.ipcw$slopes.pspec
    calib.blr.slopes.list[[i]][["IPCW.DGMSPEC"]][sim, ] <- temp.calib.blr.ipcw$slopes.DGMspec
    
    calib.blr.diff.list[[i]][["NO.WEIGHT"]][sim, ] <- temp.calib.blr$diff.pred.obs
    calib.blr.diff.list[[i]][["IPCW.MSPEC"]][sim, ] <- temp.calib.blr.ipcw$diff.pred.obs.mspec
    calib.blr.diff.list[[i]][["IPCW.PSPEC"]][sim, ] <- temp.calib.blr.ipcw$diff.pred.obs.pspec
    calib.blr.diff.list[[i]][["IPCW.DGMSPEC"]][sim, ] <- temp.calib.blr.ipcw$diff.pred.obs.DGMspec
    
    ###
    ### Calculate the intercepts and slopes using MLR
    ###
    temp.calib.mlr <- calc.calib.mlr(data.mstate.reduc,
                                     data.raw.reduc,
                                     t.eval = t.eval,
                                     p.est = p.est[[i]])
    
    temp.calib.mlr.ipcw <- calc.calib.mlr.ipcw(data.mstate.reduc,
                                               data.raw.reduc,
                                               t.eval = t.eval,
                                               p.est = p.est[[i]])
    
    ### Assign to output
    calib.mlr.int.list[[i]][["NO.WEIGHT"]][sim, 2:5] <- temp.calib.mlr$int
    calib.mlr.int.list[[i]][["IPCW.MSPEC"]][sim, 2:5] <- temp.calib.mlr.ipcw$int.mspec
    calib.mlr.int.list[[i]][["IPCW.PSPEC"]][sim, 2:5] <- temp.calib.mlr.ipcw$int.pspec
    calib.mlr.int.list[[i]][["IPCW.DGMSPEC"]][sim, 2:5] <- temp.calib.mlr.ipcw$int.DGMspec
    
    calib.mlr.slopes.list[[i]][["NO.WEIGHT"]][sim, 2:5] <- temp.calib.mlr$slopes
    calib.mlr.slopes.list[[i]][["IPCW.MSPEC"]][sim, 2:5] <- temp.calib.mlr.ipcw$slopes.mspec
    calib.mlr.slopes.list[[i]][["IPCW.PSPEC"]][sim, 2:5] <- temp.calib.mlr.ipcw$slopes.pspec
    calib.mlr.slopes.list[[i]][["IPCW.DGMSPEC"]][sim, 2:5] <- temp.calib.mlr.ipcw$slopes.DGMspec
    
    calib.mlr.diff.list[[i]][["NO.WEIGHT"]][sim, ] <- temp.calib.mlr$diff.pred.obs
    calib.mlr.diff.list[[i]][["IPCW.MSPEC"]][sim, ] <- temp.calib.mlr.ipcw$diff.pred.obs.mspec
    calib.mlr.diff.list[[i]][["IPCW.PSPEC"]][sim, ] <- temp.calib.mlr.ipcw$diff.pred.obs.pspec
    calib.mlr.diff.list[[i]][["IPCW.DGMSPEC"]][sim, ] <- temp.calib.mlr.ipcw$diff.pred.obs.DGMspec
    
    ###
    ### Aalen-Johansen estimator
    ###
    
    ### Deduct AJ estimate from mean risk
    calib.aj.int.list[[i]][sim, ] <- obs.aj.pctls.mean - colMeans(p.est[[i]])
    
    ### 
    ### Save true risks
    ###
    calib.true.list[[i]][sim, ] <- colMeans(p.true - p.est[[i]])
  }
}

###
### Save image
###
rm(list = setdiff(ls(), list("scen", "n.cohort", "n.sim", "n.pctls",
                             "calib.blr.int.list",
                             "calib.blr.slopes.list",
                             "calib.blr.diff.list",
                             "calib.mlr.int.list",
                             "calib.mlr.slopes.list",
                             "calib.mlr.diff.list",
                             "calib.aj.int.list",
                             "calib.true.list")))

save.image(paste("data/small_sample_analysis_mean_", scen, "_n", n.cohort, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")




