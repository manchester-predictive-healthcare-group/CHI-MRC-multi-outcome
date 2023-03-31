### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
n.cohort <- as.numeric(args[1])
n.pctls <- as.numeric(args[2])
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))

### Create a table for the estimates
tables.est <- vector("list", 3)
names(tables.est) <- paste("tpnum", 1:3, sep = "")
for (tpnum in 1:3){
  tables.est[[tpnum]] <- vector("list", 5)
  names(tables.est[[tpnum]]) <- paste("State", 1:5, sep = "")
  for (state in 1:5){
    tables.est[[tpnum]][[state]] <- vector("list", 3)
    names(tables.est[[tpnum]][[state]]) <- c("M1C1", "M1C2", "M1C3")
  }
}

### Create a table for the standard error
tables.est.se <- vector("list", 3)
names(tables.est.se) <- paste("tpnum", 1:3, sep = "")
for (tpnum in 1:3){
  tables.est.se[[tpnum]] <- vector("list", 5)
  names(tables.est.se[[tpnum]]) <- paste("State", 1:5, sep = "")
  for (state in 1:5){
    tables.est.se[[tpnum]][[state]] <- vector("list", 3)
    names(tables.est.se[[tpnum]][[state]]) <- c("M1C1", "M1C2", "M1C3")
  }
}

####################
### Do BLR first ###
####################

###
### Create data for plots
###

### Cycle through preicted transition probabilities
for (tpnum in 1:3){
  
  ### Cycle through scenario
  for (scen in c(do.call(paste0, expand.grid(c("M1"), c("C1", "C2", "C3"))))){
    
    ### Load worksapce
    load(paste("data/large_sample_analysis_mean_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
    
    ### Create list of true calibration (estimand)
    calib.true <- list(colMeans(p.true) - colMeans(p.est.mean[[1]]),
                       colMeans(p.true) - colMeans(p.est.mean[[2]]),
                       colMeans(p.true) - colMeans(p.est.mean[[3]]))
    
    ### Cycle through states
    for (state in 1:5){
      tables.est[[tpnum]][[state]][[scen]] <- round(c(calib.blr.ipcw.list[[tpnum]]$diff.pred.obs.pspec[[state]],
                                                      calib.blr.ipcw.list[[tpnum]]$diff.pred.obs.DGMspec[[state]],
                                                      calib.blr.ipcw.list[[tpnum]]$diff.pred.obs.mspec[[state]],
                                                      calib.blr.list[[tpnum]]$diff.pred.obs[[state]]) - calib.true[[tpnum]][state], 3)
      
      tables.est.se[[tpnum]][[state]][[scen]] <- round(c(calib.blr.ipcw.boot.se.list[[tpnum]][["se"]][paste("diff.pred.obs.pspec", state, sep = "")],
                                                         calib.blr.ipcw.boot.se.list[[tpnum]][["se"]][paste("diff.pred.obs.DGMspec", state, sep = "")],
                                                         calib.blr.ipcw.boot.se.list[[tpnum]][["se"]][paste("diff.pred.obs.mspec", state, sep = "")],
                                                         calib.blr.ipcw.boot.se.list[[tpnum]][["se"]][paste("diff.pred.obs", state, sep = "")]), 3)
    }
  }
}

###
### Create table to be fed into dot and whisker plot in correct format
###
data.dw <- vector("list", 3)
names(data.dw) <- paste("tpnum", 1:3, sep = "")
for (tpnum in 1:3){
  data.dw[[tpnum]] <- vector("list", 5)
  names(data.dw[[tpnum]]) <- paste("State", 1:5, sep = "")
  for (state in 1:5){
    data.dw[[tpnum]][[state]] <- data.frame("term" = c(rep("NIC", 4), rep("WIC", 4), rep("SIC", 4)), 
                                            "model" = rep(c("BLR-IPCW", "BLR-IPCW-DGM", "BLR-IPCW-MISS", "BLR"), 3),
                                            "estimate" = c(tables.est[[tpnum]][[state]][["M1C1"]], tables.est[[tpnum]][[state]][["M1C2"]], tables.est[[tpnum]][[state]][["M1C3"]]),
                                            "std.error" = c(tables.est.se[[tpnum]][[state]][["M1C1"]], tables.est.se[[tpnum]][[state]][["M1C2"]], tables.est.se[[tpnum]][[state]][["M1C3"]]),
                                            "tpnum" = rep(c("Perfect", "Over predict", "Under predict")[tpnum], 12),
                                            "state" = rep(paste("State", state, sep = " "), 12))
  }
}


###
### Combine all data for a single plot with facet_grid
###
data.dw.facet <- do.call("rbind", lapply(data.dw, function(x){do.call("rbind", x)}))
data.dw.facet$tpnum <- factor(data.dw.facet$tpnum, levels = c("Perfect", "Over predict", "Under predict"))
gg.dw <- dwplot(data.dw.facet) +# scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
  facet_grid(rows = vars(tpnum), cols = vars(state))  + 
  xlab("Bias") + ylab("Scenario")

### Save main plot
CairoPNG(paste("figures/gg_sens_blr_large_sample_mean_N", n.cohort, "_npctls", n.pctls, ".png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.dw)
dev.off()


###################
### Do MLR next ###
###################

###
### Create data for plots
###

### Cycle through preicted transition probabilities
for (tpnum in 1:3){
  
  ### Cycle through scenario
  for (scen in c(do.call(paste0, expand.grid(c("M1"), c("C1", "C2", "C3"))))){
    
    ### Load worksapce
    load(paste("data/large_sample_analysis_mean_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
    
    ### Create list of true calibration (estimand)
    calib.true <- list(colMeans(p.true) - colMeans(p.est.mean[[1]]),
                       colMeans(p.true) - colMeans(p.est.mean[[2]]),
                       colMeans(p.true) - colMeans(p.est.mean[[3]]))
    
    ### Cycle through states
    for (state in 1:5){
      tables.est[[tpnum]][[state]][[scen]] <- round(c(calib.mlr.ipcw.list[[tpnum]]$diff.pred.obs.pspec[[state]],
                                                      calib.mlr.ipcw.list[[tpnum]]$diff.pred.obs.DGMspec[[state]],
                                                      calib.mlr.ipcw.list[[tpnum]]$diff.pred.obs.mspec[[state]],
                                                      calib.mlr.list[[tpnum]]$diff.pred.obs[[state]]) - calib.true[[tpnum]][state], 3)
      
      tables.est.se[[tpnum]][[state]][[scen]] <- round(c(calib.mlr.ipcw.boot.se.list[[tpnum]][["se"]][paste("diff.pred.obs.pspec", state, sep = "")],
                                                         calib.mlr.ipcw.boot.se.list[[tpnum]][["se"]][paste("diff.pred.obs.DGMspec", state, sep = "")],
                                                         calib.mlr.ipcw.boot.se.list[[tpnum]][["se"]][paste("diff.pred.obs.mspec", state, sep = "")],
                                                         calib.mlr.ipcw.boot.se.list[[tpnum]][["se"]][paste("diff.pred.obs", state, sep = "")]), 3)
    }
  }
}

###
### Create table to be fed into dot and whisker plot in correct format
###
data.dw <- vector("list", 3)
names(data.dw) <- paste("tpnum", 1:3, sep = "")
for (tpnum in 1:3){
  data.dw[[tpnum]] <- vector("list", 5)
  names(data.dw[[tpnum]]) <- paste("State", 1:5, sep = "")
  for (state in 1:5){
    data.dw[[tpnum]][[state]] <- data.frame("term" = c(rep("NIC", 4), rep("WIC", 4), rep("SIC", 4)), 
                                            "model" = rep(c("MLR-IPCW", "MLR-IPCW-DGM", "MLR-IPCW-MISS", "MLR"), 3),
                                            "estimate" = c(tables.est[[tpnum]][[state]][["M1C1"]], tables.est[[tpnum]][[state]][["M1C2"]], tables.est[[tpnum]][[state]][["M1C3"]]),
                                            "std.error" = c(tables.est.se[[tpnum]][[state]][["M1C1"]], tables.est.se[[tpnum]][[state]][["M1C2"]], tables.est.se[[tpnum]][[state]][["M1C3"]]),
                                            "tpnum" = rep(c("Perfect", "Over predict", "Under predict")[tpnum], 12),
                                            "state" = rep(paste("State", state, sep = " "), 12))
  }
}


###
### Combine all data for a single plot with facet_grid
###
data.dw.facet <- do.call("rbind", lapply(data.dw, function(x){do.call("rbind", x)}))
data.dw.facet$tpnum <- factor(data.dw.facet$tpnum, levels = c("Perfect", "Over predict", "Under predict"))
gg.dw <- dwplot(data.dw.facet) +# scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
  facet_grid(rows = vars(tpnum), cols = vars(state))  + 
  xlab("Bias") + ylab("Scenario")

### Save main plot
CairoPNG(paste("figures/gg_sens_mlr_large_sample_mean_N", n.cohort, "_npctls", n.pctls, ".png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.dw)
dev.off()

