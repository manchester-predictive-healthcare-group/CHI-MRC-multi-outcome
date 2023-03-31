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
print(paste("n.pctls = ", n.cohort, sep = ""))

### Create a table for the mean estimates
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

### Create a table for the mdian estimates
tables.est.median <- vector("list", 3)
names(tables.est.median) <- paste("tpnum", 1:3, sep = "")
for (tpnum in 1:3){
  tables.est.median[[tpnum]] <- vector("list", 5)
  names(tables.est.median[[tpnum]]) <- paste("State", 1:5, sep = "")
  for (state in 1:5){
    tables.est.median[[tpnum]][[state]] <- vector("list", 3)
    names(tables.est.median[[tpnum]][[state]]) <- c("M1C1", "M1C2", "M1C3")
  }
}

### Create a table for the percentile range lower bound
tables.est.lower <- vector("list", 3)
names(tables.est.lower) <- paste("tpnum", 1:3, sep = "")
for (tpnum in 1:3){
  tables.est.lower[[tpnum]] <- vector("list", 5)
  names(tables.est.lower[[tpnum]]) <- paste("State", 1:5, sep = "")
  for (state in 1:5){
    tables.est.lower[[tpnum]][[state]] <- vector("list", 3)
    names(tables.est.lower[[tpnum]][[state]]) <- c("M1C1", "M1C2", "M1C3")
  }
}

### Create a table for the percentile range upper bound
tables.est.upper <- vector("list", 3)
names(tables.est.upper) <- paste("tpnum", 1:3, sep = "")
for (tpnum in 1:3){
  tables.est.upper[[tpnum]] <- vector("list", 5)
  names(tables.est.upper[[tpnum]]) <- paste("State", 1:5, sep = "")
  for (state in 1:5){
    tables.est.upper[[tpnum]][[state]] <- vector("list", 3)
    names(tables.est.upper[[tpnum]][[state]]) <- c("M1C1", "M1C2", "M1C3")
  }
}

###
### Create data for plots
###

### Cycle through preicted transition probabilities
for (tpnum in 1:3){
  ### Cycle through scenario
  for (scen in c(do.call(paste0, expand.grid(c("M1"), c("C1", "C2", "C3"))))){

    ### Load worksapce
    load(paste("data/small_sample_analysis_mean_", scen, "_n", n.cohort, "_npctls", n.pctls, ".RData", sep = ""))
  
    ### Cycle through states
    for (state in 1:5){
      tables.est[[tpnum]][[state]][[scen]] <- c(round(colMeans(calib.aj.int.list[[tpnum]] - calib.true.list[[tpnum]]), 3)[state],
                                                round(colMeans(calib.blr.diff.list[[tpnum]][["IPCW.PSPEC"]] - calib.true.list[[tpnum]]), 3)[state],
                                                round(colMeans(calib.mlr.diff.list[[tpnum]][["IPCW.PSPEC"]] - calib.true.list[[tpnum]]), 3)[state])
      
      tables.est.median[[tpnum]][[state]][[scen]] <- c(round(apply(calib.aj.int.list[[tpnum]] - calib.true.list[[tpnum]], 2, quantile, probs = 0.5), 3)[state],
                                                       round(apply(calib.blr.diff.list[[tpnum]][["IPCW.PSPEC"]] - calib.true.list[[tpnum]], 2, quantile, probs = 0.5), 3)[state],
                                                       round(apply(calib.mlr.diff.list[[tpnum]][["IPCW.PSPEC"]] - calib.true.list[[tpnum]], 2, quantile, probs = 0.5), 3)[state])
      
      tables.est.lower[[tpnum]][[state]][[scen]] <- c(round(apply(calib.aj.int.list[[tpnum]] - calib.true.list[[tpnum]], 2, quantile, probs = 0.025), 3)[state],
                                                   round(apply(calib.blr.diff.list[[tpnum]][["IPCW.PSPEC"]] - calib.true.list[[tpnum]], 2, quantile, probs = 0.025), 3)[state],
                                                   round(apply(calib.mlr.diff.list[[tpnum]][["IPCW.PSPEC"]] - calib.true.list[[tpnum]], 2, quantile, probs = 0.025), 3)[state])
      
      tables.est.upper[[tpnum]][[state]][[scen]] <- c(round(apply(calib.aj.int.list[[tpnum]] - calib.true.list[[tpnum]], 2, quantile, probs = 0.975), 3)[state],
                                                      round(apply(calib.blr.diff.list[[tpnum]][["IPCW.PSPEC"]] - calib.true.list[[tpnum]], 2, quantile, probs = 0.975), 3)[state],
                                                      round(apply(calib.mlr.diff.list[[tpnum]][["IPCW.PSPEC"]] - calib.true.list[[tpnum]], 2, quantile, probs = 0.975), 3)[state])
    }
  }
}

#################
### Mean plot ###
#################

###
### Create table to be fed into dot and whisker plot in correct format
###
data.dw <- vector("list", 3)
names(data.dw) <- paste("tpnum", 1:3, sep = "")
for (tpnum in 1:3){
  data.dw[[tpnum]] <- vector("list", 5)
  names(data.dw[[tpnum]]) <- paste("State", 1:5, sep = "")
  for (state in 1:5){
    data.dw[[tpnum]][[state]] <- data.frame("term" = c(rep("NIC", 3), rep("WIC", 3), rep("SIC", 3)), 
                                            "model" = rep(c("AJ", "BLR-IPCW", "MLR-IPCW"), 3),
                                            "estimate" = c(tables.est[[tpnum]][[state]][["M1C1"]], tables.est[[tpnum]][[state]][["M1C2"]], tables.est[[tpnum]][[state]][["M1C3"]]),
                                            "conf.low" = c(tables.est.lower[[tpnum]][[state]][["M1C1"]], tables.est.lower[[tpnum]][[state]][["M1C2"]], tables.est.lower[[tpnum]][[state]][["M1C3"]]),
                                            "conf.high" = c(tables.est.upper[[tpnum]][[state]][["M1C1"]], tables.est.upper[[tpnum]][[state]][["M1C2"]], tables.est.upper[[tpnum]][[state]][["M1C3"]]),
                                            "tpnum" = rep(c("Perfect", "Over predict", "Under predict")[tpnum], 9),
                                            "state" = rep(paste("State", state, sep = " "), 9))
  }
}


###
### Combine all data for a single plot with facet_grid
###
data.dw.facet <- do.call("rbind", lapply(data.dw, function(x){do.call("rbind", x)}))
data.dw.facet$tpnum <- factor(data.dw.facet$tpnum, levels = c("Perfect", "Over predict", "Under predict"))
gg.dw <- dwplot(data.dw.facet) + #scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
  facet_grid(rows = vars(tpnum), cols = vars(state))  + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  xlab("Bias") + ylab("Scenario")

### Save main plot
CairoPNG(paste("figures/gg_small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_mean.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.dw)
dev.off()
CairoTIFF(paste("figures/gg_small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_mean.tiff", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.dw)
dev.off()
CairoPDF(paste("figures/gg_small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_mean.pdf", sep = ""), 
         width = 15, height = 10)
print(gg.dw)
dev.off()

###################
### Median plot ###
###################

###
### Create table to be fed into dot and whisker plot in correct format
###
data.dw <- vector("list", 3)
names(data.dw) <- paste("tpnum", 1:3, sep = "")
for (tpnum in 1:3){
  data.dw[[tpnum]] <- vector("list", 5)
  names(data.dw[[tpnum]]) <- paste("State", 1:5, sep = "")
  for (state in 1:5){
    data.dw[[tpnum]][[state]] <- data.frame("term" = c(rep("NIC", 3), rep("WIC", 3), rep("SIC", 3)), 
                                            "model" = rep(c("AJ", "BLR-IPCW", "MLR-IPCW"), 3),
                                            "estimate" = c(tables.est.median[[tpnum]][[state]][["M1C1"]], tables.est.median[[tpnum]][[state]][["M1C2"]], tables.est.median[[tpnum]][[state]][["M1C3"]]),
                                            "conf.low" = c(tables.est.lower[[tpnum]][[state]][["M1C1"]], tables.est.lower[[tpnum]][[state]][["M1C2"]], tables.est.lower[[tpnum]][[state]][["M1C3"]]),
                                            "conf.high" = c(tables.est.upper[[tpnum]][[state]][["M1C1"]], tables.est.upper[[tpnum]][[state]][["M1C2"]], tables.est.upper[[tpnum]][[state]][["M1C3"]]),
                                            "tpnum" = rep(c("Perfect", "Over predict", "Under predict")[tpnum], 9),
                                            "state" = rep(paste("State", state, sep = " "), 9))
  }
}


###
### Combine all data for a single plot with facet_grid
###
data.dw.facet <- do.call("rbind", lapply(data.dw, function(x){do.call("rbind", x)}))
data.dw.facet$tpnum <- factor(data.dw.facet$tpnum, levels = c("Perfect", "Over predict", "Under predict"))
gg.dw <- dwplot(data.dw.facet) + #scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
  facet_grid(rows = vars(tpnum), cols = vars(state))  + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  xlab("Bias") + ylab("Scenario")

### Save main plot
CairoPNG(paste("figures/gg_small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_median.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.dw)
dev.off()
CairoTIFF(paste("figures/gg_small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_median.tiff", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.dw)
dev.off()
CairoPDF(paste("figures/gg_small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_median.pdf", sep = ""), 
         width = 15, height = 10)
print(gg.dw)
dev.off()

# ###
# ### Create plots
# ###
# gg.dw <- vector("list", 3)
# names(gg.dw) <- paste("tpnum", 1:3, sep = "")
# for (tpnum in 1:3){
#   gg.dw[[tpnum]] <- vector("list", 5)
#   names(gg.dw[[tpnum]]) <- paste("State", 1:5, sep = "")
#   for (state in 1:5){
#     gg.dw[[tpnum]][[state]] <- dwplot(data.dw[[tpnum]][[state]])
#   }
# }
