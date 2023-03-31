###
### This program will create combined plots for the large sample analyses
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
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

### Load prepped dat afor pv (this contains the appropriate datasets, and also estimate of Aalen-Johansen)
if (scen %in% c("M1C1", "M1C2", "M1C3")){
  load(paste("data/large_sample_analysis_moderate_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
} else if (!(scen %in% c("M1C1", "M1C2", "M1C3"))){
  load(paste("data/large_sample_analysis_moderate_N", n.cohort, "_", scen, ".RData", sep = ""))
}
source("code/z_functions.R")

### Run seperately for each p.est (iterate over i)
for (i in 1:3){
  
  print(paste("i = ", i, Sys.time(), sep = " "))
  
  ###
  ### Create plots to validate BLR and MLR weighting (sensitivity analyis, what happens when we don't use weighting, or weighing is miss-specified)
  ###
  
  ###
  ### BLR first 
  
  ### Put into a list
  plots.blr.sens.list <- c(calib.blr.ipcw.list[[i]][["plots.pspec.list"]],
                           calib.blr.ipcw.list[[i]][["plots.DGMspec.list"]],
                           calib.blr.ipcw.list[[i]][["plots.mspec.list"]],
                           calib.blr.list[[i]][["plots.list"]])
  
  ### Create main plots
  plots.blr.sens.gg <- ggarrange(plotlist = plots.blr.sens.list, nrow = 4, ncol = 5)
  plots.blr.sens.gg <- 
    annotate_figure(plots.blr.sens.gg, left = "       BLR                           BLR-IPCW-MISS                             BLR-IPCW-DGM                             BLR-IPCW")
  
  ### Save main plot
  CairoPNG(paste("figures/gg_large_sample_moderate_sens_blr_N", n.cohort, "_", scen, "_est", i, ".png", sep = ""), 
           dpi = 300, width = 15, height = 10, unit = "in")
  print(plots.blr.sens.gg)
  dev.off()
  
  ###
  ### MLR second 
  
  ### Put into a list
  plots.mlr.sens.list <- c(calib.mlr.ipcw.list[[i]][["plots.pspec.list"]],
                           calib.mlr.ipcw.list[[i]][["plots.DGMspec.list"]],
                           calib.mlr.ipcw.list[[i]][["plots.mspec.list"]],
                           calib.mlr.list[[i]][["plots.list"]])
  
  ### Create main plots
  plots.mlr.sens.gg <- ggarrange(plotlist = plots.mlr.sens.list, nrow = 4, ncol = 5)
  plots.mlr.sens.gg <- 
    annotate_figure(plots.mlr.sens.gg, left = "       MLR                           MLR-IPCW-MISS                             MLR-IPCW-DGM                             MLR-IPCW")
  
  ### Save main plot
  CairoPNG(paste("figures/gg_large_sample_moderate_sens_mlr_N", n.cohort, "_", scen, "_est", i, ".png", sep = ""), 
           dpi = 300, width = 15, height = 10, unit = "in")
  print(plots.mlr.sens.gg)
  dev.off()
  
}
