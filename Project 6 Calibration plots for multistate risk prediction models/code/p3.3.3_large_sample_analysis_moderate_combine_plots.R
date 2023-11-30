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

### Load prepped data for pv (this contains the appropriate datasets, and also estimate of Aalen-Johansen)
if (scen %in% c("M1C1", "M1C2", "M1C3")){
  load(paste("data/large_sample_analysis_moderate_N", n.cohort, "_", scen, "_npctls", n.pctls, ".RData", sep = ""))
} else if (!(scen %in% c("M1C1", "M1C2", "M1C3"))){
  load(paste("data/large_sample_analysis_moderate_N", n.cohort, "_", scen, ".RData", sep = ""))
}
source("code/z_functions.R")

#############################################################################################################
### First create plots where perfect, over and under predicting transition probabilities are on each row, ###
### with a different graph for each method                                                                ###
#############################################################################################################

###
### Create main plot BLR-IPCW
###

### Put into a list
plots.comb.list <- c(calib.blr.ipcw.list[[1]][["plots.pspec.list"]],
                     calib.blr.ipcw.list[[2]][["plots.pspec.list"]],
                     calib.blr.ipcw.list[[3]][["plots.pspec.list"]])
plots.comb.gg <- ggarrange(plotlist = plots.comb.list, nrow = 3, ncol = 5, common.legend = TRUE, legend = "bottom")

plots.comb.gg <- 
  annotate_figure(plots.comb.gg, left = "         Under predicting                                               Over predicting                                               Perfectly predicting        ")

### Save main plot
CairoPNG(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_BLRIPCW.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(plots.comb.gg)
dev.off()
CairoTIFF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_BLRIPCW.tiff", sep = ""), 
          dpi = 300, width = 15, height = 10, unit = "in")
print(plots.comb.gg)
dev.off()
CairoPDF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_BLRIPCW.pdf", sep = ""), 
         width = 15, height = 10)
print(plots.comb.gg)
dev.off()


###
### Create main plot MLR-IPCW
###

### Put into a list
plots.comb.list <- c(calib.mlr.ipcw.list[[1]][["plots.pspec.list"]],
                     calib.mlr.ipcw.list[[2]][["plots.pspec.list"]],
                     calib.mlr.ipcw.list[[3]][["plots.pspec.list"]])
plots.comb.gg <- ggarrange(plotlist = plots.comb.list, nrow = 3, ncol = 5, common.legend = TRUE, legend = "bottom")
plots.comb.gg <- 
  annotate_figure(plots.comb.gg, left = "         Under predicting                                               Over predicting                                               Perfectly predicting        ")

### Save main plot
CairoPNG(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_MLRIPCW.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(plots.comb.gg)
dev.off()
CairoTIFF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_MLRIPCW.tiff", sep = ""), 
          dpi = 300, width = 15, height = 10, unit = "in")
print(plots.comb.gg)
dev.off()
CairoPDF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_MLRIPCW.pdf", sep = ""), 
         width = 15, height = 10)
print(plots.comb.gg)
dev.off()

###
### Create main plot pseudo-value
###

### Put into a list
plots.comb.list <- c(calib.pv.list[[1]][["plots.list"]],
                     calib.pv.list[[2]][["plots.list"]],
                     calib.pv.list[[3]][["plots.list"]])
plots.comb.gg <- ggarrange(plotlist = plots.comb.list, nrow = 3, ncol = 5, common.legend = TRUE, legend = "bottom")
plots.comb.gg <- 
  annotate_figure(plots.comb.gg, left = "         Under predicting                                               Over predicting                                               Perfectly predicting        ")

### Save main plot
CairoPNG(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_PV.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(plots.comb.gg)
dev.off()
CairoTIFF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_PV.tiff", sep = ""), 
          dpi = 300, width = 15, height = 10, unit = "in")
print(plots.comb.gg)
dev.off()
CairoPDF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_PV.pdf", sep = ""), 
         width = 15, height = 10)
print(plots.comb.gg)
dev.off()


#############################################################################################################
### First create plots where perfect, over and under predicting transition probabilities are on each row,
### Now to do a seperate plot for each type of predicted transition probability, where the methods are compared
### row by row.
#############################################################################################################
print("START PLOTS PART 2")

### Run seperately for each p.est (iterate over i)
for (i in 1:3){
  
  print(paste("i = ", i, Sys.time(), sep = " "))
  
  ###
  ### Create main plot (BLR IPCW, MLR IPCW, AJ and PV)
  ###
  
  ### Run through each scenario
  if (scen %in% c("M1C1", "M1C2", "M1C3")){
    ### Put each analysis into seperate ggarrange object with its own legend
    plots.pv <- ggarrange(plotlist = calib.pv.list[[i]][["plots.list"]], nrow = 1, ncol = 5, common.legend = TRUE, legend = "bottom")
    plots.blr <- ggarrange(plotlist = calib.blr.ipcw.list[[i]][["plots.pspec.list"]], nrow = 1, ncol = 5, common.legend = TRUE, legend = "bottom")
    plots.mlr <- ggarrange(plotlist = calib.mlr.ipcw.list[[i]][["plots.pspec.list"]], nrow = 1, ncol = 5, common.legend = TRUE, legend = "bottom")
    plots.comb.gg <- ggarrange(plots.pv, plots.blr, plots.mlr, nrow = 3)
    plots.comb.gg <- 
      annotate_figure(plots.comb.gg, left = "         MLR-IPCW                                                     BLR-IPCW                                                       Pseudo-value ")
    
    ### Save main plot
    CairoPNG(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_est", i, ".png", sep = ""), 
             dpi = 300, width = 15, height = 10.5, unit = "in")
    print(plots.comb.gg)
    dev.off()
    CairoTIFF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_est", i, ".tiff", sep = ""), 
              dpi = 300, width = 15, height = 10.5, unit = "in")
    print(plots.comb.gg)
    dev.off()
    CairoPDF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_est", i, ".pdf", sep = ""), 
             width = 15, height = 10.5)
    print(plots.comb.gg)
    dev.off()
    
  } else if (!(scen %in% c("M1C1", "M1C2", "M1C3"))){
    ### Put each analysis into seperate ggarrange object with its own legend
    plots.aj <- ggarrange(calib.aj[[i]][["plots.list"]], nrow = 1, ncol = 5, common.legend = TRUE, legend = "bottom")
    plots.blr <- ggarrange(calib.blr.ipcw.list[[i]][["plots.pspec.list"]], nrow = 1, ncol = 5, common.legend = TRUE, legend = "bottom")
    plots.mlr <- ggarrange(calib.mlr.ipcw.list[[i]][["plots.pspec.list"]], nrow = 1, ncol = 5, common.legend = TRUE, legend = "bottom")
    plots.comb.gg <- ggarrange(plots.aj, plots.blr, plots.mlr, nrow = 3)
    plots.comb.gg <- 
      annotate_figure(plots.comb.gg, left = "            MLR-IPCW                                                                 BLR-IPCW                                                                      AJ          ")
    
    ### Save main plot
    CairoPNG(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_est", i, ".png", sep = ""), 
             dpi = 300, width = 15, height = 10, unit = "in")
    print(plots.comb.gg)
    dev.off()
    CairoTIFF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_est", i, ".tiff", sep = ""), 
              dpi = 300, width = 15, height = 10, unit = "in")
    print(plots.comb.gg)
    dev.off()
    CairoPDF(paste("figures/gg_large_sample_moderate_main_N", n.cohort, "_", scen, "_npctls", n.pctls, "_est", i, ".pdf", sep = ""), 
             width = 15, height = 10)
    print(plots.comb.gg)
    dev.off()
  }
  
}
