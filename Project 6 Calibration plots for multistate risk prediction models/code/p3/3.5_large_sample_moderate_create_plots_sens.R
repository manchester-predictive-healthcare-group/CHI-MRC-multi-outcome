###
### Create plots (sensitivity analyses).
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions.R")
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Load the true calibration data
true.calib <- readRDS(paste("data/sim/true.calib.lineplot.rcs.nk", 5, ".rds", sep = ""))

###
### For the sensitivity analysis plots, we plot the BLR-IPCW and MLR-IPCW analyses, with weights estimated
### from the data, weights taken from the DGM, and no weight
###

###
### Write function to add truedata to the calibration plot data for the calibration curves
###
add_true_data_curve <- function(calibdata, truedata){
  
  ### Add true risks to plotdata
  for (state in 1:5){
    
    ### Reduce truedata to landmarked individuals
    truedata[[state]] <- subset(truedata[[state]], id %in% calibdata[[state]]$id) 
    if (sum(calibdata[[state]]$pred != truedata[[state]]$pred)){
      stop("Mismatch between predicted risks to genereate true calibration curve, and predicted risks used to estimate calibration curve")
    }
    
    ### Assign true calibration curve values to calibdata, so we can plot together with the estimated calibration curve
    calibdata[[state]]$true <- truedata[[state]]$obs
    
    ### Reduce to 2500 observations to reduce size of plot file
    calibdata[[state]] <- calibdata[[state]][1:2500, ]
    
  }
  
  ### Return calibdata
  return(calibdata)
  
}

###
### Write function to add truedata to the calibraiton plot data for the calibration scatter plots
###
add_true_data_scatter <- function(calibdata, truedata){
  
  ### Add true risks to plotdata
  for (state in 1:5){
    
    ### Assign true calibration curve values to calibdata, so we can plot together with the estimated calibration curve
    calibdata[[state]]$true <- truedata[,paste("p.true", state, sep = "")]
    
    ### Reduce to 5,000 observations, as plotting ~ 120,000 observations is unnecesary
    calibdata[[state]] <- calibdata[[state]][1:2500, ]
    
  }
  
  ### Return calibdata
  return(calibdata)
  
}

#################################################
### Sensitivity analysis for the IPCW methods ###
#################################################

### Write a function to produce calibration plots for a given scenario (scen) and set of predicted transition probabilities (i)
### Calibration plots will be using no weights, estimating weights from data (main analysis) and using weights from the DGM
create_plots_sens_IPCW <- function(scen.tp, scen){

  ### Load the calibration plot data
  load(paste("data/sim/large_sample_analysis_moderate_blr_", scen, ".RData", sep = ""))
  load(paste("data/sim/large_sample_analysis_moderate_mlr_", scen, ".RData", sep = ""))
  load(paste("data/sim/large_sample_analysis_moderate_sens_", scen, ".RData", sep = ""))
  source("code/z_functions.R", local = TRUE)
  
  ### Create an object to store the three plots
  calib.plots <- vector("list", 3)
  
  ###
  ### Start with plots for BLR-IPCW
  ###
  
  ### Assign true calibration curve data
  truedata.curve <- true.calib[[scen.tp]]
  
  ###
  ### Create plots estimating weights from data (same as main analysis)

  ### Extract plot data and true curve data
  calibdata <- add_true_data_curve(calibdata = calib.moderate.blr[[scen.tp]][["plotdata"]],
                                   truedata = truedata.curve)

  ### Create calibration plots
  calib.plots[[1]] <- 
    plot_calib_msm(calibdata, 
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1, 
                   inclu.legend = FALSE,
                   axis.titles.x = 0,
                   axis.titles.y = 1, 
                   axis.titles.text.y = "Weights estimated from data\nObserved risk",
                   size = 16,
                   CI = FALSE)
  
  ###
  ### Create plots using DGM specified weights

  ### Extract plot data and true curve data
  calibdata <- add_true_data_curve(calibdata = calib.moderate.blr.DGMw[[scen.tp]][["plotdata"]],
                                   truedata = truedata.curve)

  ### Create calibration plots
  calib.plots[[2]] <- 
    plot_calib_msm(calibdata, 
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1, 
                   inclu.legend = FALSE,
                   inclu.title = FALSE,
                   axis.titles.x = 0,
                   axis.titles.y = 1, 
                   axis.titles.text.y = "Weights taken from DGM\nObserved risk",
                   size = 16,
                   CI = FALSE)
  
  ###
  ### Create plots using no weights

  ### Extract plot data and true curve data
  calibdata <- add_true_data_curve(calibdata = calib.moderate.blr.NOw[[scen.tp]][["plotdata"]],
                                   truedata = truedata.curve)
  
  ### Create calibration plots
  calib.plots[[3]] <- 
    plot_calib_msm(calibdata, 
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1,
                   inclu.legend = TRUE,
                   legend.seperate = TRUE,
                   inclu.title = FALSE,
                   axis.titles.x = c(1,2,3,4,5),
                   axis.titles.y = 1, 
                   axis.titles.text.y = "No weights\nObserved risk",
                   size = 16,
                   CI = FALSE)

  ### Combine into a single plot for the Figure
  plots.fig <- gridExtra::arrangeGrob(calib.plots[[1]], calib.plots[[2]], calib.plots[[3]][["plots"]], calib.plots[[3]][["legend"]], 
                                      nrow = 4, heights = c(1, 1, 1, 0.05))
  
  CairoPNG(paste("figures/sens_large_sample_moderate_", scen, "_tp", scen.tp, "_blr.png", sep = ""), width = 20, height = 12, unit = "in", res = 150)
  grid::grid.draw(plots.fig)
  dev.off()
  CairoPDF(paste("figures/sens_large_sample_moderate_", scen,  "_tp", scen.tp, "_blr.pdf", sep = ""), width = 20, height = 12)
  grid::grid.draw(plots.fig)
  dev.off()
  
  
  ###
  ### Now do plots for MLR-IPCW
  ###
  
  ### Assign true calibration curve and scatter data
  truedata.scatter <- 
    subset(data.raw, id %in% calib.moderate.mlr[[1]][["plotdata"]][[1]]$id) |> 
    dplyr::select(paste("p.true", 1:5, sep = ""))
  
  
  ###
  ### Create plots estimating weights from data (same as main analysis)

  ### Extract plot data and true curve data
  calibdata <- add_true_data_scatter(calibdata = calib.moderate.mlr[[scen.tp]][["plotdata"]],
                                   truedata = truedata.scatter)
  
  ### Create calibration plots
  calib.plots[[1]] <- 
    plot_calib_mlr(calibdata, 
                   transparency.plot = 0.1,
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1, 
                   inclu.legend = FALSE,
                   axis.titles.x = 0,
                   axis.titles.y = 1, 
                   axis.titles.text.y = "Weights estimated from data\nObserved risk",
                   size = 16)
  
  ###
  ### Create plots using DGM specified weights

  ### Extract plot data and true curve data
  calibdata <- add_true_data_scatter(calibdata = calib.moderate.mlr.DGMw[[scen.tp]][["plotdata"]],
                                   truedata = truedata.scatter)
  
  ### Create calibration plots
  calib.plots[[2]] <- 
    plot_calib_mlr(calibdata, 
                   transparency.plot = 0.1,
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1, 
                   inclu.legend = FALSE,
                   inclu.title = FALSE,
                   axis.titles.x = 0,
                   axis.titles.y = 1, 
                   axis.titles.text.y = "Weights taken from DGM\nObserved risk",
                   size = 16)
  
  ###
  ### Create plots using no weights
  
  ### Extract plot data and true curve data
  calibdata <- add_true_data_scatter(calibdata = calib.moderate.mlr.NOw[[scen.tp]][["plotdata"]],
                                   truedata = truedata.scatter)
  
  ### Create calibration plots
  calib.plots[[3]] <- 
    plot_calib_mlr(calibdata, 
                   transparency.plot = 0.1,
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1,
                   inclu.legend = TRUE,
                   legend.seperate = TRUE,
                   inclu.title = FALSE,
                   axis.titles.x = c(1,2,3,4,5),
                   axis.titles.y = 1, 
                   axis.titles.text.y = "No weights\nObserved risk",
                   size = 16)
  
  ### Combine into a single plot for the Figure
  plots.fig <- gridExtra::arrangeGrob(calib.plots[[1]], calib.plots[[2]], calib.plots[[3]][["plots"]], calib.plots[[3]][["legend"]], 
                                      nrow = 4, heights = c(1, 1, 1, 0.05))
  
  CairoPNG(paste("figures/sens_large_sample_moderate_", scen, "_tp", scen.tp, "_mlr.png", sep = ""), width = 20, height = 12, unit = "in", res = 150)
  grid::grid.draw(plots.fig)
  dev.off()
  CairoPDF(paste("figures/sens_large_sample_moderate_", scen,  "_tp", scen.tp, "_mlr.pdf", sep = ""), width = 20, height = 12)
  grid::grid.draw(plots.fig)
  dev.off()
  
}

###
### Generate the plots 
###
lapply(c(1,2,3), create_plots_sens_IPCW, scen = "C1")
lapply(c(1,2,3), create_plots_sens_IPCW, scen = "C2")
lapply(c(1,2,3), create_plots_sens_IPCW, scen = "C3")


########################################################
### Sensitivity analysis for the pseudo-value method ###
########################################################

### Write a function to produce calibration plots for a given scenario (scen) and set of predicted transition probabilities (i)
### Calibration plots will be made grouping individuals by predicted risks, and not grouping individuals

### Write a function to produce calibration plots for a given scenario (scen) and set of predicted transition probabilities (i)
### Calibration plots will be using no weights, estimating weights from data (main analysis) and using weights from the DGM
create_plots_sens_PV <- function(scen.tp, scen){
  
  ### Load the calibration plot data
  load(paste("data/sim/large_sample_analysis_moderate_pv_", scen, ".RData", sep = ""))
  load(paste("data/sim/large_sample_analysis_moderate_sens_pv_", scen, ".RData", sep = ""))
  source("code/z_functions.R", local = TRUE)
  
  ### Create an object to store the three plots
  calib.plots <- vector("list", 3)
  
  ###
  ### Start with plots for BLR-IPCW
  ###
  
  ### Assign true calibration curve data
  truedata.curve <- true.calib[[scen.tp]]
  
  ###
  ### Create plots where individuals were grouped by predicted risk before calculating pseudo-values (same as main analysis)
  
  ### Extract plot data and true curve data
  calibdata <- add_true_data_curve(calibdata = calib.moderate.pv[[scen.tp]][["plotdata"]],
                                   truedata = truedata.curve)
  
  ### Create calibration plots
  calib.plots[[1]] <- 
    plot_calib_msm(calibdata, 
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1, 
                   inclu.legend = FALSE,
                   axis.titles.x = 0,
                   axis.titles.y = 1, 
                   axis.titles.text.y = "Grouped by\npredicted risk",
                   size = 16,
                   CI = TRUE)
  
  ###
  ### Create plots where individuals are not grouped by predicted risk before calculating pseudo-values
  
  ### Extract plot data and true curve data
  calibdata <- add_true_data_curve(calibdata = calib.moderate.pv.sens[[scen.tp]][["plotdata"]],
                                   truedata = truedata.curve)
  
  ### Create calibration plots
  calib.plots[[2]] <- 
    plot_calib_msm(calibdata, 
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1, 
                   inclu.legend = TRUE,
                   legend.seperate = TRUE,
                   inclu.title = FALSE,
                   axis.titles.x = c(1,2,3,4,5),
                   axis.titles.y = 1, 
                   axis.titles.text.y = "Not grouped by\npredicted risk",
                   size = 16,
                   CI = TRUE)
  
  ### Combine into a single plot for the Figure
  plots.fig <- gridExtra::arrangeGrob(calib.plots[[1]], calib.plots[[2]][["plots"]], calib.plots[[2]][["legend"]], 
                                      nrow = 3, heights = c(1, 1, 0.05))
  
  CairoPNG(paste("figures/sens_large_sample_moderate_", scen, "_tp", scen.tp, "_pv.png", sep = ""), width = 20, height = 8, unit = "in", res = 150)
  grid::grid.draw(plots.fig)
  dev.off()
  CairoPDF(paste("figures/sens_large_sample_moderate_", scen, "_tp", scen.tp, "_pv.pdf", sep = ""), width = 20, height = 8)
  grid::grid.draw(plots.fig)
  dev.off()
  
}

###
### Generate the plots 
###
lapply(c(1,2,3), create_plots_sens_PV, scen = "C1")
lapply(c(1,2,3), create_plots_sens_PV, scen = "C2")
lapply(c(1,2,3), create_plots_sens_PV, scen = "C3")