###
### Create plots.
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Load the true calibration data
true.calib <- readRDS(paste("data/sim/true.calib.lineplot.rcs.nk", 5, ".rds", sep = ""))

###################
### Scenario C1 ###
###################

### Set scenario
scen <- "C1"

### Load the calibration plot data
load(paste("data/sim/large_sample_analysis_moderate_blr_", scen, ".RData", sep = ""))
load(paste("data/sim/large_sample_analysis_moderate_mlr_", scen, ".RData", sep = ""))
load(paste("data/sim/large_sample_analysis_moderate_pv_", scen, ".RData", sep = ""))
source("code/z_functions.R")

###
### Cycle through predicted risks and create plots for BLR-IPCW
###

### Create an object to store the three plots
calib.plots <- vector("list", 3)

### Cycle through the predicted risks and create plots
for (scen.tp in 1:3){

  ### Extract plot data and true curve data
  calibdata <- calib.moderate.blr[[scen.tp]][["plotdata"]]
  truedata <- true.calib[[scen.tp]]
  
  ### Add true risks to plotdata
  for (state in 1:5){
    
    ### Reduce truedata to landmarked individuals
    truedata[[state]] <- subset(truedata[[state]], id %in% calibdata[[state]]$id) 
    if (sum(calibdata[[state]]$pred != truedata[[state]]$pred)){
      stop("Mismatch between predicted risks to genereate true calibration curve, and predicted risks used to estimate calibration curve")
    }
    
    ### Assign true calibration curve values to calibdata, so we can plot together with the estimated calibration curve
    calibdata[[state]]$true <- truedata[[state]]$obs
    
    ### Reduce to 5000 observations to reduce size of plot file
    calibdata[[state]] <- calibdata[[state]][1:2500, ]
  }
  
  if (scen.tp == 1){
    calib.plots[[scen.tp]] <- 
      plot_calib_msm(calibdata, 
                     marg.density = TRUE, 
                     marg.density.size = 8, 
                     ncol = 5, 
                     nrow = 1, 
                     inclu.legend = FALSE,
                     axis.titles.x = 0,
                     axis.titles.y = 1, 
                     axis.titles.text.y = "Perfectly calibrated\nObserved risk",
                     size = 16)
  } else if (scen.tp == 2){
    calib.plots[[scen.tp]] <- 
      plot_calib_msm(calibdata, 
                     marg.density = TRUE, 
                     marg.density.size = 8, 
                     ncol = 5, 
                     nrow = 1, 
                     inclu.legend = FALSE,
                     inclu.title = FALSE,
                     axis.titles.x = 0,
                     axis.titles.y = 1, 
                     axis.titles.text.y = "Miscalibration 1\nObserved risk",
                     size = 16)
  } else if (scen.tp == 3){
    calib.plots[[scen.tp]] <- 
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
                     axis.titles.text.y = "Miscalibration 2\nObserved risk",
                     size = 16)
  }
  
}

### Combine into a single plot for the Figure
plots.fig <- gridExtra::arrangeGrob(calib.plots[[1]], calib.plots[[2]], calib.plots[[3]][["plots"]], calib.plots[[3]][["legend"]], 
                                      nrow = 4, heights = c(1, 1, 1, 0.05))

CairoPNG(paste("figures/large_sample_moderate_", scen, "_blr.png", sep = ""), width = 20, height = 12, unit = "in", res = 150)
grid::grid.draw(plots.fig)
dev.off()
CairoPDF(paste("figures/large_sample_moderate_", scen, "_blr.pdf", sep = ""), width = 20, height = 12)
grid::grid.draw(plots.fig)
dev.off()


###
### Cycle through predicted risks and create plots for PV
###

### Create an object to store the three plots
calib.plots <- vector("list", 3)

### Cycle through the predicted risks and create plots
for (scen.tp in 1:3){
  
  ### Extract plot data and true curve data
  calibdata <- calib.moderate.pv[[scen.tp]][["plotdata"]]
  truedata <- true.calib[[scen.tp]]
  
  ### Add true risks to plotdata
  for (state in 1:5){
    
    ### Reduce truedata to landmarked individuals
    truedata[[state]] <- subset(truedata[[state]], id %in% calibdata[[state]]$id) 
    if (sum(calibdata[[state]]$pred != truedata[[state]]$pred)){
      stop("Mismatch between predicted risks to genereate true calibration curve, and predicted risks used to estimate calibration curve")
    }
    
    ### Assign true calibration curve values to calibdata, so we can plot together with the estimated calibration curve
    calibdata[[state]]$true <- truedata[[state]]$obs
    
    ### Reduce to 5000 observations to reduce size of plot file
    calibdata[[state]] <- calibdata[[state]][1:2500, ]
  }
  
  if (scen.tp == 1){
    calib.plots[[scen.tp]] <- 
      plot_calib_msm(calibdata, 
                     marg.density = TRUE, 
                     marg.density.size = 8, 
                     ncol = 5, 
                     nrow = 1, 
                     inclu.legend = FALSE,
                     axis.titles.x = 0,
                     axis.titles.y = 1, 
                     axis.titles.text.y = "Perfectly calibrated\nObserved risk",
                     size = 16)
  } else if (scen.tp == 2){
    calib.plots[[scen.tp]] <- 
      plot_calib_msm(calibdata, 
                     marg.density = TRUE, 
                     marg.density.size = 8, 
                     ncol = 5, 
                     nrow = 1, 
                     inclu.legend = FALSE,
                     inclu.title = FALSE,
                     axis.titles.x = 0,
                     axis.titles.y = 1, 
                     axis.titles.text.y = "Miscalibration 1\nObserved risk",
                     size = 16)
  } else if (scen.tp == 3){
    calib.plots[[scen.tp]] <- 
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
                     axis.titles.text.y = "Miscalibration 2\nObserved risk",
                     size = 16)
  }
  
}

### Combine into a single plot for the Figure
plots.fig <- gridExtra::arrangeGrob(calib.plots[[1]], calib.plots[[2]], calib.plots[[3]][["plots"]], calib.plots[[3]][["legend"]], 
                                     nrow = 4, heights = c(1, 1, 1, 0.05))

CairoPNG(paste("figures/large_sample_moderate_", scen, "_pv.png", sep = ""), width = 20, height = 12, unit = "in", res = 150)
grid::grid.draw(plots.fig)
dev.off()
CairoPDF(paste("figures/large_sample_moderate_", scen, "_pv.pdf", sep = ""), width = 20, height = 12)
grid::grid.draw(plots.fig)
dev.off()


###
### Cycle through predicted risks and create plots for MLR-IPCW
###

### Create an object to store the three plots
calib.plots <- vector("list", 3)

### Note the true y-axis data (observed risk) is now the same for cal1, cal2 and cal3, 
### as we do not need to obtain observed risk as a function of predicted risk

### First Reduce truedata to landmarked individuals
truedata <- subset(data.raw, id %in% calib.moderate.mlr[[1]][["plotdata"]][[1]]$id) |> 
  dplyr::select(paste("p.true", 1:5, sep = ""))

### Cycle through predicted risks
for (scen.tp in 1:3){
  
  ### Extract plot data
  calibdata <- calib.moderate.mlr[[scen.tp]][["plotdata"]]
  
  ### Add true risks to plotdata
  for (state in 1:5){
    
    ### Assign true calibration curve values to calibdata, so we can plot together with the estimated calibration curve
    calibdata[[state]]$true <- truedata[,paste("p.true", state, sep = "")]
    
    ### Reduce to 5,000 observations, as plotting ~ 120,000 observations is unnecesary
    calibdata[[state]] <- calibdata[[state]][1:2500, ]
    
  }
  
  if (scen.tp == 1){
    calib.plots[[scen.tp]] <- 
      plot_calib_mlr(calibdata, 
                     transparency.plot = 0.1,
                     marg.density = TRUE, 
                     marg.density.size = 8, 
                     ncol = 5, 
                     nrow = 1, 
                     inclu.legend = FALSE,
                     axis.titles.x = 0,
                     axis.titles.y = 1, 
                     axis.titles.text.y = "Perfectly calibrated\nObserved risk",
                     size = 16)
  } else if (scen.tp == 2){
    calib.plots[[scen.tp]] <- 
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
                     axis.titles.text.y = "Miscalibration 1\nObserved risk",
                     size = 16)
  } else if (scen.tp == 3){
    calib.plots[[scen.tp]] <- 
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
                     axis.titles.text.y = "Miscalibration 2\nObserved risk",
                     size = 16)
  }

}

### Combine into a single plot for the Figure
plots.fig <- gridExtra::arrangeGrob(calib.plots[[1]], calib.plots[[2]], calib.plots[[3]][["plots"]], calib.plots[[3]][["legend"]], 
                                      nrow = 4, heights = c(1, 1, 1, 0.05))

CairoPNG(paste("figures/large_sample_moderate_", scen, "_mlr.png", sep = ""), width = 20, height = 12, unit = "in", res = 150)
grid::grid.draw(plots.fig)
dev.off()
CairoPDF(paste("figures/large_sample_moderate_", scen, "_mlr.pdf", sep = ""), width = 20, height = 12)
grid::grid.draw(plots.fig)
dev.off()



###########################
### Scenarios C2 and C3 ###
###########################

### For these scenarios, we want to combine plots from each method for the same set of predicted risks
### (for C1, we plotted the same method with different sets of predicted risks)

###
### Write function to add truedata to the calibraiton plot data for the calibration curves
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
    
    ### Reduce to 5000 observations to reduce size of plot file
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

###
### Write a function to cycle through calibration methods for a given set of predicted transition probabilities
###
create_plots_C2C3 <- function(scen.tp, scen){

  ### Load the calibration plot data
  load(paste("data/sim/large_sample_analysis_moderate_blr_", scen, ".RData", sep = ""))
  load(paste("data/sim/large_sample_analysis_moderate_mlr_", scen, ".RData", sep = ""))
  load(paste("data/sim/large_sample_analysis_moderate_pv_", scen, ".RData", sep = ""))
  source("code/z_functions.R", local = TRUE)
  
  ### Create an object to store the three plots
  calib.plots <- vector("list", 3)
  
  ### Assign true calibration curve and scatter data
  truedata.curve <- true.calib[[scen.tp]]
  truedata.scatter <- 
    subset(data.raw, id %in% calib.moderate.mlr[[1]][["plotdata"]][[1]]$id) |> 
    dplyr::select(paste("p.true", 1:5, sep = ""))
  
  
  ###
  ### Create plot for BLR-IPCW
  
  ### Extract plot data and true curve data
  calibdata <- add_true_data_curve(calibdata = calib.moderate.blr[[scen.tp]][["plotdata"]],
                                   truedata = truedata.curve)
  
  ### Create calibration plots for BLR-IPCW
  calib.plots[[1]] <- 
    plot_calib_msm(calibdata, 
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1, 
                   inclu.legend = FALSE,
                   axis.titles.x = 0,
                   axis.titles.y = 1, 
                   axis.titles.text.y = "BLR-IPCW\nObserved risk",
                   size = 16)
  
  ###
  ### Create plot for PV
  
  ### Combine plot data and true curve data
  calibdata <- add_true_data_curve(calibdata = calib.moderate.pv[[scen.tp]][["plotdata"]],
                                   truedata = truedata.curve)
  
  ### Create calibration plots for PV
  calib.plots[[2]] <- 
    plot_calib_msm(calibdata, 
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1, 
                   inclu.legend = TRUE,
                   legend.title = "Calibration curves:",
                   axis.titles.x = 0,
                   axis.titles.y = 1, 
                   axis.titles.text.y = "Pseudo-value\nObserved risk",
                   size = 16)
  
  ###
  ### Create plot for MLR-IPCW
  
  ### Combine plot data and true curve data
  calibdata <- add_true_data_scatter(calibdata = calib.moderate.mlr[[scen.tp]][["plotdata"]],
                                     truedata = truedata.scatter)
  
  ### Create calibration plots for MLR-IPCW
  calib.plots[[3]] <- 
    plot_calib_mlr(calibdata, 
                   transparency.plot = 0.1,
                   marg.density = TRUE, 
                   marg.density.size = 8, 
                   ncol = 5, 
                   nrow = 1,
                   inclu.legend = TRUE,
                   legend.seperate = TRUE,
                   legend.title = "Calibration scatter plot:",
                   inclu.title = FALSE,
                   axis.titles.x = c(1,2,3,4,5),
                   axis.titles.y = 1, 
                   axis.titles.text.y = "MLR-IPCW\nObserved risk",
                   size = 16)
  
  ### Create a grob of the two legends
  legends.grob <- gridExtra::arrangeGrob(calib.plots[[2]][["legend"]], calib.plots[[3]][["legend"]], 
                                         nrow = 1)
  
  ### Combine into a single plot for the Figure
  plots.fig <- gridExtra::arrangeGrob(calib.plots[[1]], calib.plots[[2]][["plots"]], calib.plots[[3]][["plots"]], legends.grob, 
                                      nrow = 4, heights = c(1, 1, 1, 0.05))
  
  CairoPNG(paste("figures/large_sample_moderate_", scen, "_tp", scen.tp, ".png", sep = ""), width = 20, height = 12, unit = "in", res = 150)
  grid::grid.draw(plots.fig)
  dev.off()
  CairoPDF(paste("figures/large_sample_moderate_", scen,  "_tp", scen.tp, ".pdf", sep = ""), width = 20, height = 12)
  grid::grid.draw(plots.fig)
  dev.off()
  
}

###
### Generate the plots 
###
lapply(c(1,2,3), create_plots_C2C3, scen = "C2")
lapply(c(1,2,3), create_plots_C2C3, scen = "C3")