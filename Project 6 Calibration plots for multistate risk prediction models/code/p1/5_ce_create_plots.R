###
### This program will create calibration plots
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### Write a function to produce plots for a given sample size
create_plots_ce <- function(n){
  
  ### Save image
  load(paste("data/ce/ce_assess_calib_N", n, "_npctls", 20, ".RData", sep = ""))
  
  ### Reduce calib.model.mlr to a smaller number of data points
  for (state in 1:6){
    calib.mod.pv[["plotdata"]][[state]] <- calib.mod.pv[["plotdata"]][[state]][1:2500, ]
    calib.mod.blr[["plotdata"]][[state]] <- calib.mod.blr[["plotdata"]][[state]][1:2500, ]
    calib.mod.mlr[["plotdata"]][[state]] <- calib.mod.mlr[["plotdata"]][[state]][1:5000, ]
  }
  
  ### Create plot objects
  plot.pv <- plot(calib.mod.pv, nrow = 1, ncol = 6, size.text = 16, size.line = 0.5,
                  marg.density = TRUE,
                  axis.titles.x = 0, axis.titles.y = 1, axis.titles.text.y = "Pseudo-value\nObserved risk", 
                  legend.include = TRUE, legend.seperate = TRUE, 
                  titles = c("Healthy", "CVD", "T2D", "CKD", "Multimorbidity", "Death"))
  plot.blr <- plot(calib.mod.blr, nrow = 1, ncol = 6, size.text = 16, size.line = 0.5,
                   marg.density = TRUE,
                   axis.titles.x = 0, axis.titles.y = 1, axis.titles.text.y = "BLR-IPCW\nObserved risk", 
                   legend.include = FALSE, titles.include = FALSE)
  plot.mlr <- plot(calib.mod.mlr, nrow = 1, ncol = 6, transparency.plot = 0.05, size.text = 16, size.point = 0.75,
                   marg.density = TRUE,
                   axis.titles.x = 0, axis.titles.y = 1, axis.titles.text.y = "MLR-IPCW\nObserved risk", 
                   titles.include = FALSE)
  
  ### AJ
  ### Some extra manipulation required for AJ given plots were not produced by calibmsm
  calib.plots.aj[[1]] <- calib.plots.aj[[1]] + ylab("AJ\nObserved risk") + ggtitle(NULL) + theme(text = element_text(size = 16)) +
    xlim(c(0,0.8)) + ylim(c(0,0.8))
  
  for (i in 2:6){
    calib.plots.aj[[i]] <- calib.plots.aj[[i]] + ylab(NULL) + ggtitle(NULL) + theme(text = element_text(size = 16)) +
      xlim(list(c(0,0.15), c(0,0.2), c(0,0.2), c(0,0.3), c(0,0.8))[[(i-1)]]) +
      ylim(list(c(0,0.15), c(0,0.2), c(0,0.2), c(0,0.3), c(0,0.8))[[(i-1)]])
  }
  
  plot.aj <- gridExtra::arrangeGrob(calib.plots.aj[[1]],  calib.plots.aj[[2]],  calib.plots.aj[[3]],
                                           calib.plots.aj[[4]],  calib.plots.aj[[5]],  calib.plots.aj[[6]],
                                           ncol = 6)
  
  ### Combine into a single plot for the Figure
  plots.fig <- gridExtra::arrangeGrob(plot.pv[["plots"]],
                                      plot.blr,
                                      plot.mlr,
                                      plot.aj,
                                      plot.pv[["legend"]], nrow = 5, heights = c(1, 1, 1, 1, 0.05))
  
  CairoPNG(paste("figures/ce_N", n, ".png", sep = ""), width = 20, height = 13.3, unit = "in", res = 300)
  grid::grid.draw(plots.fig)
  dev.off()
  CairoPDF(paste("figures/ce_N", n, ".pdf", sep = ""), width = 20, height = 13.3)
  grid::grid.draw(plots.fig)
  dev.off()
  
}

create_plots_ce(5000)
create_plots_ce(100000)
