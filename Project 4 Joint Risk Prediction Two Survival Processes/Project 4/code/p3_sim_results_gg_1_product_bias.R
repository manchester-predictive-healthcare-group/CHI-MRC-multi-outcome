###########
### This program produces the ggplot for calibraiton of the product method, panelled by level of residual correlation.
### Mean is superimposed over the 1000 calibration curves
###########

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Source packages
source("code/sim_function_load_packages.R")

#######################################################
#######################################################
### Create plots panelled by scearios, for a given sample size and DGM
### To highlight the impact of increasing residual correlation
#######################################################
#######################################################

### Must write a seperate loop for scenarios s1 and s2

### Create vector to loop through
n.vec <- c(1000, 2500, 5000)

### Loop through vectors and create plots
for (j in 1:length(n.vec)){
  
  print(n.vec[j])
  
  ### Create vectors to store individual plots
  gg.product.panel.s.list <- vector("list", 6)
  names(gg.product.panel.s.list) <- c("LN (no corr)", "LL (low corr)", "LH (high corr)", 
                                     "HN (no corr)", "HL (low corr)", "HH (high corr)")

  
  ### Nocorr scenario s1
  ## load workspace
  load(paste("data/sim_results_s1_n", n.vec[j], "v1000.RData", sep = ""))
  ### Load results functions
  source("code/sim_function_results.R")
  ### Source the most updated input parameters for plots
  source("code/sim_function_calibration_input_parameters_s1.R")
  
  ### Create ggplot
  ## Define limits
  x.lim <- c(0, 0.075)
  y.lim <- c(0, 0.15)
  ## Create plot
  gg.product.panel.s.list[[1]] <- create.individual.ggplots.product(analysed.data = res.DGM.nocorr, x.lim.in = x.lim, y.lim.in = y.lim, 
                                                                   font.size.in = font.size, linecolours.in = linecolours, 
                                                                   linewidths.in = linewidths, anal.methods.sim.inn = anal.methods.sim)
 
  ### Low corr scenario s1
  ## load workspace
  load(paste("data/sim_results_s1.1_n", n.vec[j], "v1000.RData", sep = ""))
  ### Load results functions
  source("code/sim_function_results.R")
  ### Source the most updated input parameters for plots
  source("code/sim_function_calibration_input_parameters_s1.1.R")

  ### Create ggplot
  ## Define limits
  x.lim <- c(0, 0.075)
  y.lim <- c(0, 0.15)
  ## Create plot
  gg.product.panel.s.list[[2]] <- create.individual.ggplots.product(analysed.data = res.DGM.clay, x.lim.in = x.lim, y.lim.in = y.lim, 
                                                                   font.size.in = font.size, linecolours.in = linecolours, 
                                                                   linewidths.in = linewidths, anal.methods.sim.inn = anal.methods.sim)
  
  ### High corr scenario s1
  ## load workspace
  load(paste("data/sim_results_s1.2_n", n.vec[j], "v1000.RData", sep = ""))
  ### Load results functions
  source("code/sim_function_results.R")
  ### Source the most updated input parameters for plots
  source("code/sim_function_calibration_input_parameters_s1.2.R")

  ### Create ggplot
  ## Define limits
  x.lim <- c(0, 0.075)
  y.lim <- c(0, 0.15)
  ## Create plot
  gg.product.panel.s.list[[3]] <- create.individual.ggplots.product(analysed.data = res.DGM.clay, x.lim.in = x.lim, y.lim.in = y.lim, 
                                                                   font.size.in = font.size, linecolours.in = linecolours, 
                                                                   linewidths.in = linewidths, anal.methods.sim.inn = anal.methods.sim)
  
  ### Nocorr scenario s2
  ## load workspace
  load(paste("data/sim_results_s2_n", n.vec[j], "v1000.RData", sep = ""))
  ### Load results functions
  source("code/sim_function_results.R")
  ### Source the most updated input parameters for plots
  source("code/sim_function_calibration_input_parameters_s2.R")

  ### Create ggplot
  ## Define limits
  x.lim <- c(0, 0.62)
  y.lim <- c(0, 0.8)
  ## Create plot
  gg.product.panel.s.list[[4]] <- create.individual.ggplots.product(analysed.data = res.DGM.nocorr, x.lim.in = x.lim, y.lim.in = y.lim, 
                                                                   font.size.in = font.size, linecolours.in = linecolours, 
                                                                   linewidths.in = linewidths, anal.methods.sim.inn = anal.methods.sim)
  
  ### Low corr scenario s2
  ## load workspace
  load(paste("data/sim_results_s2.1_n", n.vec[j], "v1000.RData", sep = ""))
  ### Load results functions
  source("code/sim_function_results.R")
  ### Source the most updated input parameters for plots
  source("code/sim_function_calibration_input_parameters_s2.1.R")

  ### Create ggplot
  ## Define limits
  x.lim <- c(0, 0.62)
  y.lim <- c(0, 0.8)
  ## Create plot
  gg.product.panel.s.list[[5]] <- create.individual.ggplots.product(analysed.data = res.DGM.clay, x.lim.in = x.lim, y.lim.in = y.lim, 
                                                                   font.size.in = font.size, linecolours.in = linecolours, 
                                                                   linewidths.in = linewidths, anal.methods.sim.inn = anal.methods.sim)
  
  ### High corr scenario s2
  ## load workspace
  load(paste("data/sim_results_s2.2_n", n.vec[j], "v1000.RData", sep = ""))
  ### Load results functions
  source("code/sim_function_results.R")
  ### Source the most updated input parameters for plots
  source("code/sim_function_calibration_input_parameters_s2.2.R")
  
  ### Create ggplot
  ## Define limits
  x.lim <- c(0, 0.62)
  y.lim <- c(0, 0.8)
  ## Create plot
  gg.product.panel.s.list[[6]] <- create.individual.ggplots.product(analysed.data = res.DGM.clay, x.lim.in = x.lim, y.lim.in = y.lim, 
                                                                   font.size.in = font.size, linecolours.in = linecolours, 
                                                                   linewidths.in = linewidths, anal.methods.sim.inn = anal.methods.sim)
  
  #######
  ## Combine into single ggplot
  #######
  ## Create grid labels based off this vector also
  grid.labels.in  <- names(gg.product.panel.s.list)
  
  gg.product.panel.s.comb <- ggarrange(gg.product.panel.s.list[[1]], gg.product.panel.s.list[[2]], gg.product.panel.s.list[[3]], 
                                      gg.product.panel.s.list[[4]], gg.product.panel.s.list[[5]], gg.product.panel.s.list[[6]],
                                      labels = grid.labels.in, nrow = 2, ncol = 3, 
                                      font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in, 
                                      common.legend = TRUE, legend = "right")
  

  #######
  ### Save combined plots
  #######
  CairoTIFF(filename = paste("figures/gg_1.summary.product.", n.vec[j], ".tiff", sep = ""),
            dpi = 1000, width = 10, height = 4.62, units = "in")
  plot(gg.product.panel.s.comb)
  dev.off()
  
  CairoPNG(filename = paste("figures/gg_1.summary.product.", n.vec[j], ".png", sep = ""),
            dpi = 1200, width = 10, height = 4.62, units = "in")
  plot(gg.product.panel.s.comb)
  dev.off()
  
  CairoPDF(file = paste("figures/gg_1.summary.product.", n.vec[j], ".pdf", sep = ""),
            width = 10, height = 4.62)
  plot(gg.product.panel.s.comb)
  dev.off()
  
  ## clear workspace
  rm(gg.median.panel.s.list, gg.p.range.panel.s.list, gg.median.panel.s.comb, gg.p.range.panel.s.comb)
}
