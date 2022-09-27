###
### This program will produce plots for the small sample analysis
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.4/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions_results.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")

### Run through DGMs/scenarios and create plots
for (K.sim in c(3, 5)){
  for (DGM.in in c("mult", "seqlog")){
    for (seed.coef.sim in c(1, 2, 3)){
      
      print(paste("K.sim = ", K.sim, ", DGM.in = ", DGM.in, ", seed.coef.sim = ", seed.coef.sim, sep = ""))
      
      ### Create plot data for median
      plot.data.p50 <- create.small.sample.ggplot.data.specific(summary.metrics = c("p50"), N.vec = c(100, 250, 500), scen.vec = c(1, 2, 3, 4), 
         
                                                                DGM.in = DGM.in, K.sim = K.sim, seed.coef.sim = seed.coef.sim)
      ###
      ### Create ggplot
      ###
      plot.p50 <- create.small.sample.ggplots(data.in = plot.data.p50, 
                                               x.lim.in  = c(0,1),
                                               y.lim.in = c(0,1),
                                               xlab.in = "Predicted risk",
                                               ylab.in = "Observed risk", 
                                               font.size.in = 10)
      
      ### Save ggplot
      ggsave(paste("figures/gg_small_samp_p50_DGM", DGM.in, "_K", K.sim, "_P", P.sim, "_seed", seed.coef.sim, ".png", sep =""), 
             plot.p50, dpi = 300)
      print("SAVED P50")
      
      
      ###
      ### Create plot data for mean
      ###
      plot.data.mean <- create.small.sample.ggplot.data.specific(summary.metrics = c("mean"), N.vec = c(100, 250, 500), scen.vec = c(1, 2, 3, 4), 
                                                                DGM.in = DGM.in, K.sim = K.sim, seed.coef.sim = seed.coef.sim)
      
      ### Create ggplot
      plot.mean <- create.small.sample.ggplots(data.in = plot.data.mean, 
                                              x.lim.in  = c(0,1),
                                              y.lim.in = c(0,1),
                                              xlab.in = "Predicted risk",
                                              ylab.in = "Observed risk", 
                                              font.size.in = 10)
      
      ### Save ggplot
      ggsave(paste("figures/gg_small_samp_mean_DGM", DGM.in, "_K", K.sim, "_P", P.sim, "_seed", seed.coef.sim, ".png", sep =""), 
             plot.mean, dpi = 300)
      print("SAVED MEAN")
      
      
      ###
      ### Create plot data for range
      ###
      plot.data.range <- create.small.sample.ggplot.data.specific(summary.metrics = c("range"), N.vec = c(100, 250, 500), scen.vec = c(1, 2, 3, 4), 
                                                                DGM.in = DGM.in, K.sim = K.sim, seed.coef.sim = seed.coef.sim)
      
      ### Create ggplot
      plot.range <- create.small.sample.ggplots.range(data.in = plot.data.range, 
                                                      x.lim.in  = c(0,1),
                                                      y.lim.in = c(0,1),
                                              xlab.in = "Predicted risk",
                                              ylab.in = "5th - 95th percentile range in observed risk", 
                                              font.size.in = 10)
      
      ### Save ggplot
      ggsave(paste("figures/gg_small_samp_range_DGM", DGM.in, "_K", K.sim, "_P", P.sim, "_seed", seed.coef.sim, ".png", sep =""), 
             plot.range, dpi = 300)
      print("SAVED RANGE")
      
    }
  }
}


