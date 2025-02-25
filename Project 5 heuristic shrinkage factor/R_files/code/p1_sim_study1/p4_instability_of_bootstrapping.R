###
### This program will highlight that bootstrapping is not a panacea, and that while bootstrapping works well on average,
### in any one given scenario, it may not.
###

### Set working directory
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_8.2/")

### Source functions
R.func.sources <- list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load functions
source("R/sim_functions.R")

### Call libraries
library(ggplot2)


### Write a function to produce instability plot
create_instability_plot <- function(n, correlated = TRUE){
# n <- 1000
# correlated <- TRUE
# metadata[248,]
  ### Read in data
  metadata <- readRDS("data/heuristic_shrinkage_sim_metadata.rds")
  metadata <- subset(metadata, C.pop != 0 & corr == correlated & cat == FALSE & P.meas == 10)
  
  ### Identify some scenarios where we bootstrapping does well
  metadata <- subset(metadata, abs(S.boot.diff.mean) < 0.01)
  
  ### Get set and iter
  set_iter <- c(metadata$set[which.min(abs(metadata$n - n))], metadata$sim.iter[which.min(abs(metadata$n - n))])

  ### Assign
  set <- set_iter[1]
  sim.iter <- set_iter[2]
  
  ### Load
  if (correlated == TRUE){
    load(paste("data/run_sim_corr", set, ".RData", sep = ""))
  } else {
    load(paste("data/run_sim_nocorr", set, ".RData", sep = ""))
  }

  ### Get range of S.boot values
  plotdat <- data.frame("S.boot" = output.S.boot[[sim.iter]], "S.pop" = output.S.pop[[sim.iter]])

  ### Plot them
  myplot <- ggplot(data = plotdat) + 
    geom_point(aes(x = S.boot, y = S.pop), alpha = 0.25) + 
    geom_abline(intercept = 0, slope = 1) + 
    coord_fixed(ratio = 1, 
                xlim = c(min(plotdat$S.boot, plotdat$S.pop), max(plotdat$S.boot, plotdat$S.pop)),
                ylim = c(min(plotdat$S.boot, plotdat$S.pop), max(plotdat$S.boot, plotdat$S.pop))) + 
    ggtitle(paste("N = ", input.data[[sim.iter]]$nreq)) +
    xlab(expression(S[boot])) + ylab(expression(S[pop]))
  
 
  ### Create output
  output <- list("plot" = myplot, "input" = input.data[[sim.iter]])
  return(output)
}

###
### Function to create a combined plot
###
create_combined_plot <- function(correlated){
  
  ### Create plots
  plot250 <- create_instability_plot(250, correlated)
  plot500 <- create_instability_plot(500, correlated)
  plot750 <- create_instability_plot(750, correlated)
  plot1000 <- create_instability_plot(1000, correlated)
  plot2000 <- create_instability_plot(2000, correlated)
  plot4000 <- create_instability_plot(4000, correlated)
  
  ### Combine
  plotlist <- list(plot250$plot, plot500$plot, plot750$plot, plot1000$plot, plot2000$plot, plot4000$plot)
  plots.comb <- ggpubr::ggarrange(plotlist = plotlist, nrow = 2, ncol = 3)
  
  # Save to disk
  Cairo::CairoPNG(paste("figures/gg.instability.plot.corr", as.numeric(correlated), ".png", sep = ""),
                  width = 10.5, height = 7, unit = "in", dpi = 300)
  plot(plots.comb)
  dev.off()
}

create_combined_plot(TRUE)
create_combined_plot(FALSE)