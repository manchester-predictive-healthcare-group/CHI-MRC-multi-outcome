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

### Write a function to extract set and iter for a given n
get_set_iter <- function(n, correlated = TRUE, P = NULL){
  
  ### Read in data
  metadata <- readRDS("data/heuristic_shrinkage_sim_metadata_ssd.rds")
  if (is.null(P)){
    metadata <- subset(metadata, corr == correlated & !is.na(S.pop.mean.pavlou))
    
  } else {
    metadata <- subset(metadata, corr == correlated & P.meas == P & !is.na(S.pop.mean.pavlou))
    
  }
  
  ### Get set and iter
  set_iter <- c(metadata$set[which.min(abs(as.numeric(metadata$nreq.pavlou) - n))], metadata$sim.iter[which.min(abs(as.numeric(metadata$nreq.pavlou) - n))])
  
  return(as.numeric(set_iter))
  
}

### Get values of set and iter for the plot
set_iter_list_true <- lapply(c(250, 500, 750, 1000, 2000, 4000), function(x) {get_set_iter(x, TRUE, P = 10)})
set_iter_list_false <- lapply(c(250, 500, 750, 1000, 2000, 4000), function(x) {get_set_iter(x, FALSE, P = 10)})

set_iter_list_true_noP <- lapply(c(250, 500, 750, 1000, 2000, 4000), function(x) {get_set_iter(x, TRUE)})
set_iter_list_false_noP <- lapply(c(250, 500, 750, 1000, 2000, 4000), function(x) {get_set_iter(x, FALSE)})


### Write a function to produce instability plot
create_instability_plot <- function(set_iter_list, correlated = TRUE, output.text = NULL){
  
  ### Create a data frame
  plotdat <- data.frame("scenario" = numeric(0), "N" = numeric(0), "S.pop"= numeric(0))

  ### Get plot data
  for (loopvar in 1:length(set_iter_list)){
    print(loopvar)
    ### Assign set_iter
    set_iter <- set_iter_list[[loopvar]]
    
    ### Load
    if (correlated == TRUE){
      load(paste("data/run_sim_validate_ssd_corr", set_iter[1], ".RData", sep = ""))
    } else {
      load(paste("data/run_sim_validate_ssd_nocorr", set_iter[1], ".RData", sep = ""))
    }
    
    ### Create data frame to add to plotdat
    addition <- data.frame("scenario" = loopvar, 
                           "N" = input.data[[set_iter[2]]]$nreq.pavlou, 
                           "S.pop" = output.pavlou.S.pop[[set_iter[2]]])
    plotdat <- rbind(plotdat, addition)
    
  }
  
  ### Create scenario variable
  plotdat$scenario <- as.factor(plotdat$scenario)
  levels(plotdat$scenario) <- paste("N = ", unique(plotdat$N))
  
  ### Plot them
  myplot <- ggplot(data = plotdat, aes(scenario, S.pop)) + 
    geom_violin(draw_quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) + 
    ylab(expression(S[pop])) + ylim(c(0.6, 1.3))
    
  # Save to disk
  Cairo::CairoPNG(paste("figures/gg.instability.ssd.S.pop.corr", as.numeric(correlated), output.text,".png", sep = ""),
                  width = 8, height = 6, unit = "in", dpi = 300)
  plot(myplot)
  dev.off()
  
}

### Apply function to create plot
create_instability_plot(set_iter_list_true, TRUE)
# create_instability_plot(set_iter_list_false, FALSE)