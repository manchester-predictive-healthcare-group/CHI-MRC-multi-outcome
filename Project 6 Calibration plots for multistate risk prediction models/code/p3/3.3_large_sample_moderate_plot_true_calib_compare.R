###
### Compare the true calibration curves
###

###
### This program will estimate the true calibration curve based off the true risks and the predicted risks
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_functions_true_transition_probs.R")
source("code/z_load_packages.R")

### True calibration is the same irrespective of the censoring that has been applied, 
### so just need to load data for one scenario
load(paste("data/sim/large_sample_prep_data_C1.RData", sep = ""))
source("code/z_functions.R")


rcs.nk3 <- readRDS(paste("data/sim/true.calib.lineplot.rcs.nk3.rds", sep = ""))
rcs.nk4 <- readRDS(paste("data/sim/true.calib.lineplot.rcs.nk4.rds", sep = ""))
rcs.nk5 <- readRDS(paste("data/sim/true.calib.lineplot.rcs.nk5.rds", sep = ""))
loess.span0.1 <- readRDS(paste("data/sim/true.calib.lineplot.loess.span0.1.rds", sep = ""))
loess.span0.1.nostat <- readRDS(paste("data/sim/true.calib.lineplot.loess.span0.1.nostat.rds", sep = ""))

### Combine plot data
plotdata.rcs.nk3 <- do.call("rbind", 
                            lapply(c(1,2,3), 
                                   function(y) {
                                     data.frame("tp" = paste("tp", y, sep = ""), 
                                                do.call("rbind", 
                                                        lapply(c(1,2,3,4,5), 
                                                               function(x) {data.frame("state" = x, 
                                                                                       "model" = "rcs.nk3", 
                                                                                       rcs.nk3[[y]][[x]][1:5000, ])})))}))

plotdata.rcs.nk4 <- do.call("rbind", 
                            lapply(c(1,2,3), 
                                   function(y) {
                                     data.frame("tp" = paste("tp", y, sep = ""), 
                                                do.call("rbind", 
                                                        lapply(c(1,2,3,4,5), 
                                                               function(x) {data.frame("state" = x, 
                                                                                       "model" = "rcs.nk4", 
                                                                                       rcs.nk4[[y]][[x]][1:5000, ])})))}))


plotdata.rcs.nk5 <- do.call("rbind", 
                            lapply(c(1,2,3), 
                                   function(y) {
                                     data.frame("tp" = paste("tp", y, sep = ""), 
                                                do.call("rbind", 
                                                        lapply(c(1,2,3,4,5), 
                                                               function(x) {data.frame("state" = x, 
                                                                                       "model" = "rcs.nk5", 
                                                                                       rcs.nk5[[y]][[x]][1:5000, ])})))}))

plotdata.loess.span0.1 <- do.call("rbind", 
                            lapply(c(1,2,3), 
                                   function(y) {
                                     data.frame("tp" = paste("tp", y, sep = ""), 
                                                do.call("rbind", 
                                                        lapply(c(1,2,3,4,5), 
                                                               function(x) {data.frame("state" = x, 
                                                                                       "model" = "loess.span0.1", 
                                                                                       loess.span0.1[[y]][[x]][1:5000, ])})))}))


plotdata.loess.span0.1.nostat <- do.call("rbind", 
                                  lapply(c(1,2,3), 
                                         function(y) {
                                           data.frame("tp" = paste("tp", y, sep = ""), 
                                                      do.call("rbind", 
                                                              lapply(c(1,2,3,4,5), 
                                                                     function(x) {data.frame("state" = x, 
                                                                                             "model" = "loess.span0.1.nostat", 
                                                                                             loess.span0.1[[y]][[x]][1:5000, ])})))}))

### Combine into one dataset
plotdata <- rbind(plotdata.rcs.nk3,
                  plotdata.rcs.nk4,
                  plotdata.rcs.nk5,
                  plotdata.loess.span0.1,
                  plotdata.loess.span0.1.nostat)

### Make plot
truecalib.plot <- ggplot2::ggplot(data = plotdata) +
  ggplot2::geom_line(ggplot2::aes(x = pred, y = obs, color = model)) +
  ggplot2::facet_wrap(~ tp + state, nrow = 3, ncol = 5) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggplot2::xlab("Predicted risk") +
  ggplot2::ylab("Observed risk") +
  ggplot2::theme(legend.position = "bottom")

CairoPNG("figures/compare_true_calib.png", width = 20, height = 12, res = 300, unit = "in")
truecalib.plot
dev.off()
