###
### Create plots for mean calibration.
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

# ### Extract scenario from command line
# args <- commandArgs(trailingOnly = T)
# scen <- args[1]
# print(paste("scen = ", scen, sep = ""))
### Load calibration data for BLR and MLR


### Load packages
source("code/z_load_packages.R")
source("code/z_functions.R")

### This function will get the mean calibration data from the large sample analysis and get in format for plotting, for a given scenario
### and Method
combine_mean_plotdata <- function(scen, method, sens = FALSE){

  ### Assign vector for naming scenario in terms of censoring
  scen.names <- c("RC", "WAC", "SAC")
  names(scen.names) <- c("C1", "C2", "C3")
  
  ### Load calibration data
  if (sens == FALSE){
    load(paste("data/sim/large_sample_mean_", scen, ".RData", sep = ""))
  } else if (sens == TRUE){
    load(paste("data/sim/large_sample_mean_sens_", scen, ".RData", sep = ""))
  }

  
  ### Save BLR data in required format
  if (method == "blr"){
    plotdata <- do.call("rbind", 
                        lapply(c(1,2,3), 
                               function(x) 
                               {out <- data.frame(
                                 "model" = "BLR-IPCW",
                                 "tp" = x, 
                                 "scen" = scen.names[scen], 
                                 "state" = paste("State", 1:5),
                                 do.call("rbind", calib.mean.blr[[x]][["mean"]]))
                                rownames(out) <- NULL
                                return(out)}))
  } else if (method == "mlr"){
    plotdata <- do.call("rbind", 
                        lapply(c(1,2,3), 
                               function(x) 
                               {out <- data.frame(
                                 "model" = "MLR-IPCW",
                                 "tp" = x, 
                                 "scen" = scen.names[scen], 
                                 "state" = paste("State", 1:5),
                                 do.call("rbind", calib.mean.mlr[[x]][["mean"]]))
                                rownames(out) <- NULL
                                return(out)}))
  } else if (method == "aj"){
    plotdata <- do.call("rbind", 
                        lapply(c(1,2,3), 
                               function(x) 
                               {out <- data.frame(
                                 "model" = "AJ",
                                 "tp" = x, 
                                 "scen" = scen.names[scen], 
                                 "state" = paste("State", 1:5),
                                 do.call("rbind", calib.mean.aj[[x]][["mean"]]))
                                rownames(out) <- NULL
                                return(out)}))
  }

  return(plotdata)
  
}

###
### Get true calibration for each state in the three prediction scenarios

###
### Write a function to produce the plots
###
create_mean_plot <- function(sens){
  
  ### Load calibration data
  load(paste("data/sim/large_sample_mean_C1.RData", sep = ""))
  calib.true.mean.list <- lapply(c(1,2,3), function(x) {colMeans(tp.true - tp.pred[[x]])})
  
  ### Extract and format plotdata
  if (sens == FALSE){
    plotdata.aj <- do.call("rbind", lapply(c("C1", "C2", "C3"), function(x) {combine_mean_plotdata(x, method = "aj")}))
    plotdata.blr <- do.call("rbind", lapply(c("C1", "C2", "C3"), function(x) {combine_mean_plotdata(x, method = "blr")}))
    plotdata.mlr <- do.call("rbind", lapply(c("C1", "C2", "C3"), function(x) {combine_mean_plotdata(x, method = "mlr")}))
  } else if (sens == TRUE){
    plotdata.aj <- do.call("rbind", lapply(c("C1", "C2", "C3"), function(x) {combine_mean_plotdata(x, method = "aj", sens = TRUE)}))
    plotdata.blr <- do.call("rbind", lapply(c("C1", "C2", "C3"), function(x) {combine_mean_plotdata(x, method = "blr", sens = TRUE)}))
    plotdata.mlr <- do.call("rbind", lapply(c("C1", "C2", "C3"), function(x) {combine_mean_plotdata(x, method = "mlr", sens = TRUE)}))
  }
  
  ### Combine into one dataset
  plotdata.mean <- rbind(plotdata.aj, plotdata.blr, plotdata.mlr)
  
  ### Assign variable formats
  plotdata.mean$tp <- factor(plotdata.mean$tp, labels = c("Perfect", "Miscalibrated 1", "Miscalibrated 2"))
  
  ### Give variables correct names for dwplot
  plotdata.mean <- dplyr::rename(plotdata.mean, term = scen, estimate = mean, conf.low = mean.lower, conf.high = mean.upper)
  
  ### Put these into a dataframe
  int.dat <- data.frame(int = do.call("c", calib.true.mean.list), 
                        state = paste("State", 1:5), 
                        tp = c(rep(1,5), rep(2,5), rep(3,5)))
  int.dat$tp <- factor(int.dat$tp, labels = c("Perfect", "Miscalibrated 1", "Miscalibrated 2"))
  
  ### Put into dotwhisker plots
  gg.dw.mean <- dwplot(plotdata.mean) + 
    geom_vline(data = int.dat, aes(xintercept = int, color = "True\ncalibration"), show.legend = TRUE) + #scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
    facet_grid(rows = vars(tp), cols = vars(state))  + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    xlab("Mean calibration") + ylab("Scenario") + scale_x_continuous(breaks = seq(-0.04, 0.04, 0.02), labels = function(x) ifelse(x == 0, "0", x)) +
    guides(color = guide_legend(title = "Method", override.aes = list(linetype = c(0, 0, 0, 1), shape = c(19, 19, 19, NA)))) + 
    theme(legend.position = "bottom")
  
  return(gg.dw.mean)

}


###
### Create and save images
###

### Main analysis
CairoPNG(paste("figures/large_sample_mean.png", sep = ""), 
         dpi = 150, width = 15, height = 10, unit = "in")
print(create_mean_plot(sens = FALSE))
dev.off()
CairoPDF(paste("figures/large_sample_mean.pdf", sep = ""), 
         width = 15, height = 10)
print(create_mean_plot(sens = FALSE))
dev.off()


### Sensitivity analyses
CairoPNG(paste("figures/large_sample_mean_sens.png", sep = ""), 
         dpi = 150, width = 15, height = 10, unit = "in")
print(create_mean_plot(sens = TRUE))
dev.off()
CairoPDF(paste("figures/large_sample_mean_sens.pdf", sep = ""), 
         width = 15, height = 10)
print(create_mean_plot(sens = TRUE))
dev.off()
