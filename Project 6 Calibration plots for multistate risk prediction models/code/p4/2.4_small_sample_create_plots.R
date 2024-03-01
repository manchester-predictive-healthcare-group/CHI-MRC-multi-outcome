###
### This program will create plots for the small sample analysis, mean calibration.
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")
source("code/z_functions.R")

### This function will get the mean calibration data from the small sample analysis, calculate the bias in the estimtea of mean calibration,
### and get in format for plotting, for a given scenario, n.cohort and n.pctls
combine_mean_plotdata <- function(scen, n.cohort, n.pctls, type, sens = FALSE){
  
  ### Load data
  if (sens == FALSE){
    load(paste("data/sim/small_sample_mean_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, ".RData", sep = ""))
  } else if (sens == TRUE){
    load(paste("data/sim/small_sample_mean_analysis_sens_", scen, "_n", n.cohort, ".RData", sep = ""))
  }
  
  
  ### First calculate the difference between each calibration and the true risk
  bias.aj <- lapply(c(1,2,3), function(x){calib.aj.mean.list.comb[[x]] - calib.true.mean.list.comb[[x]]})
  bias.blr <- lapply(c(1,2,3), function(x){calib.blr.mean.list.comb[[x]] - calib.true.mean.list.comb[[x]]})
  bias.mlr <- lapply(c(1,2,3), function(x){calib.mlr.mean.list.comb[[x]] - calib.true.mean.list.comb[[x]]})
  
  ### Put into a list
  bias.all <- list(bias.aj, bias.blr, bias.mlr)
  names(bias.all) <- c("AJ", "BLR-IPCW", "MLR-IPCW")

  ### Create a table of mean bias and median bias, with CI or perncentile range for each
  bias.mean <- lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {lapply(bias.all[[x]], colMeans)})
  names(bias.mean) <- c("AJ", "BLR-IPCW", "MLR-IPCW")
  
  bias.sd <- lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {lapply(bias.all[[x]], function(y) {apply(y, 2, sd)})})
  names(bias.sd) <- c("AJ", "BLR-IPCW", "MLR-IPCW")
  
  bias.mean.lower <- lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {lapply(c(1,2,3), function(y) {bias.mean[[x]][[y]]} - 1.96*bias.sd[[x]][[y]])})
  names(bias.mean.lower) <- c("AJ", "BLR-IPCW", "MLR-IPCW")
  
  bias.mean.upper <- lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {lapply(c(1,2,3), function(y) {bias.mean[[x]][[y]]} + 1.96*bias.sd[[x]][[y]])})
  names(bias.mean.upper) <- c("AJ", "BLR-IPCW", "MLR-IPCW")
  
  bias.median <- lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {lapply(bias.all[[x]], function(y) {apply(y, 2, function(z) 
  {quantile(z, probs = 0.5)})})})
  names(bias.median) <- c("AJ", "BLR-IPCW", "MLR-IPCW")
  
  bias.p025 <- lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {lapply(bias.all[[x]], function(y) {apply(y, 2, function(z) 
  {quantile(z, probs = 0.025)})})})
  names(bias.p025) <- c("AJ", "BLR-IPCW", "MLR-IPCW")
  
  bias.p975 <- lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {lapply(bias.all[[x]], function(y) {apply(y, 2, function(z) 
  {quantile(z, probs = 0.975)})})})
  names(bias.p975) <- c("AJ", "BLR-IPCW", "MLR-IPCW")
  
  ### Assign vector for naming scenario in terms of censoring
  scen.names <- c("RC", "WAC", "SAC")
  names(scen.names) <- c("C1", "C2", "C3")
  
  ### Creata  data frame
  if (type == "mean"){
    plotdata <- data.frame(estimate = do.call("c", lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {do.call("c", bias.mean[[x]])})), 
                                conf.low = do.call("c", lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {do.call("c", bias.mean.lower[[x]])})),
                                conf.high = do.call("c", lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {do.call("c", bias.mean.upper[[x]])})),
                                state = rep(paste("State", 1:5, sep = ""), 9), 
                                tp = rep(c(rep(1, 5), rep(2, 5), rep(3, 5)), 3),
                                model = c(rep("AJ", 15), rep("BLR-IPCW", 15), rep("MLR-IPCW", 15)),
                                term = rep(scen.names[scen], 45))
  } else if (type == "median"){
    plotdata <- data.frame(estimate = do.call("c", lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {do.call("c", bias.median[[x]])})), 
                                conf.low = do.call("c", lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {do.call("c", bias.p025[[x]])})),
                                conf.high = do.call("c", lapply(c("AJ", "BLR-IPCW", "MLR-IPCW"), function(x) {do.call("c", bias.p975[[x]])})),
                                state = rep(paste("State", 1:5, sep = ""), 9), 
                                tp = rep(c(rep(1, 5), rep(2, 5), rep(3, 5)), 3),
                                model = c(rep("AJ", 15), rep("BLR-IPCW", 15), rep("MLR-IPCW", 15)),
                                term = rep(scen.names[scen], 45))
  }

  ### Assign variable formats
  plotdata$tp <- factor(plotdata$tp, labels = c("Perfect", "Miscalibrated 1", "Miscalibrated 2"))
  
  ### Return
  return(plotdata)
  
}

### Define n.cohort and n.pctls
for (n.cohort in c(1500, 3000)){
  
  ###
  ### Create plots for sensitivity analyses
  ###
  
  ### Apply this functio to each scenario and combine into a single dataset
  plotdata.mean <- do.call("rbind", lapply(c("C1", "C2", "C3"), 
                                           function(x) {combine_mean_plotdata(scen = x, 
                                                                              n.cohort = n.cohort, 
                                                                              n.pctls = n.pctls, 
                                                                              type = "mean", 
                                                                              sens = TRUE)}))
  
  plotdata.median <- do.call("rbind", lapply(c("C1", "C2", "C3"), 
                                             function(x) {combine_mean_plotdata(scen = x, 
                                                                                n.cohort = n.cohort, 
                                                                                n.pctls = n.pctls, 
                                                                                type = "median", 
                                                                                sens = TRUE)}))
  
  
  ### Put into dotwhisker plots
  gg.dw.mean <- dwplot(plotdata.mean) + #scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
    facet_grid(rows = vars(tp), cols = vars(state))  + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    xlab("Bias") + ylab("Scenario") + scale_x_continuous(breaks = seq(-0.05, 0.05, 0.025), labels = function(x) ifelse(x == 0, "0", x))
  
  gg.dw.median <- dwplot(plotdata.median) + #scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
    facet_grid(rows = vars(tp), cols = vars(state))  + 
    geom_vline(xintercept = 0, linetype = "dashed") + 
    xlab("Bias") + ylab("Scenario") + scale_x_continuous(breaks = seq(-0.05, 0.05, 0.025), labels = function(x) ifelse(x == 0, "0", x))
  
  
  ### Save images
  CairoPNG(paste("figures/sens_small_sample_mean_N", n.cohort,"_mean.png", sep = ""), 
           dpi = 150, width = 15, height = 10, unit = "in")
  print(gg.dw.mean)
  dev.off()
  CairoPDF(paste("figures/sens_small_sample_mean_N", n.cohort, "_mean.pdf", sep = ""), 
           width = 15, height = 10)
  print(gg.dw.mean)
  dev.off()
  
  CairoPNG(paste("figures/sens_small_sample_mean_N", n.cohort, "_median.png", sep = ""), 
           dpi = 150, width = 15, height = 10, unit = "in")
  print(gg.dw.median)
  dev.off()
  CairoPDF(paste("figures/sens_small_sample_mean_N", n.cohort, "_median.pdf", sep = ""), 
           width = 15, height = 10)
  print(gg.dw.median)
  dev.off()
  
  for (n.pctls in c(5, 10, 20)){

    ###
    ### Create plots for main analyses, for each number of percentiles when doing the AJ analysis 
    ### (would it actually be better to put the different number of percentiles in the same plot?)
    ###
    
    ### Apply this functio to each scenario and combine into a single dataset
    plotdata.mean <- do.call("rbind", lapply(c("C1", "C2", "C3"), 
                                             function(x) {combine_mean_plotdata(scen = x, n.cohort = n.cohort, n.pctls = n.pctls, type = "mean")}))
    
    plotdata.median <- do.call("rbind", lapply(c("C1", "C2", "C3"), 
                                               function(x) {combine_mean_plotdata(scen = x, n.cohort = n.cohort, n.pctls = n.pctls, type = "median")}))
    
    
    ### Put into dotwhisker plots
    gg.dw.mean <- dwplot(plotdata.mean) + #scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
      facet_grid(rows = vars(tp), cols = vars(state))  + 
      geom_vline(xintercept = 0, linetype = "dashed") + 
      xlab("Bias") + ylab("Scenario") + scale_x_continuous(breaks = seq(-0.025, 0.025, 0.025), labels = function(x) ifelse(x == 0, "0", x))
    
    gg.dw.median <- dwplot(plotdata.median) + #scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
      facet_grid(rows = vars(tp), cols = vars(state))  + 
      geom_vline(xintercept = 0, linetype = "dashed") + 
      xlab("Bias") + ylab("Scenario") + scale_x_continuous(breaks = seq(-0.025, 0.025, 0.025), labels = function(x) ifelse(x == 0, "0", x))
    
    
    ### Save images
    CairoPNG(paste("figures/small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_mean.png", sep = ""), 
             dpi = 150, width = 15, height = 10, unit = "in")
    print(gg.dw.mean)
    dev.off()
    CairoPDF(paste("figures/small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_mean.pdf", sep = ""), 
             width = 15, height = 10)
    print(gg.dw.mean)
    dev.off()
    
    CairoPNG(paste("figures/small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_median.png", sep = ""), 
             dpi = 150, width = 15, height = 10, unit = "in")
    print(gg.dw.median)
    dev.off()
    CairoPDF(paste("figures/small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_median.pdf", sep = ""), 
             width = 15, height = 10)
    print(gg.dw.median)
    dev.off()
    
    ### Save table
    plotdata.median$estimate <- round(plotdata.median$estimate, 3)
    plotdata.median$conf.high <- round(plotdata.median$conf.high, 3)
    plotdata.median$conf.low <- round(plotdata.median$conf.low, 3)
    write.table(plotdata.median, paste("figures/table_small_sample_mean_N", n.cohort, "_npctls", n.pctls, "_median.csv", sep = ""))
    
  }
}



