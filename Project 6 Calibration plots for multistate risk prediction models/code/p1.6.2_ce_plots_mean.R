### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")

### Read in n.cohort
n.devel <- 5000
n.pctls <- 20

### Load the workspace containing results
load(paste("data/ce_assess_calib_N", n.devel, "_npctls", n.pctls, ".RData", sep = ""))

###
### Create data for plots
###
tables.est <- vector("list", 6)
tables.est.se <- vector("list", 6)

### Cycle through states
for (state in 1:6){
  tables.est[[state]] <- round(c(
    calib.mean.blr.ipcw$diff.pred.obs[[state]],
    calib.mean.mlr.ipcw$diff.pred.obs[[state]],
    calib.mean.pv), 3)
  
  tables.est.se[[state]] <- round(c(calib.mean.blr.ipcw.se[["se"]][paste("diff.pred.obs.pspec", state, sep = "")],
                                    calib.mean.mlr.ipcw.se[["se"]][paste("diff.pred.obs.pspec", state, sep = "")],
                                    0), 3)
}


###
### Create table to be fed into dot and whisker plot in correct format
###
data.dw <- vector("list", 6)
names(data.dw) <- paste("state", 1:6, sep = "")
  for (state in 1:6){
    data.dw[[state]] <- data.frame("term" = rep("Mean calibration", 3), 
                                            "model" = c("BLR-IPCW", "MLR-IPCW", "PV"),
                                            "estimate" = c(tables.est[[state]], tables.est[[state]], tables.est[[state]]),
                                            "std.error" = c(tables.est.se[[state]], tables.est.se[[state]], tables.est.se[[state]]),
                                            "state" = rep(paste("State", state, sep = " "), 3))
  }
}


###
### Combine all data for a single plot with facet_grid
###
data.dw.facet <- do.call("rbind", lapply(data.dw, function(x){do.call("rbind", x)}))
gg.dw <- dwplot(data.dw.facet) +# scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
  facet_grid(cols = vars(state))  + 
  #xlab("Bias") + 
  ylab("Scenario")

### Save main plot
CairoPNG(paste("figures/gg_ce_mean_N", n.devel, "_npctls", n.pctls, ".png", sep = ""), 
         dpi = 300, width = 15, height = 3, unit = "in")
print(gg.dw)
dev.off()
