### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Load packages
library(knitr)
library(ggpubr)

### Create vectors to loop through
scen.vec <- c("s1.1", "s1.2", "s2.1", "s2.2")
n.vec <- c(1000, 2500, 5000)

### Run a loop to calculate table for each
for (i in 1:length(scen.vec)){
  for (j in 1:length(n.vec)){
    scen <- scen.vec[i]
    n.devel <- n.vec[j]
    
    load(paste("data/sim_results_", scen, "_n", n.devel, "v1000.RData", sep = ""))
    
    C.Har.tab.all <- rbind(res.DGM.msm["C.Har.tab"]$C.Har.tab, res.DGM.clay["C.Har.tab"]$C.Har.tab, res.DGM.gumb["C.Har.tab"]$C.Har.tab, 
                           res.DGM.frank["C.Har.tab"]$C.Har.tab, res.DGM.normal["C.Har.tab"]$C.Har.tab, res.DGM.gamma["C.Har.tab"]$C.Har.tab)
    C.Uno.tab.all <- rbind(res.DGM.msm["C.Uno.tab"]$C.Uno.tab, res.DGM.clay["C.Uno.tab"]$C.Uno.tab, res.DGM.gumb["C.Uno.tab"]$C.Uno.tab, 
                           res.DGM.frank["C.Uno.tab"]$C.Uno.tab, res.DGM.normal["C.Uno.tab"]$C.Uno.tab, res.DGM.gamma["C.Uno.tab"]$C.Uno.tab)
    
    print(paste("scenario = ", scen, "n.devel = ", n.devel, sep = ""))
    print("HARREL's C")
    kable(C.Har.tab.all)
    print("UNO'S C")
    kable(C.Uno.tab.all)
  }
}


### Need to do no correlation scenario separately
scen.nocorr.vec <- c("s1", "s2")

for (i in 1:length(scen.nocorr.vec))
  for (j in 1:length(n.vec)){
    scen <- scen.nocorr.vec[i]
    n.devel <- n.vec[j]
    
    load(paste("data/sim_results_", scen, "_n", n.devel, "v1000.RData", sep = ""))
    
    print(paste("scenario = ", scen, "n.devel = ", n.devel, sep = ""))
    print("HARREL's C")
    kable(res.DGM.nocorr["C.Har.tab"]$C.Har.tab)
    print("UNO'S C")
    kable(res.DGM.nocorr["C.Uno.tab"]$C.Uno.tab)
  }


