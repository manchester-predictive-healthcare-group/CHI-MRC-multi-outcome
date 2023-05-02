###
### This program will...
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_5")
getwd()

### Load packages
source("code/z_load_packages.R")
source("code/z_functions.R")

### Assign j, s and t.eval from comand line
args <- commandArgs(trailingOnly = T)
j <- as.numeric(args[1])
s <- as.numeric(args[2])
t.eval <- as.numeric(args[3])
print(paste("j = ", j, Sys.time()))
print(paste("s = ", s, Sys.time()))
print(paste("t.eval = ", t.eval, Sys.time()))

### Load predicted risks and prepped data
load("data/prep_ebmt.RData")
load(paste("data/combine_tp_cmprsk_leaveoo_j", j, "_s", s, "_t.eval", t.eval, ".RData", sep = ""))

### Calculate ipcw weights
weights <- calc.weights(data.mstate = msebmt, 
                        data.raw = ebmt, 
                        #covs = covs, 
                        covs = c("match", "proph", "agecl", "year"), 
                        j = j,
                        landmark.type = "state", 
                        s = s, 
                        t.eval = t.eval)

### Combineweights with ebmt dataset
ebmt.j.weights <- cbind(ebmt.j, weights %>% subset(id %in% ebmt.j$id) %>% select(c(ipcw,pcw)))

### Assign nk
nk <- 3

### 
### BLR
###

### Produce calibration plots
calib.blr.rcs <- calc.calib.blr.rcs(data.mstate = msebmt.j,
                             data.raw = ebmt.j.weights,
                             tmat, 
                             j=j,
                             s=s,
                             t.eval = t.eval,
                             p.est = tp.all %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                             nk = nk,
                             CI = 95,
                             R.boot = 500)

### Arrange into one using ggarrange
if (length(calib.blr.rcs[["plot.data"]]) > 2){
  gg.plots.blr.rcs <- plot.calib.blr(calib.blr.rcs, combine = TRUE, nrow = 2, ncol = ceiling(length(calib.blr.rcs[["plot.data"]])/2))
  
  ### Save plot
  CairoPNG(paste("figures/gg_tp_cmprsk_blrval_j", j, "_s", s, "_t", t.eval, "_nk", nk, ".png", sep = ""), 
           dpi = 300, width = 15, height = 10, unit = "in")
  print(gg.plots.blr.rcs)
  dev.off()
  
} else if (length(plots.list) == 2){
  gg.plots.blr.rcs <- plot.calib.blr(calib.blr.rcs, combine = TRUE, nrow = 2, ncol = 1)

  ### Save plot
  CairoPNG(paste("figures/gg_tp_cmprsk_blrval_j", j, "_s", s, "_t", t.eval, "_nk", nk, ".png", sep = ""), 
           dpi = 300, width = 15, height = 7.5, unit = "in")
  print(gg.plots.blr.rcs)
  dev.off()
} else if (length(plots.list) == 1){
  gg.plots.blr.rcs <- plot.calib.blr(calib.blr.rcs, combine = TRUE, nrow = 1, ncol = 1)

  ### Save plot
  CairoPNG(paste("figures/gg_tp_cmprsk_blrval_j", j, "_s", s, "_t", t.eval, "_nk", nk, ".png", sep = ""), 
           dpi = 300, width = 7.5, height = 7.5, unit = "in")
  print(gg.plots.blr.rcs)
  dev.off()
}




