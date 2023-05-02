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
load(paste("data/combine_tp_cmprsk_leaveoo_j", j, "_s", s, "_t.eval", t.eval, ".RData", sep = ""))
load("data/prep_ebmt.RData")

### Define number of knots
nk <- 3

### An individual cannot go to state 4 from state 1, so remove these columns
plots.list <- calc.calib.cmprsk.mod(data.mstate = msebmt, 
                                    data.raw = ebmt, 
                                    j = j, 
                                    s = s, 
                                    t.eval = t.eval,
                                    p.est = tp.all[,paste("pstate", 1:6, sep = "")], 
                                    nk = nk)

### Arrange into one using ggarrange
if (length(plots.list) > 2){
  gg.plots.arrange <- ggarrange(plotlist = plots.list, nrow = 2, ncol = ceiling(length(plots.list)/2))
  
  ### Save plot
  CairoPNG(paste("figures/gg_tp_cmprsk_j", j, "_s", s, "_t", t.eval, "_nk", nk, ".png", sep = ""), 
           dpi = 300, width = 15, height = 10, unit = "in")
  print(gg.plots.arrange)
  dev.off()
  
} else if (length(plots.list) == 2){
  gg.plots.arrange <- ggarrange(plotlist = plots.list, nrow = 1, ncol = 2)
  
  ### Save plot
  CairoPNG(paste("figures/gg_tp_cmprsk_j", j, "_s", s, "_t", t.eval, "_nk", nk, ".png", sep = ""), 
           dpi = 300, width = 15, height = 7.5, unit = "in")
  print(gg.plots.arrange)
  dev.off()
} else if (length(plots.list) == 1){
  gg.plots.arrange <- ggarrange(plotlist = plots.list, nrow = 1, ncol = 1)
  
  ### Save plot
  CairoPNG(paste("figures/gg_tp_cmprsk_j", j, "_s", s, "_t", t.eval, "_nk", nk, ".png", sep = ""), 
           dpi = 300, width = 7.5, height = 7.5, unit = "in")
  print(gg.plots.arrange)
  dev.off()
}




