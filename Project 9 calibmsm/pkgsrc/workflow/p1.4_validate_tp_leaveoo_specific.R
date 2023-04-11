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

### Assign s and t.eval from comand line
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])
s <- as.numeric(args[2])
t.eval <- as.numeric(args[3])
print(paste("model = ", model, Sys.time()))
print(paste("s = ", s, Sys.time()))
print(paste("t.eval = ", t.eval, Sys.time()))
model <- 1
s <- 0
t.eval <- 1826
### Load predicted risks and prepped data
load(paste("data/combine_tp_leaveoo_model", model, "_s", s, "_t", t.eval, ".RData", sep = ""))
load("data/prep_ebmt.RData")

### Choose span and degree for loess plots
span <- 1
degree.blr <- 2

### Choose ps.int and degree for mlr plots
ps.int <- 4
degree.mlr <- 3

### Choose number of knots for rcs
nk <- 3

### Calibration of transition probabilities out of each state (note there is no point in doing this for states 5 and 6, which are absorbing)
### Also for s = 0, only need to do j = 1, as everyone starts in initial state (cant be j > 1 for s= 0)

if (s == 0){
  ### Set j = 1
  j <- 1
  
  ### Calculate ipcw weights
  weights.fixed <- calc.weights(data.mstate = msebmt, 
                          data.raw = ebmt, 
                          covs = covs, 
                          j = j,
                          landmark.type = "all", 
                          s = s, 
                          t.eval = t.eval)
  
  ### Combineweights with ebmt dataset
  ebmt.weights <- cbind(ebmt, weights.fixed)
  
  ### 
  ### Produce calibration data
  ###
  
  ### BLR using loess smoothers
  calib.blr.loess <- calc.calib.blr.loess(data.mstate = msebmt,
                                          data.raw = ebmt.weights,
                                          tmat, 
                                          j=j,
                                          s=s,
                                          t.eval = t.eval,
                                          p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                          span = span,
                                          degree = degree.blr,
                                          CI = FALSE,
                                          R.boot = 250,
                                          stabilised = FALSE)
  
  calib.blr.loess.ihw <- calc.calib.blr.loess.ihw(data.mstate = msebmt,
                                          data.raw = ebmt,
                                          tmat, 
                                          j=j,
                                          s=s,
                                          t.eval = t.eval,
                                          p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                          span = span,
                                          degree = degree.blr,
                                          CI = FALSE,
                                          R.boot = 250,
                                          stabilised = FALSE)
  
  ### BLR using rcs
  calib.blr.rcs <- calc.calib.blr.rcs(data.mstate = msebmt,
                                      data.raw = ebmt.weights,
                                      tmat, 
                                      j=j,
                                      s=s,
                                      t.eval = t.eval,
                                      p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                      nk = 3,
                                      CI = 95,
                                      R.boot = 500)
  
  ### MLR using nominal recalibration framework
  calib.mlr <- calc.calib.mlr(data.mstate = msebmt,
                              data.raw = ebmt.weights,
                              tmat, 
                              j=j,
                              s=s,
                              t.eval = t.eval,
                              p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                              ps.int = 4, 
                              degree = 3)
  
  ### 
  ### Produce plots
  ###
  gg.plots.blr.loess <- plot.calib.blr(calib.blr.loess, combine = TRUE, nrow = 2, ncol = 3)
  gg.plots.blr.loess.ihw <- plot.calib.blr(calib.blr.loess.ihw, combine = TRUE, nrow = 2, ncol = 3)
  gg.plots.blr.rcs <- plot.calib.blr(calib.blr.rcs, combine = TRUE, nrow = 2, ncol = 3)
  gg.plots.mlr <- plot.calib.mlr(calib.mlr, combine = TRUE, nrow = 2, ncol = 3)
  
  ### Save plots
  CairoPNG(paste("figures/gg_tp_j", j, "_s", s, "_t", t.eval, "_blr_loess_span", span, "_degree", degree.blr, ".png", sep = ""), 
           dpi = 300, width = 15, height = 10, unit = "in")
  print(gg.plots.blr.loess)
  dev.off()
  
  CairoPNG(paste("figures/gg_tp_j", j, "_s", s, "_t", t.eval, "_blr_loess_iwh_span", span, "_degree", degree.blr, ".png", sep = ""), 
           dpi = 300, width = 15, height = 10, unit = "in")
  print(gg.plots.blr.loess.ihw)
  dev.off()
  source("code/z_functions.R")
  head(calib.blr.loess[[1]][[1]])
  head(calib.blr.loess.ihw[[1]][[1]])
  
  sum(calib.blr.loess[[1]][[1]]$obs != calib.blr.loess.ihw[[1]][[1]]$obs)
  sum(calib.blr.loess[[1]][[1]]$pred != calib.blr.loess.ihw[[1]][[1]]$pred)
  CairoPNG(paste("figures/gg_tp_j", j, "_s", s, "_t", t.eval, "_blr_rcs_nk", nk, ".png", sep = ""), 
           dpi = 300, width = 15, height = 10, unit = "in")
  print(gg.plots.blr.rcs)
  dev.off()
  
  CairoPNG(paste("figures/gg_tp_j", j, "_s", s, "_t", t.eval, "_mlr_psint", ps.int, "_degree", degree.mlr, ".png", sep = ""), 
           dpi = 300, width = 15, height = 10, unit = "in")
  print(gg.plots.mlr)
  dev.off()

} else if (s > 0){
  ### Cycle through j's
  for (j in 1:4){
    
    ### Calculate ipcw weights
    weights <- calc.weights(data.mstate = msebmt, 
                            data.raw = ebmt, 
                            covs = covs, 
                            j = j,
                            landmark.type = "all", 
                            s = s, 
                            t.eval = t.eval)
    
    ### Combineweights with ebmt dataset
    ebmt.weights <- cbind(ebmt, weights)
    
    ### 
    ### Produce calibration data
    ###
    
    ### BLR using loess smoothers
    calib.blr.loess <- calc.calib.blr.loess(data.mstate = msebmt,
                                            data.raw = ebmt.weights,
                                            tmat, 
                                            j=j,
                                            s=s,
                                            t.eval = t.eval,
                                            p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                            span = span,
                                            degree = degree.blr,
                                            CI = 95,
                                            R.boot = 500)
    
    ### BLR using rcs
    calib.blr.rcs <- calc.calib.blr.rcs(data.mstate = msebmt,
                                        data.raw = ebmt.weights,
                                        tmat, 
                                        j=j,
                                        s=s,
                                        t.eval = t.eval,
                                        p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                        nk = 3,
                                        CI = 95,
                                        R.boot = 500)
    
    ### MLR using nominal recalibration framework
    calib.mlr <- calc.calib.mlr(data.mstate = msebmt,
                                data.raw = ebmt.weights,
                                tmat, 
                                j=j,
                                s=s,
                                t.eval = t.eval,
                                p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                ps.int = 4, 
                                degree = 3)
    
    ### 
    ### Produce plots
    ###
    gg.plots.blr.loess <- plot.calib.blr(calib.blr.loess, combine = TRUE, nrow = 2, ncol = ceiling(length(calib.blr.loess[["metadata"]][["valid.transitions"]])/2))
    gg.plots.blr.rcs <- plot.calib.blr(calib.blr.rcs, combine = TRUE, nrow = 2, ncol = ceiling(length(calib.blr.loess[["metadata"]][["valid.transitions"]])/2))
    gg.plots.mlr <- plot.calib.mlr(calib.mlr, combine = TRUE, nrow = 2, ncol = ceiling(length(calib.blr.loess[["metadata"]][["valid.transitions"]])/2))
    
    ### Save plots
    CairoPNG(paste("figures/gg_tp_j", j, "_s", s, "_t", t.eval, "_blr_loess_span", span, "_degree", degree.blr, ".png", sep = ""), 
             dpi = 300, width = 15, height = 10, unit = "in")
    print(gg.plots.blr.loess)
    dev.off()
    
    CairoPNG(paste("figures/gg_tp_j", j, "_s", s, "_t", t.eval, "_blr_rcs_nk", nk, ".png", sep = ""), 
             dpi = 300, width = 15, height = 10, unit = "in")
    print(gg.plots.blr.rcs)
    dev.off()
    
    CairoPNG(paste("figures/gg_tp_j", j, "_s", s, "_t", t.eval, "_mlr_psint", ps.int, "_degree", degree.mlr, ".png", sep = ""), 
             dpi = 300, width = 15, height = 10, unit = "in")
    print(gg.plots.mlr)
    dev.off()
    
    ### Clear workspace
    rm(calib.blr.ipcw.mod, gg.plots.arrange)
    
  }
}

