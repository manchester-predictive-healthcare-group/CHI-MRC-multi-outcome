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

#######################
### DO LOESS FIRST
#######################

###
### Create data for plots
###
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
                                        R.boot = 100,
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
                                                CI = 95,
                                                R.boot = 100,
                                                stabilised = FALSE)

calib.blr.loess.gen <- calc.calib.blr(data.mstate = msebmt,
                                       data.raw = ebmt,
                                       tmat, 
                                       j=j,
                                       s=s,
                                       t.eval = t.eval,
                                       p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                       curve.type = "loess",
                                       loess.span = span,
                                       loess.degree = degree.blr,
                                       weights = weights.fixed$ipcw,
                                       w.covs = covs,
                                       w.landmark.type = "all",
                                       w.max = 10,
                                       w.stabilised = FALSE,
                                       CI = 95,
                                       CI.R.boot = 100)

calib.blr.loess.ihw.gen <- calc.calib.blr(data.mstate = msebmt,
                                          data.raw = ebmt,
                                          tmat, 
                                          j=j,
                                          s=s,
                                          t.eval = t.eval,
                                          p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                          curve.type = "loess",
                                          loess.span = span,
                                          loess.degree = degree.blr,
                                          weights = NULL,
                                          w.covs = covs,
                                          w.landmark.type = "all",
                                          w.max = 10,
                                          w.stabilised = FALSE,
                                          CI = 95,
                                          CI.R.boot = 100)

###
### Create plots themselves
###
gg.plots.blr.loess <- plot.calib.blr(calib.blr.loess, combine = TRUE, nrow = 2, ncol = 3)
gg.plots.blr.loess.ihw <- plot.calib.blr(calib.blr.loess.ihw, combine = TRUE, nrow = 2, ncol = 3)
gg.plots.blr.loess.gen <- plot.calib.blr(calib.blr.loess.gen, combine = TRUE, nrow = 2, ncol = 3)
gg.plots.blr.loess.ihw.gen <- plot.calib.blr(calib.blr.loess.ihw.gen, combine = TRUE, nrow = 2, ncol = 3)

###
### Save plots
###
CairoPNG(paste("figures/gg_TEST_loess.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.loess)
dev.off()
CairoPNG(paste("figures/gg_TEST_loess.ihw.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.loess.ihw)
dev.off()
CairoPNG(paste("figures/gg_TEST_loess.gen.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.loess.gen)
dev.off()
CairoPNG(paste("figures/gg_TEST_loess.ihw.gen.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.loess.ihw.gen)
dev.off()


#######################
### NOW DO RCS
#######################

###
### Create data for plots
###
calib.blr.rcs <- calc.calib.blr.rcs(data.mstate = msebmt,
                                        data.raw = ebmt.weights,
                                        tmat, 
                                        j=j,
                                        s=s,
                                        t.eval = t.eval,
                                        p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                        nk = nk,
                                        CI = 95,
                                        R.boot = 100,
                                        stabilised = FALSE)

calib.blr.rcs.ihw <- calc.calib.blr.rcs.ihw(data.mstate = msebmt,
                                                data.raw = ebmt,
                                                tmat, 
                                                j=j,
                                                s=s,
                                                t.eval = t.eval,
                                                p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                                curve.nk = nk,
                                                CI = 95,
                                                R.boot = 500,
                                                stabilised = FALSE)

calib.blr.rcs.gen <- calc.calib.blr(data.mstate = msebmt,
                                      data.raw = ebmt,
                                      tmat, 
                                      j=j,
                                      s=s,
                                      t.eval = t.eval,
                                      p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                      curve.type = "rcs",
                                      rcs.nk = nk,
                                      weights = weights.fixed$ipcw,
                                      w.covs = covs,
                                      w.landmark.type = "all",
                                      w.max = 10,
                                      w.stabilised = FALSE,
                                      CI = 95,
                                      CI.R.boot = 100)

calib.blr.rcs.ihw.gen <- calc.calib.blr(data.mstate = msebmt,
                                          data.raw = ebmt,
                                          tmat, 
                                          j=j,
                                          s=s,
                                          t.eval = t.eval,
                                          p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                          curve.type = "rcs",
                                          rcs.nk = nk,
                                          weights = NULL,
                                          w.covs = covs,
                                          w.landmark.type = "all",
                                          w.max = 10,
                                          w.stabilised = FALSE,
                                          CI = 95,
                                          CI.R.boot = 100)

###
### Create plots themselves
###
gg.plots.blr.rcs <- plot.calib.blr(calib.blr.rcs, combine = TRUE, nrow = 2, ncol = 3)
gg.plots.blr.rcs.ihw <- plot.calib.blr(calib.blr.rcs.ihw, combine = TRUE, nrow = 2, ncol = 3)
gg.plots.blr.rcs.gen <- plot.calib.blr(calib.blr.rcs.gen, combine = TRUE, nrow = 2, ncol = 3)
gg.plots.blr.rcs.ihw.gen <- plot.calib.blr(calib.blr.rcs.ihw.gen, combine = TRUE, nrow = 2, ncol = 3)

###
### Save plots
###
CairoPNG(paste("figures/gg_TEST_rcs.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.rcs)
dev.off()
CairoPNG(paste("figures/gg_TEST_rcs.ihw.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.rcs.ihw)
dev.off()
CairoPNG(paste("figures/gg_TEST_rcs.gen.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.rcs.gen)
dev.off()
CairoPNG(paste("figures/gg_TEST_rcs.ihw.gen.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.rcs.ihw.gen)
dev.off()


###
### Test different weights
###
calib.blr.rcs.ihw.gen <- calc.calib.blr(data.mstate = msebmt,
                                        data.raw = ebmt,
                                        tmat, 
                                        j=j,
                                        s=s,
                                        t.eval = t.eval,
                                        p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                        curve.type = "rcs",
                                        rcs.nk = nk,
                                        weights = NULL,
                                        w.covs = covs,
                                        w.landmark.type = "all",
                                        w.max = 10,
                                        w.stabilised = FALSE,
                                        CI = 95,
                                        CI.R.boot = 500)


calib.blr.rcs.ihw.gen.wstab <- calc.calib.blr(data.mstate = msebmt,
                                        data.raw = ebmt,
                                        tmat, 
                                        j=j,
                                        s=s,
                                        t.eval = t.eval,
                                        p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                        curve.type = "rcs",
                                        rcs.nk = nk,
                                        weights = NULL,
                                        w.covs = covs,
                                        w.landmark.type = "all",
                                        w.max = 10,
                                        w.stabilised = TRUE,
                                        CI = 95,
                                        CI.R.boot = 500)

calib.blr.rcs.ihw.gen.wnull <- calc.calib.blr(data.mstate = msebmt,
                                           data.raw = ebmt,
                                           tmat, 
                                           j=j,
                                           s=s,
                                           t.eval = t.eval,
                                           p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                           curve.type = "rcs",
                                           rcs.nk = nk,
                                           weights = NULL,
                                           w.covs = NULL,
                                           w.landmark.type = "all",
                                           w.max = 10,
                                           w.stabilised = FALSE,
                                           CI = 95,
                                           CI.R.boot = 500)


calib.blr.rcs.ihw.gen.wnull.wstab <- calc.calib.blr(data.mstate = msebmt,
                                              data.raw = ebmt,
                                              tmat, 
                                              j=j,
                                              s=s,
                                              t.eval = t.eval,
                                              p.est = tp.all[[paste("j", j, sep = "")]] %>% select(any_of(paste("pstate", 1:6, sep = ""))),
                                              curve.type = "rcs",
                                              rcs.nk = nk,
                                              weights = NULL,
                                              w.covs = NULL,
                                              w.landmark.type = "all",
                                              w.max = 10,
                                              w.stabilised = TRUE,
                                              CI = 95,
                                              CI.R.boot = 500)

gg.plots.blr.rcs.ihw.gen <- plot.calib.blr(calib.blr.rcs.ihw.gen, combine = TRUE, nrow = 2, ncol = 3)
gg.plots.blr.rcs.ihw.gen.wstab <- plot.calib.blr(calib.blr.rcs.ihw.gen.wstab, combine = TRUE, nrow = 2, ncol = 3)
gg.plots.blr.rcs.ihw.gen.wnull <- plot.calib.blr(calib.blr.rcs.ihw.gen.wnull, combine = TRUE, nrow = 2, ncol = 3)
gg.plots.blr.rcs.ihw.gen.wnull.wstab <- plot.calib.blr(calib.blr.rcs.ihw.gen.wnull.wstab, combine = TRUE, nrow = 2, ncol = 3)

CairoPNG(paste("figures/gg_TEST_rcs.ihw.gen.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.rcs.ihw.gen)
dev.off()
CairoPNG(paste("figures/gg_TEST_rcs.ihw.gen.wstab.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.rcs.ihw.gen.wstab)
dev.off()
CairoPNG(paste("figures/gg_TEST_rcs.ihw.gen.wnull.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.rcs.ihw.gen.wnull)
dev.off()
CairoPNG(paste("figures/gg_TEST_rcs.ihw.gen.wnull.wstab.png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(gg.plots.blr.rcs.ihw.gen.wnull.wstab)
dev.off()