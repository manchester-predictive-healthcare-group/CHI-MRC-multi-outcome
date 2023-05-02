###
### This program will...
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_5")
getwd()

### Load packages
source("code/z_load_packages.R")

### Set the seed
set.seed(505)

###
### Prepare the data for analysis
###

### Load the data
data("ebmt4")
ebmt <- ebmt4
rm(ebmt4)

### Define tmat
tmat <- transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
                          c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
tmat

### Create mstate format
msebmt <- msprep(data = ebmt, trans = tmat, time = c(NA, "rec", "ae","recae", "rel", "srv"), 
                 status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"), 
                 keep = c("match", "proph", "year", "agecl"))
msebmt[msebmt$id == 1, c(1:8, 10:12)]

### Define covariates for model
covs <- c("match", "proph", "year", "agecl")
msebmt <- expand.covs(msebmt, covs, longnames = FALSE)
msebmt[msebmt$id == 1, -c(9, 10, 12:48, 61:84)]

### Change time scale of model into years
msebmt[, c("Tstart", "Tstop", "time")] <- msebmt[, c("Tstart","Tstop", "time")]

### Create a variable which is maximum observed follow up time for all individuals, this is when they were either censored, relapsed or died
ebmt$dtcens <- pmin(ebmt$rel, ebmt$srv)
ebmt$dtcens.s <- 1 - pmax(ebmt$rel.s, ebmt$srv.s)

### Save image
save.image("data/prep_ebmt.RData")
