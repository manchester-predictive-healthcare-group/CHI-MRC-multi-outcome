###
### Calculate Monte Carlo Standard Error (MCE)
###

### Set working directory
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_8.2/")

### Source functions
R.func.sources <- list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load functions
source("R/sim_functions.R")

### Call libraries
library(ggplot2)

### Read in data
metadata <- readRDS("data/heuristic_shrinkage_sim_metadata.rds")
metadata <- subset(metadata, C.pop != 0)
metadata <- subset(metadata, cat == FALSE)

### Add SE for the diff
metadata$S.VH.diff.se <- metadata$S.VH.diff.sd/sqrt(250)

### Get quantiles
quantile(metadata$S.VH.diff.se, p = c(0.5,0.9,0.99,1))

### Take a look at those scenarios with high SE
metadata.high.MCE <- subset(metadata, S.VH.diff.se > 0.1)
metadata.high.MCE[1:30, c("C.pop", "P.meas", "P.unmeas", "nreq")]

### Write a function to get the quantiles for a given range of C-statistics
get_quantiles <- function(C1,C2){
  df <- subset(metadata, C.pop > C1 & C.pop <= C2)
  quantiles.MCE <- quantile(as.numeric(df$S.VH.diff.se), p = c(0.5,0.9,0.99,1))
  return(quantiles.MCE)
}

### Put MCE by quantiles ranges into a list

### Define cut-offs for C-statistic
quantiles.list <- list(c(0.6, 0.65),
                       c(0.65, 0.7),
                       c(0.7, 0.75),
                       c(0.75, 0.8),
                       c(0.8, 0.85),
                       c(0.85, 0.9))

### Get the quantiles of MCE across range of C-statistics
MCE.by.quantiles <- lapply(quantiles.list, function(x) {get_quantiles(x[1],x[2])})
MCE.by.quantiles <- do.call("rbind", MCE.by.quantiles)

### Apply this function to all (i.e. no subsettings on C-statistics)
MCE.all <- get_quantiles(0,1)

### Combine
MCE.by.quantiles <- rbind(MCE.all, S.VH.by.quantiles)


