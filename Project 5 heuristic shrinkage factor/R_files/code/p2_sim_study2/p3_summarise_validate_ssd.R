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
metadata <- readRDS("data/heuristic_shrinkage_sim_metadata_ssd.rds")

### Split metadata into boservations both are recorded, and obserations just Riley recorded
metadata.both <- subset(metadata, !is.na(S.pop.mean.riley) & !is.na(S.pop.mean.pavlou))
metadata.all.riley <- subset(metadata, !is.na(S.pop.mean.riley))
metadata.just.riley <- subset(metadata, !is.na(S.pop.mean.riley) & is.na(S.pop.mean.pavlou))
str(metadata.both)
### Write a function to get the quantiles for a given range of C-statistics
get_quantiles <- function(C1,C2){
  df <- subset(metadata.both, C.pop > C1 & C.pop <= C2)
  quantiles.riley <- quantile(as.numeric(df$S.pop.mean.riley), p = c(0.025,0.25,0.5,0.75,0.975))
  quantiles.pavlou <- quantile(as.numeric(df$S.pop.mean.pavlou), p = c(0.025,0.25,0.5,0.75,0.975))
  return(list("riley" = quantiles.riley,
              "pavlou" = quantiles.pavlou))
}

### Put S_VH by quantiles ranges into a list
quantiles.list <- list(c(0.6, 0.65),
                       c(0.65, 0.7),
                       c(0.7, 0.75),
                       c(0.75, 0.8),
                       c(0.8, 0.85),
                       c(0.85, 0.9))
S.VH.by.quantiles <- lapply(quantiles.list, function(x) {get_quantiles(x[1],x[2])})
S.VH.all <- get_quantiles(0,1)

### Put into a table
S.VH.by.quantiles <- lapply(S.VH.by.quantiles, function(x) {
  out <- rbind(x[[1]], x[[2]])
  rownames(out) <- c("RvS", "Rvs-adj")
  return(out)
})
