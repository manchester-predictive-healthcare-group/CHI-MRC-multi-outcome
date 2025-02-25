### Set working directory
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_8.2/")

### Source functions
R.func.sources <- list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load functions
source("R/sim_functions.R")

### Write a function to extract the metadata
extract_metadata <- function(set, corr, cat, noNA){
  
  # set <- 221
  # corr <- TRUE
  # cat <- TRUE
  # noNA <- TRUE 
  # 
  ### Create object to store plot data
  metadata <- NULL
  
  print(paste("SET", set))
  ### Load workspace
  if (corr == TRUE){
    if (cat == FALSE){
      load(paste("data/run_sim_corr", set, ".RData", sep = ""))
    } else {
      load(paste("data/run_sim_cat_corr", set, ".RData", sep = ""))
    }
  } else {
    if (cat == FALSE){
      load(paste("data/run_sim_nocorr", set, ".RData", sep = ""))
    } else {
      load(paste("data/run_sim_cat_nocorr", set, ".RData", sep = ""))
    }
  }
  
  ### First start by getting the number of complete simulations
  n.comp <- max(which(sapply(output.S.VH, is.null) == FALSE))
  
  ### Create a list to store the vectors in which we will combine
  vec.list <- vector("list", n.comp)
  
  ### Create each vector
  for (sim.iter in 1:n.comp){
    ### Create vector
    if (noNA == TRUE){
      vec.list[[sim.iter]] <- c("S.VH.mean" = mean(output.S.VH[[sim.iter]], na.rm = TRUE),
                                "S.VH.sd" = sd(output.S.VH[[sim.iter]], na.rm = TRUE),
                                "S.VH.median" = median(output.S.VH[[sim.iter]], na.rm = TRUE),
                                "S.boot.mean" = mean(output.S.boot[[sim.iter]], na.rm = TRUE),
                                "S.boot.sd" = sd(output.S.boot[[sim.iter]], na.rm = TRUE),
                                "S.boot.median" = median(output.S.boot[[sim.iter]], na.rm = TRUE),
                                "S.pop.mean" = mean(output.S.pop[[sim.iter]], na.rm = TRUE),
                                "S.pop.sd" = sd(output.S.pop[[sim.iter]], na.rm = TRUE),
                                "S.pop.median" = median(output.S.pop[[sim.iter]], na.rm = TRUE),
                                "S.VH.diff.mean" = mean(output.S.VH[[sim.iter]] - output.S.pop[[sim.iter]], na.rm = TRUE),
                                "S.VH.diff.sd" = sd(output.S.VH[[sim.iter]] - output.S.pop[[sim.iter]], na.rm = TRUE),
                                "S.VH.diff.median" = median(output.S.VH[[sim.iter]] - output.S.pop[[sim.iter]], na.rm = TRUE),
                                "S.boot.diff.mean" = mean(output.S.boot[[sim.iter]] - output.S.pop[[sim.iter]], na.rm = TRUE),
                                "S.boot.diff.sd" = sd(output.S.boot[[sim.iter]] - output.S.pop[[sim.iter]], na.rm = TRUE),
                                "S.boot.diff.median" = median(output.S.boot[[sim.iter]] - output.S.pop[[sim.iter]], na.rm = TRUE),
                                "LR.mean" = mean(output.LR[[sim.iter]], na.rm = TRUE),
                                "LR.sd" = sd(output.LR[[sim.iter]], na.rm = TRUE),
                                "LR.median" = median(output.LR[[sim.iter]], na.rm = TRUE),
                                "C.mean" = mean(output.C[[sim.iter]], na.rm = TRUE),
                                "C.sd" = sd(output.C[[sim.iter]], na.rm = TRUE),
                                "C.median" = median(output.C[[sim.iter]], na.rm = TRUE),
                                "D.mean" = mean(output.D[[sim.iter]], na.rm = TRUE),
                                "D.sd" = sd(output.D[[sim.iter]], na.rm = TRUE),
                                "D.median" = median(output.D[[sim.iter]], na.rm = TRUE),
                                "R2.CS.app.mean" = mean(output.R2.CS.app[[sim.iter]], na.rm = TRUE),
                                "R2.CS.app.sd" = sd(output.R2.CS.app[[sim.iter]], na.rm = TRUE),
                                "R2.CS.app.median" = median(output.R2.CS.app[[sim.iter]], na.rm = TRUE),
                                "C.pop" = input.data[[sim.iter]][["C.pop"]],
                                "R2.CS.pop" = input.data[[sim.iter]][["R2.CS.pop"]],
                                "P.meas" = input.data[[sim.iter]][["P"]],
                                "P.unmeas" = input.data[[sim.iter]][["P.total"]] - input.data[[sim.iter]][["P"]],
                                "nreq" = input.data[[sim.iter]][["nreq"]],
                                "prop" = input.data[[sim.iter]][["prop"]], 
                                set = set, 
                                sim.iter = sim.iter)
    } else {
      vec.list[[sim.iter]] <- c("S.VH.mean" = mean(output.S.VH[[sim.iter]]),
                                "S.VH.sd" = sd(output.S.VH[[sim.iter]]),
                                "S.VH.median" = median(output.S.VH[[sim.iter]]),
                                "S.boot.mean" = mean(output.S.boot[[sim.iter]]),
                                "S.boot.sd" = sd(output.S.boot[[sim.iter]]),
                                "S.boot.median" = median(output.S.boot[[sim.iter]]),
                                "S.pop.mean" = mean(output.S.pop[[sim.iter]]),
                                "S.pop.sd" = sd(output.S.pop[[sim.iter]]),
                                "S.pop.median" = median(output.S.pop[[sim.iter]]),
                                "S.VH.diff.mean" = mean(output.S.VH[[sim.iter]] - output.S.pop[[sim.iter]]),
                                "S.VH.diff.sd" = sd(output.S.VH[[sim.iter]] - output.S.pop[[sim.iter]]),
                                "S.VH.diff.median" = median(output.S.VH[[sim.iter]] - output.S.pop[[sim.iter]]),
                                "S.boot.diff.mean" = mean(output.S.boot[[sim.iter]] - output.S.pop[[sim.iter]]),
                                "S.boot.diff.sd" = sd(output.S.boot[[sim.iter]] - output.S.pop[[sim.iter]]),
                                "S.boot.diff.median" = median(output.S.boot[[sim.iter]] - output.S.pop[[sim.iter]]),
                                "LR.mean" = mean(output.LR[[sim.iter]]),
                                "LR.sd" = sd(output.LR[[sim.iter]]),
                                "LR.median" = median(output.LR[[sim.iter]]),
                                "C.mean" = mean(output.C[[sim.iter]]),
                                "C.sd" = sd(output.C[[sim.iter]]),
                                "C.median" = median(output.C[[sim.iter]]),
                                "D.mean" = mean(output.D[[sim.iter]]),
                                "D.sd" = sd(output.D[[sim.iter]]),
                                "D.median" = median(output.D[[sim.iter]]),
                                "R2.CS.app.mean" = mean(output.R2.CS.app[[sim.iter]]),
                                "R2.CS.app.sd" = sd(output.R2.CS.app[[sim.iter]]),
                                "R2.CS.app.median" = median(output.R2.CS.app[[sim.iter]]),
                                "C.pop" = input.data[[sim.iter]][["C.pop"]],
                                "R2.CS.pop" = input.data[[sim.iter]][["R2.CS.pop"]],
                                "P.meas" = input.data[[sim.iter]][["P"]],
                                "P.unmeas" = input.data[[sim.iter]][["P.total"]] - input.data[[sim.iter]][["P"]],
                                "nreq" = input.data[[sim.iter]][["nreq"]],
                                "prop" = input.data[[sim.iter]][["prop"]], 
                                set = set, 
                                sim.iter = sim.iter)
    }
  }
  
  ### Create temp plot data and turn into a list, so can be combined with existing list
  metadata <- data.frame(do.call("rbind", vec.list))
  
  ### Remove rows with negative LR, driven by a small number of obscure results
  metadata <- metadata[metadata$R2.CS.app.mean > 0, ]
  metadata <- metadata[metadata$S.VH.mean > 0, ]
  
  ## Add variable for correlation
  if (corr == TRUE){
    metadata$corr <- TRUE
  } else {
    metadata$corr <- FALSE
  }
  
  ## Add variable for categorical
  if (cat == TRUE){
    metadata$cat <- TRUE
  } else {
    metadata$cat <- FALSE
  }
  
  return(metadata)
  
}

### Find numbers for which simulations have run
library(stringr)
valid.files.corr <- as.numeric(str_extract(list.files("data"), "(?<=run_sim_corr)[0-9]*")[!is.na(str_extract(list.files("data"), "(?<=run_sim_corr)[0-9]*"))])
valid.files.corr <- valid.files.corr[!is.na(valid.files.corr)]

valid.files.nocorr <- as.numeric(str_extract(list.files("data"), "(?<=run_sim_nocorr)[0-9]*")[!is.na(str_extract(list.files("data"), "(?<=run_sim_nocorr)[0-9]*"))])
valid.files.nocorr <- valid.files.nocorr[!is.na(valid.files.nocorr)]

### Previously did scenarios with set > 300, but these used an old DGM, so remove
valid.files.corr <- valid.files.corr[valid.files.corr<=300]
valid.files.nocorr <- valid.files.nocorr[valid.files.nocorr<=300]

# valid.files.cat.corr <- as.numeric(str_extract(list.files("data"), "(?<=run_sim_cat_corr)[0-9]*")[!is.na(str_extract(list.files("data"), "(?<=run_sim_cat_corr)[0-9]*"))])
# valid.files.cat.corr <- valid.files.cat.corr[!is.na(valid.files.cat.corr)]
# ## Sim run 63 didn't produce any scenarios, so remove
# valid.files.cat.corr <- valid.files.cat.corr[!(valid.files.cat.corr==63)]
# 
# valid.files.cat.nocorr <- as.numeric(str_extract(list.files("data"), "(?<=run_sim_cat_nocorr)[0-9]*")[!is.na(str_extract(list.files("data"), "(?<=run_sim_cat_nocorr)[0-9]*"))])
# valid.files.cat.nocorr <- valid.files.cat.nocorr[!is.na(valid.files.cat.nocorr)]


### Extract metadata
metadata.corr <- lapply(valid.files.corr, extract_metadata, corr = TRUE, cat = FALSE, noNA = TRUE)
metadata.nocorr <- lapply(valid.files.nocorr, extract_metadata, corr = FALSE, cat = FALSE, noNA = TRUE)

# metadata.cat.corr <- lapply(valid.files.cat.corr, extract_metadata, corr = TRUE, cat = TRUE, noNA = TRUE)
# metadata.cat.nocorr <- lapply(valid.files.cat.nocorr, extract_metadata, corr = FALSE, cat = TRUE, noNA = TRUE)

### Combine metadata
metadata.corr <- do.call("rbind", metadata.corr)
metadata.nocorr <- do.call("rbind", metadata.nocorr)

# metadata.cat.corr <- do.call("rbind", metadata.cat.corr)
# metadata.cat.nocorr <- do.call("rbind", metadata.cat.nocorr)

### Combine the metadata
# metadata.comb <- rbind(metadata.corr, metadata.nocorr, metadata.cat.corr, metadata.cat.nocorr)
metadata.comb <- rbind(metadata.corr, metadata.nocorr)

### Save into a few formats
saveRDS(metadata.comb, "data/heuristic_shrinkage_sim_metadata.rds")
write.table(metadata.comb, "data/heuristic_shrinkage_sim_metadata.txt", sep = "\t")
write.table(metadata.comb, "data/heuristic_shrinkage_sim_metadata.csv", sep = ",")

### The categorical scenarios are problematic. 
### Most fail after a small number of simulation runs with the following error message:
### 'contrasts can be applied only to factors with 2 or more levels'
### I think when the sample size is small, and we apply bootstrapping, 
### we end up with a dataset with no individuals with one of the levels of one of the predictors
### So maybe increase the minimum sample size?

### Fixed! Need to re-run

