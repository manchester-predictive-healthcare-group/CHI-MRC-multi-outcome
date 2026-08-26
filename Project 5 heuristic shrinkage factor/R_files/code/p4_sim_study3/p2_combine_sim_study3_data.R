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
extract_metadata <- function(file, corr, cat, noNA){
  
  # file <- filenames[10]
  # file <- "run_sim_study3_corr_set10.RData"
  # corr <- TRUE
  # cat <- TRUE
  # noNA <- TRUE
  print(file)
  
  ### Create object to store plot data
  metadata <- NULL
  
  ### Load workspace
  load(paste("data/", file, sep = ""))
  
  ### First start by getting the number of complete simulations
  n.comp <- max(which(sapply(output.S.VH, is.null) == FALSE))
  
  ### Only do the following if n.comp is not = -Inf, which occurs when there were no successful simulation runs in the .RData file
  ### This occurs when the first scenario in a given set results in an error (e.g. during the bootstrapping process).
  if (n.comp != -Inf){
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
                                  "C.target" = input.data[[sim.iter]][["C.target"]],
                                  "R2.CS.pop" = input.data[[sim.iter]][["R2.CS.pop"]],
                                  "P" = input.data[[sim.iter]][["P"]],
                                  "nreq" = input.data[[sim.iter]][["nreq"]],
                                  "prop" = input.data[[sim.iter]][["prop"]], 
                                  "scen.possible" = input.data[[sim.iter]][["scen.possible"]], 
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
                                  "C.target" = input.data[[sim.iter]][["C.target"]],
                                  "R2.CS.pop" = input.data[[sim.iter]][["R2.CS.pop"]],
                                  "P" = input.data[[sim.iter]][["P"]],
                                  "nreq" = input.data[[sim.iter]][["nreq"]],
                                  "prop" = input.data[[sim.iter]][["prop"]], 
                                  "scen.possible" = input.data[[sim.iter]][["scen.possible"]], 
                                  set = set, 
                                  sim.iter = sim.iter)
      }
    }
    
    ### Create temp plot data and turn into a list, so can be combined with existing list
    metadata <- data.frame(do.call("rbind", vec.list))
    
    ### Remove rows where we could not find coefficients to match the desired C-statistic within 5 minutes of waiting (i.e., randomly simulated
    ### scenario parameters were not possible)
    metadata <- metadata[metadata$scen.possible == 1, ]
    
    ### Remove rows with negative LR, driven by a small number of obscure results
    metadata <- metadata[metadata$R2.CS.app.mean > 0, ]
    metadata <- metadata[metadata$S.VH.mean > 0, ]
    
    ## Add variable for correlation
    if (stringr::str_detect(file, "nocorr") == FALSE){
      metadata$corr <- TRUE
    } else {
      metadata$corr <- FALSE
    }
    
  } else {
    
    metadata <- NULL
    
  }
  
  return(metadata)
  
}

### Find names of files for which simulations have run
filenames <- list.files("data")
filenames <- filenames[stringr::str_detect(filenames, "run_sim_study3")]

### Remove set = 26, because this started with a P = 28 and n = 102 scenario, which failed during boostrapping, and therefore contained
### no results, and is causing an error.
# filenames <- filenames[!grepl("run_sim_study3_corr_set26.RData", filenames)]
# filenames <- filenames[!grepl("run_sim_study3_corr_set426.RData", filenames)]

### Extract metadata
metadata <- lapply(filenames, extract_metadata, noNA = TRUE)

### Remove NULL elements (small number of sets of simulation runs where the first scenario resulted in an error during bootstrapping, resulting in
### no output in the .RData file)
which(sapply(metadata, is.null))

### Combine metadata
metadata.comb <- do.call("rbind", metadata)

### Calculate EPV
metadata.comb$EPP <- metadata.comb$prop*metadata.comb$nreq/metadata.comb$P

### Save into a few formats
saveRDS(metadata.comb, "data/sim_study3_metadata.rds")
write.table(metadata.comb, "data/sim_study3_metadata.txt", sep = "\t")
write.table(metadata.comb, "data/sim_study3_metadata.csv", sep = ",")



