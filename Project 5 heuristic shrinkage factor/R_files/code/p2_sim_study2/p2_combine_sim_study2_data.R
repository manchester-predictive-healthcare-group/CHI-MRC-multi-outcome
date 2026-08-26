### Clear workspace
rm(list=ls())

### Set working directory
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_8.2/")

### Source functions
R.func.sources <- list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load functions
source("R/sim_functions.R")

### Write a function to extract the metadata
### The corr argument defines whether data was generated with correlation or not
### noNA = TRUE will remove all NA's before calculating mean or median
### Note, the Pavlou formula may have been unsuccessful, in which case we assign NA, as the simulation was not ran.
### Note, if the Riley sample size was bigger than 10,000, we assign NA, as the simulation was not ran.
extract_metadata <- function(set, corr, noNA){

  ### Create object to store plot data
  metadata <- NULL
  
  print(paste("SET", set))
  ### Load workspace
  if (corr == TRUE){
    load(paste("data/run_sim_validate_ssd_corr", set, ".RData", sep = ""))
  } else {
    load(paste("data/run_sim_validate_ssd_nocorr", set, ".RData", sep = ""))
  }
  
  
  # class(input.data[[8]]$nreq.pavlou)
  # class(input.data[[7]]$nreq.pavlou)
  # input.data[[7]]
  # input.data[[8]]
  # output.riley.S.VH[[8]]
  
  ### First start by getting the number of complete simulations
  n.comp <- max(which(sapply(output.riley.S.pop, is.null) == FALSE))
  
  ### Create a list to store the vectors in which we will combine
  vec.list <- vector("list", n.comp)

  ### Create each vector
  for (sim.iter in 1:n.comp){
    
    ### Vector for input data
    vec.list.input <- c("P.meas" = input.data[[sim.iter]][["P"]],
                        "P.unmeas" = input.data[[sim.iter]][["P.total"]] - input.data[[sim.iter]][["P"]],
                        "nreq.pavlou" = ifelse(class(input.data[[sim.iter]]$nreq.pavlou) != "try-error", input.data[[sim.iter]][["nreq.pavlou"]], NA),
                        "nreq.riley" = input.data[[sim.iter]][["nreq.riley"]],
                        "prop" = input.data[[sim.iter]][["prop"]], 
                        "C.pop" = input.data[[sim.iter]][["C.pop"]], 
                        "R2.CS.pop" = input.data[[sim.iter]][["R2.CS.pop"]], 
                        set = set, 
                        sim.iter = sim.iter)
    
    ### Create vector of outputs from the simulation
    if (noNA == TRUE){
      ### Create vector for output from Pavlou sample size
      if (is.null(output.pavlou.S.pop[[sim.iter]])){
        vec.list.pavlou <- c("S.pop.mean.pavlou" = NA,
                             "S.pop.sd.pavlou" = NA,
                             "S.pop.median.pavlou" = NA)
      } else {
        vec.list.pavlou <- c("S.pop.mean.pavlou" = mean(output.pavlou.S.pop[[sim.iter]], na.rm = TRUE),
                             "S.pop.sd.pavlou" = sd(output.pavlou.S.pop[[sim.iter]], na.rm = TRUE),
                             "S.pop.median.pavlou" = median(output.pavlou.S.pop[[sim.iter]], na.rm = TRUE))
      }
      
      ### Vector for output from Riley sample size
      if (is.null(output.riley.S.pop[[sim.iter]])){
        vec.list.riley <- c("S.pop.mean.riley" = NA,
                            "S.pop.sd.riley" = NA,
                            "S.pop.median.riley" = NA)
      } else {
        vec.list.riley <- c("S.pop.mean.riley" = mean(output.riley.S.pop[[sim.iter]], na.rm = TRUE),
                            "S.pop.sd.riley" = sd(output.riley.S.pop[[sim.iter]], na.rm = TRUE),
                            "S.pop.median.riley" = median(output.riley.S.pop[[sim.iter]], na.rm = TRUE))
      }
      
    } else {
      ### Create vector for output from Pavlou sample size
      if (class(input.data[[sim.iter]]$nreq.pavlou) == "try-error"){
        vec.list.pavlou <- c("S.pop.mean.pavlou" = NA,
                             "S.pop.sd.pavlou" = NA,
                             "S.pop.median.pavlou" = NA)
      } else {
        vec.list.pavlou <- c( "S.pop.mean.pavlou" = mean(output.pavlou.S.pop[[sim.iter]]),
                             "S.pop.sd.pavlou" = sd(output.pavlou.S.pop[[sim.iter]]),
                             "S.pop.median.pavlou" = median(output.pavlou.S.pop[[sim.iter]]))
      }
      
      ### Vector for output from Riley sample size
      if (is.null(output.riley.S.VH[[sim.iter]])){
        vec.list.riley <- c("S.pop.mean.riley" = NA,
                            "S.pop.sd.riley" = NA,
                            "S.pop.median.riley" = NA)
      } else {
        vec.list.riley <- c("S.pop.mean.riley" = mean(output.riley.S.pop[[sim.iter]]),
                            "S.pop.sd.riley" = sd(output.riley.S.pop[[sim.iter]]),
                            "S.pop.median.riley" = median(output.riley.S.pop[[sim.iter]]))
      }
      
    }
    
    ### Assign vector of outputs to be saved
    vec.list[[sim.iter]] <- c(vec.list.pavlou, vec.list.riley, vec.list.input)
    
  }
  
  ### Create temp plot data and turn into a list, so can be combined with existing list
  metadata <- data.frame(do.call("rbind", vec.list))
  
  ## Add variable for correlation
  if (corr == TRUE){
    metadata$corr <- TRUE
  } else {
    metadata$corr <- FALSE
  }
  
  return(metadata)
}

### Find numbers for which simulations have run
library(stringr)
valid.files.corr <- as.numeric(str_extract(list.files("data"), "(?<=run_sim_validate_ssd_corr)[0-9]*")[!is.na(str_extract(list.files("data"), "(?<=run_sim_validate_ssd_corr)[0-9]*"))])
valid.files.corr <- valid.files.corr[!is.na(valid.files.corr)]

valid.files.nocorr <- as.numeric(str_extract(list.files("data"), "(?<=run_sim_validate_ssd_nocorr)[0-9]*")[!is.na(str_extract(list.files("data"), "(?<=run_sim_validate_ssd_nocorr)[0-9]*"))])
valid.files.nocorr <- valid.files.nocorr[!is.na(valid.files.nocorr)]

### Extract metadata
metadata.corr <- lapply(valid.files.corr, extract_metadata, corr = TRUE, noNA = TRUE)
metadata.nocorr <- lapply(valid.files.nocorr, extract_metadata, corr = FALSE, noNA = TRUE)

### Combine metadata
metadata.corr <- do.call("rbind", metadata.corr[unlist(lapply(metadata.corr, function(x) {!is.null(x)}))])
metadata.nocorr <- do.call("rbind", metadata.nocorr[unlist(lapply(metadata.nocorr, function(x) {!is.null(x)}))])

### Combine the metadata
metadata.comb <- rbind(metadata.corr, metadata.nocorr)
# metadata.comb <- metadata.corr

### Save into a few formats
saveRDS(metadata.comb, "data/heuristic_shrinkage_sim_metadata_ssd.rds")
write.table(metadata.comb, "data/heuristic_shrinkage_sim_metadata_ssd.txt", sep = "\t")
write.table(metadata.comb, "data/heuristic_shrinkage_sim_metadata_ssd.csv", sep = ",")

###
### After initially running this simulation, the samplesizedev packae was updated and convergence was improved.
### I therefore re-ran all the scenarios where Pavlou approach failed to converge. Read these in now.
###

### Write a function to extract the metadata
### The corr argument defines whether data was generated with correlation or not
### noNA = TRUE will remove all NA's before calculating mean or median
### Note, the Pavlou formula may have been unsuccessful, in which case we assign NA, as the simulation was not ran.
### Note, if the Riley sample size was bigger than 10,000, we assign NA, as the simulation was not ran.
extract_metadata_rerun <- function(row, noNA){

  ### Create object to store plot data
  metadata <- NULL
  
  print(paste("SET", row))
  ### Load workspace
  load(paste("data/run_sim_validate_ssd_", row, "_rerun_failed_scenarios.RData", sep = ""))

  # class(input.data[[8]]$nreq.pavlou)
  # class(input.data[[7]]$nreq.pavlou)
  # input.data[[7]]
  # input.data[[8]]
  # output.riley.S.VH[[8]]
  
  # ### First start by getting the number of complete simulations
  # n.comp <- max(which(sapply(output.riley.S.pop, is.null) == FALSE))
  # 
  # ### Create a list to store the vectors in which we will combine
  # vec.list <- vector("list", n.comp)

    
    ### Vector for input data
    vec.list.input <- c("P.meas" = input.data[[1]][["P"]],
                        "P.unmeas" = input.data[[1]][["P.total"]] - input.data[[1]][["P"]],
                        "nreq.pavlou" = ifelse(class(input.data[[1]]$nreq.pavlou) != "try-error" & input.data[[1]]$nreq.pavlou != "Too many iterations", input.data[[1]][["nreq.pavlou"]], NA),
                        "nreq.riley" = input.data[[1]][["nreq.riley"]],
                        "prop" = input.data[[1]][["prop"]], 
                        "C.pop" = input.data[[1]][["C.pop"]], 
                        "R2.CS.pop" = input.data[[1]][["R2.CS.pop"]], 
                        set = set, 
                        sim.iter = sim)
    
    ### Create vector of outputs from the simulation
    if (noNA == TRUE){
      ### Create vector for output from Pavlou sample size
      if (is.null(output.pavlou.S.pop[[1]])){
        vec.list.pavlou <- c("S.pop.mean.pavlou" = NA,
                             "S.pop.sd.pavlou" = NA,
                             "S.pop.median.pavlou" = NA)
      } else {
        vec.list.pavlou <- c("S.pop.mean.pavlou" = mean(output.pavlou.S.pop[[1]], na.rm = TRUE),
                             "S.pop.sd.pavlou" = sd(output.pavlou.S.pop[[1]], na.rm = TRUE),
                             "S.pop.median.pavlou" = median(output.pavlou.S.pop[[1]], na.rm = TRUE))
      }
      
      ### Vector for output from Riley sample size
      if (is.null(output.riley.S.pop[[1]])){
        vec.list.riley <- c("S.pop.mean.riley" = NA,
                            "S.pop.sd.riley" = NA,
                            "S.pop.median.riley" = NA)
      } else {
        vec.list.riley <- c("S.pop.mean.riley" = mean(output.riley.S.pop[[1]], na.rm = TRUE),
                            "S.pop.sd.riley" = sd(output.riley.S.pop[[1]], na.rm = TRUE),
                            "S.pop.median.riley" = median(output.riley.S.pop[[1]], na.rm = TRUE))
      }
      
    } else {
      ### Create vector for output from Pavlou sample size
      if (class(input.data[[1]]$nreq.pavlou) == "try-error"){
        vec.list.pavlou <- c("S.pop.mean.pavlou" = NA,
                             "S.pop.sd.pavlou" = NA,
                             "S.pop.median.pavlou" = NA)
      } else {
        vec.list.pavlou <- c( "S.pop.mean.pavlou" = mean(output.pavlou.S.pop[[1]]),
                              "S.pop.sd.pavlou" = sd(output.pavlou.S.pop[[1]]),
                              "S.pop.median.pavlou" = median(output.pavlou.S.pop[[1]]))
      }
      
      ### Vector for output from Riley sample size
      if (is.null(output.riley.S.VH[[1]])){
        vec.list.riley <- c("S.pop.mean.riley" = NA,
                            "S.pop.sd.riley" = NA,
                            "S.pop.median.riley" = NA)
      } else {
        vec.list.riley <- c("S.pop.mean.riley" = mean(output.riley.S.pop[[1]]),
                            "S.pop.sd.riley" = sd(output.riley.S.pop[[1]]),
                            "S.pop.median.riley" = median(output.riley.S.pop[[1]]))
      }
      
    }
    
    ### Assign vector of outputs to be saved
    vec.list <- c(vec.list.pavlou, vec.list.riley, vec.list.input)
  
    # ### Format
    # for (col in c(1:7, 9:15)){
    #   metadata[col] <- as.numeric(metadata[col])
    # }
    
  ## Add variable for correlation
  if (corr == TRUE){
    vec.list[["corr"]] <- TRUE
  } else {
    vec.list[["corr"]] <- FALSE
  }
    
  return(vec.list)
    
}

### Extract metadata
metadata.rerun <- lapply(1:3686, extract_metadata_rerun, noNA = TRUE)

### Combine metadata
metadata.rerun <- do.call("rbind", metadata.rerun)
metadata.rerun <- data.frame(metadata.rerun)
metadata.rerun <- dplyr::mutate(metadata.rerun, corr = dplyr::case_when(corr == 1 ~ TRUE,
                                                                        corr == 0 ~ FALSE))

### Read in metadata from original run of the simulation 
### (note this was created above, but I am re-reading this in, so code can be easily seperated)
### Remove scenarios that were re-run (where Pavlou approach failed and Riley was < 10000)
metadata.original <- readRDS("data/heuristic_shrinkage_sim_metadata_ssd.rds")
metadata.original <- subset(metadata.original, !(is.na(S.pop.mean.pavlou) & !is.na(S.pop.mean.riley)))

### Combine new scenarios that were re-run with original metadata
metadata.both <- rbind(metadata.original, metadata.rerun)

### Check they have same number of rows
metadata.original <- readRDS("data/heuristic_shrinkage_sim_metadata_ssd.rds")
testthat::expect_equal(nrow(metadata.both), nrow(metadata.original))

### Count rows with negative meanR2 and meanS.VH or meanS.boot, driven by a small number of obscure results over the simulation iterations
metadata.remove.pavlou <- subset(metadata.both, S.pop.mean.pavlou < 0 | S.pop.mean.pavlou > 2)
nrow(metadata.remove.pavlou)
metadata.remove.riley <- subset(metadata.both, S.pop.mean.riley < 0 | S.pop.mean.riley > 2)
nrow(metadata.remove.riley)

### Remove these prior to saving
metadata.both <- subset(metadata.both, !(!is.na(S.pop.mean.riley) & (S.pop.mean.riley < 0 | S.pop.mean.riley > 2)))
metadata.both <- subset(metadata.both, !(!is.na(S.pop.mean.pavlou) & (S.pop.mean.pavlou < 0 | S.pop.mean.pavlou > 2)))

### Save into a few formats
saveRDS(metadata.both, "data/sim_study2_metadata.rds")
write.table(metadata.both, "data/sim_study2_metadata.txt", sep = "\t")
write.table(metadata.both, "data/sim_study2_metadata.csv", sep = ",")
