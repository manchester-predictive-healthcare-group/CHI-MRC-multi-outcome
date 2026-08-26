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
extract_metadata <- function(file, noNA){
 
  # file <- testfile
  # noNA <- TRUE

  ### Create object to store plot data
  metadata <- NULL
  
  ### Load workspace
  print(file)
  load(paste("data/", file, sep = ""))
  
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
    vec.list.input <- c("P" = input.data[[sim.iter]][["P"]],
                        "pavlou.no.converge" = ifelse(class(input.data[[sim.iter]]$nreq.pavlou) != "try-error", 
                                                      ifelse(input.data[[sim.iter]]$nreq.pavlou != "Too many iterations", "converged", "too many iterations"), 
                                                      "try_error"),
                        "nreq.pavlou" = ifelse(class(input.data[[sim.iter]]$nreq.pavlou) != "try-error" & input.data[[sim.iter]]$nreq.pavlou != "Too many iterations", input.data[[sim.iter]][["nreq.pavlou"]], NA),
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
  for (col in c(1:7, 9:15)){
    metadata[,col] <- as.numeric(metadata[,col])
  }
  
  ## Add variable for correlation
  if (stringr::str_detect(file, "nocorr") == FALSE){
    metadata$corr <- TRUE
  } else {
    metadata$corr <- FALSE
  }
  
  return(metadata)
}

## Find names of files for which simulations have run
filenames <- list.files("data")
filenames <- filenames[stringr::str_detect(filenames, "run_sim_study4")]

# filenames <- filenames[!stringr::str_detect(filenames, "run_sim_study4_corr1746")]
# filenames <- filenames[!stringr::str_detect(filenames, "run_sim_study4_corr1748")]
# filenames <- filenames[!stringr::str_detect(filenames, "run_sim_study4_corr1559")]
# filenames <- filenames[!stringr::str_detect(filenames, "run_sim_study4_nocorr1581")]
# filenames <- filenames[!stringr::str_detect(filenames, "run_sim_study4_nocorr1591")]

### Extract metadata
metadata <- lapply(filenames, extract_metadata, noNA = TRUE)

### Combine metadata
metadata.comb <- do.call("rbind", metadata)
str(metadata.comb)

### Calculate EPV
metadata.comb$EPP.riley <- metadata.comb$prop*metadata.comb$nreq.riley/metadata.comb$P
metadata.comb$EPP.pavlou <- metadata.comb$prop*metadata.comb$nreq.pavlou/metadata.comb$P

### Count rows with negative meanR2 and meanS.VH or meanS.boot, driven by a small number of obscure results over the simulation iterations
metadata.remove.pavlou <- subset(metadata.comb, S.pop.mean.pavlou < 0 | S.pop.mean.pavlou > 2)
nrow(metadata.remove.pavlou)
metadata.remove.riley <- subset(metadata.comb, S.pop.mean.riley < 0 | S.pop.mean.riley > 2)
nrow(metadata.remove.riley)

### Remove these prior to saving
metadata.comb <- subset(metadata.comb, !(!is.na(S.pop.mean.riley) & (S.pop.mean.riley < 0 | S.pop.mean.riley > 2)))
metadata.comb <- subset(metadata.comb, !(!is.na(S.pop.mean.pavlou) & (S.pop.mean.pavlou < 0 | S.pop.mean.pavlou > 2)))

### Save into a few formats
saveRDS(metadata.comb, "data/sim_study4_metadata.rds")
write.table(metadata.comb, "data/sim_study4_metadata.txt", sep = "\t")
write.table(metadata.comb, "data/sim_study4_metadata.csv", sep = ",")
