### This code is for the large sample validation, so the models are developed and calibrated on the same dataset

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.5/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions.R")
source("code/sim_functions_results.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")


### Load and save figures
for (seed.coef.sim in c(1)){
  for (K.sim in c(5)){
    for (scenario.sim in c(1)){
    Sys.time()
    print(paste("K.sim = ", K.sim, " scenario = ", scenario.sim, " seed.coef.sim = ", seed.coef.sim, sep = ""))
      
    ## Load data
    load(paste("data/sim_run_large_sample_s", scenario.sim, ".", seed.coef.sim, "_K", K.sim, ".RData", sep = ""))

    ## Load updated results functions
    source("code/sim_functions_results.R")
    
    ## Create plots for each DGM using mlr recalibration framework
    sim.res.plot.flex.mlr <- vector("list", 2)
    names(sim.res.plot.flex.mlr) <- c("DGM.multinomial", "DGM.seqlog")
    for (i in 1:2){
      sim.res.plot.flex.mlr[[i]] <- create.plot.flex.mlr(sim.out.list[[i]])
      print(paste("plot mlr ", i))
    }
    
    
    ## Create plots for each DGM using binary logistic recalibration framework
    sim.res.plot.flex.po <- vector("list", 2)
    names(sim.res.plot.flex.po) <- c("DGM.multinomial", "DGM.seqlog")
    for (i in 1:2){
      sim.res.plot.flex.po[[i]] <- create.plot.flex.po(sim.out.list[[i]])
      print(paste("plot po ", i))
    }
    
    ## Save plots
    ggsave(paste("figures/gg_large_sample_flex_mlr_s", c("A", "B", "C", "D")[as.numeric(scenario.sim)], ".", seed.coef.sim, "_K", K.sim, "_DGMmult.png", sep = ""), sim.res.plot.flex.mlr[[1]], dpi = 300)
    ggsave(paste("figures/gg_large_sample_flex_mlr_s", c("A", "B", "C", "D")[as.numeric(scenario.sim)], ".", seed.coef.sim, "_K", K.sim, "_DGMseqlog.png", sep = ""), sim.res.plot.flex.mlr[[2]], dpi = 300)
    
    ggsave(paste("figures/gg_large_sample_flex_po_s", c("A", "B", "C", "D")[as.numeric(scenario.sim)], ".", seed.coef.sim, "_K", K.sim, "_DGMmult.png", sep = ""), sim.res.plot.flex.po[[1]], dpi = 300)
    ggsave(paste("figures/gg_large_sample_flex_po_s", c("A", "B", "C", "D")[as.numeric(scenario.sim)], ".", seed.coef.sim, "_K", K.sim, "_DGMseqlog.png", sep = ""), sim.res.plot.flex.po[[2]], dpi = 300)
    }
  }
}


