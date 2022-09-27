### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_1.4/")

### Load packages
source("code/sim_load_packages.R")

### Load functions
source("code/sim_functions.R")
source("code/sim_functions_results.R")

### Load generic input parameters
source("code/sim_load_generic_input_parameters.R")

### Size of development datasets
N.devel.sim <- 500000

### Cycle through scenario, seed.coef.sim and K.sim
for (K.sim in c(3, 5)){
  for (seed.coef.sim in c(1, 2, 3)){
    
    ### Create lists to store data in for each scenario, which will then be concatenated
    ggdata.DGMmult.list <- vector("list", 4)
    ggdata.DGMseqlog.list <- vector("list", 4)
    
    ### Loop through scenarios and create datasets
    for (scenario.sim in c(1,2,3,4)){
       
      ### Load .RData file for DGM mult
      load(paste("data/s", scenario.sim, "/sim_run_large_sample_DGMmult_scen" , scenario.sim, ".", seed.coef.sim, "_K", 
                 K.sim, "_P", P.sim, "_N", N.devel.sim, ".RData", sep =""))
      
      ### Load functions
      source("code/sim_functions_results.R")
      
      ### Create data for ggplot
      ggdata.DGMmult.list[[scenario.sim]] <- create.large.sample.ggplot.data.DGMmult(sim.out.binary = sim.out.binary, 
                                                                                     sim.out.multinomial = sim.out.multinomial, 
                                                                                     sim.out.binary.rcs3 = sim.out.binary.rcs3, 
                                                                                     sim.out.binary.rcs5 = sim.out.binary.rcs5,
                                                                                     scenario.sim = scenario.sim,
                                                                                     seed.coef.sim = seed.coef.sim)
      
      ### Load .RData file for DGM seqlog
      load(paste("data/s", scenario.sim, "/sim_run_large_sample_DGMseqlog_scen" , scenario.sim, ".", seed.coef.sim, "_K", 
                 K.sim, "_P", P.sim, "_N", N.devel.sim, ".RData", sep =""))
      
      ### Load functions
      source("code/sim_functions_results.R")
      
      ### Create data for ggplot
      ggdata.DGMseqlog.list[[scenario.sim]] <- create.large.sample.ggplot.data.DGMseqlog(sim.out.multinomial = sim.out.multinomial, 
                                                                                         sim.out.binary = sim.out.binary, 
                                                                                         sim.out.multinomial.rcs3 = sim.out.multinomial.rcs3, 
                                                                                         sim.out.multinomial.rcs5 = sim.out.multinomial.rcs5, 
                                                                                         scenario.sim = scenario.sim,
                                                                                         seed.coef.sim = seed.coef.sim)
      
    }
    
    ### Combine the data into a single dataset (seperately for DGMmult and DGMseqlog)
    gg.data.DGMmult <- do.call("rbind", ggdata.DGMmult.list)
    gg.data.DGMseqlog <- do.call("rbind", ggdata.DGMseqlog.list)
    
    ### Create ggplot
    ggplot.DGMmult <- create.large.sample.ggplot(data.in = gg.data.DGMmult, font.size = 10, 
                                                 x.lim.in = c(0,1), y.lim.in = c(0,1), 
                                                 xlab.in = "Observed risk", ylab.in = "Predicted risk")
    
    ### Save ggplot
    ggsave(paste("figures/gg_large_samp_DGMmult_K", K.sim, "_P", P.sim, "_seed", seed.coef.sim, ".png", sep =""), ggplot.DGMmult,
           dpi = 300)
    
    ### Create ggplot
    ggplot.DGMseqlog <- create.large.sample.ggplot(data.in = gg.data.DGMseqlog, font.size = 10, 
                                                   x.lim.in = c(0,1), y.lim.in = c(0,1), 
                                                   xlab.in = "Observed risk", ylab.in = "Predicted risk")
    
    ### Save ggplot
    ggsave(paste("figures/gg_large_samp_DGMseqlog_K", K.sim, "_P", P.sim, "_seed", seed.coef.sim, ".png", sep =""), ggplot.DGMseqlog,
           dpi = 300)
  }
}



