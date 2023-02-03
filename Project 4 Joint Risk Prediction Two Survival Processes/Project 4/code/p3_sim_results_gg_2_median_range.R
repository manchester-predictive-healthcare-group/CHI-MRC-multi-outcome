###########
### This program produces the horizontal plots of median calibration of each method, and percentile range of each method,
### panelled by DGM for a given scenario and sample size
###########

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Source packages
source("code/sim_function_load_packages.R")


### Create vectors to loop through
# scen.vec <- c("s2.2")
# n.vec <- c(5000)

scen.vec <- c("s1.1", "s1.2", "s2.1", "s2.2")
n.vec <- c(1000, 2500, 5000)

print(scen.vec)
print(n.vec)

### Loop through vectors and create plots for each scenario
for (i in 1:length(scen.vec)){
  for (j in 1:length(n.vec)){

    ## load workspace
    load(paste("data/sim_results_", scen.vec[i], "_n", n.vec[j], "v1000.RData", sep = ""))
    ### Load results functions
    source("code/sim_function_results.R")
    ### Source the most updated input parameters for plots
    source(paste("code/sim_function_calibration_input_parameters_", scen.vec[i], ".R", sep =""))
    
    ########################
    ### Create  median plots
    ########################
    
    ### Define extra required parameters for median
    sum.vec <- "obs.median.diff"
    abline <- 0.15
    
    ### Create vector to store individual plots
    gg.median.list <- vector("list", 6)
    names(gg.median.list) <- c("DGM-1:MSM", "DGM-2:Clayton", "DGM-3:Gumbel", "DGM-4:Frank", "DGM-5:Normal", "DGM-6:Gamma")
    ## Create grid labels based off this vector also
    grid.labels.in  <- names(gg.median.list)
    
    ### Create each of the plots
    gg.median.list[["DGM-1:MSM"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.msm, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                     xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                     linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                     anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.median.list[["DGM-2:Clayton"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.clay, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                      xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                      linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                      anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.median.list[["DGM-3:Gumbel"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.gumb, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                      xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                      linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                      anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.median.list[["DGM-4:Frank"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.frank, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                       xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                       linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                       anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.median.list[["DGM-5:Normal"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.normal, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                        xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                        linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                        anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.median.list[["DGM-6:Gamma"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.gamma, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                       xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                       linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                       anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    ### Create the combined plot
    create.combined.calib.summary.ggplot.horiz(ggplots.summary.in = gg.median.list, 
                                               scenario = scen.vec[i], n.size = n.vec[j], suffix.to.save = "median")
    
    
    ########################
    ### Create  mean plots
    ########################
    
    ### Define extra required parameters for mean
    sum.vec <- "obs.mean.diff"
    abline <- 0.15
    
    ### Create vector to store individual plots
    gg.mean.list <- vector("list", 6)
    names(gg.mean.list) <- c("DGM-1:MSM", "DGM-2:Clayton", "DGM-3:Gumbel", "DGM-4:Frank", "DGM-5:Normal", "DGM-6:Gamma")
    ## Create grid labels based off this vector also
    grid.labels.in  <- names(gg.mean.list)
    
    ### Create each of the plots
    gg.mean.list[["DGM-1:MSM"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.msm, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                   xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                   linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                   anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.mean.list[["DGM-2:Clayton"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.clay, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                    xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                    linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                    anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.mean.list[["DGM-3:Gumbel"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.gumb, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                    xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                    linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                    anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.mean.list[["DGM-4:Frank"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.frank, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                     xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                     linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                     anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.mean.list[["DGM-5:Normal"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.normal, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                      xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                      linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                      anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.mean.list[["DGM-6:Gamma"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.gamma, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                     xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                     linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                     anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    ### Create the combined plot
    create.combined.calib.summary.ggplot.horiz(ggplots.summary.in = gg.mean.list, 
                                               scenario = scen.vec[i], n.size = n.vec[j], suffix.to.save = "mean")
    
    
    ########################
    ### Create p.range plots
    ########################
    
    ### Define extra required parameters for p.range
    sum.vec <- "obs.p.range"
    abline <- 0
    
    ### Create vector to store individual plots
    gg.p.range.list <- vector("list", 6)
    names(gg.p.range.list) <- c("DGM-1:MSM", "DGM-2:Clayton", "DGM-3:Gumbel", "DGM-4:Frank", "DGM-5:Normal", "DGM-6:Gamma")
    ## Create grid labels based off this vector also
    grid.labels.in  <- names(gg.median.list)
    
    
    ### Create each of the plots
    gg.p.range.list[["DGM-1:MSM"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.msm, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.p.range, 
                                                                      xlab.in = xlab, ylab.in = ylab.p.range, font.size.in = font.size, 
                                                                      linecolours.in = linecolours.horiz, lwd.in = lwd.horiz, abline.size.in = abline,
                                                                      anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.p.range.list[["DGM-2:Clayton"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.clay, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.p.range, 
                                                                       xlab.in = xlab, ylab.in = ylab.p.range, font.size.in = font.size, 
                                                                       linecolours.in = linecolours.horiz, lwd.in = lwd.horiz, abline.size.in = abline,
                                                                       anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.p.range.list[["DGM-3:Gumbel"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.gumb, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.p.range, 
                                                                       xlab.in = xlab, ylab.in = ylab.p.range, font.size.in = font.size, 
                                                                       linecolours.in = linecolours.horiz, lwd.in = lwd.horiz, abline.size.in = abline,
                                                                       anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.p.range.list[["DGM-4:Frank"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.frank, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.p.range, 
                                                                        xlab.in = xlab, ylab.in = ylab.p.range, font.size.in = font.size, 
                                                                        linecolours.in = linecolours.horiz, lwd.in = lwd.horiz, abline.size.in = abline,
                                                                        anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.p.range.list[["DGM-5:Normal"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.normal, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.p.range, 
                                                                         xlab.in = xlab, ylab.in = ylab.p.range, font.size.in = font.size, 
                                                                         linecolours.in = linecolours.horiz, lwd.in = lwd.horiz, abline.size.in = abline,
                                                                         anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    gg.p.range.list[["DGM-6:Gamma"]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.gamma, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.p.range, 
                                                                        xlab.in = xlab, ylab.in = ylab.p.range, font.size.in = font.size, 
                                                                        linecolours.in = linecolours.horiz, lwd.in = lwd.horiz, abline.size.in = abline,
                                                                        anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    ### Create the combined plot
    create.combined.calib.summary.ggplot.horiz(ggplots.summary.in = gg.p.range.list, 
                                               scenario = scen.vec[i], n.size = n.vec[j], suffix.to.save = "p.range")
    
    
    ## clear workspace
    rm(gg.p.range.list, gg.median.list, res.DGM.clay, res.DGM.gumb, res.DGM.frank, res.DGM.msm, res.DGM.normal, res.DGM.gamma)
  }
}

warnings()
