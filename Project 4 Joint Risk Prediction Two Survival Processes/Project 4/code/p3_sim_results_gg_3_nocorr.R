###########
### This program produces the horizontal plots of median calibration of each method, and percentile range of each method,
### but for the no correlation scenario, panelled by sample size
###########

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project 4/")

### Source packages
source("code/sim_function_load_packages.R")


#######################################################
#######################################################
### Create a plot for the no correlation scenarios
#######################################################
#######################################################



### Must write a seperate loop for scenarios s1 and s2

### Create vectors to store individual plots
gg.median.nocorr.list <- vector("list", 6)
names(gg.median.nocorr.list) <- c("LN, n = 1000", "LN, n = 2500", "LN, n = 5000", "HN, n = 1000", "HN, n = 2500", "HN, n = 5000")
## Create grid labels based off this vector also
grid.labels.in  <- names(gg.median.nocorr.list)


gg.p.range.nocorr.list <- vector("list", 6)
names(gg.p.range.nocorr.list) <- c("LN, n = 1000", "LN, n = 2500", "LN, n = 5000", "HN, n = 1000", "HN, n = 2500", "HN, n = 5000")
## Create grid labels based off this vector also
grid.labels.in  <- names(gg.p.range.nocorr.list)

### Create vectors to loop through
#scen.vec <- c("s2")
#n.vec <- c(1000, 2500, 5000)

scen.vec <- c("s1", "s2")
n.vec <- c(1000, 2500, 5000)

### Loop through vectors and create plots
k <- 1
for (i in 1:length(scen.vec)){
  for (j in 1:length(n.vec)){
    print(scen.vec[i])
    print(n.vec[j])
    ## load workspace
    load(paste("data/sim_results_", scen.vec[i], "_n", n.vec[j], "v1000.RData", sep = ""))
    ### Load results functions
    source("code/sim_function_results.R")
    ### Source the most updated input parameters for plots
    source(paste("code/sim_function_calibration_input_parameters_", scen.vec[i], ".R", sep =""))
    
    ## Create grid labels
    grid.labels.in  <- names(gg.median.nocorr.list)
    
    
    ########################
    ### Create  median plots
    ########################
    
    ### Define extra required parameters for median
    sum.vec <- "obs.median.diff"
    abline <- 0.15
    
    ### Create ggplot
    gg.median.nocorr.list[[k]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.nocorr, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.med, 
                                                                    xlab.in = xlab, ylab.in = ylab.med, font.size.in = font.size, 
                                                                    linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                    anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    ########################
    ### Createp.range plots
    ########################
    
    ### Define extra required parameters for median
    sum.vec <- "obs.p.range"
    abline <- 0
    
    ### Create ggplot
    gg.p.range.nocorr.list[[k]] <- create.summary.horizontal.ggplots(analysed.data = res.DGM.nocorr, x.lim.in = x.lim.horiz, y.lim.in = y.lim.horiz.p.range, 
                                                                     xlab.in = xlab, ylab.in = ylab.p.range, font.size.in = font.size, 
                                                                     linecolours.in = linecolours.horiz, lwd.in = lwd.horiz,  abline.size.in = abline,
                                                                     anal.methods.sim.inn = anal.methods.sim, sum.vec.in = sum.vec)
    
    ## Add one to counter
    k <- k+1
    
    ## clear workspace
    rm(res.DGM.nocorr)
  }
}
rm(k)

## Combine into single ggplot
gg.median.nocorr.comb <- ggarrange(gg.median.nocorr.list[[1]], gg.median.nocorr.list[[2]], gg.median.nocorr.list[[3]], 
                                   gg.median.nocorr.list[[4]], gg.median.nocorr.list[[5]], gg.median.nocorr.list[[6]],
                                   labels = grid.labels.in, nrow = 2, ncol = 3, 
                                   font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in, 
                                   common.legend = TRUE, legend = "right")

gg.p.range.nocorr.comb <- ggarrange(gg.p.range.nocorr.list[[1]], gg.p.range.nocorr.list[[2]], gg.p.range.nocorr.list[[3]], 
                                    gg.p.range.nocorr.list[[4]], gg.p.range.nocorr.list[[5]], gg.p.range.nocorr.list[[6]],
                                    labels = grid.labels.in, nrow = 2, ncol = 3, 
                                    font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in, 
                                    common.legend = TRUE, legend = "right")

## Save combined plots
ggsave(paste("figures/gg_3.summary.nocorr.median.tiff", sep = ""), 
       gg.median.nocorr.comb, device = "tiff", 
       dpi = 800, width = 10, height = 4.62, units = "in", compression = 'lzw')

ggsave(paste("figures/gg_3.summary.nocorr.p.range.tiff", sep = ""), 
       gg.p.range.nocorr.comb, device = "tiff", 
       dpi = 800, width = 10, height = 4.62, units = "in", compression = 'lzw')

ggsave(paste("figures/gg_3.summary.nocorr.median.png", sep = ""), 
       gg.median.nocorr.comb, device = "png", 
       dpi = 600, width = 10, height = 4.62, units = "in")

ggsave(paste("figures/gg_3.summary.nocorr.p.range.png", sep = ""), 
       gg.p.range.nocorr.comb, device = "png", 
       dpi = 600, width = 10, height = 4.62, units = "in")

ggsave(paste("figures/gg_3.summary.nocorr.median.pdf", sep = ""), 
       gg.median.nocorr.comb, device = "pdf", 
       dpi = 600, width = 10, height = 4.62)

ggsave(paste("figures/gg_3.summary.nocorr.p.range.pdf", sep = ""), 
       gg.p.range.nocorr.comb, device = "pdf", 
       dpi = 600, width = 10, height = 4.62)

