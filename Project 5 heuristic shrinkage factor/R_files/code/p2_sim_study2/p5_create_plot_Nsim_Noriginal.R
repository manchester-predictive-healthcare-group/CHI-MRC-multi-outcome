###
### Explore sample sizes
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
metadata <- readRDS("data/sim_study2_metadata.rds")
metadata <- subset(metadata, nreq.riley < 10000 & !is.na(nreq.pavlou) & corr == TRUE)

### Create plot S.VH
# fun_color_range <- colorRampPalette(c("green", "yellow", "red", "magenta", "blue", "cyan"))
# my_colors <- fun_color_range(60)
my_colors <- colorspace::rainbow_hcl(50, c = 70, l = 70)

gg.out <- ggplot(aes(x = nreq.riley, y = nreq.pavlou, color = C.pop), data = metadata) +
  scale_colour_gradientn(colors = my_colors) + xlim(c(0,2500)) + ylim(c(0,2500))
gg.out <- gg.out + geom_point(size = 0.5) + geom_abline(intercept = 0, slope = 1) +
  xlab(expression(N[original])) + ylab(expression(N[sim])) +
  labs(color = expression(C[pop]))

Cairo::CairoPNG("figures/ss2/ss2_sensitivity_compare_Nsim_Noriginal.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(gg.out)
dev.off()
gg.out
