###
### This file will produce plots to help establish why EPV is a confusing metric and we must explain this to the reader
###

### HOLD C-STATISTIC
### VARY EPP

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
metadata <- readRDS("data/sim_study1_metadata.rds")
metadata <- subset(metadata, C.pop != 0)



### Write a function which will combine data for a specific version
create_plot <- function(df, C1 = NULL, C2 = NULL, N1 = NULL, N2 = NULL, P = NULL, correlated = TRUE, cat = FALSE, xvar, yvar, byvar, bygroups = NULL, xrange = NULL, yrange = NULL, alpha = NULL){
  
  ### Subset to C-statistics
  if (!is.null(C1) & !is.null(C2)){
    df <- subset(df, C.pop < C2 & C.pop > C1)
  }
  
  if (!is.null(N1) & !is.null(N2)){
    df <- subset(df, nreq < N2 & nreq > N1)
  }
  
  if (!is.null(P)){
    df <- subset(df, P.meas == P)
  }
  
  ### Reduce to correlation or uncorrelated
  df <- subset(df, corr == as.numeric(correlated))
  
  ### Reduce to categorical or not categorical
  df <- subset(df, cat == as.numeric(cat))
  
  ### Subset if required
  if (!is.null(xrange)){
    df <- df[df[[xvar]] > xrange[1] & df[[xvar]] < xrange[2], ]
  }
  
  if (!is.null(yrange)){
    df <- df[df[[yvar]] > yrange[1] & df[[yvar]] < yrange[2], ]
  }
  
  ### If more than 5000 rows, reduce to 5000
  if (nrow(df) > 5000){
    df <- df[1:5000, ]
  }
  
  ### Create colour range for ggplot
  # fun_color_range <- colorRampPalette(c("green", "blue"))
  # my_colors <- fun_color_range(20)
  my_colors <- colorspace::rainbow_hcl(50, c = 70, l = 70)
  
  ### Assign label for axes and legend
  x.label <- NULL
  y.label <- NULL
  by.label <- NULL
  
  # x label
  if (xvar == "S.VH.mean"){
    x.label <- expression(paste("mean(", S[VH], ")", sep = ""))
  } else if (xvar == "S.boot.mean"){
    x.label <- expression(paste("mean(", S[boot], ")", sep = ""))
  } else if (xvar == "S.VH.median"){
    x.label <- expression(paste("median(", S[VH], ")", sep = ""))
  } else if (xvar == "S.boot.median"){
    x.label <- expression(paste("median(", S[boot], ")", sep = ""))
  }
  
  # y label
  if (yvar == "S.pop.mean"){
    y.label <- expression(paste("mean(", S[opt], ")", sep = ""))
  } else if (yvar == "S.pop.median"){
    y.label <- expression(paste("median(", S[opt], ")", sep = ""))
  }
  
  # by label
  if (byvar == "R2.CS.app.mean"){
    by.label <- expression(paste("mean(", R^2*phantom()[paste(CS,",",app, sep = "")], ")", sep = ""))
  } else if (byvar == "R2.CS.app.median"){
    by.label <- expression(paste("median(", R^2*phantom()[paste(CS,",",app, sep = "")], ")", sep = ""))
  } else if (byvar == "R2.CS.pop"){
    by.label <- expression(paste(R^2*phantom()[paste(CS,",",pop,sep = "")], sep = ""))
  } else if (byvar == "C.mean"){
    by.label <- expression(paste("mean(", C[app], ")", sep = ""))
  } else if (byvar == "C.median"){
    by.label <- expression(paste("median(", C[app], ")", sep = ""))
  } else if (byvar == "C.pop"){
    by.label <- expression(C[pop])
  } else if (byvar == "P.meas"){
    by.label <- expression(Q[meas])
  } else if (byvar == "EPP"){
    by.label <- expression(EPP)
  } else if (byvar == "nreq"){
    by.label <- "N"
  }
  
  ### Formatting
  xvar <- rlang::sym(xvar)
  yvar <- rlang::sym(yvar)
  byvar <- rlang::sym(byvar)
  
  ### Create base for plot
  if (is.null(bygroups)){
    gg.out <- ggplot(aes(x = !! xvar, y = !! yvar, color = !! byvar), data = df) +
      scale_colour_gradientn(colors = my_colors)
  } else {
    gg.out <- ggplot(aes(x = !! xvar, y = !! yvar, color = cut(!! byvar, bygroups, right = FALSE)), data = df)
  }
  
  ### Add data points
  if (is.null(alpha)){
    gg.out <- gg.out + geom_point(size = 0.5) + geom_abline(intercept = 0, slope = 1) + xlim(c(0, 1.25)) + ylim(c(0, 1.25))
  } else {
    gg.out <- gg.out + geom_point(size = 0.5, alpha = alpha) + geom_abline(intercept = 0, slope = 1)
  }
  
  ### Add extra bits
  gg.out <- gg.out + 
    tune::coord_obs_pred() +
    theme(legend.position = "bottom", 
          text = element_text(size = 8), 
          legend.text = element_text(size = 8), 
          legend.key.height = unit(0.3, "cm"),
          legend.margin = unit(c(-0.2,0,0,0), "cm"),
          legend.box.margin = unit(c(-0.2,0,0,0), "cm"),
          plot.margin = grid::unit(c(0.2,0.2,0.2,0.2), "cm")) +
    labs(color = by.label, x = x.label, y = y.label) + 
    ggtitle(paste(C1, " < ", expression(C[pop]), " < ", C2, sep = ""))
  
  ### If bygroups specified, split over two rows
  if (!is.null(bygroups)){
    gg.out <- gg.out + guides(color = guide_legend(nrow = 2))
  }
  
  return(gg.out)
  
}
####################
### Create plots ###
####################

###
### EPP
###
myplot1 <- create_plot(df = metadata, C1 = 0.64, C2 = 0.66,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "P.meas")
myplot2 <- create_plot(df = metadata, C1 = 0.69, C2 = 0.71,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "P.meas")
myplot3 <- create_plot(df = metadata, C1 = 0.74, C2 = 0.76,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "P.meas")
myplot4 <- create_plot(df = metadata, C1 = 0.79, C2 = 0.8,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "P.meas")


### Arrange
plotlist <- list(myplot1, myplot2, myplot3, myplot4)
plots.comb <- ggpubr::ggarrange(plotlist = plotlist, nrow = 2, ncol = 2)

# Save to disk
Cairo::CairoPNG("figures/gg.explore.Q.meas.by.Cstatistic.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plots.comb)
dev.off()


###
### EPP
###
myplot1 <- create_plot(df = metadata, C1 = 0.64, C2 = 0.66,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "EPP",
                       bygroups = c(0, 5, 10, 20, 40, 100, Inf))
myplot2 <- create_plot(df = metadata, C1 = 0.69, C2 = 0.71,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "EPP",
                       bygroups = c(0, 5, 10, 20, 40, 100, Inf))
myplot3 <- create_plot(df = metadata, C1 = 0.74, C2 = 0.76,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "EPP",
                       bygroups = c(0, 5, 10, 20, 40, 100, Inf))
myplot4 <- create_plot(df = metadata, C1 = 0.79, C2 = 0.8,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "EPP",
                       bygroups = c(0, 5, 10, 20, 40, 100, Inf))


### Arrange
plotlist <- list(myplot1, myplot2, myplot3, myplot4)
plots.comb <- ggpubr::ggarrange(plotlist = plotlist, nrow = 2, ncol = 2)

# Save to disk
Cairo::CairoPNG("figures/gg.explore.EPP.subset.Cstatistic.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plots.comb)
dev.off()

###
### N
###
myplot1 <- create_plot(df = metadata, C1 = 0.64, C2 = 0.66,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "nreq", 
                       bygroups = c(100, 150, 300, 500, Inf))
myplot2 <- create_plot(df = metadata, C1 = 0.69, C2 = 0.71,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "nreq", 
                       bygroups = c(100, 150, 300, 500, Inf))
myplot3 <- create_plot(df = metadata, C1 = 0.74, C2 = 0.76,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "nreq", 
                       bygroups = c(100, 150, 300, 500, Inf))
myplot4 <- create_plot(df = metadata, C1 = 0.79, C2 = 0.8,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "nreq", 
                       bygroups = c(100, 150, 300, 500, Inf))


### Arrange
plotlist <- list(myplot1, myplot2, myplot3, myplot4)
plots.comb <- ggpubr::ggarrange(plotlist = plotlist, nrow = 2, ncol = 2)

# Save to disk
Cairo::CairoPNG("figures/gg.explore.nreq.subset.Cstatistic.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plots.comb)
dev.off()


###
### C.pop
###
myplot1 <- create_plot(df = metadata, C1 = 0.64, C2 = 0.66,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "C.pop")
myplot2 <- create_plot(df = metadata, C1 = 0.69, C2 = 0.71,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "C.pop")
myplot3 <- create_plot(df = metadata, C1 = 0.74, C2 = 0.76,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "C.pop")
myplot4 <- create_plot(df = metadata, C1 = 0.79, C2 = 0.8,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "C.pop")


### Arrange
plotlist <- list(myplot1, myplot2, myplot3, myplot4)
plots.comb <- ggpubr::ggarrange(plotlist = plotlist, nrow = 2, ncol = 2)

# Save to disk
Cairo::CairoPNG("figures/gg.explore.C.pop.subset.Cstatistic.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plots.comb)
dev.off()


###
### C.pop
###
myplot1 <- create_plot(df = metadata, P = 4,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "C.pop")
myplot2 <- create_plot(df = metadata, P = 8,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "C.pop")
myplot3 <- create_plot(df = metadata, P = 12,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "C.pop")
myplot4 <- create_plot(df = metadata, P = 16,
                       xvar = "S.VH.mean", yvar = "S.pop.mean", 
                       byvar = "C.pop")


### Arrange
plotlist <- list(myplot1, myplot2, myplot3, myplot4)
plots.comb <- ggpubr::ggarrange(plotlist = plotlist, nrow = 2, ncol = 2)

# Save to disk
Cairo::CairoPNG("figures/gg.explore.C.pop.subset.Q.meas.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plots.comb)
dev.off()
colnames(metadata)
