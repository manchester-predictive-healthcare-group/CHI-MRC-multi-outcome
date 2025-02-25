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
metadata <- readRDS("data/heuristic_shrinkage_sim_metadata.rds")
metadata <- subset(metadata, C.pop != 0)

### Write a function which will combine data for a specific version
create_plot <- function(df, correlated = TRUE, cat = FALSE, xvar, yvar, byvar, bygroups = NULL, xrange = NULL, yrange = NULL, alpha = NULL){
  
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
  fun_color_range <- colorRampPalette(c("green", "blue"))
  my_colors <- fun_color_range(20)
  
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
    y.label <- expression(paste("mean(", S[pop], ")", sep = ""))
  } else if (yvar == "S.pop.median"){
    y.label <- expression(paste("median(", S[pop], ")", sep = ""))
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
    gg.out <- gg.out + geom_point(size = 0.5) + geom_abline(intercept = 0, slope = 1)
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
          plot.margin = grid::unit(c(0,0,0,0), "cm")) +
    labs(color = by.label, x = x.label, y = y.label)
  
  ### If bygroups specified, split over two rows
  if (!is.null(bygroups)){
    gg.out <- gg.out + guides(color = guide_legend(nrow = 2))
  }
  
  return(gg.out)
  
}

#########################
### Create mean plots ###
#########################

### Write a little function to do this
save_plots <- function(xvar2, type, correlated, cat, alpha = NULL){
  
  # Create plots
  plot.corr.R2 <- create_plot(df = metadata, 
                              correlated = correlated,
                              cat = cat,
                              xvar = paste(xvar2, ".", type, sep = ""),
                              yvar = paste("S.pop", ".", type, sep = ""),
                              byvar = paste("R2.CS.app.", type, sep = ""), 
                              yrange = c(0,1), 
                              alpha = alpha)
  plot.corr.R2.pop <- create_plot(df = metadata, 
                              correlated = correlated,
                              cat = cat,
                              xvar = paste(xvar2, ".", type, sep = ""),
                              yvar = paste("S.pop", ".", type, sep = ""),
                              byvar = "R2.CS.pop", 
                              yrange = c(0,1), 
                              alpha = alpha)
  plot.corr.C <- create_plot(df = metadata, 
                             correlated = correlated,
                             cat = cat,
                             xvar = paste(xvar2, ".", type, sep = ""),
                             yvar = paste("S.pop", ".", type, sep = ""),
                             byvar = paste("C.", type, sep = ""),  
                             yrange = c(0,1), 
                             alpha = alpha)
  plot.corr.C.pop <- create_plot(df = metadata, 
                             correlated = correlated,
                             cat = cat,
                             xvar = paste(xvar2, ".", type, sep = ""),
                             yvar = paste("S.pop", ".", type, sep = ""),
                             byvar = "C.pop",  
                             yrange = c(0,1), 
                             alpha = alpha)
  plot.corr.P <- create_plot(df = metadata, 
                             correlated = correlated,
                             cat = cat,
                             xvar = paste(xvar2, ".", type, sep = ""),
                             yvar = paste("S.pop", ".", type, sep = ""),
                             byvar = "P.meas", 
                             yrange = c(0,1), 
                             alpha = alpha)
  plot.corr.n <- create_plot(df = metadata, 
                             correlated = correlated,
                             cat = cat,
                             xvar = paste(xvar2, ".", type, sep = ""),
                             yvar = paste("S.pop", ".", type, sep = ""),
                             byvar = "nreq", 
                             bygroups = c(100, 150, 300, 500, Inf), 
                             yrange = c(0,1), 
                             alpha = alpha)
  
  # Combine plots
  plotlist <- list(plot.corr.R2.pop, plot.corr.C.pop, plot.corr.P, plot.corr.n)
  plotlist.app <- list(plot.corr.R2, plot.corr.C)
  plots.corr.comb <- ggpubr::ggarrange(plotlist = plotlist, nrow = 2, ncol = 2)
  plots.corr.comb.app <- ggpubr::ggarrange(plotlist = plotlist.app, nrow = 1, ncol = 2)
  
  # Save to disk
  Cairo::CairoPNG(paste("figures/gg.corr", as.numeric(correlated), ".cat", as.numeric(cat), ".", xvar2, ".", type, ".alpha", ifelse(is.null(alpha), 0, alpha), ".png", sep = ""),
                  width = 7, height = 7, unit = "in", dpi = 300)
  plot(plots.corr.comb)
  dev.off()
  Cairo::CairoPNG(paste("figures/gg.app.corr", as.numeric(correlated), ".cat", as.numeric(cat), ".", xvar2, ".", type, ".alpha", ifelse(is.null(alpha), 0, alpha), ".png", sep = ""),
                  width = 7, height = 3.5, unit = "in", dpi = 300)
  plot(plots.corr.comb.app)
  dev.off()
  
}

###
### Create plots
###
for (xvar2 in c("S.VH", "S.boot")){
  for (type in c("mean", "median")){
    for (correlated in c(TRUE, FALSE)){
      for (cat in c(FALSE)){
        for (alpha in c(0.25)){
          print(paste(xvar2, type, correlated, cat, alpha))
          save_plots(xvar2 = xvar2, type = type, correlated = correlated, cat = cat, alpha = alpha)
        }
      }
    }
  }
}

###
### Reproduce C-statistic plot highlighting models with C-statistic around 0.8
###
plot.corr.C <- create_plot(df = metadata, xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "C.mean", 
                           bygroups = c(0, 0.775, 0.825, 1), yrange = c(0,1))

plot.corr.C.pop <- create_plot(df = metadata, xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "C.pop", 
                           bygroups = c(0, 0.73, 0.8, 1), yrange = c(0,1))

plot.nocorr.C <- create_plot(df = metadata, correlated = FALSE, xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "C.mean", 
                           bygroups = c(0, 0.775, 0.825, 1), yrange = c(0,1))

plot.nocorr.C.pop <- create_plot(df = metadata, correlated = FALSE, xvar = "S.VH.mean", yvar = "S.pop.mean", byvar = "C.pop", 
                               bygroups = c(0, 0.7, 0.8, 1), yrange = c(0,1))


# Save to disk
Cairo::CairoPNG("figures/gg.corr1.S.VH.mean.C.group.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plot.corr.C)
dev.off()

Cairo::CairoPNG("figures/gg.corr1.S.VH.mean.C.pop.group.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plot.corr.C.pop)
dev.off()

Cairo::CairoPNG("figures/gg.corr0.S.VH.mean.C.group.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plot.nocorr.C)
dev.off()

Cairo::CairoPNG("figures/gg.corr0.S.VH.mean.C.pop.group.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plot.nocorr.C.pop)
dev.off()
