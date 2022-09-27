###
### Functions to analyse results are in a seperate source file, so then when funcitons are loaded to run the simulation,
### these functions are not loaded and stored in the resulting .RData files. That makes it difficult to alter analysis functions
### at a later date, when older versions of the functions are saved in the .RData files
###

################################################################################################
### 5.1) Functions to create combine output plot data for ggplots, for large sample analysis ###
################################################################################################

### For data generated using DGM seqlog (only consider rcs for multinomial models)
create.large.sample.ggplot.data.DGMseqlog <- function(sim.out.binary, sim.out.multinomial, sim.out.multinomial.rcs3, sim.out.multinomial.rcs5, scenario.sim, seed.coef.sim){
  dat.ggplot <- rbind(data.frame(sim.out.binary[[1]][seq(1, nrow(sim.out.binary[[1]]), 100), ], "anal.method" = "BLR"),
                      data.frame(sim.out.multinomial[[1]][seq(1, nrow(sim.out.multinomial[[1]]), 100), ], "anal.method" = "MLR"),
                      data.frame(sim.out.multinomial.rcs3[[1]][seq(1, nrow(sim.out.multinomial.rcs3[[1]]), 100), ], "anal.method" = "MLR, rcs, 3 knots"),
                      data.frame(sim.out.multinomial.rcs5[[1]][seq(1, nrow(sim.out.multinomial.rcs5[[1]]), 100), ], "anal.method" = "MLR, rcs, 5 knots"))
  dat.ggplot$anal.method <- factor(dat.ggplot$anal.method, levels = c("BLR", "MLR", "MLR, rcs, 3 knots", "MLR, rcs, 5 knots"))
  dat.ggplot$scenario.sim <- paste("scenario ", c("A", "B", "C", "D")[as.numeric(scenario.sim)], ".", seed.coef.sim, sep = "")
  return(dat.ggplot)}

### For data generated using DGM mult (only consider rcs for binary models)
create.large.sample.ggplot.data.DGMmult <- function(sim.out.binary, sim.out.multinomial, sim.out.binary.rcs3, sim.out.binary.rcs5, scenario.sim, seed.coef.sim){
  dat.ggplot <- rbind(data.frame(sim.out.binary[[1]][seq(1, nrow(sim.out.binary[[1]]), 100), ], "anal.method" = "BLR"),
                      data.frame(sim.out.multinomial[[1]][seq(1, nrow(sim.out.multinomial[[1]]), 100), ], "anal.method" = "MLR"),
                      data.frame(sim.out.binary.rcs3[[1]][seq(1, nrow(sim.out.binary.rcs3[[1]]), 100), ], "anal.method" = "BLR, rcs, 3 knots"),
                      data.frame(sim.out.binary.rcs5[[1]][seq(1, nrow(sim.out.binary.rcs5[[1]]), 100), ], "anal.method" = "BLR, rcs, 5 knots"))
  dat.ggplot$anal.method <- factor(dat.ggplot$anal.method, levels =  c("BLR", "MLR", "BLR, rcs, 3 knots", "BLR, rcs, 5 knots"))
  dat.ggplot$scenario.sim <- paste("scenario ", c("A", "B", "C", "D")[as.numeric(scenario.sim)], ".", seed.coef.sim, sep = "")
  return(dat.ggplot)}


##########################################################################
### 5.2) Define function to generate ggplots for large sample analysis ###
##########################################################################
create.large.sample.ggplot <- function(data.in, font.size, 
                                       x.lim.in, y.lim.in, xlab.in, ylab.in){
  env <- new.env(parent = globalenv())
  env$subset <- data.in
  env$font.size <- font.size
  env$x.lim <- x.lim.in
  env$y.lim <- y.lim.in
  env$xlab.in <- xlab.in
  env$ylab.in <- ylab.in
  
  ggplot.out <- with(env, {
    ggplot.out <- ggplot(subset) +
      geom_abline(intercept = 0, slope = 1, color = "black", size = 0.75, linetype = "dotted") +
      geom_line(aes(x = pred.Y1, y = pred.obs.Y1.loess), size = 0.75, color = "red", show.legend = FALSE) +  
      facet_grid(anal.method ~ scenario.sim) +
      #xlim(x.lim) + 
      xlab(xlab.in) +
      #ylim(y.lim) + 
      ylab(ylab.in) + theme_bw(base_size = font.size)
    return(ggplot.out)
  })
  
  return(ggplot.out)
}


#############################################################################################################################
### 6.1) This function takes the calibration plot output from the simulation, and gets in a format to produce the ggplots ###
### This gets the data from a specific scenario and get into long format, and calculates mean/median/p5/p90/range. Then   ###
### data from all these scenarios then needs to be combined into a single dataset for plotting, which is done with the next function.###
#############################################################################################################################
create.small.sample.ggplot.data.general <- function(sim.out.calib.data, pred.eval){
  
  ### sim.out.calib.data should be a list of all the calibration data outputs from a given scenario
  
  ### First re-assign all the sim.run variables
  for (i in 1:length(sim.out.calib.data)){
    sim.out.calib.data[[i]]$sim.run <- rep(i, nrow(sim.out.calib.data[[i]]))
  }
  
  ### Want to combine the calibration data into one dataset for ggplot
  calib.data.comb <- do.call("rbind", sim.out.calib.data)
  calib.data.comb$map.var <- rep("group1", nrow(calib.data.comb))
  
  
  ### Create temporary datasets in which we can calculate mean, p5 and p5 of each observed risk across all simulations
  ### for each method
  calib.data.comb.cbind <- do.call("cbind", sim.out.calib.data)
  
  calib.data.lin <- calib.data.comb.cbind[ , colnames(calib.data.comb.cbind) == "pred.obs.Y1.lin"]
  
  #   calib.data.rcs <- calib.data.comb.cbind[ , colnames(calib.data.comb.cbind) == "pred.obs.Y1.rcs"]
  
  calib.data.loess <- calib.data.comb.cbind[ , colnames(calib.data.comb.cbind) == "pred.obs.Y1.loess"]
  
  
  ### Calculate the mean, p5 and p95 for each analysis method (and smoother type)
  ## Linear recalibration model
  obs.mean.lin <- rowMeans(calib.data.lin, na.rm = TRUE)
  obs.sd.lin <- apply( calib.data.lin, 1, var, na.rm = TRUE)
  obs.p5.lin <- apply( calib.data.lin, 1, quantile, probs = 0.05, na.rm = TRUE)
  obs.p50.lin <- apply( calib.data.lin, 1, quantile, probs = 0.5, na.rm = TRUE)
  obs.p95.lin <- apply( calib.data.lin, 1, quantile, probs = 0.95, na.rm = TRUE)
  obs.range.lin <- obs.p95.lin - obs.p5.lin
  
  #   ## rcs recalibration model
  #   obs.mean.rcs <- rowMeans(calib.data.rcs, na.rm = TRUE)
  #   obs.p5.rcs <- apply( calib.data.rcs, 1, quantile, probs = 0.05, na.rm = TRUE)
  #   obs.p95.rcs <- apply( calib.data.rcs, 1, quantile, probs = 0.95, na.rm = TRUE)
  
  
  ## Loess smoother recalibraiton model
  obs.mean.loess <- rowMeans(calib.data.loess, na.rm = TRUE)
  obs.sd.loess <- apply( calib.data.loess, 1, var, na.rm = TRUE)
  obs.p5.loess <- apply( calib.data.loess, 1, quantile, probs = 0.05, na.rm = TRUE)
  obs.p50.loess <- apply( calib.data.loess, 1, quantile, probs = 0.5, na.rm = TRUE)
  obs.p95.loess <- apply( calib.data.loess, 1, quantile, probs = 0.95, na.rm = TRUE)
  obs.range.loess <- obs.p95.loess - obs.p5.loess
  
  
  ### Put these into a dataset to rbind with the main data going into the ggplot
  summary.data <- data.frame("pred.Y1" = rep(pred.eval, 6),
                             "pred.obs.Y1.lin" = c(obs.mean.lin, obs.sd.lin, obs.p50.lin, obs.p5.lin, obs.p95.lin, obs.range.lin),
                             #"pred.obs.Y1.rcs" = c(obs.mean.rcs, obs.p5.rcs, obs.p95.rcs),
                             "pred.obs.Y1.loess" = c(obs.mean.loess, obs.sd.loess, obs.p50.loess, obs.p5.loess, obs.p95.loess, obs.range.loess),
                             "sim.run" = c(rep(length(sim.out.calib.data) + 1, length(pred.eval)), 
                                           rep(length(sim.out.calib.data) + 2, length(pred.eval)), 
                                           rep(length(sim.out.calib.data) + 3, length(pred.eval)),
                                           rep(length(sim.out.calib.data) + 4, length(pred.eval)),
                                           rep(length(sim.out.calib.data) + 5, length(pred.eval)),
                                           rep(length(sim.out.calib.data) + 6, length(pred.eval))),
                             "map.var" = c(rep("mean", length(pred.eval)), rep("sd", length(pred.eval)), rep("p50", length(pred.eval)),
                                           rep("p5", length(pred.eval)), rep("p95", length(pred.eval)), rep("range", length(pred.eval))))
  
  ### Add to the main datasets
  calib.data.comb <- rbind(calib.data.comb, summary.data)
  
  ### Return the output
  return(calib.data.comb)
}

#############################################################################################################
### 6.2) This function loads scenarios of interest, creates the data for plots from different scenarios,  ###
### and combines it into a single dataset for plotting                                                    ###
#############################################################################################################
# summary.metrics = c("p50")
# N.vec = c(1000, 2500)
# scen.vec = c(1, 2)
# DGM.in = "mult"
# K.sim = 5
# seed.coef.sim = 1
# N.devel.sim = 1000
# scenario.sim = 1
create.small.sample.ggplot.data.specific <- function(summary.metrics, N.vec, scen.vec, DGM.in, K.sim, seed.coef.sim){
  
  ### Start a counter for the loops
  counter <- 1
  
  ### Load output for each scenario of interest, create get into appropriate format, and rbind into a dataset
  for (N.devel.sim in N.vec){
    print(N.devel.sim)
    for (scenario.sim in scen.vec){
      print(scenario.sim)
      ### Load combined simulation results
      load(paste("data/results_small_samp_DGM", DGM.in, "_scen" , scenario.sim, ".", seed.coef.sim, "_K", 
                 K.sim, "_P", P.sim, "_N", N.devel.sim, ".RData", sep =""))
      
      ### First create plot data for multinomial and binary analyses
      print("do multinomial")
      data.plot.multinomial <- create.small.sample.ggplot.data.general(comb.sim.out.calib.data.multinomial, pred.eval = pred.eval)
      print("do binary")
      data.plot.binary <- create.small.sample.ggplot.data.general(comb.sim.out.calib.data.binary, pred.eval = pred.eval)
      
      ### These will be ploted together in a single graph, so combine into a single dataset
      data.both <- rbind(data.frame(data.plot.multinomial, "anal.method" = "multinomial"), 
                         data.frame(data.plot.binary, "anal.method" = "binary"))
      
      ### Add a variable to indicate the scenario
      data.both$scen <- paste("scenario ", scenario.sim, sep = "")
      
      ### Add a variable to indicate sample size
      data.both$samp.size <- paste("N = ", N.devel.sim, sep = "")
      
      ### Add to dataset
      if (counter == 1){
        data.out <- data.both
      } else if(counter > 1){
        data.out <- rbind(data.out, data.both)
      }
      
      ### Add 1 to counter
      counter <- counter + 1
      
      ### Remove items to avoid crossover
      rm(data.plot.multinomial, data.plot.binary, data.both, comb.sim.out.calib.data.multinomial, comb.sim.out.calib.data.binary)
    }
  }
  
  ### Finally reduce output dataset to observations which we want to plot
  data.out <- data.out[data.out$map.var %in% summary.metrics, ]
  
  ### Also turn samp.size into a factor variable with appropriate levels
  data.out$samp.size <- factor(data.out$samp.size, levels = paste("N = ", N.vec, sep = ""))
  
  ### Also give the variable scen the appropriate names
  data.out$scen[data.out$scen == "scenario 1"] <- paste("scenario A", ".", seed.coef.sim, sep = "")
  data.out$scen[data.out$scen == "scenario 2"] <- paste("scenario B", ".", seed.coef.sim, sep = "")
  data.out$scen[data.out$scen == "scenario 3"] <- paste("scenario C", ".", seed.coef.sim, sep = "")
  data.out$scen[data.out$scen == "scenario 4"] <- paste("scenario D", ".", seed.coef.sim, sep = "")
  
  ### Return the dataset
  return(data.out)
}


##########################################################################
### 6.3) This function takes ggplot data, and creates a ggplot we want ###
##########################################################################
create.small.sample.ggplots <- function(data.in, x.lim.in, y.lim.in, xlab.in, ylab.in, font.size.in){
  
  #   data.ggplot.comb.in = data.ggplot.comb
  #   x.lim.in = x.lim
  #   y.lim.in = y.lim
  #   xlab.in = xlab
  #   ylab.in = ylab
  #   font.size.in = font.size
  #   linecolours.in = linecolours
  #   lwd.in = lwd
  #   abline.size.in = abline.size
  #   sum.vec.in = sum.vec
  
  env <- new.env(parent = globalenv())
  env$subset <- data.in
  #env$group.colours <- group.colours.in
  #env$lwd.fix <- lwd.fix.in
  env$font.size <- font.size.in
  env$x.lim <- x.lim.in
  env$y.lim <- y.lim.in
  env$xlab.in <- xlab.in
  env$ylab.in <- ylab.in
  #env$abline.size <- abline.size.in
  
  ggplot.out <- with(env, {
    ggplot.out <- ggplot(subset) +
      geom_abline(intercept = 0, slope = 1, color = "black", size = 0.75, linetype = "dotted") + 
      geom_line(aes(x = pred.Y1, y = pred.obs.Y1.loess, color = anal.method, 
                    linetype = anal.method), 
                size = 0.75, show.legend = TRUE) +
      facet_grid(samp.size ~ scen) +
      #scale_colour_brewer(palette = group.colours) +
      #scale_linetype_manual(values = group.lty) +
      #scale_size_manual(values = group.lwd) +
      #xlim(x.lim) + 
      xlab(xlab.in) +
      #ylim(y.lim) + 
      ylab(ylab.in) + 
      theme_bw(base_size = font.size) +
      theme(legend.position = "bottom")
    return(ggplot.out)
  })
  
  return(ggplot.out)
}


##########################################################################
### 6.4) This function takes ggplot data, and creates a ggplot that will be used to plot calibration variation. ###
### The code is very similar to 6.4, except there is no abline           ###
##########################################################################
create.small.sample.ggplots.range <- function(data.in, x.lim.in, y.lim.in, xlab.in, ylab.in, font.size.in){
  
  #   data.ggplot.comb.in = data.ggplot.comb
  #   x.lim.in = x.lim
  #   y.lim.in = y.lim
  #   xlab.in = xlab
  #   ylab.in = ylab
  #   font.size.in = font.size
  #   linecolours.in = linecolours
  #   lwd.in = lwd
  #   abline.size.in = abline.size
  #   sum.vec.in = sum.vec
  
  env <- new.env(parent = globalenv())
  env$subset <- data.in
  #env$group.colours <- group.colours.in
  #env$lwd.fix <- lwd.fix.in
  env$font.size <- font.size.in
  env$x.lim <- x.lim.in
  env$y.lim <- y.lim.in
  env$xlab.in <- xlab.in
  env$ylab.in <- ylab.in
  #env$abline.size <- abline.size.in
  
  ggplot.out <- with(env, {
    ggplot.out <- ggplot(subset) +
      geom_line(aes(x = pred.Y1, y = pred.obs.Y1.loess, color = anal.method, 
                    linetype = anal.method), 
                size = 0.75, show.legend = TRUE) +
      facet_grid(samp.size ~ scen) +
      #scale_colour_brewer(palette = group.colours) +
      #scale_linetype_manual(values = group.lty) +
      #scale_size_manual(values = group.lwd) +
      #xlim(x.lim) + 
      xlab(xlab.in) +
      #ylim(y.lim) + 
      ylab(ylab.in) + 
      theme_bw(base_size = font.size) +
      theme(legend.position = "bottom")
    return(ggplot.out)
  })
  
  return(ggplot.out)
}