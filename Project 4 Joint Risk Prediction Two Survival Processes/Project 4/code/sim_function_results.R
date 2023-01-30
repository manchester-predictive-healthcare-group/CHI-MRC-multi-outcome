########################
########################
### Define functions ###
########################
########################


###############################################
### 1) Function to summarise discrimination ###
###############################################

summarise.discrimination.func <- function(data.in){
  discrim.tab <- c("product" = paste(round(mean(data.in$product), 3), " (", round(sd(data.in$product),3), ")", sep = ""),
                   "joint" = paste(round(mean(data.in$joint), 3), " (", round(sd(data.in$joint),3), ")", sep = ""),
                   "msm" = paste(round(mean(data.in$msm), 3), " (", round(sd(data.in$msm),3), ")", sep = ""),
                   "c.clay" = paste(round(mean(data.in$c.clay), 3), " (", round(sd(data.in$c.clay),3), ")", sep = ""),
                   "c.gumb" = paste(round(mean(data.in$c.gumb), 3), " (", round(sd(data.in$c.gumb),3), ")", sep = ""),
                   "c.frank" = paste(round(mean(data.in$c.frank), 3), " (", round(sd(data.in$c.frank),3), ")", sep = ""),
                   "f.normal_weib" = paste(round(mean(data.in$f.normal_weib), 3), " (", round(sd(data.in$f.normal_weib),3), ")", sep = ""),
                   "f.gamma_weib" = paste(round(mean(data.in$f.gamma_weib), 3), " (", round(sd(data.in$f.gamma_weib),3), ")", sep = ""))
  return(discrim.tab)
}



########################################################################################
### 2.1) Function to take the parallelised output stored in all the .RData files,      ###
### and combine into a single object which can be analysed using the above functions ###
########################################################################################

### Create a function which will combine all the analysis results from the 1000 .RData files
create.combined.output.objects <- function(dgm.in){
  
  ### Load simulation output from first run of simulation
  load(paste("data/", scen, "/sim_run_", scen, "_dgm", dgm.in, "_n", n.devel.fix, "v", n.valid.fix, "s", 1, ".RData", sep = ""))
  
  ### Assign objects of interest to a new name
  C.Har.tab.comb <- sim.out.C.Har
  C.Uno.tab.comb <- sim.out.C.Uno
  sim.out.calib.plot.comb <- sim.out.calib.plot
  failed.sim <- 0
  
  ### Remove everything except what we want to retain. In workspaces where the simulation failed, 
  ### we don't want leave the results from the previous iteration in the workspace
  rm(list=setdiff(ls(), list("C.Har.tab.comb", "C.Uno.tab.comb", "sim.out.calib.plot.comb", "failed.sim",
                             "dgm.in", "scen", "n.devel.fix", "n.valid.fix", 
                             "n.knot", "x.lim", "y.lim", "pred.eval", 
                             "create.dat.calib.plot", "create.dat.calib.plot.run","summarise.discrimination.func", 
                             "create.combined.output.objects", "analyse.results", "anal.methods.sim")))
  
  
  ### For the next 999 simulation runs, want to load the workspace, check if simulation was successful, then append the results to the combined object               
  for (i in 2:1000){
    
    print(i)
    
    if (file.exists(paste("data/", scen, "/sim_run_", scen, "_dgm", dgm.in, "_n", n.devel.fix, "v", n.valid.fix, "s", i, ".RData", sep = "")) == TRUE){
      
      ## Load data
      load(paste("data/", scen, "/sim_run_", scen, "_dgm", dgm.in, "_n", n.devel.fix, "v", n.valid.fix, "s", i, ".RData", sep = ""))
      
      ## Add objects to combined
      C.Har.tab.comb <- rbind(C.Har.tab.comb, sim.out.C.Har)
      C.Uno.tab.comb <- rbind(C.Uno.tab.comb, sim.out.C.Uno)
      sim.out.calib.plot.comb <- append(sim.out.calib.plot.comb, sim.out.calib.plot)
    } else {
      if (!(failed.sim[1] == 0)){
        failed.sim <- c(failed.sim, i)
      } else if (failed.sim[1] == 0){
        failed.sim <- i
      }
    }
    
    ### Remove excess
    rm(list=setdiff(ls(), list("C.Har.tab.comb", "C.Uno.tab.comb", "sim.out.calib.plot.comb", "failed.sim",
                               "dgm.in", "scen", "n.devel.fix", "n.valid.fix", 
                               "n.knot", "x.lim", "y.lim", "pred.eval", 
                               "create.dat.calib.plot", "create.dat.calib.plot.run","summarise.discrimination.func", 
                               "create.combined.output.objects", "analyse.results", "anal.methods.sim")))
  }
  
  return(list("C.Har.tab.comb" = C.Har.tab.comb, 
              "C.Uno.tab.comb" = C.Uno.tab.comb, 
              "sim.out.calib.plot.comb" = sim.out.calib.plot.comb,
              "failed.sim" = failed.sim))
}


########################################################################################
### 2.2) Function to take the parallelised output stored in all the .RData files,    ###
### and combine into a single object which can be analysed using the above functions ###
### but specifically for the nocorr scenarios, as the .RData files have slightly different names ###
########################################################################################

### Create a function which will combine all the analysis results from the 1000 .RData files
create.combined.output.objects.nocorr <- function(){
  
  ### Load simulation output from first run of simulation
  load(paste("data/", scen, ".1/sim_run_", scen, "_dgmnocorr_n", n.devel.fix, "v", n.valid.fix, "s", 1, ".RData", sep = ""))
  
  ### Assign objects of interest to a new name
  C.Har.tab.comb <- sim.out.C.Har
  C.Uno.tab.comb <- sim.out.C.Uno
  sim.out.calib.plot.comb <- sim.out.calib.plot
  failed.sim <- 0
  
  ### Remove everything except what we want to retain. In workspaces where the simulation failed, 
  ### we don't want leave the results from the previous iteration in the workspace
  rm(list=setdiff(ls(), list("C.Har.tab.comb", "C.Uno.tab.comb", "sim.out.calib.plot.comb", "failed.sim",
                             "dgm.in", "scen", "n.devel.fix", "n.valid.fix", 
                             "n.knot", "x.lim", "y.lim", "pred.eval", 
                             "create.dat.calib.plot", "create.dat.calib.plot.run","summarise.discrimination.func", 
                             "create.combined.output.objects", "analyse.results", "anal.methods.sim")))
  
  
  ### For the next 999 simulation runs, want to load the workspace, check if simulation was successful, then append the results to the combined object               
  for (i in 2:1000){
    
    print(i)
    
    if (file.exists(paste("data/", scen, ".1/sim_run_", scen, "_dgmnocorr_n", n.devel.fix, "v", n.valid.fix, "s", i, ".RData", sep = "")) == TRUE){
      
      ### Load simulation output from first run of simulation
      load(paste("data/", scen, ".1/sim_run_", scen, "_dgmnocorr_n", n.devel.fix, "v", n.valid.fix, "s", i, ".RData", sep = ""))
      
      ## Add objects to combined
      C.Har.tab.comb <- rbind(C.Har.tab.comb, sim.out.C.Har)
      C.Uno.tab.comb <- rbind(C.Uno.tab.comb, sim.out.C.Uno)
      sim.out.calib.plot.comb <- append(sim.out.calib.plot.comb, sim.out.calib.plot)
    } else {
      if (!(failed.sim[1] == 0)){
        failed.sim <- c(failed.sim, i)
      } else if (failed.sim[1] == 0){
        failed.sim <- i
      }
    }
    
    ### Remove excess
    rm(list=setdiff(ls(), list("C.Har.tab.comb", "C.Uno.tab.comb", "sim.out.calib.plot.comb", "failed.sim",
                               "dgm.in", "scen", "n.devel.fix", "n.valid.fix", 
                               "n.knot", "x.lim", "y.lim", "pred.eval", 
                               "create.dat.calib.plot", "create.dat.calib.plot.run","summarise.discrimination.func", 
                               "create.combined.output.objects", "analyse.results", "anal.methods.sim")))
  }
  
  return(list("C.Har.tab.comb" = C.Har.tab.comb, 
              "C.Uno.tab.comb" = C.Uno.tab.comb, 
              "sim.out.calib.plot.comb" = sim.out.calib.plot.comb,
              "failed.sim" = failed.sim))
}



#######################################################################################################################################
### 3) Function to create a dataset which contains smoothed observed risks, corresponding to predicted risks, for a single dataset. ###
### This function is called on in function 3 ###
#######################################################################################################################################

### Define a function that will create data for a calibration plot
### Input is a table with containing predicted risks for each method from the simulation, the true risks for each individual,
### and a seperate vector indicating what values of predicted risk we want to do the plot over

### We calculate the corresponding observed risk, for a given predicted risk, according to each model
create.dat.calib.plot <- function(data.in, pred.eval.in, n.knot.in, anal.methods.sim.in){
  
  ### For each analysis method, create a model by regressing true/observed risk on predicted risk using rcs
  model.calib.list <- vector("list", length(anal.methods.sim.in))
  for (i in 1:length(anal.methods.sim.in)){
    # Create a dataset with which to fit model
    temp.dat <- data.frame("truerisk" = data.in[ , "truerisk"], "predrisk" = data.in[ , paste("predrisk.", anal.methods.sim.in[[i]], sep = "")])
    # Arrange by predrisk
    temp.dat <- arrange(temp.dat, predrisk)
    # Remove observations where the predrisk risk is clearly erroneous (above 1 or below 0)
    ## Note that this is rare, and occurs solely when outcome prevalence are small (scenarios 1.1 and 1.2), the sample size is small (n = 1000),
    ## and the analysis method is the multistate model. We will make a note of this at publication.
    temp.dat <- temp.dat[temp.dat$predrisk < 1 & temp.dat$predrisk > 0 , ]
    # Fit the model
    model.calib.list[[i]] <- mfp(truerisk ~ fp(predrisk), data = temp.dat, family = "gaussian")
  }
  
  
  ### Now make observed risk as a function of predicted risk by making predictions from these models
  ## Create object to store data
  obs.smooth.list <- vector("list", length(anal.methods.sim.in))
  
  ## Create a data frame containing vector of predicted risks, which we map to observed risks using each model
  newdata.pred <- data.frame("predrisk" = pred.eval.in)
  
  ## Generate the observed risks
  for (i in 1:length(anal.methods.sim.in)){
    newdata.pred[, paste("obs.smooth.", anal.methods.sim.in[i], sep = "")] <- predict(model.calib.list[[i]], newdata = newdata.pred, type = 'response')
  }
  
  return(newdata.pred)
} 

#############################################################################
### 4) Function to create the data for calibration plots (uses function 3, and output from function 2), ###
### running through every iteration.                                                                     ###
### Note that data.inn is the combined output from the simulation (comes from function 2)                ###
### It's a list of all the simulation outputs (sim.out.calib.plot)                                      ###
#############################################################################
create.dat.calib.plot.run <- function(data.inn, n.knot.in, pred.eval.in, anal.methods.sim.inn){
  
  #   data.inn <- dgm.msm.out.comb[[3]][1:25]
  #   x.lim.in <- x.lim
  #   y.lim.in <- y.lim
  #   font.size.in <- font.size
  #   pred.eval.in <- pred.eval
  #   anal.methods.sim.inn <- anal.methods.sim
  #   pred.eval.in <- pred.eval
  
  ### Let k.sim.suc.in be the number of successful simulations
  k.sim.suc.in <- length(data.inn)
  #k.sim.suc.in <- 50
  
  
  ### Create an object to store all the converted datasets
  dat.calib.plot.list <- vector("list", k.sim.suc.in)
  
  ### Apply the function defined above to each piece of simulations output
  for (i in 1:k.sim.suc.in){
    print(paste("analyse", i, sep = ""))
    dat.calib.plot.list[[i]] <- create.dat.calib.plot(data.in = data.inn[[i]], pred.eval.in = pred.eval.in, n.knot.in = n.knot.in, 
                                                      anal.methods.sim.in = anal.methods.sim.inn)
  }
  
  ### Combine all the datasets into one dataset, so all the observed risks are in one place (each row is a predicted risk)
  print("memA")
  dat.calib.plot <- do.call(cbind, dat.calib.plot.list)
  
  ### Seperate out the observed risks for each analysis method into different datasets
  print("memB")
  dat.calib.plot.seperate <- vector("list", length(anal.methods.sim.inn))
  names(dat.calib.plot.seperate) <- anal.methods.sim.inn
  
  for (i in 1:length(anal.methods.sim.inn)){
    ### Seperate them out
    print("memC")
    temp.data <- dat.calib.plot[ , colnames(dat.calib.plot) == paste("obs.smooth.", anal.methods.sim.inn[i], sep= "")]
    colnames(temp.data) <- paste("obs.smooth.", 1:k.sim.suc.in, sep = "")
    print("memD")
    ### Now calculate row means and percentiles from each dataset, and also store the vector of predicted risks
    ### They must be calculated and put in an object, then added in after, to not impact the calculation of percentiles and mean
    
    # Calculate
    print("memE")
    obs.mean <- rowMeans(temp.data, na.rm = TRUE)
    print("memF")
    obs.p5 <- apply( temp.data, 1, quantile, probs = 0.05, na.rm = TRUE)
    print("memG")
    obs.p95 <- apply( temp.data, 1, quantile, probs = 0.95, na.rm = TRUE)
    obs.p50 <- apply( temp.data, 1, quantile, probs = 0.5, na.rm = TRUE)
    
    # Add to object
    temp.data$predrisk <- pred.eval.in
    temp.data$obs.mean <- obs.mean
    temp.data$obs.p5 <- obs.p5
    temp.data$obs.p50 <- obs.p50
    temp.data$obs.p95 <- obs.p95
    temp.data$obs.mean.diff <- obs.mean - pred.eval.in
    temp.data$obs.median.diff <- obs.p50 - pred.eval.in
    temp.data$obs.p.range <- obs.p95 - obs.p5
    
    print("memH")
    ### Create long data for plot
    temp.data <- melt(temp.data, id = "predrisk")
    
    ### Want to create a new variable to map to colours, line thickness, etc, which will be same for all individual plots, then different for the mean and percentiles
    print("memI")
    map.var <- rep(0, nrow(temp.data))
    print("memJ")
    map.var[grep("obs.smooth", temp.data$variable, fixed = TRUE)] <- "ind"
    print("memK")
    map.var[temp.data$variable == "obs.p5"] <- "p5"
    map.var[temp.data$variable == "obs.p50"] <- "p50"
    map.var[temp.data$variable == "obs.p95"] <- "p95"
    map.var[temp.data$variable == "obs.mean"] <- "mean"
    map.var[temp.data$variable == "obs.mean.diff"] <- "mean.diff"
    map.var[temp.data$variable == "obs.median.diff"] <- "median.diff"
    map.var[temp.data$variable == "obs.p.range"] <- "p.range"
    temp.data$map.var <- map.var
    
    dat.calib.plot.seperate[[anal.methods.sim.inn[i]]] <- temp.data
  }
  
  return(dat.calib.plot.seperate)
}



#######################################################################################################
### 5) A function to analyse the results as outlined in methods section (utilises function 4) ###
#######################################################################################################

### Write a function to do this
### Input must be the output from the previous function create.combined.output.object
analyse.results <- function(data.in){
  
  ## Create output object
  output.obj <- vector("list", 12)
  names(output.obj) <- c("C.Har.tab", "C.Uno.tab", "data.calib.product", "data.calib.joint", "data.calib.msm", "data.calib.c.clay", 
                         "data.calib.c.gumb", "data.calib.c.frank", "data.calib.f.normal_weib", "data.calib.f.gamma_weib", "failed.sim", "n.sim.suc")
  
  ## Create discrimination tables
  output.obj[["C.Har.tab"]] <- summarise.discrimination.func(data.in[["C.Har.tab.comb"]])
  output.obj[["C.Uno.tab"]] <- summarise.discrimination.func(data.in[["C.Uno.tab.comb"]])
  
  ## Create calibration plots and assign to object
  data.calib <- create.dat.calib.plot.run(data.inn = data.in[["sim.out.calib.plot.comb"]], n.knot.in = n.knot, 
                                          pred.eval.in = pred.eval, anal.methods.sim.inn = anal.methods.sim)
  
  output.obj[["data.calib.msm"]] <- data.calib[["msm"]]
  output.obj[["data.calib.product"]] <- data.calib[["product"]]
  output.obj[["data.calib.joint"]] <- data.calib[["joint"]]
  output.obj[["data.calib.c.clay"]] <- data.calib[["c.clay"]]
  output.obj[["data.calib.c.gumb"]] <- data.calib[["c.gumb"]]
  output.obj[["data.calib.c.frank"]] <- data.calib[["c.frank"]]
  output.obj[["data.calib.f.normal_weib"]] <- data.calib[["f.normal_weib"]]
  output.obj[["data.calib.f.gamma_weib"]] <- data.calib[["f.gamma_weib"]]
  
  ## Add vector of failed simulation runs
  output.obj[["failed.sim"]] <- data.in[["failed.sim"]]
  
  ## Add a number of succesful simulations
  output.obj[["n.sim.suc"]] <- nrow(data.in[["C.Har.tab.comb"]])
  
  ## Return output object
  return(output.obj)
}


#####################################################################
### 6) Create ggplots from analysed data (output from function 5) ###
#####################################################################
create.individual.ggplots <- function(analysed.data, x.lim.in, y.lim.in, font.size.in, linecolours.in, linewidths.in, anal.methods.sim.inn){
  
  ### Cycle through and create the dataset with only the individual curves, the p5, p95, p50 and mean curves
  for (i in 1:length(anal.methods.sim.inn)){
    ## Create a temporary dataset
    temp.data <- analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]]
    ## Extract the rows which contain summary data
    temp.data <- temp.data[temp.data$map.var %in% c("ind", "mean", "p5", "p95", "p50"), ]
    ## Assign back
    analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]] <- temp.data
  }
  
  ### Define groups colours and line widths
  group.colours <- c("ind" = linecolours.in[1], "mean" = linecolours.in[2], "p5" = linecolours.in[3], "p95"  = linecolours.in[3], "p50" = linecolours.in[4])
  group.lwd <-c("ind" = linewidths.in[1], "mean" = linewidths.in[2], "p5" = linewidths.in[2], "p95" = linewidths.in[2], "mean" = linewidths.in[2])
  
  ### Create object for output to be stored in
  plot.calib.out <- vector("list", length(anal.methods.sim.inn))
  names(plot.calib.out) <- anal.methods.sim.inn
  
  ### Define function to create the plots (it's done in a seperate environment, to reduce size of sized .RData file)
  create.ggplot <- function(data.in, group.colours.in, group.lwd.in, font.size.in, x.lim.in, y.lim.in){
    env <- new.env(parent = globalenv())
    env$subset <- data.in
    env$group.colours <- group.colours.in
    env$group.lwd <- group.lwd.in
    env$font.size <- font.size.in
    env$x.lim <- x.lim.in
    env$y.lim <- y.lim.in
    
    ggplot.out <- with(env, {
      ggplot.out <- ggplot(subset) +
        geom_line(aes(x = predrisk, y = value, group = variable, color = map.var, size = map.var), show.legend = FALSE) +  
        scale_colour_manual(values = group.colours) +
        scale_size_manual(values = group.lwd) +
        geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "magenta", size = 0.9) + 
        xlim(x.lim) + xlab("Predicted probability") +
        ylim(y.lim) + ylab("Observed probability") + theme_bw(base_size = font.size)
      return(ggplot.out)
    })
    
    return(ggplot.out)
  }
  
  for (i in 1:length(anal.methods.sim.inn)){
    print(anal.methods.sim.inn[i])
    plot.calib.out[[anal.methods.sim.inn[i]]] <- create.ggplot(analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]], 
                                                               group.colours, group.lwd, font.size.in, x.lim.in, y.lim.in)
  }
  
  return(plot.calib.out)
}


##########################################################################################################################################
### 6.2) Create ggplots from analysed data (output from function 5), but only have n simulations runs of the calibration curves showing ###
##########################################################################################################################################
create.individual.ggplots.ncurves <- function(analysed.data, x.lim.in, y.lim.in, font.size.in, linecolours.in, linewidths.in, anal.methods.sim.inn, n){

  ### Start by remove sims from data so that only n of the individual calibration curves remain
  ## Create a new object to store the new datasets in
  analysed.data.temp <- vector("list", length(anal.methods.sim.inn))
  names(analysed.data.temp) <- paste("data.calib.",anal.methods.sim.inn, sep = "")

  ### Cycle through and create the dataset with n caibration curves
  for (i in 1:length(anal.methods.sim.inn)){
    ## Create a temporary dataset
    temp.data <- analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]]
    ## Extract the rows which contain summary data
    temp.data.summary <- temp.data[temp.data$variable %in% c("obs.mean", "obs.p5", "obs.p95"), ]
    ## Extract the length of each calibration curve (required for subsetting exactly n of these)
    length.curve <- as.numeric(sum(temp.data$variable == "obs.mean"))
    ## Extract n*length.curve rows
    temp.data.curves <- temp.data[1:(n*length.curve), ]
    ## Combine into a single dataset and assign to the storage object
    analysed.data.temp[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]] <- rbind(temp.data.curves, temp.data.summary)
  }
  
  ### Rename storage object to fit in with code, and remove temporary file
  analysed.data <- analysed.data.temp
  rm(analysed.data.temp, temp.data, temp.data.summary, length.curve, temp.data.curves)

  ### Define groups colours and line widths
  group.colours <- c("ind" = linecolours.in[1], "mean" = linecolours.in[2], "p5" = linecolours.in[3], "p95"  = linecolours.in[3])
  group.lwd <-c("ind" = linewidths.in[1], "mean" = linewidths.in[2], "p5" = linewidths.in[2], "p95" = linewidths.in[2])
  
  ### Create object for output to be stored in
  plot.calib.out <- vector("list", length(anal.methods.sim.inn))
  names(plot.calib.out) <- anal.methods.sim.inn
  
  ### Define function to create the plots (it's done in a seperate environment, to reduce size of sized .RData file)
  create.ggplot <- function(data.in, group.colours.in, group.lwd.in, font.size.in, x.lim.in, y.lim.in){
    env <- new.env(parent = globalenv())
    env$subset <- data.in
    env$group.colours <- group.colours.in
    env$group.lwd <- group.lwd.in
    env$font.size <- font.size.in
    env$x.lim <- x.lim.in
    env$y.lim <- y.lim.in
    
    ggplot.out <- with(env, {
      ggplot.out <- ggplot(subset) +
        geom_line(aes(x = predrisk, y = value, group = variable, color = map.var, size = map.var), show.legend = FALSE) +  
        scale_colour_manual(values = group.colours) +
        scale_size_manual(values = group.lwd) +
        geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "magenta", size = 0.9) + 
        xlim(x.lim) + xlab("Predicted probability") +
        ylim(y.lim) + ylab("Observed probability") + theme_bw(base_size = font.size)
      return(ggplot.out)
    })
    
    return(ggplot.out)
  }
  
  for (i in 1:length(anal.methods.sim.inn)){
    print(anal.methods.sim.inn[i])
    plot.calib.out[[anal.methods.sim.inn[i]]] <- create.ggplot(analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]], 
                                                               group.colours, group.lwd, font.size.in, x.lim.in, y.lim.in)
  }
  
  return(plot.calib.out)
}


###########################################################################
### 6.3) Create ggplots from analysed data (output from function 5),    ###
### that only show the mean/median/p5/p95 summary lines, but for each   ###
### analysis method superimposed on the same plot                       ###
###########################################################################
create.summary.ggplots <- function(analysed.data, x.lim.in, y.lim.in, xlab.in, ylab.in, font.size.in, linecolours.in, lwd.in, 
                                   abline.size.in, anal.methods.sim.inn, sum.vec.in){
  
  
  ## Create empty objects for median difference from y=x, and for p90 - p5
  ## Reduce each analysed dataset to only contian summary data and add a variable which indicates which analysis method it was
  for (i in 1:length(anal.methods.sim.inn)){
    print(anal.methods.sim.inn[i])
    analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]] <- analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]] %>%
      ## Reduce
      filter(variable == sum.vec.in) %>%
      ## Add variable
      mutate(anal.method = anal.methods.sim.inn[i])
  }
  
  ### Now combine into a single dataset
  temp.summary.data <- analysed.data[paste("data.calib.",anal.methods.sim.inn, sep = "")]
  temp.summary.data <- do.call("rbind", temp.summary.data)
  
  ### Define function to create the plots (it's done in a seperate environment, to reduce size of .RData file)
  create.ggplot <- function(data.in, group.colours.in, lwd.fix.in, font.size.in, abline.size.in, x.lim.in, y.lim.in, xlab.in, ylab.in){
    env <- new.env(parent = globalenv())
    env$subset <- data.in
    env$group.colours <- group.colours.in
    env$lwd.fix <- lwd.fix.in
    env$font.size <- font.size.in
    env$x.lim <- x.lim.in
    env$y.lim <- y.lim.in
    env$xlab.in <- xlab.in
    env$ylab.in <- ylab.in
    env$abline.size <- abline.size.in
    
    ggplot.out <- with(env, {
      ggplot.out <- ggplot(subset) +
        geom_abline(intercept = 0, slope = 1, color = "black", size = abline.size) + 
        geom_line(aes(x = predrisk, y = value, group = anal.method, color = anal.method, linetype = anal.method), 
                  size = lwd.fix, show.legend = TRUE) +  
        #scale_colour_brewer(palette = group.colours) +
        #scale_linetype_manual(values = group.lty) +
        #scale_size_manual(values = group.lwd) +
        xlim(x.lim) + xlab(xlab.in) +
        ylim(y.lim) + ylab(ylab.in) + theme_bw(base_size = font.size)
      return(ggplot.out)
    })
    
    return(ggplot.out)
  }
  
  
  ### Create the ggplot
  plot.calib.out <- create.ggplot(data.in = temp.summary.data, 
                                  group.colours.in = linecolours.in, 
                                  lwd.fix.in = lwd.in,
                                  font.size.in = font.size.in, 
                                  abline.size.in = abline.size.in, 
                                  x.lim.in = x.lim.in, 
                                  y.lim.in = y.lim.in,
                                  xlab.in = xlab.in,
                                  ylab.in = ylab.in)
  
  return(plot.calib.out)
}



###########################################################################
### 6.4) Create ggplots from analysed data (output from function 5),    ###
### that only show the mean/median/p5/p95 summary lines, but for each   ###
### analysis method superimposed on the same plot, but horizontally     ###
###########################################################################
create.summary.horizontal.ggplots <- function(analysed.data, x.lim.in, y.lim.in, xlab.in, ylab.in, font.size.in, linecolours.in, lwd.in, 
                                              abline.size.in, anal.methods.sim.inn, sum.vec.in){
  
  #     analysed.data <- res.DGM.clay
  #     x.lim.in <- x.lim
  #     y.lim.in <- y.lim
  #     font.size.in <- font.size
  #     linecolours.in <- 1:8
  #     linetyps.in <- linetyps
  #     linewidths.in <- lwds
  #     anal.methods.sim.inn <- anal.methods.sim
  #     sum.vec.in <- sum.vec
  
  ## Create empty objects for median difference from y=x, and for p90 - p5
  ## Reduce each analysed dataset to only contian summary data and add a variable which indicates which analysis method it was
  for (i in 1:length(anal.methods.sim.inn)){
    print(anal.methods.sim.inn[i])
    analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]] <- analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]] %>%
      ## Reduce
      filter(variable == sum.vec.in) %>%
      ## Add variable
      mutate(anal.method = anal.methods.sim.inn[i])
    
    ## Create an extra step to change names of analysis methods to be in line with manuscript
    if (anal.methods.sim.inn[i] == "joint"){
      analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]]$anal.method <- "dual-o"
    } else if (anal.methods.sim.inn[i] == "f.gamma_weib"){
      analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]]$anal.method <- "f-gam"
    } else if (anal.methods.sim.inn[i] == "f.normal_weib"){
      analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]]$anal.method <- "f-norm"
    } else if (anal.methods.sim.inn[i] == "c.clay"){
      analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]]$anal.method <- "c-clay"
    } else if (anal.methods.sim.inn[i] == "c.gumb"){
      analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]]$anal.method <- "c-gumb"
    } else if (anal.methods.sim.inn[i] == "c.frank"){
      analysed.data[[paste("data.calib.",anal.methods.sim.inn[i], sep = "")]]$anal.method <- "c-frank"
    }
  }

  ### Now combine into a single dataset
  temp.summary.data <- analysed.data[paste("data.calib.",anal.methods.sim.inn, sep = "")]
  temp.summary.data <- do.call("rbind", temp.summary.data)
  
  ### Define function to create the plots (it's done in a seperate environment, to reduce size of .RData file)
  create.ggplot <- function(data.in, group.colours.in, lwd.fix.in, font.size.in, abline.size.in, x.lim.in, y.lim.in, xlab.in, ylab.in){
    env <- new.env(parent = globalenv())
    env$subset <- data.in
    env$group.colours <- group.colours.in
    env$lwd.fix <- lwd.fix.in
    env$font.size <- font.size.in
    env$x.lim <- x.lim.in
    env$y.lim <- y.lim.in
    env$xlab.in <- xlab.in
    env$ylab.in <- ylab.in
    env$abline.size <- abline.size.in
    
    ggplot.out <- with(env, {
      ggplot.out <- ggplot(subset) +
        geom_abline(intercept = 0, slope = 0, color = "black", size = abline.size) + 
        geom_line(aes(x = predrisk, y = value, group = anal.method, color = anal.method, linetype = anal.method), 
                  size = lwd.fix, show.legend = TRUE) +  
        #scale_colour_brewer(palette = group.colours) +
        #scale_linetype_manual(values = group.lty) +
        #scale_size_manual(values = group.lwd) +
        xlim(x.lim) + xlab(xlab.in) +
        ylim(y.lim) + ylab(ylab.in) + theme_bw(base_size = font.size) + 
        guides(color = guide_legend(title = "Analysis\nMethod"),
               linetype = guide_legend(title = "Analysis\nMethod"))
      return(ggplot.out)
    })
    
    return(ggplot.out)
  }
  
  
  ### Create the ggplot
  plot.calib.out <- create.ggplot(data.in = temp.summary.data, 
                                  group.colours.in = linecolours.in, 
                                  lwd.fix.in = lwd.in,
                                  font.size.in = font.size.in, 
                                  abline.size.in = abline.size.in, 
                                  x.lim.in = x.lim.in, 
                                  y.lim.in = y.lim.in,
                                  xlab.in = xlab.in,
                                  ylab.in = ylab.in)
  
  return(plot.calib.out)
}

###################################################################################################
### 6.5) Define function to calculate ggplot, wiht superimposed lines, but just for product method 
###################################################################################################

create.individual.ggplots.product <- function(analysed.data, x.lim.in, y.lim.in, font.size.in, linecolours.in, linewidths.in, anal.methods.sim.inn){
  
#   analysed.data <- res.DGM.nocorr
#   x.lim.in <- x.lim
#   y.lim.in <- y.lim
#   font.size.in <- font.size
#   linecolours.in <- linecolours
#   linewidths.in <- linewidths
#   anal.methods.sim.inn <- anal.methods.sim
  
  ## Extract n.sim.suc
  n.sim.suc <- analysed.data$n.sim.suc

  ## Extra relevant data (product method)
  analysed.data <- analysed.data[["data.calib.product"]]
  ## Extract the rows which contain summary data
  analysed.data <- analysed.data[analysed.data$map.var %in% c("ind", "p5", "p95", "p50"), ]
  ## Re-order the data so that p5, p95 and p50 are plotted before ind
  
  ### Define groups colours and line widths
  group.colours <- c("ind" = linecolours.in[1], "mean" = linecolours.in[2], "p5" = linecolours.in[3], "p95"  = linecolours.in[3], "p50" = linecolours.in[2], "perfect" = "purple")
  group.lwd <-c("ind" = linewidths.in[1], "mean" = linewidths.in[2], "p5" = linewidths.in[2], "p95" = linewidths.in[2], "p50" = linewidths.in[2], "perfect" = 0.7)

  ### Remove mean from group.colours as it isnt being plotted
  group.colours <- group.colours[names(group.colours) != "mean"]

  ### Add a straight line to dataset
  straight.line.data <- data.frame("predrisk" = seq(0, 1,0.01), 
                                   "value" = seq(0, 1, 0.01), 
                                   "variable" = rep("perfect", length(seq(0, 1, 0.01))),
                                   "map.var" = rep("perfect", length(seq(0, 1, 0.01)))) 
  analysed.data <- rbind(straight.line.data, analysed.data)
  
  ### Create a factor variable for the different groups being plotted and give it the correct order
  analysed.data$variable <- factor(analysed.data$variable, levels = c(paste("obs.smooth.", 1:n.sim.suc, sep = ""), "perfect", "obs.p50", "obs.p5", "obs.p95"))
  ### Define function to create the plots (it's done in a seperate environment, to reduce size of sized .RData file)
  create.ggplot <- function(data.in, group.colours.in, group.lwd.in, font.size.in, x.lim.in, y.lim.in){
    env <- new.env(parent = globalenv())
    env$subset <- data.in
    env$group.colours <- group.colours.in
    env$group.lwd <- group.lwd.in
    env$font.size <- font.size.in
    env$x.lim <- x.lim.in
    env$y.lim <- y.lim.in
    
    ggplot.out <- with(env, {
      ggplot.out <- ggplot(subset) +
        geom_line(aes(x = predrisk, y = value, group = variable, color = map.var, size = map.var)) +  
        scale_colour_manual(values = group.colours, labels = c("Single\niteration", "p5", "p95", "Median", "1:1")) +
        scale_size_manual(values = group.lwd, guide = "none") +
        #geom_abline(aes(intercept = 0, slope = 1, color = "purple"), size = 0.9) + 
        xlim(x.lim) + xlab("Predicted probability") +
        ylim(y.lim) + ylab("Observed probability") + theme_bw(base_size = font.size) + labs(color = "Calibraiton\ncurve")
      return(ggplot.out)
    })
    
    return(ggplot.out)
  }
  
  ## Create ggplot
  plot.calib.out <- create.ggplot(analysed.data, group.colours, group.lwd, font.size.in, x.lim.in, y.lim.in)
   
  return(plot.calib.out)
}



#########################################################################################################################
### 7.1) Create a function to load the plots of interest for each scenario and sample size, and combine into one ggplot ###
#########################################################################################################################
create.combined.calib.ggplot <- function(scenario, n.size, suffix.to.save){
  
  ### Load worksapce
  load(paste("data/sim_results_", scenario, "_n", n.size, "v1000.RData", sep = ""))
  
  ### Create the plots and save them
  ## DGM MSM
  DGM.msm.plots <- ggarrange(ggplots.DGM.msm[["product"]], ggplots.DGM.msm[["joint"]], ggplots.DGM.msm[["msm"]], 
                             ggplots.DGM.msm[["c.clay"]],ggplots.DGM.msm[["c.gumb"]],ggplots.DGM.msm[["c.frank"]],
                             ggplots.DGM.msm[["f.normal_weib"]],ggplots.DGM.msm[["f.gamma_weib"]],
                             labels = grid.labels.in, nrow = 2, ncol = 4, 
                             font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.msm.calib.plots.png", sep = ""), DGM.msm.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM clay
  DGM.clay.plots <- ggarrange(ggplots.DGM.clay[["product"]], ggplots.DGM.clay[["joint"]], ggplots.DGM.clay[["msm"]], 
                              ggplots.DGM.clay[["c.clay"]],ggplots.DGM.clay[["c.gumb"]],ggplots.DGM.clay[["c.frank"]],
                              ggplots.DGM.clay[["f.normal_weib"]],ggplots.DGM.clay[["f.gamma_weib"]],
                              labels = grid.labels.in, nrow = 2, ncol = 4, 
                              font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.clay.calib.plots.png", sep = ""), DGM.clay.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM gumb
  DGM.gumb.plots <- ggarrange(ggplots.DGM.gumb[["product"]], ggplots.DGM.gumb[["joint"]], ggplots.DGM.gumb[["msm"]], 
                              ggplots.DGM.gumb[["c.clay"]],ggplots.DGM.gumb[["c.gumb"]],ggplots.DGM.gumb[["c.frank"]],
                              ggplots.DGM.gumb[["f.normal_weib"]],ggplots.DGM.gumb[["f.gamma_weib"]],
                              labels = grid.labels.in, nrow = 2, ncol = 4, 
                              font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.gumb.calib.plots.png", sep = ""), DGM.gumb.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM frank
  DGM.frank.plots <- ggarrange(ggplots.DGM.frank[["product"]], ggplots.DGM.frank[["joint"]], ggplots.DGM.frank[["msm"]], 
                               ggplots.DGM.frank[["c.clay"]],ggplots.DGM.frank[["c.gumb"]],ggplots.DGM.frank[["c.frank"]],
                               ggplots.DGM.frank[["f.normal_weib"]],ggplots.DGM.frank[["f.gamma_weib"]],
                               labels = grid.labels.in, nrow = 2, ncol = 4, 
                               font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.frank.calib.plots.png", sep = ""), DGM.frank.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM normal
  DGM.normal.plots <- ggarrange(ggplots.DGM.normal[["product"]], ggplots.DGM.normal[["joint"]], ggplots.DGM.normal[["msm"]], 
                                ggplots.DGM.normal[["c.clay"]],ggplots.DGM.normal[["c.gumb"]],ggplots.DGM.normal[["c.frank"]],
                                ggplots.DGM.normal[["f.normal_weib"]],ggplots.DGM.normal[["f.gamma_weib"]],
                                labels = grid.labels.in, nrow = 2, ncol = 4, 
                                font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.normal.calib.plots.png", sep = ""), DGM.normal.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM gamma
  DGM.gamma.plots <- ggarrange(ggplots.DGM.gamma[["product"]], ggplots.DGM.gamma[["joint"]], ggplots.DGM.gamma[["msm"]], 
                               ggplots.DGM.gamma[["c.clay"]],ggplots.DGM.gamma[["c.gumb"]],ggplots.DGM.gamma[["c.frank"]],
                               ggplots.DGM.gamma[["f.normal_weib"]],ggplots.DGM.gamma[["f.gamma_weib"]],
                               labels = grid.labels.in, nrow = 2, ncol = 4, 
                               font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.gamma.calib.plots.png", sep = ""), DGM.gamma.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
}


#########################################################################################################################
### 7.2) Create a function to load the plots of interest for the no correlation scenario, and combine into one ggplot ###
#########################################################################################################################
create.combined.calib.ggplot.nocorr <- function(scenario, n.size, suffix.to.save){
  
  ### Load worksapce
  load(paste("data/sim_results_", scenario, "_n", n.size, "v1000.RData", sep = ""))
  
  ### Create the plots and save them
  ## DGM MSM
  DGM.nocorr.plots <- ggarrange(ggplots.DGM.nocorr[["product"]], ggplots.DGM.nocorr[["joint"]], ggplots.DGM.nocorr[["msm"]], 
                                ggplots.DGM.nocorr[["c.clay"]],ggplots.DGM.nocorr[["c.gumb"]],ggplots.DGM.nocorr[["c.frank"]],
                                ggplots.DGM.nocorr[["f.normal_weib"]],ggplots.DGM.nocorr[["f.gamma_weib"]],
                                labels = grid.labels.in, nrow = 2, ncol = 4, 
                                font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.nocorr.calib.plots.png", sep = ""), DGM.nocorr.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
}


####################################################################################################
### 7.3) Same as 7.1, but doesn't load in data, assuming relevant files are already in workspace ###
####################################################################################################
create.combined.calib.ggplot.noload <- function(scenario, n.size, suffix.to.save){
  

  
  ### Create the plots and save them
  ## DGM MSM
  DGM.msm.plots <- ggarrange(ggplots.DGM.msm[["product"]], ggplots.DGM.msm[["joint"]], ggplots.DGM.msm[["msm"]], 
                             ggplots.DGM.msm[["c.clay"]],ggplots.DGM.msm[["c.gumb"]],ggplots.DGM.msm[["c.frank"]],
                             ggplots.DGM.msm[["f.normal_weib"]],ggplots.DGM.msm[["f.gamma_weib"]],
                             labels = grid.labels.in, nrow = 2, ncol = 4, 
                             font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.msm.calib.plots.png", sep = ""), DGM.msm.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM clay
  DGM.clay.plots <- ggarrange(ggplots.DGM.clay[["product"]], ggplots.DGM.clay[["joint"]], ggplots.DGM.clay[["msm"]], 
                              ggplots.DGM.clay[["c.clay"]],ggplots.DGM.clay[["c.gumb"]],ggplots.DGM.clay[["c.frank"]],
                              ggplots.DGM.clay[["f.normal_weib"]],ggplots.DGM.clay[["f.gamma_weib"]],
                              labels = grid.labels.in, nrow = 2, ncol = 4, 
                              font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.clay.calib.plots.png", sep = ""), DGM.clay.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM gumb
  DGM.gumb.plots <- ggarrange(ggplots.DGM.gumb[["product"]], ggplots.DGM.gumb[["joint"]], ggplots.DGM.gumb[["msm"]], 
                              ggplots.DGM.gumb[["c.clay"]],ggplots.DGM.gumb[["c.gumb"]],ggplots.DGM.gumb[["c.frank"]],
                              ggplots.DGM.gumb[["f.normal_weib"]],ggplots.DGM.gumb[["f.gamma_weib"]],
                              labels = grid.labels.in, nrow = 2, ncol = 4, 
                              font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.gumb.calib.plots.png", sep = ""), DGM.gumb.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM frank
  DGM.frank.plots <- ggarrange(ggplots.DGM.frank[["product"]], ggplots.DGM.frank[["joint"]], ggplots.DGM.frank[["msm"]], 
                               ggplots.DGM.frank[["c.clay"]],ggplots.DGM.frank[["c.gumb"]],ggplots.DGM.frank[["c.frank"]],
                               ggplots.DGM.frank[["f.normal_weib"]],ggplots.DGM.frank[["f.gamma_weib"]],
                               labels = grid.labels.in, nrow = 2, ncol = 4, 
                               font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.frank.calib.plots.png", sep = ""), DGM.frank.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM normal
  DGM.normal.plots <- ggarrange(ggplots.DGM.normal[["product"]], ggplots.DGM.normal[["joint"]], ggplots.DGM.normal[["msm"]], 
                                ggplots.DGM.normal[["c.clay"]],ggplots.DGM.normal[["c.gumb"]],ggplots.DGM.normal[["c.frank"]],
                                ggplots.DGM.normal[["f.normal_weib"]],ggplots.DGM.normal[["f.gamma_weib"]],
                                labels = grid.labels.in, nrow = 2, ncol = 4, 
                                font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.normal.calib.plots.png", sep = ""), DGM.normal.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
  ## DGM gamma
  DGM.gamma.plots <- ggarrange(ggplots.DGM.gamma[["product"]], ggplots.DGM.gamma[["joint"]], ggplots.DGM.gamma[["msm"]], 
                               ggplots.DGM.gamma[["c.clay"]],ggplots.DGM.gamma[["c.gumb"]],ggplots.DGM.gamma[["c.frank"]],
                               ggplots.DGM.gamma[["f.normal_weib"]],ggplots.DGM.gamma[["f.gamma_weib"]],
                               labels = grid.labels.in, nrow = 2, ncol = 4, 
                               font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.gamma.calib.plots.png", sep = ""), DGM.gamma.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
}


####################################################################################################
### 7.4) Same as 7.2, but doesn't load in data, assuming relevant files are already in workspace ###
####################################################################################################
create.combined.calib.ggplot.nocorr.noload <- function(scenario, n.size, suffix.to.save){
  
  ### Create the plots and save them
  ## DGM MSM
  DGM.nocorr.plots <- ggarrange(ggplots.DGM.nocorr[["product"]], ggplots.DGM.nocorr[["joint"]], ggplots.DGM.nocorr[["msm"]], 
                                ggplots.DGM.nocorr[["c.clay"]],ggplots.DGM.nocorr[["c.gumb"]],ggplots.DGM.nocorr[["c.frank"]],
                                ggplots.DGM.nocorr[["f.normal_weib"]],ggplots.DGM.nocorr[["f.gamma_weib"]],
                                labels = grid.labels.in, nrow = 2, ncol = 4, 
                                font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in)
  
  ggsave(paste("figures/", scenario, "_n", n.size, "v1000.curves", suffix.to.save, ".DGM.nocorr.calib.plots.png", sep = ""), DGM.nocorr.plots, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
  
}


####################################################################################################
### 7.5) Combines the summary ggplots, stratified by DGM  ###
####################################################################################################
create.combined.calib.summary.ggplot <- function(ggplots.summary.in, scenario, n.size, suffix.to.save){
  
  ### Create the plots and save them
  ## DGM MSM
  ggplots.summary.arranged <- ggarrange(ggplots.summary.in[["DGM.msm"]], ggplots.summary.in[["DGM.clay"]], ggplots.summary.in[["DGM.gumb"]], 
                                        ggplots.summary.in[["DGM.frank"]], ggplots.summary.in[["DGM.normal"]], ggplots.summary.in[["DGM.gamma"]],
                                        labels = grid.labels.in, nrow = 2, ncol = 3, 
                                        font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in, 
                                        common.legend = TRUE, legend = "right")
  
  ggsave(paste("figures/gg.summary.", scenario, "_n", n.size, "v1000.", suffix.to.save, ".png", sep = ""), 
         ggplots.summary.arranged, type = "cairo-png", 
         dpi = 300,  width = 1143, height = 529, units = "mm")
}


create.combined.calib.summary.ggplot.horiz <- function(ggplots.summary.in, scenario, n.size, suffix.to.save){
  
  ### Create the plots and save them
  ## DGM MSM
  ggplots.summary.arranged <- ggarrange(ggplots.summary.in[["DGM-1:MSM"]], ggplots.summary.in[["DGM-2:Clayton"]], ggplots.summary.in[["DGM-3:Gumbel"]], 
                                        ggplots.summary.in[["DGM-4:Frank"]], ggplots.summary.in[["DGM-5:Normal"]], ggplots.summary.in[["DGM-6:Gamma"]],
                                        labels = grid.labels.in, nrow = 2, ncol = 3, 
                                        font.label = list(size = font.label.size.in), label.x = label.x.in, label.y = label.y.in, 
                                        common.legend = TRUE, legend = "right")
  
  ggsave(paste("figures/gg_2.summary.horiz.", scenario, "_n", n.size, "v1000.", suffix.to.save, ".tiff", sep = ""), 
         ggplots.summary.arranged, device = "tiff", 
         dpi = 800, width = 10, height = 4.62, units = "in", compression = 'lzw')
  
  ggsave(paste("figures/gg_2.summary.horiz.", scenario, "_n", n.size, "v1000.", suffix.to.save, ".png", sep = ""), 
         ggplots.summary.arranged, device = "png", 
         dpi = 600, width = 10, height = 4.62, units = "in")
  
  ggsave(paste("figures/gg_2.summary.horiz.", scenario, "_n", n.size, "v1000.", suffix.to.save, ".pdf", sep = ""), 
         ggplots.summary.arranged, device = "pdf", 
         width = 10, height = 4.62, units = "in")
}


# ##########################################################################
# ### OLD version of create.dat.calib.plot.run, which didn't do ggplots in seperate environment, making .RData file massive ###
# ##########################################################################
# create.dat.calib.plot.run_OLD2 <- function(data.inn, n.knot.in, x.lim.in, y.lim.in, font.size.in, linecolours.in, linewidths.in, pred.eval.in, anal.methods.sim.inn){
#   
#   #   data.inn <- dgm.msm.out.comb[[3]][1:25]
#   #   x.lim.in <- x.lim
#   #   y.lim.in <- y.lim
#   #   font.size.in <- font.size
#   #   pred.eval.in <- pred.eval
#   #   anal.methods.sim.inn <- anal.methods.sim
#   #   pred.eval.in <- pred.eval
#   
#   ### Let k.sim.suc.in be the number of successful simulations
#   k.sim.suc.in <- length(data.inn)
#   #k.sim.suc.in <- 50
#   
#   
#   ### Create an object to store all the converted datasets
#   dat.calib.plot.list <- vector("list", k.sim.suc.in)
#   
#   ### Apply the function defined above to each piece of simulations output
#   for (i in 1:k.sim.suc.in){
#     print(paste("analyse", i, sep = ""))
#     dat.calib.plot.list[[i]] <- create.dat.calib.plot(data.in = data.inn[[i]], pred.eval.in = pred.eval.in, n.knot.in = n.knot.in, 
#                                                       anal.methods.sim.in = anal.methods.sim.inn)
#   }
#   
#   ### Combine all the datasets into one dataset, so all the observed risks are in one place (each row is a predicted risk)
#   print("memA")
#   dat.calib.plot <- do.call(cbind, dat.calib.plot.list)
#   
#   ### Seperate out the observed risks for each analysis method into different datasets
#   print("memB")
#   dat.calib.plot.seperate <- vector("list", length(anal.methods.sim.inn))
#   names(dat.calib.plot.seperate) <- anal.methods.sim.inn
#   
#   for (i in 1:length(anal.methods.sim.inn)){
#     ### Seperate them out
#     print("memC")
#     temp.data <- dat.calib.plot[ , colnames(dat.calib.plot) == paste("obs.smooth.", anal.methods.sim.inn[i], sep= "")]
#     colnames(temp.data) <- paste("obs.smooth.", 1:k.sim.suc.in, sep = "")
#     print("memD")
#     ### Now calculate row means and percentiles from each dataset, and also store the vector of predicted risks
#     ### They must be calculated and put in an object, then added in after, to not impact the calculation of percentiles and mean
#     
#     # Calculate
#     print("memE")
#     obs.mean <- rowMeans(temp.data, na.rm = TRUE)
#     print("memF")
#     obs.p5 <- apply( temp.data, 1, quantile, probs = 0.05, na.rm = TRUE)
#     print("memG")
#     obs.p95 <- apply( temp.data, 1, quantile, probs = 0.95, na.rm = TRUE)
#     
#     # Add to object
#     temp.data$predrisk <- pred.eval.in
#     temp.data$obs.mean <- obs.mean
#     temp.data$obs.p5 <- obs.p5
#     temp.data$obs.p95 <- obs.p95
#     print("memH")
#     ### Create long data for plot
#     temp.data <- melt(temp.data, id = "predrisk")
#     
#     ### Want to create a new variable to map to colours, line thickness, etc, which will be same for all individual plots, then different for the mean and percentiles
#     print("memI")
#     map.var <- rep(0, nrow(temp.data))
#     print("memJ")
#     map.var[grep("obs.smooth", temp.data$variable, fixed = TRUE)] <- "ind"
#     print("memK")
#     map.var[temp.data$variable == "obs.p5"] <- "p5"
#     map.var[temp.data$variable == "obs.p95"] <- "p95"
#     map.var[temp.data$variable == "obs.mean"] <- "mean"
#     temp.data$map.var <- map.var
#     
#     dat.calib.plot.seperate[[i]] <- temp.data
#   }
#   
#   
#   ####################################
#   ### Fit the ggplots to this data ###
#   ####################################
#   
#   ### Define groups colours and line widths
#   group.colours <- c("ind" = linecolours.in[1], "mean" = linecolours.in[2], "p5" = linecolours.in[3], "p95"  = linecolours.in[3])
#   group.lwd <-c("ind" = linewidths.in[1], "mean" = linewidths.in[2], "p5" = linewidths.in[2], "p95" = linewidths.in[2])
#   
#   ### Create object for output to be stored in
#   plot.calib.out <- vector("list", length(anal.methods.sim.inn))
#   names(plot.calib.out) <- paste("plot.calib.", anal.methods.sim.inn[i], sep= "")
#   
#   ### Define function to create the plots (probably shouldn't do this within another function)
#   create.ggplot <- function(data.in){
#     ggplot(data.in) +
#       geom_line(aes(x = predrisk, y = value, group = variable, color = map.var, size = map.var), show.legend = FALSE) +  
#       scale_colour_manual(values = group.colours) +
#       scale_size_manual(values = group.lwd) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "magenta", size = 0.9) + 
#       xlim(x.lim.in) + xlab("Predicted probability") +
#       ylim(y.lim.in) + ylab("Observed probability") + theme_bw(base_size = font.size.in)
#   }
#   
#   plot.calib.out[["plot.calib.product"]] <- create.ggplot(dat.calib.plot.seperate[["product"]])
#   plot.calib.out[["plot.calib.joint"]] <- create.ggplot(dat.calib.plot.seperate[["joint"]])
#   plot.calib.out[["plot.calib.msm"]] <- create.ggplot(dat.calib.plot.seperate[["msm"]])
#   plot.calib.out[["plot.calib.c.clay"]] <- create.ggplot(dat.calib.plot.seperate[["c.clay"]])
#   plot.calib.out[["plot.calib.c.gumb"]] <- create.ggplot(dat.calib.plot.seperate[["c.gumb"]])
#   plot.calib.out[["plot.calib.c.frank"]] <- create.ggplot(dat.calib.plot.seperate[["c.frank"]])
#   plot.calib.out[["plot.calib.f.normal"]] <- create.ggplot(dat.calib.plot.seperate[["f.normal"]])
#   plot.calib.out[["plot.calib.f.gamma"]] <- create.ggplot(dat.calib.plot.seperate[["f.gamma"]])
#   
#   return(plot.calib.out)
# }
# 
# 
# ##########################################################################
# ### OLD version of create.dat.calib.plot.run, which just plots mean, p5 and p95 ###
# ##########################################################################
# create.dat.calib.plot.run_OLD <- function(data.inn, n.knot.in, x.lim.in, y.lim.in, font.size.in, pred.eval.in, anal.methods.sim.inn){
#   
#   ### Let k.sim.suc.in be the number of successful simulations
#   k.sim.suc.in <- length(data.inn)
#   
#   ### Create an object to store all the converted datasets
#   dat.calib.plot.list <- vector("list", k.sim.suc.in)
#   
#   ### Apply the function defined above to each piece of simulations output
#   for (i in 1:k.sim.suc.in){
#     print(paste("analyse", i, sep = ""))
#     dat.calib.plot.list[[i]] <- create.dat.calib.plot(data.in = data.inn[[i]], pred.eval.in = pred.eval.in, n.knot.in = n.knot.in, 
#                                                       anal.methods.sim.in = anal.methods.sim.inn)
#   }
#   
#   ### Combine all the datasets into one dataset, so all the observed risks are in one place (each row is a predicted risk)
#   dat.calib.plot <- do.call(cbind, dat.calib.plot.list)
#   
#   ### Seperate out the observed risks for each analysis method into different datasets
#   dat.calib.plot.seperate <- vector("list", length(anal.methods.sim.inn))
#   names(dat.calib.plot.seperate) <- anal.methods.sim.inn
#   
#   for (i in 1:length(anal.methods.sim.inn)){
#     dat.calib.plot.seperate[[i]] <- dat.calib.plot[ , colnames(dat.calib.plot) == paste("obs.smooth.", anal.methods.sim.inn[i], sep= "")]
#   }
#   
#   ### Calculate row means and percentiles from each dataset, and store into another dataframe, which will be used for calibration plot
#   dat.calib.plot.av <- vector("list", length(anal.methods.sim.inn))
#   names(dat.calib.plot.av) <- anal.methods.sim.inn
#   
#   for (i in 1:length(anal.methods.sim.inn)){
#     # Create data frame to store averages/percentiles
#     dat.calib.plot.av[[i]] <- data.frame("predrisk" = pred.eval.in)
#     
#     # Calculate and assign
#     dat.calib.plot.av[[i]][ , "obs.mean"] <- rowMeans(dat.calib.plot.seperate[[i]])
#     dat.calib.plot.av[[i]][ , "obs.p5"] <- apply( dat.calib.plot.seperate[[i]], 1, quantile, probs = 0.05, na.rm = TRUE)
#     dat.calib.plot.av[[i]][ , "obs.p95"] <- apply( dat.calib.plot.seperate[[i]], 1, quantile, probs = 0.95, na.rm = TRUE)
#   }
#   
#   
#   ####################################
#   ### Fit the ggplots to this data ###
#   ####################################
#   
#   plot.calib.out <- vector("list", length(anal.methods.sim.inn))
#   names(plot.calib.out) <- paste("plot.calib.", anal.methods.sim.inn[i], sep= "")
#   
#   
#   plot.calib.out[["plot.calib.msm"]] <- ggplot(dat.calib.plot.av[["msm"]]) +
#     geom_line(aes(x = predrisk, y = obs.mean)) + 
#     geom_line(aes(x = predrisk, y = obs.p95)) + 
#     geom_line(aes(x = predrisk, y = obs.p5)) + 
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
#     xlim(x.lim.in) + xlab("Predicted probability") +
#     ylim(y.lim.in) + ylab("Observed probability") + theme_bw(base_size = font.size.in)
#   
#   
#   plot.calib.out[["plot.calib.product"]] <- ggplot(dat.calib.plot.av[["product"]]) +
#     geom_line(aes(x = predrisk, y = obs.mean)) + 
#     geom_line(aes(x = predrisk, y = obs.p95)) + 
#     geom_line(aes(x = predrisk, y = obs.p5)) + 
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
#     xlim(x.lim.in) + xlab("Predicted probability") +
#     ylim(y.lim.in) + ylab("Observed probability") + theme_bw(base_size = font.size.in)
#   
#   
#   plot.calib.out[["plot.calib.joint"]] <- ggplot(dat.calib.plot.av[["joint"]]) +
#     geom_line(aes(x = predrisk, y = obs.mean)) + 
#     geom_line(aes(x = predrisk, y = obs.p95)) + 
#     geom_line(aes(x = predrisk, y = obs.p5)) + 
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
#     xlim(x.lim.in) + xlab("Predicted probability") +
#     ylim(y.lim.in) + ylab("Observed probability") + theme_bw(base_size = font.size.in)
#   
#   plot.calib.out[["plot.calib.c.clay"]] <- ggplot(dat.calib.plot.av[["c.clay"]]) +
#     geom_line(aes(x = predrisk, y = obs.mean)) + 
#     geom_line(aes(x = predrisk, y = obs.p95)) + 
#     geom_line(aes(x = predrisk, y = obs.p5)) + 
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
#     xlim(x.lim.in) + xlab("Predicted probability") +
#     ylim(y.lim.in) + ylab("Observed probability") + theme_bw(base_size = font.size.in)
#   
#   plot.calib.out[["plot.calib.c.gumb"]] <- ggplot(dat.calib.plot.av[["c.gumb"]]) +
#     geom_line(aes(x = predrisk, y = obs.mean)) + 
#     geom_line(aes(x = predrisk, y = obs.p95)) + 
#     geom_line(aes(x = predrisk, y = obs.p5)) + 
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
#     xlim(x.lim.in) + xlab("Predicted probability") +
#     ylim(y.lim.in) + ylab("Observed probability") + theme_bw(base_size = font.size.in)
#   
#   plot.calib.out[["plot.calib.c.frank"]] <- ggplot(dat.calib.plot.av[["c.frank"]]) +
#     geom_line(aes(x = predrisk, y = obs.mean)) + 
#     geom_line(aes(x = predrisk, y = obs.p95)) + 
#     geom_line(aes(x = predrisk, y = obs.p5)) + 
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
#     xlim(x.lim.in) + xlab("Predicted probability") +
#     ylim(y.lim.in) + ylab("Observed probability") + theme_bw(base_size = font.size.in)
#   
#   plot.calib.out[["plot.calib.f.normal_weib"]] <- ggplot(dat.calib.plot.av[["f.normal_weib"]]) +
#     geom_line(aes(x = predrisk, y = obs.mean)) + 
#     geom_line(aes(x = predrisk, y = obs.p95)) + 
#     geom_line(aes(x = predrisk, y = obs.p5)) + 
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
#     xlim(x.lim.in) + xlab("Predicted probability") +
#     ylim(y.lim.in) + ylab("Observed probability") + theme_bw(base_size = font.size.in)
#   
#   plot.calib.out[["plot.calib.f.gamma_weib"]] <- ggplot(dat.calib.plot.av[["f.gamma_weib"]]) +
#     geom_line(aes(x = predrisk, y = obs.mean)) + 
#     geom_line(aes(x = predrisk, y = obs.p95)) + 
#     geom_line(aes(x = predrisk, y = obs.p5)) + 
#     geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
#     xlim(x.lim.in) + xlab("Predicted probability") +
#     ylim(y.lim.in) + ylab("Observed probability") + theme_bw(base_size = font.size.in)
#   
#   return(plot.calib.out)
# }