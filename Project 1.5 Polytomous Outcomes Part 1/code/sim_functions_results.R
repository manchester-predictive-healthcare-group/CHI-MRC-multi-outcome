#################################################################
### 4.1A) Create table for results from large sample analysis ###
#################################################################

create.table.large.sample <- function(data.in, dp){
  table1 <- cbind(
    #rbind(format(round(data.in$OvE.multinomial, dp), nsmall = dp),
          #format(round(data.in$OvE.seqlog, dp), nsmall = dp),
          #format(round(data.in$OvE.OvA, dp), nsmall = dp),
          #format(round(data.in$OvE.OvO.PC, dp), nsmall = dp)),
    
#     c(paste(format(round(data.in$wc.po.multinomial$alpha[1], dp), nsmall = dp), "/", format(round(data.in$wc.po.multinomial$beta[1], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.po.seqlog$alpha[1], dp), nsmall = dp), "/", format(round(data.in$wc.po.seqlog$beta[1], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.po.OvA$alpha[1], dp), nsmall = dp), "/", format(round(data.in$wc.po.OvA$beta[1], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.po.OvO.PC$alpha[1], dp), nsmall = dp), "/", format(round(data.in$wc.po.OvO.PC$beta[1], dp), nsmall = dp), sep = "")),
#     
#     c(paste(format(round(data.in$wc.po.multinomial$alpha[2], dp), nsmall = dp), "/", format(round(data.in$wc.po.multinomial$beta[2], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.po.seqlog$alpha[2], dp), nsmall = dp), "/", format(round(data.in$wc.po.seqlog$beta[2], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.po.OvA$alpha[2], dp), nsmall = dp), "/", format(round(data.in$wc.po.OvA$beta[2], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.po.OvO.PC$alpha[2], dp), nsmall = dp), "/", format(round(data.in$wc.po.OvO.PC$beta[2], dp), nsmall = dp), sep = "")),
#     
#     c(paste(format(round(data.in$wc.po.multinomial$alpha[3], dp), nsmall = dp), "/", format(round(data.in$wc.po.multinomial$beta[3], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.po.seqlog$alpha[3], dp), nsmall = dp), "/", format(round(data.in$wc.po.seqlog$beta[3], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.po.OvA$alpha[3], dp), nsmall = dp), "/", format(round(data.in$wc.po.OvA$beta[3], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.po.OvO.PC$alpha[3], dp), nsmall = dp), "/", format(round(data.in$wc.po.OvO.PC$beta[3], dp), nsmall = dp), sep = "")),
    
#     c(paste(format(round(data.in$wc.ms.multinomial$alpha[1], dp), nsmall = dp), "/", format(round(data.in$wc.ms.multinomial$beta[1], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.ms.seqlog$alpha[1], dp), nsmall = dp), "/", format(round(data.in$wc.ms.seqlog$beta[1], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.ms.OvA$alpha[1], dp), nsmall = dp), "/", format(round(data.in$wc.ms.OvA$beta[1], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.ms.OvO.PC$alpha[1], dp), nsmall = dp), "/", format(round(data.in$wc.ms.OvO.PC$beta[1], dp), nsmall = dp), sep = "")),
#     
#     c(paste(format(round(data.in$wc.ms.multinomial$alpha[2], dp), nsmall = dp), "/", format(round(data.in$wc.ms.multinomial$beta[2], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.ms.seqlog$alpha[2], dp), nsmall = dp), "/", format(round(data.in$wc.ms.seqlog$beta[2], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.ms.OvA$alpha[2], dp), nsmall = dp), "/", format(round(data.in$wc.ms.OvA$beta[2], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.ms.OvO.PC$alpha[2], dp), nsmall = dp), "/", format(round(data.in$wc.ms.OvO.PC$beta[2], dp), nsmall = dp), sep = "")),
#     
#     c(paste(format(round(data.in$wc.ms.multinomial$alpha[3], dp), nsmall = dp), "/", format(round(data.in$wc.ms.multinomial$beta[3], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.ms.seqlog$alpha[3], dp), nsmall = dp), "/", format(round(data.in$wc.ms.seqlog$beta[3], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.ms.OvA$alpha[3], dp), nsmall = dp), "/", format(round(data.in$wc.ms.OvA$beta[3], dp), nsmall = dp), sep = ""),
#       paste(format(round(data.in$wc.ms.OvO.PC$alpha[3], dp), nsmall = dp), "/", format(round(data.in$wc.ms.OvO.PC$beta[3], dp), nsmall = dp), sep = "")),
    
    c(format(round(data.in$calib.flex.mlr.multinomial$ECI, dp), nsmall = dp), 
      format(round(data.in$calib.flex.mlr.seqlog$ECI, dp), nsmall = dp), 
      format(round(data.in$calib.flex.mlr.OvA$ECI, dp), nsmall = dp), 
      format(round(data.in$calib.flex.mlr.OvO.PC$ECI, dp), nsmall = dp)),
    
    rbind(format(round(data.in$Cstat.multinomial, dp), nsmall = dp),
          format(round(data.in$Cstat.seqlog, dp), nsmall = dp),
          format(round(data.in$Cstat.OvA, dp), nsmall = dp),
          format(round(data.in$Cstat.OvO.PC, dp), nsmall = dp)),
    
    c(format(round(data.in$PDI.multinomial, dp), nsmall = dp), 
      format(round(data.in$PDI.seqlog, dp), nsmall = dp), 
      format(round(data.in$PDI.OvA, dp), nsmall = dp), 
      format(round(data.in$PDI.OvO.PC, dp), nsmall = dp))
  )
  
  colnames(table1) <- c("ECI", "AUC.Y=1", "AUC.Y=2", "AUC.Y=3", "AUC.Y=4", "AUC.Y=5", "PDI")
  
  return(table1)
}




##################################################################################################################################
### 4.2) Function to take flexible calibration plot data from the simulation, and turn it into long format ready to be plotted ###
##################################################################################################################################

turn.data.long <- function(data.in){
  ### Create a variable with number of outcome categories
  K <- as.numeric(max(data.in$Y))
  
  ## First first 1000 rows of dataset
  data.in <- data.in[1:1000, ]
  ## Create id variable
  data.in$id <- 1:nrow(data.in)
  ## Do the melting and return the variable
  if (K == 3){
    data.out <- data.frame(melt(data.in, id.vars = c("id","Y"), measure.vars = c("p1", "p2", "p3"), variable.name = "outcome.cat", value.name = "pred"),
                           "obs" = melt(data.in, id.vars = c("id"), measure.vars = c("p.obs1", "p.obs2", "p.obs3"), value.name = "obs")[, c("obs")])
  } else if (K == 5){
    data.out <- data.frame(melt(data.in, id.vars = c("id","Y"), measure.vars = c("p1", "p2", "p3", "p4", "p5"), variable.name = "outcome.cat", value.name = "pred"),
                           "obs" = melt(data.in, id.vars = c("id"), measure.vars = c("p.obs1", "p.obs2", "p.obs3", "p.obs4", "p.obs5"), value.name = "obs")[, c("obs")])
  }
  return(data.out)
}



#####################################################################################################################
### 4.3) Create plot for flexible calibration plots using polytomous recalibration framework of van Hoorde (2014) ###
#####################################################################################################################

create.plot.flex.mlr <- function(data.in){
  
  ### Create long dataset for calibration data from each analysis method
  dat.calib.flex.mlr.multinomial <- turn.data.long(data.in$calib.flex.mlr.multinomial$calib.data)
  dat.calib.flex.mlr.seqlog <- turn.data.long(data.in$calib.flex.mlr.seqlog$calib.data)
  dat.calib.flex.mlr.OvA <- turn.data.long(data.in$calib.flex.mlr.OvA$calib.data)
  dat.calib.flex.mlr.OvO.PC <- turn.data.long(data.in$calib.flex.mlr.OvO.PC$calib.data)
  
  ### Create new combined dataset so can plot all methods in one graph/grid
  dat.calib.flex.mlr <- rbind(data.frame(dat.calib.flex.mlr.multinomial, "anal.meth" = rep("MLR", nrow(dat.calib.flex.mlr.multinomial))),
                              data.frame(dat.calib.flex.mlr.seqlog, "anal.meth" = rep("c-ratio", nrow(dat.calib.flex.mlr.seqlog))),
                              data.frame(dat.calib.flex.mlr.OvA, "anal.meth" = rep("OvA-N", nrow(dat.calib.flex.mlr.OvA))),
                              data.frame(dat.calib.flex.mlr.OvO.PC, "anal.meth" = rep("OvO-PC", nrow(dat.calib.flex.mlr.OvO.PC))))
  dat.calib.flex.mlr$anal.meth <- factor(dat.calib.flex.mlr$anal.meth, levels = c("MLR", "c-ratio", "OvA-N", "OvO-PC"))
  
  ### Now lets make a ggplot
  plot.calib.flex.mlr <- ggplot(dat.calib.flex.mlr, aes(x = pred, y = obs)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(size = 0.5, color = "red") +
    xlab("Predicted") + ylab("Observed") +
    facet_grid(anal.meth ~ outcome.cat)
  
  return(plot.calib.flex.mlr)
}



################################################################################################################
### 4.4) Create plot for flexible calibration plots using binary logistic recalibraiton per outcome category ###
################################################################################################################

create.plot.flex.po <- function(data.in){
  
  ### Create long dataset for calibration data from each analysis method
  dat.calib.flex.po.multinomial <- turn.data.long(data.in$calib.flex.po.multinomial)
  dat.calib.flex.po.seqlog <- turn.data.long(data.in$calib.flex.po.seqlog)
  dat.calib.flex.po.OvA <- turn.data.long(data.in$calib.flex.po.OvA)
  dat.calib.flex.po.OvO.PC <- turn.data.long(data.in$calib.flex.po.OvO.PC)
  
  ### Create new combined dataset so can plot all methods in one graph/grid
  dat.calib.flex.po <- rbind(data.frame(dat.calib.flex.po.multinomial, "anal.meth" = rep("MLR", nrow(dat.calib.flex.po.multinomial))),
                              data.frame(dat.calib.flex.po.seqlog, "anal.meth" = rep("c-ratio", nrow(dat.calib.flex.po.seqlog))),
                              data.frame(dat.calib.flex.po.OvA, "anal.meth" = rep("OvA-N", nrow(dat.calib.flex.po.OvA))),
                              data.frame(dat.calib.flex.po.OvO.PC, "anal.meth" = rep("OvO-PC", nrow(dat.calib.flex.po.OvO.PC))))
  dat.calib.flex.po$anal.meth <- factor(dat.calib.flex.po$anal.meth, levels = c("MLR", "c-ratio", "OvA-N", "OvO-PC"))
  
  ### Now lets make a ggplot
  plot.calib.flex.po <- ggplot(dat.calib.flex.po, aes(x = pred, y = obs)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_line(color = "red") +
    xlab("Predicted") + ylab("Observed") +
    facet_grid(anal.meth ~ outcome.cat)
  
  return(plot.calib.flex.po)
}




#########################################################################################################################
### 5.1) Function to turn list output from small sample simulations, into a dataframe, so I can rbind all the results ###
### from different simulation runs. Used as part of function 5.3.                                                     ###
#########################################################################################################################
create.data.frame.per.row <- function(data.in){
  
  dat <- data.frame(
    ## OvE metrics
    "OvE.multinomial.p1" = as.numeric(data.in$OvE.multinomial["p1"]),
    "OvE.multinomial.p2" = as.numeric(data.in$OvE.multinomial["p2"]),
    "OvE.multinomial.p3" = as.numeric(data.in$OvE.multinomial["p3"]),
    "OvE.seqlog.p1" = as.numeric(data.in$OvE.seqlog["p1"]),
    "OvE.seqlog.p2" = as.numeric(data.in$OvE.seqlog["p2"]),
    "OvE.seqlog.p3" = as.numeric(data.in$OvE.seqlog["p3"]),
    "OvE.OvA.p1" = as.numeric(data.in$OvE.OvA["p1"]),
    "OvE.OvA.p2" = as.numeric(data.in$OvE.OvA["p2"]),
    "OvE.OvA.p3" = as.numeric(data.in$OvE.OvA["p3"]),
    "OvE.OvO.PC.p1" = as.numeric(data.in$OvE.OvO.PC["p1"]),
    "OvE.OvO.PC.p2" = as.numeric(data.in$OvE.OvO.PC["p2"]),
    "OvE.OvO.PC.p3" = as.numeric(data.in$OvE.OvO.PC["p3"]),
    
    ## Weak calibration per outcome metrics
    "wc.po.multinomial.p1.alpha" = as.numeric(data.in$wc.po.multinomial$alpha["p1"]),
    "wc.po.multinomial.p1.beta" = as.numeric(data.in$wc.po.multinomial$beta["p1"]),
    "wc.po.multinomial.p2.alpha" = as.numeric(data.in$wc.po.multinomial$alpha["p2"]),
    "wc.po.multinomial.p2.beta" = as.numeric(data.in$wc.po.multinomial$beta["p2"]),
    "wc.po.multinomial.p3.alpha" = as.numeric(data.in$wc.po.multinomial$alpha["p3"]),
    "wc.po.multinomial.p3.beta" = as.numeric(data.in$wc.po.multinomial$beta["p3"]),
    
    "wc.po.seqlog.p1.alpha" = -as.numeric(data.in$wc.po.seqlog$alpha["p1"]),
    "wc.po.seqlog.p1.beta" = as.numeric(data.in$wc.po.seqlog$beta["p1"]),
    "wc.po.seqlog.p2.alpha" = -as.numeric(data.in$wc.po.seqlog$alpha["p2"]),
    "wc.po.seqlog.p2.beta" = as.numeric(data.in$wc.po.seqlog$beta["p2"]),
    "wc.po.seqlog.p3.alpha" = -as.numeric(data.in$wc.po.seqlog$alpha["p3"]),
    "wc.po.seqlog.p3.beta" = as.numeric(data.in$wc.po.seqlog$beta["p3"]),
    
    "wc.po.OvA.p1.alpha" = as.numeric(data.in$wc.po.OvA$alpha["p1"]),
    "wc.po.OvA.p1.beta" = as.numeric(data.in$wc.po.OvA$beta["p1"]),
    "wc.po.OvA.p2.alpha" = as.numeric(data.in$wc.po.OvA$alpha["p2"]),
    "wc.po.OvA.p2.beta" = as.numeric(data.in$wc.po.OvA$beta["p2"]),
    "wc.po.OvA.p3.alpha" = as.numeric(data.in$wc.po.OvA$alpha["p3"]),
    "wc.po.OvA.p3.beta" = as.numeric(data.in$wc.po.OvA$beta["p3"]),
    
    "wc.po.OvO.PC.p1.alpha" = as.numeric(data.in$wc.po.OvO.PC$alpha["p1"]),
    "wc.po.OvO.PC.p1.beta" = as.numeric(data.in$wc.po.OvO.PC$beta["p1"]),
    "wc.po.OvO.PC.p2.alpha" = as.numeric(data.in$wc.po.OvO.PC$alpha["p2"]),
    "wc.po.OvO.PC.p2.beta" = as.numeric(data.in$wc.po.OvO.PC$beta["p2"]),
    "wc.po.OvO.PC.p3.alpha" = as.numeric(data.in$wc.po.OvO.PC$alpha["p3"]),
    "wc.po.OvO.PC.p3.beta" = as.numeric(data.in$wc.po.OvO.PC$beta["p3"]),
    
    ## Weak calibration model specific metrics
    "wc.ms.multinomial.lp1.alpha" = as.numeric(data.in$wc.ms.multinomial$alpha["lp1"]),
    "wc.ms.multinomial.lp1.beta" = as.numeric(data.in$wc.ms.multinomial$beta["lp1"]),
    "wc.ms.multinomial.lp2.alpha" = as.numeric(data.in$wc.ms.multinomial$alpha["lp2"]),
    "wc.ms.multinomial.lp2.beta" = as.numeric(data.in$wc.ms.multinomial$beta["lp2"]),
    "wc.ms.multinomial.lp3.alpha" = as.numeric(data.in$wc.ms.multinomial$alpha["lp3"]),
    "wc.ms.multinomial.lp3.beta" = as.numeric(data.in$wc.ms.multinomial$beta["lp3"]),
    
    "wc.ms.seqlog.lp1.alpha" = as.numeric(data.in$wc.ms.seqlog$alpha["lp1"]),
    "wc.ms.seqlog.lp1.beta" = as.numeric(data.in$wc.ms.seqlog$beta["lp1"]),
    "wc.ms.seqlog.lp2.alpha" = as.numeric(data.in$wc.ms.seqlog$alpha["lp2"]),
    "wc.ms.seqlog.lp2.beta" = as.numeric(data.in$wc.ms.seqlog$beta["lp2"]),
    "wc.ms.seqlog.lp3.alpha" = as.numeric(data.in$wc.ms.seqlog$alpha["lp3"]),
    "wc.ms.seqlog.lp3.beta" = as.numeric(data.in$wc.ms.seqlog$beta["lp3"]),
    
    "wc.ms.OvA.lp1.alpha" = as.numeric(data.in$wc.ms.OvA$alpha["lp1"]),
    "wc.ms.OvA.lp1.beta" = as.numeric(data.in$wc.ms.OvA$beta["lp1"]),
    "wc.ms.OvA.lp2.alpha" = as.numeric(data.in$wc.ms.OvA$alpha["lp2"]),
    "wc.ms.OvA.lp2.beta" = as.numeric(data.in$wc.ms.OvA$beta["lp2"]),
    "wc.ms.OvA.lp3.alpha" = as.numeric(data.in$wc.ms.OvA$alpha["lp3"]),
    "wc.ms.OvA.lp3.beta" = as.numeric(data.in$wc.ms.OvA$beta["lp3"]),
    
    "wc.ms.OvO.PC.lp1.alpha" = as.numeric(data.in$wc.ms.OvO.PC$alpha["lp1"]),
    "wc.ms.OvO.PC.lp1.beta" = as.numeric(data.in$wc.ms.OvO.PC$beta["lp1"]),
    "wc.ms.OvO.PC.lp2.alpha" = as.numeric(data.in$wc.ms.OvO.PC$alpha["lp2"]),
    "wc.ms.OvO.PC.lp2.beta" = as.numeric(data.in$wc.ms.OvO.PC$beta["lp2"]),
    "wc.ms.OvO.PC.lp3.alpha" = as.numeric(data.in$wc.ms.OvO.PC$alpha["lp3"]),
    "wc.ms.OvO.PC.lp3.beta" = as.numeric(data.in$wc.ms.OvO.PC$beta["lp3"]),
    
    ## C-statistic metrics
    "Cstat.multinomial.p1" = as.numeric(data.in$Cstat.multinomial["p1"]),
    "Cstat.multinomial.p2" = as.numeric(data.in$Cstat.multinomial["p2"]),
    "Cstat.multinomial.p3" = as.numeric(data.in$Cstat.multinomial["p3"]),
    "Cstat.seqlog.p1" = as.numeric(data.in$Cstat.seqlog["p1"]),
    "Cstat.seqlog.p2" = as.numeric(data.in$Cstat.seqlog["p2"]),
    "Cstat.seqlog.p3" = as.numeric(data.in$Cstat.seqlog["p3"]),
    "Cstat.OvA.p1" = as.numeric(data.in$Cstat.OvA["p1"]),
    "Cstat.OvA.p2" = as.numeric(data.in$Cstat.OvA["p2"]),
    "Cstat.OvA.p3" = as.numeric(data.in$Cstat.OvA["p3"]),
    "Cstat.OvO.PC.p1" = as.numeric(data.in$Cstat.OvO.PC["p1"]),
    "Cstat.OvO.PC.p2" = as.numeric(data.in$Cstat.OvO.PC["p2"]),
    "Cstat.OvO.PC.p3" = as.numeric(data.in$Cstat.OvO.PC["p3"]),
    
    ## PDI metrics
    "PDI.multinomial" = as.numeric(data.in$PDI.multinomial),
    "PDI.seqlog" = as.numeric(data.in$PDI.seqlog),
    "PDI.OvA" = as.numeric(data.in$PDI.OvA),
    "PDI.OvO.PC" = as.numeric(data.in$PDI.OvO.PC))
  
  return(dat)
}


################################################################################################################################
### 5.2) Function to take processed output from small sample simulation and create an output table. Input must be means      ###
### Processing is done using function 5.1, and is done within 5.3. This function is used as part of 5.3                      ###
################################################################################################################################
create.table.small.sample <- function(sim.out.means, dp){
  
  output.table <- cbind(
    
    ## OvE metrics
    c(format(round(sim.out.means["OvE.multinomial.p1"], dp), nsmall = dp),
      format(round(sim.out.means["OvE.seqlog.p1"], dp), nsmall = dp),
      format(round(sim.out.means["OvE.OvA.p1"], dp), nsmall = dp),
      format(round(sim.out.means["OvE.OvO.PC.p1"], dp), nsmall = dp)),
    
    c(format(round(sim.out.means["OvE.multinomial.p2"], dp), nsmall = dp),
      format(round(sim.out.means["OvE.seqlog.p2"], dp), nsmall = dp),
      format(round(sim.out.means["OvE.OvA.p2"], dp), nsmall = dp),
      format(round(sim.out.means["OvE.OvO.PC.p2"], dp), nsmall = dp)),
    
    c(format(round(sim.out.means["OvE.multinomial.p3"], dp), nsmall = dp),
      format(round(sim.out.means["OvE.seqlog.p3"], dp), nsmall = dp),
      format(round(sim.out.means["OvE.OvA.p3"], dp), nsmall = dp),
      format(round(sim.out.means["OvE.OvO.PC.p3"], dp), nsmall = dp)),
    
    ## Weak calibration per outcome metrics
    c(paste(format(round(sim.out.means["wc.po.multinomial.p1.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.multinomial.p1.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.po.seqlog.p1.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.seqlog.p1.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.po.OvA.p1.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.OvA.p1.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.po.OvO.PC.p1.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.OvO.PC.p1.beta"], dp), nsmall = dp), sep = "")),
    
    c(paste(format(round(sim.out.means["wc.po.multinomial.p2.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.multinomial.p2.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.po.seqlog.p2.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.seqlog.p2.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.po.OvA.p2.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.OvA.p2.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.po.OvO.PC.p2.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.OvO.PC.p2.beta"], dp), nsmall = dp), sep = "")),
    
    c(paste(format(round(sim.out.means["wc.po.multinomial.p3.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.multinomial.p3.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.po.seqlog.p3.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.seqlog.p3.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.po.OvA.p3.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.OvA.p3.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.po.OvO.PC.p3.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.po.OvO.PC.p3.beta"], dp), nsmall = dp), sep = "")),
    
    ## Weak calibration model specific metrics
    c(paste(format(round(sim.out.means["wc.ms.multinomial.lp1.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.multinomial.lp1.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.ms.seqlog.lp1.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.seqlog.lp1.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.ms.OvA.lp1.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.OvA.lp1.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.ms.OvO.PC.lp1.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.OvO.PC.lp1.beta"], dp), nsmall = dp), sep = "")),
    
    c(paste(format(round(sim.out.means["wc.ms.multinomial.lp2.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.multinomial.lp2.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.ms.seqlog.lp2.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.seqlog.lp2.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.ms.OvA.lp2.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.OvA.lp2.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.ms.OvO.PC.lp2.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.OvO.PC.lp2.beta"], dp), nsmall = dp), sep = "")),
    
    c(paste(format(round(sim.out.means["wc.ms.multinomial.lp3.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.multinomial.lp3.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.ms.seqlog.lp3.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.seqlog.lp3.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.ms.OvA.lp3.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.OvA.lp3.beta"], dp), nsmall = dp), sep = ""),
      paste(format(round(sim.out.means["wc.ms.OvO.PC.lp3.alpha"], dp), nsmall = dp), "/", 
            format(round(sim.out.means["wc.ms.OvO.PC.lp3.beta"], dp), nsmall = dp), sep = "")),
    
    ## Cstat per outcome metrics
    c(format(round(sim.out.means["Cstat.multinomial.p1"], dp), nsmall = dp),
      format(round(sim.out.means["Cstat.seqlog.p1"], dp), nsmall = dp),
      format(round(sim.out.means["Cstat.OvA.p1"], dp), nsmall = dp),
      format(round(sim.out.means["Cstat.OvO.PC.p1"], dp), nsmall = dp)),
    
    c(format(round(sim.out.means["Cstat.multinomial.p2"], dp), nsmall = dp),
      format(round(sim.out.means["Cstat.seqlog.p2"], dp), nsmall = dp),
      format(round(sim.out.means["Cstat.OvA.p2"], dp), nsmall = dp),
      format(round(sim.out.means["Cstat.OvO.PC.p2"], dp), nsmall = dp)),
    
    c(format(round(sim.out.means["Cstat.multinomial.p3"], dp), nsmall = dp),
      format(round(sim.out.means["Cstat.seqlog.p3"], dp), nsmall = dp),
      format(round(sim.out.means["Cstat.OvA.p3"], dp), nsmall = dp),
      format(round(sim.out.means["Cstat.OvO.PC.p3"], dp), nsmall = dp)),
    
    ## PDI metrics
    c(format(round(sim.out.means["PDI.multinomial"], dp), nsmall = dp),
      format(round(sim.out.means["PDI.seqlog"], dp), nsmall = dp),
      format(round(sim.out.means["PDI.OvA"], dp), nsmall = dp),
      format(round(sim.out.means["PDI.OvO.PC"], dp), nsmall = dp))
  )
  
  rownames(output.table) <- c("Mult", "c-ratio", "OvA-N", "OvO-PC")
  colnames(output.table) <- c("OvE.Y=1", "OvE.Y=2", "OvE.Y=3", "wc.po.Y=1", "wc.po.Y=2", "wc.po.Y=3", 
                              "wc.ms.LP1", "wc.ms.LP2", "wc.ms.LP3", "AUC.Y=1", "AUC.Y=2", "AUC.Y=3", "PDI")
  
  return(output.table)
}



##############################################################################################################################
### 5.3) This function is a wrapper for a couple of other functions, which in conjunction are used to create the table with
### means of each simulation output, for the small sample simulations. It is applied directly to the simulation output list
### after the different parallelised output lists have been concatenated
##############################################################################################################################
create.table.small.sample.wrapper <- function(sim.out, dp){

  ### Seperate out the results into each DGM, by extracting the 1st, 2nd, 3rd or 4th list elements
  sim.out.DGM.multinomial <- lapply(sim.out, FUN = function(x) x[[1]])
  sim.out.DGM.seqlog <- lapply(sim.out, FUN = function(x) x[[2]])
#   sim.out.DGM.OvA <- lapply(sim.out, FUN = function(x) x[[3]])
#   sim.out.DGM.OvO.PC <- lapply(sim.out, FUN = function(x) x[[4]])

  ### Turn output from simulation into a dataframe, so we can calculate means of each piece of output
  ## First turn the individual output from each simulation into a dataframe of one row
  sim.out.DGM.multinomial.dat <- lapply(sim.out.DGM.multinomial, create.data.frame.per.row)
  sim.out.DGM.seqlog.dat <- lapply(sim.out.DGM.seqlog, create.data.frame.per.row)
#   sim.out.DGM.OvA.dat <- lapply(sim.out.DGM.OvA, create.data.frame.per.row)
#   sim.out.DGM.OvO.PC.dat <- lapply(sim.out.DGM.OvO.PC, create.data.frame.per.row)
  
  ## Then rbind these rows into a single dataset
  sim.out.DGM.multinomial.dat <- do.call("rbind", sim.out.DGM.multinomial.dat)
  sim.out.DGM.seqlog.dat <- do.call("rbind", sim.out.DGM.seqlog.dat)
#   sim.out.DGM.OvA.dat <- do.call("rbind", sim.out.DGM.OvA.dat)
#   sim.out.DGM.OvO.PC.dat <- do.call("rbind", sim.out.DGM.OvO.PC.dat)
  
  ### Calculate the means for each of these variables across all simulations
  sim.out.DGM.multinomial.means <- colMeans(sim.out.DGM.multinomial.dat)
  sim.out.DGM.seqlog.means <- colMeans(sim.out.DGM.seqlog.dat)
#   sim.out.DGM.OvA.means <- colMeans(sim.out.DGM.OvA.dat)
#   sim.out.DGM.OvO.PC.means <- colMeans(sim.out.DGM.OvO.PC.dat)
  
  ### Create the final output tables
  sim.out.DGM.multinomial.table <- create.table.small.sample(sim.out.DGM.multinomial.means, dp)
  sim.out.DGM.seqlog.table <- create.table.small.sample(sim.out.DGM.seqlog.means, dp)
#   sim.out.DGM.OvA.table <- create.table.small.sample(sim.out.DGM.OvA.means, dp)
#   sim.out.DGM.OvO.PC.table <- create.table.small.sample(sim.out.DGM.OvO.PC.means, dp)
  
  ### Output this as a list with appropriate names
  output.list <- list(sim.out.DGM.multinomial.table,
                      sim.out.DGM.seqlog.table)
                      #sim.out.DGM.OvA.table,
                      #sim.out.DGM.OvO.PC.table)
  names(output.list) <- c("DGM.multinomial", "DGM.seqlog")
  
  return(output.list)
}