###
### Assess calibration using each approach
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_6")
getwd()

### Load packages
source("code/z_load_packages.R")

### Read in arguments
args <- commandArgs(trailingOnly = T)
n.devel <- as.numeric(args[1])
n.pctls <- as.numeric(args[2])
print(paste("n.devel = ", n.devel, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))

### Load predicted risks
load(paste("data/ce_combine_probtrans_N", n.devel, ".RData", sep = ""))
load(paste("data/ce_pv_combine_N", n.devel, "_npctls", n.pctls, ".RData", sep = ""))
source("code/z_functions_ce.R")

###
### For N = 3000, there are some negative predictions
### Remove these
###
data.raw <- subset(data.raw, pstate1 > 0 & pstate2 > 0 & pstate3 > 0 & pstate4 > 0 & pstate5 > 0 & pstate6 > 0)
pv.comb <- subset(pv.comb, person_id %in% data.raw$person_id)
data.mstate <- subset(data.mstate, person_id %in% data.raw$person_id)
str(data.raw)
str(pv.comb)
str(data.mstate)

############################
### Moderate calibration ###
############################

###
### Calc calibration using Aalen-Johansen estimator
###

###
### Do it seperately for each transition
###

### Create objects to store plot data
aj.decile.pred <- vector("list", 6)
aj.decile.obs <- vector("list", 6)

### Create object to store sub-cohorts
data.mstate.deciles <- vector("list", 6)
data.raw.deciles <- vector("list", 6)

### Assign number of groups
num_groups <- 10

### Loop through states
for (i in 1:6){
  
  print(paste("i = ", i, Sys.time()))
  
  ### Create output vector for observed and predicted
  aj.decile.pred[[i]] <- vector(mode = "numeric", num_groups)
  aj.decile.obs[[i]] <- vector(mode = "numeric", num_groups)
  
  ### Arrange data by the risk for each transition
  data.temp <- arrange(data.raw, !! rlang::sym(paste("pstate", i, sep = "")))
  
  ### Create num_groups subcohorts
  data.mstate.deciles[[i]] <- vector("list", num_groups)
  data.raw.deciles[[i]] <- data.temp %>%
    group_by((row_number() - 1) %/% (n()/num_groups)) %>%
    nest %>% pull(data)
  
  ### Calculate AJ within each cohort, and the mean predicted risk
  for (j in 1:num_groups){
    
    print(paste("j = ", j, Sys.time()))
    
    ### Turn into dataframe from tibble
    data.raw.deciles[[i]][[j]] <- data.frame(data.raw.deciles[[i]][[j]])
    
    ### Create mstate dataframes, based on who is in the data.raw.deciles dataframes
    data.mstate.deciles[[i]][[j]] <- subset(data.mstate, person_id %in% data.raw.deciles[[i]][[j]]$person_id)
    
    ### Calc AJ
    obs.aj.temp <- calc.calib.aj.ce(data.mstate = data.mstate.deciles[[i]][[j]],
                                    tmat = tmat, 
                                    t.eval = t.eval)
    obs.aj.temp <- obs.aj.temp[["obs.aj"]]
    
    ### Assign predicted and observed risk
    aj.decile.obs[[i]][j] <- as.numeric(obs.aj.temp[paste("pstate", i, sep = "")])
    aj.decile.pred[[i]][j] <- mean(select(data.raw.deciles[[i]][[j]], paste("pstate", i, sep = ""))[,1])
    
  }
  
}


###
### Create list of plots for AJ
###

### Create empty lists to store data and plots in
calib.plots.data.aj <- vector("list", 6)
calib.plots.aj <- vector("list", 6)

### Run loop and create plots
for (i in 1:6){
  ## Create dataset for plot
  calib.plots.data.aj[[i]] <- data.frame("pred" = unlist(aj.decile.pred[[i]]), 
                                   "obs" = unlist(aj.decile.obs[[i]]))
  
  ## Create plot
  calib.plots.aj[[i]] <- ggplot(data = calib.plots.data.aj[[i]], aes(x = pred, y = obs, color = "red")) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    xlab("Predicted risk") + ylab("Observed risk") + ggtitle(paste("State", i, sep = " ")) + 
    theme(legend.position = "none")
}

print(paste("AJ DONE", Sys.time()))
   

######################################################################
### Calc calibration using binary logistic recalibration framework
######################################################################

calib.mod.blr.ipcw <- calc.calib.blr.ipcw.mod.ce(data.mstate = data.mstate, 
                                                 data.raw = data.raw, 
                                                 t.eval = t.eval, 
                                                 p.est = select(data.raw, paste("pstate", 1:6, sep = "")), 
                                                 max.weight = 10)
print(paste("BLR DONE", Sys.time()))


######################################################################
### Calc calibration using multinomial logistic recalibration framework
######################################################################

calib.mod.mlr.ipcw <- calc.calib.mlr.ipcw.mod.ce(data.mstate = data.mstate, 
                                                 data.raw = data.raw, 
                                                 t.eval = t.eval, 
                                                 p.est = select(data.raw, paste("pstate", 1:6, sep = "")), 
                                                 ps.int = 4, degree = 3, max.weight = 10)
print(paste("MLR DONE", Sys.time()))


######################################################################
### Calc calibration using pseudo-values
######################################################################

### Make sure data.raw is arranged by person_id
data.raw <- arrange(data.raw, person_id)

### Calculate calibration
calib.mod.pv <- calc.calib.pv.moderate.ce(data.raw = data.raw, 
                                       pv.comb = pv.comb, 
                                       p.est = select(data.raw, paste("pstate", 1:6, sep = "")))
print(paste("PV DONE", Sys.time()))

### Save image
print("IMAGED SAVED")
save.image(paste("data/ce_assess_calib_N", n.devel, "_npctls", n.pctls, ".RData", sep = ""))


#################
### Create plot
#################

### Create a list containing the plots 
plots.comb.list <- c(calib.plots.aj,
                     calib.mod.pv[["plots.list"]],
                     calib.mod.blr.ipcw[["plots.list"]],
                     calib.mod.mlr.ipcw[["plots.list"]])

### Arrange into one using ggarrange
plots.comb.gg <- ggarrange(plotlist = plots.comb.list, nrow = 4, ncol = 6)
plots.comb.gg <- 
  annotate_figure(plots.comb.gg, left = "           MLR-IPCW                                    BLR-IPCW                                    PV                                    AJ             ")

### Save plot
CairoPNG(paste("figures/gg_ce_moderate_N", n.devel, ".png", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(plots.comb.gg)
dev.off()
CairoTIFF(paste("figures/gg_ce_moderate_N", n.devel, ".tiff", sep = ""), 
         dpi = 300, width = 15, height = 10, unit = "in")
print(plots.comb.gg)
dev.off()
CairoPDF(paste("figures/gg_ce_moderate_N", n.devel, ".pdf", sep = ""), 
         width = 15, height = 10)
print(plots.comb.gg)
dev.off()


### Save image
print("IMAGED SAVED")
save.image(paste("data/ce_assess_calib_N", n.devel, "_npctls", n.pctls, ".RData", sep = ""))


########################
### Mean calibration ###
########################
print("START MEAN")

######################################################################
### Calc cailbration using Aalen-Johansen estimator
######################################################################

### Calculate mean observed risk using Aalen-Johansen
obs.aj <- sapply(aj.decile.obs, mean)

### Haven't bootstrapped standard errors for computational reasons. Focus is not on mean calibration.
obs.aj.se <- rep(0, 6)

### Calculate mean calibration
calib.mean.aj <- obs.aj - colMeans(select(data.raw, paste("pstate", 1:6, sep = "")))

print(paste("AJ DONE", Sys.time()))

######################################################################
### Calc calibration using binary logistic recalibration framework
######################################################################

### Calculate the mean difference in predicted and observed
calib.mean.blr.ipcw <- calc.calib.blr.ipcw.ce(data.mstate = data.mstate, 
                                                 data.raw = data.raw, 
                                                 t.eval = t.eval, 
                                                 p.est = select(data.raw, paste("pstate", 1:6, sep = "")), 
                                                 max.weight = 10)
      
### Calculate the standard error of this quantity
calib.mean.blr.ipcw.se <- calc.calib.blr.ipcw.boot.ce(data.mstate = data.mstate, 
                                              data.raw = data.raw, 
                                              t.eval = t.eval, 
                                              p.est = select(data.raw, paste("pstate", 1:6, sep = "")), 
                                              max.weight = 10,
                                              n.boot = 200)
print(paste("BLR DONE", Sys.time()))


######################################################################
### Calc calibration using multinomial logistic recalibration framework
######################################################################

### Calculate the mean difference in predicted and observed
calib.mean.mlr.ipcw <- calc.calib.mlr.ipcw.ce(data.mstate = data.mstate, 
                                                 data.raw = data.raw, 
                                                 t.eval = t.eval, 
                                                 p.est = select(data.raw, paste("pstate", 1:6, sep = "")), 
                                                 max.weight = 10)

### Calculate the standard error of this quantity
calib.mean.mlr.ipcw.se <- calc.calib.mlr.ipcw.boot.ce(data.mstate = data.mstate, 
                                                      data.raw = data.raw, 
                                                      t.eval = t.eval, 
                                                      p.est = select(data.raw, paste("pstate", 1:6, sep = "")), 
                                                      max.weight = 10,
                                                      n.boot = 200)
print(paste("MLR DONE", Sys.time()))

######################################################################
### Calibration using pseudo values
######################################################################

calib.mean.pv <- colMeans(pv.comb[, -1]) - colMeans(select(data.raw, paste("pstate", 1:6, sep = "")))
print(paste("PV DONE", Sys.time()))


### Save image
save.image(paste("data/ce_assess_calib_N", n.devel, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGED SAVED")


#####################################
### Create plot for mean calibration
#####################################

###
### Create lists to store data for plots
###
tables.est <- vector("list", 6)
tables.est.se <- vector("list", 6)

###
### Create data for plots
###
for (state in 1:6){
  tables.est[[state]] <- round(c(
    calib.mean.blr.ipcw$diff.pred.obs[[state]],
    calib.mean.mlr.ipcw$diff.pred.obs[[state]],
    calib.mean.pv[[state]]), 3)
  
  tables.est.se[[state]] <- round(c(calib.mean.blr.ipcw.se[["se"]][paste("diff.pred.obs", state, sep = "")],
                                    calib.mean.mlr.ipcw.se[["se"]][paste("diff.pred.obs", state, sep = "")],
                                    0), 3)
}


###
### Create table to be fed into dot and whisker plot in correct format
###
data.dw <- vector("list", 6)
names(data.dw) <- paste("state", 1:6, sep = "")
for (state in 1:6){
  data.dw[[state]] <- data.frame("term" = rep("Mean calibration", 3), 
                                 "model" = c("BLR-IPCW", "MLR-IPCW", "PV"),
                                 "estimate" = c(tables.est[[state]], tables.est[[state]], tables.est[[state]]),
                                 "std.error" = c(tables.est.se[[state]], tables.est.se[[state]], tables.est.se[[state]]),
                                 "state" = rep(paste("State", state, sep = " "), 3))
}


###
### Combine all data for a single plot with facet_grid
###
data.dw.facet <- do.call("rbind", data.dw)
gg.dw <- dwplot(data.dw.facet) +# scale_x_continuous(breaks = seq(-0.15, 0, 0.05)) 
  facet_wrap(vars(state), nrow = 2, ncol = 3)  + 
  #xlab("Bias") + 
  ylab("Scenario")

### Save main plot
CairoPNG(paste("figures/gg_ce_mean_N", n.devel, "_npctls", n.pctls, ".png", sep = ""), 
         dpi = 300, width = 15, height = 3, unit = "in")
print(gg.dw)
dev.off()
CairoPDF(paste("figures/gg_ce_mean_N", n.devel, "_npctls", n.pctls, ".pdf", sep = ""), 
         width = 15, height = 3)
print(gg.dw)
dev.off()
CairoTIFF(paste("figures/gg_ce_mean_N", n.devel, "_npctls", n.pctls, ".tiff", sep = ""), 
         dpi = 300, width = 15, height = 3, unit = "in")
print(gg.dw)
dev.off()

