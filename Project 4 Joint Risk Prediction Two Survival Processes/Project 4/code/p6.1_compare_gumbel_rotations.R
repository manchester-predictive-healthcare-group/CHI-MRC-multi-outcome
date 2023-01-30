###
### Check which rotation of copula provides best fit
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Source packages
source("Project 4/code/sim_function_load_packages.R")
### Load functions
source("Project 4/code/sim_function_clinical_example.R")


### Assessment will be based on calibration plots with 5 cubic splines

###
### Assess model C, 200000
###

### Load predicted risks
load(paste("Project 4/data/clinical_example_fit_models", "C", "_n", 200000,".RData"))

### Need to turn data.valid into wide format, so we can extract time_C and status_C
data.valid.plot <- tidyr::pivot_wider(data.valid, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status)) 

### Check for NA's
for (i in 1:4){
  print(sum(is.na(predrisk.gumbel.list[[i]][["risk.joint.est"]])))
}

### Create calibration plots for the predicted risks, for ones which fitted properly
dat.validplot.list <- vector("list", 4)
for (i in c(1,3)){
  dat.validplot.list[[i]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                             status.in = data.valid.plot$status_C,
                                             risk.pred.in = predrisk.gumbel.list[[i]][["risk.joint.est"]],
                                             n.knots.in = 0,
                                             t.eval = t.eval,
                                             model.name = names(predrisk.gumbel.list)[i])
}

### I should then be combining these into a single dataset for ggplot!
data.validplot <- rbind(dat.validplot.list[[1]], dat.validplot.list[[3]])
str(data.validplot)
### Create ggplot
validplot <- ggplot(aes(x = risk.pred, y = risk.obs, colour = model), data = data.validplot) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed")
validplot
