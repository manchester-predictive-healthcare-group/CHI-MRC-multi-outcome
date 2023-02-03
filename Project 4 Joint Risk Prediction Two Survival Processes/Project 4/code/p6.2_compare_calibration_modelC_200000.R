### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome")
getwd()

### Source packages
source("Project 4/code/sim_function_load_packages.R")
### Load functions
source("Project 4/code/sim_function_clinical_example.R")

### Set xlim and ylim for plots
xlim <- c(0, 0.15)
ylim <- c(0, 0.15)

xlim.zoom <- c(0, 0.06)
ylim.zoom <- c(0, 0.06)

### Assessment will be based on calibration plots with 5 cubic splines

###
### Define model type and sample size
###
model.type <- "C"
data.size <- 200000

### Load predicted risks
load(paste("Project 4/data/clinical_example_fit_models", model.type, "_n", data.size, ".RData", sep = ""))
### Load functions
source("Project 4/code/p4_clinical_example_functions.R")

#########################################################################################################################
#########################################################################################################################
### Step 1: Compare Clayton and Gumbel rotated 0 and 180 degrees, will use whichever is best for the final comparison ###
#########################################################################################################################
#########################################################################################################################

###
### DO CLAYTON FIRST
###

### Need to turn data.valid into wide format, so we can extract time_C and status_C
data.valid.plot <- tidyr::pivot_wider(data.valid, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status)) 


### Create calibration plots for the predicted risks, for ones which fitted properly (0 and 180 degrees)
dat.validplot.list <- vector("list", 4)
for (i in c(1,3)){
  dat.validplot.list[[i]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                        status.in = data.valid.plot$status_C,
                                        risk.pred.in = predrisk.clayton.list[[i]][["risk.joint.est"]],
                                        n.knots.in = 0,
                                        t.eval = t.eval,
                                        model.name = names(predrisk.clayton.list)[i])
}

### I should then be combining these into a single dataset for ggplot!
data.validplot <- rbind(dat.validplot.list[[1]], dat.validplot.list[[3]])

### Create ggplot
validplot <- ggplot(aes(x = risk.pred, y = risk.obs, colour = model), data = data.validplot) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed")
validplot

### ROTATED 180 IS BETTER

###
### NOW DO GUMBEL
###

###
### Gumbel 0 was the only one that would fit
###

# ### Create calibration plots for the predicted risks, for ones which fitted properly (0 and 180 degrees)
# dat.validplot.list <- vector("list", 4)
# for (i in c(1,3)){
#   dat.validplot.list[[i]] <- create.pred.obs(time.in = data.valid.plot$time_C,
#                                              status.in = data.valid.plot$status_C,
#                                              risk.pred.in = predrisk.gumbel.list[[i]][["risk.joint.est"]],
#                                              n.knots.in = 0,
#                                              t.eval = t.eval,
#                                              model.name = names(predrisk.gumbel.list)[i])
# }
# 
# ### I should then be combining these into a single dataset for ggplot!
# data.validplot <- rbind(dat.validplot.list[[1]], dat.validplot.list[[3]])
# 
# ### Create ggplot
# validplot <- ggplot(aes(x = risk.pred, y = risk.obs, colour = model), data = data.validplot) +
#   geom_line() +
#   geom_abline(intercept = 0, slope = 1, lty = "dashed")
# validplot

### ROTATED 0 DEGREES IS BETTER


########################################
########################################
### Step 2: Compare final comparison ###
########################################
########################################

### Want to produce calibraiton plot data for each method. I will then combine into a single dataset. I will do this for 0, 3 and 5 knots

###
### Product
###
dat.validplot.product.list <- vector("list", 4)

dat.validplot.product.list[[1]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                   status.in = data.valid.plot$status_C,
                                                   risk.pred.in = predrisk.product[["risk.joint.est"]],
                                                   n.knots.in = 0,
                                                   t.eval = t.eval,
                                                   model.name = "product")

dat.validplot.product.list[[2]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                   status.in = data.valid.plot$status_C,
                                                   risk.pred.in = predrisk.product[["risk.joint.est"]],
                                                   n.knots.in = 3,
                                                   t.eval = t.eval,
                                                   model.name = "product")

dat.validplot.product.list[[3]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                   status.in = data.valid.plot$status_C,
                                                   risk.pred.in = predrisk.product[["risk.joint.est"]],
                                                   n.knots.in = 5,
                                                   t.eval = t.eval,
                                                   model.name = "product")

dat.validplot.product.list[[4]] <- create.pred.obs.by.decile(time.in = data.valid.plot$time_C,
                                                         status.in = data.valid.plot$status_C,
                                                         risk.pred.in = predrisk.product[["risk.joint.est"]], 
                                                         t.eval = t.eval, model.name = "product")

###
### Dual
###
dat.validplot.dual.list <- vector("list", 4)

dat.validplot.dual.list[[1]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                status.in = data.valid.plot$status_C,
                                                risk.pred.in = predrisk.dual[["risk.joint.est"]],
                                                n.knots.in = 0,
                                                t.eval = t.eval,
                                                model.name = "dual-o")

dat.validplot.dual.list[[2]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                status.in = data.valid.plot$status_C,
                                                risk.pred.in = predrisk.dual[["risk.joint.est"]],
                                                n.knots.in = 3,
                                                t.eval = t.eval,
                                                model.name = "dual-o")

dat.validplot.dual.list[[3]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                status.in = data.valid.plot$status_C,
                                                risk.pred.in = predrisk.dual[["risk.joint.est"]],
                                                n.knots.in = 5,
                                                t.eval = t.eval,
                                                model.name = "dual-o")

dat.validplot.dual.list[[4]] <- create.pred.obs.by.decile(time.in = data.valid.plot$time_C,
                                                             status.in = data.valid.plot$status_C,
                                                             risk.pred.in = predrisk.dual[["risk.joint.est"]], 
                                                             t.eval = t.eval, model.name = "dual-o")


###
### Clayton
###

### Pick the Clayton model that fitted best
dat.validplot.clayton.list <- vector("list", 4)

dat.validplot.clayton.list[[1]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                   status.in = data.valid.plot$status_C,
                                                   risk.pred.in = predrisk.clayton.list[["r180"]][["risk.joint.est"]],
                                                   n.knots.in = 0,
                                                   t.eval = t.eval,
                                                   model.name = "c-clay")

dat.validplot.clayton.list[[2]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                   status.in = data.valid.plot$status_C,
                                                   risk.pred.in = predrisk.clayton.list[["r180"]][["risk.joint.est"]],
                                                   n.knots.in = 3,
                                                   t.eval = t.eval,
                                                   model.name = "c-clay")

dat.validplot.clayton.list[[3]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                   status.in = data.valid.plot$status_C,
                                                   risk.pred.in = predrisk.clayton.list[["r180"]][["risk.joint.est"]],
                                                   n.knots.in = 5,
                                                   t.eval = t.eval,
                                                   model.name = "c-clay")

dat.validplot.clayton.list[[4]] <- create.pred.obs.by.decile(time.in = data.valid.plot$time_C,
                                                             status.in = data.valid.plot$status_C,
                                                             risk.pred.in = predrisk.clayton.list[["r180"]][["risk.joint.est"]], 
                                                             t.eval = t.eval, model.name = "c-clay")

###
### Gumbel
###

### Pick the Gumbel model that fitted best
dat.validplot.gumbel.list <- vector("list", 4)

dat.validplot.gumbel.list[[1]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                  status.in = data.valid.plot$status_C,
                                                  risk.pred.in = predrisk.gumbel.list[["r0"]][["risk.joint.est"]],
                                                  n.knots.in = 0,
                                                  t.eval = t.eval,
                                                  model.name = "c-gumb")

dat.validplot.gumbel.list[[2]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                  status.in = data.valid.plot$status_C,
                                                  risk.pred.in = predrisk.gumbel.list[["r0"]][["risk.joint.est"]],
                                                  n.knots.in = 3,
                                                  t.eval = t.eval,
                                                  model.name = "c-gumb")

dat.validplot.gumbel.list[[3]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                  status.in = data.valid.plot$status_C,
                                                  risk.pred.in = predrisk.gumbel.list[["r0"]][["risk.joint.est"]],
                                                  n.knots.in = 5,
                                                  t.eval = t.eval,
                                                  model.name = "c-gumb")

dat.validplot.gumbel.list[[4]] <- create.pred.obs.by.decile(time.in = data.valid.plot$time_C,
                                                             status.in = data.valid.plot$status_C,
                                                             risk.pred.in = predrisk.gumbel.list[["r0"]][["risk.joint.est"]], 
                                                             t.eval = t.eval, model.name = "c-gumb")

###
### Frank
###
dat.validplot.frank.list <- vector("list", 4)

dat.validplot.frank.list[[1]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                 status.in = data.valid.plot$status_C,
                                                 risk.pred.in = predrisk.frank[["risk.joint.est"]],
                                                 n.knots.in = 0,
                                                 t.eval = t.eval,
                                                 model.name = "c-frank")

dat.validplot.frank.list[[2]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                 status.in = data.valid.plot$status_C,
                                                 risk.pred.in = predrisk.frank[["risk.joint.est"]],
                                                 n.knots.in = 3,
                                                 t.eval = t.eval,
                                                 model.name = "c-frank")

dat.validplot.frank.list[[3]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                 status.in = data.valid.plot$status_C,
                                                 risk.pred.in = predrisk.frank[["risk.joint.est"]],
                                                 n.knots.in = 5,
                                                 t.eval = t.eval,
                                                 model.name = "c-frank")

dat.validplot.frank.list[[4]] <- create.pred.obs.by.decile(time.in = data.valid.plot$time_C,
                                                             status.in = data.valid.plot$status_C,
                                                             risk.pred.in = predrisk.frank[["risk.joint.est"]], 
                                                             t.eval = t.eval, model.name = "c-frank")

###
### Frailty normal
###

### Load the required workspace
load(paste("Project 4/data/clinical_example_fit_frailty_normal_model", model.type, "_n", data.size, "_risksonly.RData", sep = ""))
## Load functions
source("Project 4/code/p4_clinical_example_functions.R")

### Calculate the calibration plot data
dat.validplot.frailty.normal.list <- vector("list", 4)

dat.validplot.frailty.normal.list[[1]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                 status.in = data.valid.plot$status_C,
                                                 risk.pred.in = predrisk.frailty.normal[["risk.joint.est"]],
                                                 n.knots.in = 0,
                                                 t.eval = t.eval,
                                                 model.name = "f-norm")

dat.validplot.frailty.normal.list[[2]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                 status.in = data.valid.plot$status_C,
                                                 risk.pred.in = predrisk.frailty.normal[["risk.joint.est"]],
                                                 n.knots.in = 3,
                                                 t.eval = t.eval,
                                                 model.name = "f-norm")

dat.validplot.frailty.normal.list[[3]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                 status.in = data.valid.plot$status_C,
                                                 risk.pred.in = predrisk.frailty.normal[["risk.joint.est"]],
                                                 n.knots.in = 5,
                                                 t.eval = t.eval,
                                                 model.name = "f-norm")

dat.validplot.frailty.normal.list[[4]] <- create.pred.obs.by.decile(time.in = data.valid.plot$time_C,
                                                             status.in = data.valid.plot$status_C,
                                                             risk.pred.in = predrisk.frailty.normal[["risk.joint.est"]], 
                                                             t.eval = t.eval, model.name = "f-norm")

###
### Frailty gamma
###

### Load the required workspace
load(paste("Project 4/data/clinical_example_fit_frailty_gamma_model", model.type, "_n", data.size, "_risksonly.RData", sep = ""))
## Load functions
source("Project 4/code/p4_clinical_example_functions.R")

### Calculate the calibration plot data
dat.validplot.frailty.gamma.list <- vector("list", 4)

dat.validplot.frailty.gamma.list[[1]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                          status.in = data.valid.plot$status_C,
                                                          risk.pred.in = predrisk.frailty.gamma[["risk.joint.est"]],
                                                          n.knots.in = 0,
                                                          t.eval = t.eval,
                                                          model.name = "f-gam")

dat.validplot.frailty.gamma.list[[2]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                          status.in = data.valid.plot$status_C,
                                                          risk.pred.in = predrisk.frailty.gamma[["risk.joint.est"]],
                                                          n.knots.in = 3,
                                                          t.eval = t.eval,
                                                          model.name = "f-gam")

dat.validplot.frailty.gamma.list[[3]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                                          status.in = data.valid.plot$status_C,
                                                          risk.pred.in = predrisk.frailty.gamma[["risk.joint.est"]],
                                                          n.knots.in = 5,
                                                          t.eval = t.eval,
                                                          model.name = "f-gam")

dat.validplot.frailty.gamma.list[[4]] <- create.pred.obs.by.decile(time.in = data.valid.plot$time_C,
                                                                    status.in = data.valid.plot$status_C,
                                                                    risk.pred.in = predrisk.frailty.gamma[["risk.joint.est"]], 
                                                                    t.eval = t.eval, model.name = "f-gam")

###
### MSM
###

### Load the required workspace
load(paste("Project 4/data/clinical_example_msm_fit_model", model.type, "_n", data.size, "_combined.RData", sep = ""))
### Load functions
source("Project 4/code/p4_clinical_example_functions.R")

### Calculate the calibration plot data
dat.validplot.msm.list <- vector("list", 4)

dat.validplot.msm.list[[1]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                               status.in = data.valid.plot$status_C,
                                               risk.pred.in = predrisk.msm,
                                               n.knots.in = 0,
                                               t.eval = t.eval,
                                               model.name = "msm")

dat.validplot.msm.list[[2]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                               status.in = data.valid.plot$status_C,
                                               risk.pred.in = predrisk.msm,
                                               n.knots.in = 3,
                                               t.eval = t.eval,
                                               model.name = "msm")

dat.validplot.msm.list[[3]] <- create.pred.obs(time.in = data.valid.plot$time_C,
                                               status.in = data.valid.plot$status_C,
                                               risk.pred.in = predrisk.msm,
                                               n.knots.in = 5,
                                               t.eval = t.eval,
                                               model.name = "msm")

dat.validplot.msm.list[[4]] <- create.pred.obs.by.decile(time.in = data.valid.plot$time_C,
                                                         status.in = data.valid.plot$status_C,
                                                         risk.pred.in = predrisk.msm, 
                                                         t.eval = t.eval, model.name = "msm")

#########################################################################################################################################
### Step 3: Combine into a single dataset, and reduce to 1000 obsverations per model (pointless plotting 100000 obs per line in plot) ###
#########################################################################################################################################
data.calib.knot0 <- rbind(dat.validplot.product.list[[1]],
                          dat.validplot.dual.list[[1]],
                          dat.validplot.clayton.list[[1]],
                          dat.validplot.gumbel.list[[1]],
                          dat.validplot.frank.list[[1]],
                          dat.validplot.frailty.gamma.list[[1]],
                          dat.validplot.frailty.normal.list[[1]],
                          dat.validplot.msm.list[[1]])
data.calib.knot0$model <- factor(data.calib.knot0$model, levels = c("c-clay", "c-frank","c-gumb", "dual-o", "f-gam", "f-norm", "msm", "product"))
data.calib.knot0 <- data.calib.knot0[seq(1, nrow(data.calib.knot0), nrow(data.calib.knot0)/8000), ]

data.calib.knot3 <- rbind(dat.validplot.product.list[[2]],
                          dat.validplot.dual.list[[2]],
                          dat.validplot.clayton.list[[2]],
                          dat.validplot.gumbel.list[[2]],
                          dat.validplot.frank.list[[2]],
                          dat.validplot.frailty.gamma.list[[2]],
                          dat.validplot.frailty.normal.list[[2]],
                          dat.validplot.msm.list[[2]])
data.calib.knot3$model <- factor(data.calib.knot3$model, levels = c("c-clay", "c-frank","c-gumb", "dual-o", "f-gam", "f-norm", "msm", "product"))
data.calib.knot3 <- data.calib.knot3[seq(1, nrow(data.calib.knot3), nrow(data.calib.knot3)/8000), ]

data.calib.knot5 <- rbind(dat.validplot.product.list[[3]],
                          dat.validplot.dual.list[[3]],
                          dat.validplot.clayton.list[[3]],
                          dat.validplot.gumbel.list[[3]],
                          dat.validplot.frank.list[[3]],
                          dat.validplot.frailty.gamma.list[[3]],
                          dat.validplot.frailty.normal.list[[3]] ,
                          dat.validplot.msm.list[[3]])
data.calib.knot5$model <- factor(data.calib.knot5$model, levels = c("c-clay", "c-frank","c-gumb", "dual-o", "f-gam", "f-norm", "msm", "product"))
data.calib.knot5 <- data.calib.knot5[seq(1, nrow(data.calib.knot5), nrow(data.calib.knot5)/8000), ]

data.calib.decile <- rbind(dat.validplot.product.list[[4]],
                           dat.validplot.dual.list[[4]],
                           dat.validplot.clayton.list[[4]],
                           dat.validplot.gumbel.list[[4]],
                           dat.validplot.frank.list[[4]],
                           dat.validplot.frailty.gamma.list[[4]],
                           dat.validplot.frailty.normal.list[[4]],
                           dat.validplot.msm.list[[4]])

########################
### Plot with ggplot ###
########################


### Create the plots
plot.calib.knot0 <- ggplot(aes(x = risk.pred, y = risk.obs, color = model, lty = model), data = data.calib.knot0) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  geom_rug(col = rgb(.5, 0, 0, alpha = .05)) + xlim(xlim) + ylim(ylim) + 
  xlab("Predicted risk") + ylab("Observed risk") +
  guides(color = guide_legend(title = "Analysis\nMethod"),
         linetype = guide_legend(title = "Analysis\nMethod")) 

plot.calib.knot3 <- ggplot(aes(x = risk.pred, y = risk.obs, color = model, lty = model), data = data.calib.knot3) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  geom_rug(col = rgb(.5, 0, 0, alpha = .05)) + xlim(xlim) + ylim(ylim) + 
  xlab("Predicted risk") + ylab("Observed risk") +
  guides(color = guide_legend(title = "Analysis\nMethod"),
         linetype = guide_legend(title = "Analysis\nMethod"))

plot.calib.knot5 <- ggplot(aes(x = risk.pred, y = risk.obs, color = model, lty = model), data = data.calib.knot5) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  geom_rug(col = rgb(.5, 0, 0, alpha = .05)) + xlim(xlim) + ylim(ylim) + 
  xlab("Predicted risk") + ylab("Observed risk") +
  guides(color = guide_legend(title = "Analysis\nMethod"),
         linetype = guide_legend(title = "Analysis\nMethod"))


plot.calib.knot0.zoom <- ggplot(aes(x = risk.pred, y = risk.obs, color = model, lty = model), data = data.calib.knot0) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  geom_rug(col = rgb(.5, 0, 0, alpha = .05)) + xlim(xlim.zoom) + ylim(ylim.zoom) + 
  xlab("Predicted risk") + ylab("Observed risk") +
  guides(color = guide_legend(title = "Analysis\nMethod"),
         linetype = guide_legend(title = "Analysis\nMethod"))

plot.calib.knot3.zoom <- ggplot(aes(x = risk.pred, y = risk.obs, color = model, lty = model), data = data.calib.knot3) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  geom_rug(col = rgb(.5, 0, 0, alpha = .05)) + xlim(xlim.zoom) + ylim(ylim.zoom) + 
  xlab("Predicted risk") + ylab("Observed risk") +
  guides(color = guide_legend(title = "Analysis\nMethod"),
         linetype = guide_legend(title = "Analysis\nMethod"))

plot.calib.knot5.zoom <- ggplot(aes(x = risk.pred, y = risk.obs, color = model, lty = model), data = data.calib.knot5) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  geom_rug(col = rgb(.5, 0, 0, alpha = .05)) + xlim(xlim.zoom) + ylim(ylim.zoom) + 
  xlab("Predicted risk") + ylab("Observed risk") +
  guides(color = guide_legend(title = "Analysis\nMethod"),
         linetype = guide_legend(title = "Analysis\nMethod"))

plot.calib.decile <- ggplot(aes(x = pred.risk.decile, y = obs.risk.decile), data = data.calib.decile) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = "dashed") +
  facet_wrap(~ model, nrow = 2) + theme(aspect.ratio = 1/1) + xlim(c(0, 0.08)) + ylim(c(0, 0.08)) +
  xlab("Average predicted risk (by decile)") + ylab("Observed risk (Kaplan Meier, by decile)")

# plot.calib.knot0
# plot.calib.knot3
# plot.calib.knot5
# plot.calib.decile

### Save the plots in TIFF format
CairoTIFF(filename = paste("Project 4/figures/clin.example.calib.knot0.model", model.type, ".n", data.size, ".tiff", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot0
dev.off()

CairoTIFF(filename = paste("Project 4/figures/clin.example.calib.knot3.model", model.type, ".n", data.size, ".tiff", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot3
dev.off()

CairoTIFF(filename = paste("Project 4/figures/clin.example.calib.knot5.model", model.type, ".n", data.size, ".tiff", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot5
dev.off()

CairoTIFF(filename = paste("Project 4/figures/clin.example.calib.knot0.zoom.model", model.type, ".n", data.size, ".tiff", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot0.zoom
dev.off()

CairoTIFF(filename = paste("Project 4/figures/clin.example.calib.knot3.zoom.model", model.type, ".n", data.size, ".tiff", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot3.zoom
dev.off()

CairoTIFF(filename = paste("Project 4/figures/clin.example.calib.knot5.zoom.model", model.type, ".n", data.size, ".tiff", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot5.zoom
dev.off()

CairoTIFF(filename = paste("Project 4/figures/clin.example.calib.decile.model", model.type, ".n", data.size, ".tiff", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.decile
dev.off()


### Save the plots in PNG format
CairoPNG(filename = paste("Project 4/figures/clin.example.calib.knot0.model", model.type, ".n", data.size, ".png", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot0
dev.off()

CairoPNG(filename = paste("Project 4/figures/clin.example.calib.knot3.model", model.type, ".n", data.size, ".png", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot3
dev.off()

CairoPNG(filename = paste("Project 4/figures/clin.example.calib.knot5.model", model.type, ".n", data.size, ".png", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot5
dev.off()

CairoPNG(filename = paste("Project 4/figures/clin.example.calib.knot0.zoom.model", model.type, ".n", data.size, ".png", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot0.zoom
dev.off()

CairoPNG(filename = paste("Project 4/figures/clin.example.calib.knot3.zoom.model", model.type, ".n", data.size, ".png", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot3.zoom
dev.off()

CairoPNG(filename = paste("Project 4/figures/clin.example.calib.knot5.zoom.model", model.type, ".n", data.size, ".png", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.knot5.zoom
dev.off()

CairoPNG(filename = paste("Project 4/figures/clin.example.calib.decile.model", model.type, ".n", data.size, ".png", sep = ""), 
          width = 6, height = 5, units = "in", dpi = 800)
plot.calib.decile
dev.off()



### Save the plots in PDF format
CairoPDF(file = paste("Project 4/figures/clin.example.calib.knot0.model", model.type, ".n", data.size, ".pdf", sep = ""), 
          width = 6, height = 5)
plot.calib.knot0
dev.off()

CairoPDF(file = paste("Project 4/figures/clin.example.calib.knot3.model", model.type, ".n", data.size, ".pdf", sep = ""), 
          width = 6, height = 5)
plot.calib.knot3
dev.off()

CairoPDF(file = paste("Project 4/figures/clin.example.calib.knot5.model", model.type, ".n", data.size, ".pdf", sep = ""), 
          width = 6, height = 5)
plot.calib.knot5
dev.off()

CairoPDF(file = paste("Project 4/figures/clin.example.calib.knot0.zoom.model", model.type, ".n", data.size, ".pdf", sep = ""), 
          width = 6, height = 5)
plot.calib.knot0.zoom
dev.off()

CairoPDF(file = paste("Project 4/figures/clin.example.calib.knot3.zoom.model", model.type, ".n", data.size, ".pdf", sep = ""), 
          width = 6, height = 5)
plot.calib.knot3.zoom
dev.off()

CairoPDF(file = paste("Project 4/figures/clin.example.calib.knot5.zoom.model", model.type, ".n", data.size, ".pdf", sep = ""), 
          width = 6, height = 5)
plot.calib.knot5.zoom
dev.off()

CairoPDF(file = paste("Project 4/figures/clin.example.calib.decile.model", model.type, ".n", data.size, ".pdf", sep = ""), 
          width = 6, height = 5)
plot.calib.decile
dev.off()






############################################
############################################
### Calculate and compare discrimination ###
############################################
############################################

### Create a dataset with all the predicted risks in for each method
dat.discrim <- data.frame(data.valid.plot, 
                          "pred.product" = arrange(dat.validplot.product.list[[1]], id)$risk.pred,
                          "pred.dual" = arrange(dat.validplot.dual.list[[1]], id)$risk.pred,
                          "pred.clayton" = arrange(dat.validplot.clayton.list[[1]], id)$risk.pred,
                          "pred.gumbel" = arrange(dat.validplot.gumbel.list[[1]], id)$risk.pred,
                          "pred.frank" = arrange(dat.validplot.frank.list[[1]], id)$risk.pred,
                          "pred.normal" = arrange(dat.validplot.frailty.normal.list[[1]], id)$risk.pred,
                          "pred.gamma" = arrange(dat.validplot.frailty.gamma.list[[1]], id)$risk.pred,
                          "pred.msm" = arrange(dat.validplot.msm.list[[1]], id)$risk.pred)

###
### Harrell's C
###

HarC.table <- c("product" = rcorr.cens(1-dat.discrim[ , "pred.product"], 
                                       Surv(dat.discrim$time_C, dat.discrim$status_C))["C Index"],
                "dual" = rcorr.cens(1-dat.discrim[ , "pred.dual"], 
                                       Surv(dat.discrim$time_C, dat.discrim$status_C))["C Index"],
                "clayton" = rcorr.cens(1-dat.discrim[ , "pred.clayton"], 
                                       Surv(dat.discrim$time_C, dat.discrim$status_C))["C Index"],
                "gumbel" = rcorr.cens(1-dat.discrim[ , "pred.gumbel"], 
                                       Surv(dat.discrim$time_C, dat.discrim$status_C))["C Index"],
                "frank" = rcorr.cens(1-dat.discrim[ , "pred.frank"], 
                                       Surv(dat.discrim$time_C, dat.discrim$status_C))["C Index"],
                "normal" = rcorr.cens(1-dat.discrim[ , "pred.normal"], 
                                       Surv(dat.discrim$time_C, dat.discrim$status_C))["C Index"],
                "gamma" = rcorr.cens(1-dat.discrim[ , "pred.gamma"], 
                                       Surv(dat.discrim$time_C, dat.discrim$status_C))["C Index"],
                "msm" = rcorr.cens(1-dat.discrim[ , "pred.msm"], 
                                       Surv(dat.discrim$time_C, dat.discrim$status_C))["C Index"]
                )

str(rcorr.cens(1-dat.discrim[ , "pred.product"], 
               Surv(dat.discrim$time_C, dat.discrim$status_C)))
print("HarC")
HarC.table

###
### Uno's C
###
data.devel.wide <- tidyr::pivot_wider(data.devel, id_cols = c(id, all_of(variables.vec)), 
                                      names_from = outcome_char, values_from = c(time, status)) 

### UnoC struggles to calculate (spits out NaN) if too many observations, so have to calculate on a subset of individuals
UnoC.nrow <- 15000

UnoC.table <- c("product" = UnoC(Surv.rsp = Surv(data.devel.wide$time_C[1:UnoC.nrow], data.devel.wide$status_C[1:UnoC.nrow]), 
                                 Surv.rsp.new = Surv(dat.discrim$time_C[1:UnoC.nrow], dat.discrim$status_C[1:UnoC.nrow]), 
                                 lpnew = 1-dat.discrim[ , "pred.product"][1:UnoC.nrow]),
                "dual" = UnoC(Surv.rsp = Surv(data.devel.wide$time_C[1:UnoC.nrow], data.devel.wide$status_C[1:UnoC.nrow]), 
                              Surv.rsp.new = Surv(dat.discrim$time_C[1:UnoC.nrow], dat.discrim$status_C[1:UnoC.nrow]), 
                              lpnew = 1-dat.discrim[ , "pred.dual"][1:UnoC.nrow]),
                "clayton" = UnoC(Surv.rsp = Surv(data.devel.wide$time_C[1:UnoC.nrow], data.devel.wide$status_C[1:UnoC.nrow]), 
                                 Surv.rsp.new = Surv(dat.discrim$time_C[1:UnoC.nrow], dat.discrim$status_C[1:UnoC.nrow]), 
                                 lpnew = 1-dat.discrim[ , "pred.clayton"][1:UnoC.nrow]),
                "gumbel" = UnoC(Surv.rsp = Surv(data.devel.wide$time_C[1:UnoC.nrow], data.devel.wide$status_C[1:UnoC.nrow]), 
                                Surv.rsp.new = Surv(dat.discrim$time_C[1:UnoC.nrow], dat.discrim$status_C[1:UnoC.nrow]), 
                                lpnew = 1-dat.discrim[ , "pred.gumbel"][1:UnoC.nrow]),
                "frank" = UnoC(Surv.rsp = Surv(data.devel.wide$time_C[1:UnoC.nrow], data.devel.wide$status_C[1:UnoC.nrow]), 
                               Surv.rsp.new = Surv(dat.discrim$time_C[1:UnoC.nrow], dat.discrim$status_C[1:UnoC.nrow]), 
                               lpnew = 1-dat.discrim[ , "pred.frank"][1:UnoC.nrow]),
                "normal" = UnoC(Surv.rsp = Surv(data.devel.wide$time_C[1:UnoC.nrow], data.devel.wide$status_C[1:UnoC.nrow]), 
                                Surv.rsp.new = Surv(dat.discrim$time_C[1:UnoC.nrow], dat.discrim$status_C[1:UnoC.nrow]), 
                                lpnew = 1-dat.discrim[ , "pred.normal"][1:UnoC.nrow]),
                "gamma" = UnoC(Surv.rsp = Surv(data.devel.wide$time_C[1:UnoC.nrow], data.devel.wide$status_C[1:UnoC.nrow]), 
                               Surv.rsp.new = Surv(dat.discrim$time_C[1:UnoC.nrow], dat.discrim$status_C[1:UnoC.nrow]), 
                               lpnew = 1-dat.discrim[ , "pred.gamma"][1:UnoC.nrow]),
                "msm" = UnoC(Surv.rsp = Surv(data.devel.wide$time_C[1:UnoC.nrow], data.devel.wide$status_C[1:UnoC.nrow]), 
                             Surv.rsp.new = Surv(dat.discrim$time_C[1:UnoC.nrow], dat.discrim$status_C[1:UnoC.nrow]), 
                             lpnew = 1-dat.discrim[ , "pred.msm"][1:UnoC.nrow])
)


print("UnoC")
UnoC.table


rm(list=setdiff(ls(), list("model.type", "data.size", "UnoC.table", "HarC.table")))
save.image(paste("Project 4/data/compare_performance_models", model.type, "_n", data.size, ".RData", sep = ""))
