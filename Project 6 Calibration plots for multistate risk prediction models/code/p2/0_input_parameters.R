### Define key data for simulation 

### Calculate the scale for each p, taken from Putter et al 2006, Estimation and prediction in a multistate model for breast cancer, Figure 2
csh.survival.probs.desired <- c("12" = 0.9, "13" = 0.8, "15" = 0.99, "24" = 0.55, "25" = 0.95, "34" = 0.7, "35" = 0.15, "45" = 0.05)

### Apply function to calculate the desired scales
scales.sim <- sapply(csh.survival.probs.desired, calc.scale)
scales.sim

### Define time at which risks will be gneerated and  calibration will be assessed
t.eval <- round(365.25*7)

### Defne the max follow up to be after t.eval
max.follow <- t.eval + 2

### Define beta's
beta.x12 = 0.5 #covariate effects for transiion 12
beta.x13 = 0.5 #covariate effects for transiion 13
beta.x15 = 0.5 #covariate effects for transiion 15
beta.x24 = 0.5 #covariate effects for transiion 24
beta.x25 = 0.5 #covariate effects for transiion 25
beta.x34 = 0.5 #covariate effects for transiion 34
beta.x35 = 0.5 #covariate effects for transiion 35
beta.x45 = 0.5 #covariate effects for transiion 45

### Set the coefficients for the censoring mechanism
if (scen == "C1"){
  cens_beta_x12 <- 0
  cens_beta_x13 <- 0
  cens_beta_x15 <- 0
  cens_beta_x24 <- 0
  cens_beta_x25 <- 0
  cens_beta_x34 <- 0
  cens_beta_x35 <- 0
  cens_beta_x45 <- 0
} else if (scen == "C2"){
  cens_beta_x12 <- 0.125
  cens_beta_x13 <- 0.125
  cens_beta_x15 <- 0.125
  cens_beta_x24 <- 0.125
  cens_beta_x25 <- 0.125
  cens_beta_x34 <- 0.125
  cens_beta_x35 <- 0.125
  cens_beta_x45 <- 0.125
} else if (scen == "C3"){
  cens_beta_x12 <- 0.25
  cens_beta_x13 <- 0.25
  cens_beta_x15 <- 0.25
  cens_beta_x24 <- 0.25
  cens_beta_x25 <- 0.25
  cens_beta_x34 <- 0.25
  cens_beta_x35 <- 0.25
  cens_beta_x45 <- 0.25
}

### Censoring is applied at a scale resulting in 20% censored by 7 years follow up
cens_scale <- calc.scale(0.6)