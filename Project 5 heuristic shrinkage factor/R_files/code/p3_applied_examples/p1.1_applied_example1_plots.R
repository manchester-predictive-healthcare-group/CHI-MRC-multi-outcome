### setwd
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_8.2")

### Load libraries
library(rms)
library(boot)
library(ggplot2)

### Set seed
set.seed(222)

### Rename
df <- readRDS("data/ce_df_mimic.rds")

###
colnames(df)

###
str(df$AKI)

### Build a model to predict AKI

### Set n.iter
n.iter <- 500

### set n.vec
n.vec <- c(100, 250, 500, 1000, 2000, 5000, 10000)

### Set number of iterations for bootstrapped
n.S.boot <- 200

### Set predictors to adjust for in models 
predictors <- c("age", "gender", "bicarbonate_mean", "creatinine_mean", "chloride_mean", "hemoglobin_mean", "platelet_mean", "potassium_mean", "ptt_mean", "inr_mean", 
                "pt_mean" , "bun_mean", "wbc_mean", "heartrate_mean", "sysbp_mean", "diasbp_mean", "resprate_mean", "tempc_mean", "spo2_mean", "glucose_mean")

### Create output
data.output <- vector("list", length(n.vec))
names(data.output) <- n.vec
for (i in 1:length(n.vec)){
  data.output[[i]] <- data.frame("S.VH" = rep(NA, n.iter), "S.boot" = rep(NA, n.iter), "R2.CS.app" = rep(NA, n.iter), "C.app" = rep(NA, n.iter))
}

###
### Define a function to calculate S using bootstrapping
###
function_calc_S_boot <- function(data, i){
  
  ### Take the bootstrap sample
  boot.dat <- data[i, ]
  
  ### Create model in bootstrapped dataset
  boot.model <- glm(as.formula(paste("AKI ~ ", paste(predictors, collapse = "+"), sep = "")), data = boot.dat, family = binomial(link = "logit")) 
  
  ### Generate predictions using this model using the new dataset
  boot.lp <- predict(boot.model, newdata = data.temp)
  
  ### Create a temporary dataset with both these things
  boot.dat.devel.temp <- data.frame(data.temp, "boot.lp" = boot.lp)
  
  ### Calculate calibration slope
  boot.calib.model <- glm(AKI ~ boot.lp, data = boot.dat.devel.temp, family = binomial(link = "logit")) 
  
  ### Save slope
  return(boot.calib.model$coefficients["boot.lp"])
  
}

### Let n.samp <- 250
n.samp <- 250

### Sample df at random
data.temp <- df[sample(1:nrow(df), n.samp, replace = FALSE), ]

### Build model
model.temp <- lrm(as.formula(paste("AKI ~ ", paste(predictors, collapse = "+"), sep = "")), data = data.temp)

### Record required statistics for prediction of S.pop
LR <- model.temp$deviance[1] - model.temp$deviance[2]
C.app <- model.temp$stats["C"]
R2.CS.app <- as.numeric(1 - exp(-(LR)/nrow(data.temp)))
P <- length(model.temp$coefficients)
S.VH <- as.numeric(1 - P/LR)

###
### Calculate S.boot

### Run the bootstrapping
boot.out <- boot(data.temp, function_calc_S_boot, R = n.S.boot)

### Assign S.boot
S.boot <- mean(boot.out$t)

### Compare
S.VH
S.boot

### Get lp
data.temp$lp <- predict(model.temp, response = "lp")

### Apply shrinkage
data.temp$lp.shrunk.VH <- data.temp$lp*S.VH
data.temp$lp.shrunk.boot <- data.temp$lp*S.boot

### Need to re-estimate intercept to ensure appropriate mean predicted risk
model.adjusted.VH <- glm(AKI ~ offset(lp.shrunk.VH), family = binomial(link = "logit"), data = data.temp)
model.adjusted.boot <- glm(AKI ~ offset(lp.shrunk.boot), family = binomial(link = "logit"), data = data.temp)

### Add to the shrunk LP
data.temp$lp.shrunk.VH <- data.temp$lp.shrunk.VH + coef(model.adjusted.VH)
data.temp$lp.shrunk.boot <- data.temp$lp.shrunk.boot + coef(model.adjusted.boot)

### Convert to probabilities
data.temp$p.shrunk.none <- exp(data.temp$lp)/(1+exp(data.temp$lp))
data.temp$p.shrunk.VH <- exp(data.temp$lp.shrunk.VH)/(1+exp(data.temp$lp.shrunk.VH))
data.temp$p.shrunk.boot <- exp(data.temp$lp.shrunk.boot)/(1+exp(data.temp$lp.shrunk.boot))

### Check mean predicted probabilities are the same after shrinkage
mean(data.temp$p.shrunk.none)
mean(data.temp$p.shrunk.VH)
mean(data.temp$p.shrunk.boot)

### Create long data
data.temp.long <- data.temp |>
  dplyr::rename(S_VH = p.shrunk.VH, S_boot = p.shrunk.boot) |> 
  tidyr::pivot_longer(cols = c(S_VH, S_boot)) |>
  dplyr::rename(p = value, shrinkage = name)

### Create plot object
plot.object <- ggplot(data = data.temp.long) + 
  geom_density(aes(x = p, color = shrinkage)) +
  scale_color_manual(labels = c(expression(S[boot]), expression(S[VH])), values = c("red", "blue")) +
  xlab("Predicted risk")

### Save to disk
Cairo::CairoPNG("figures/gg.applied.example1.png", width = 7, height = 7, unit = "in", dpi = 300)
plot(plot.object)
dev.off()
