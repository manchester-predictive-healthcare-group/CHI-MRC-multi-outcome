###
### This program will...
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_5")
getwd()

### Load packages
source("code/z_load_packages.R")

### Set the seed
set.seed(505)

### Load prepped data
load("data/prep_ebmt.RData")

### Assign id.iter from command line
args <- commandArgs(trailingOnly = T)
model <- as.numeric(args[1])
s <- as.numeric(args[2])
t.eval <- as.numeric(args[3])
id.iter <- as.numeric(args[4])
print(paste("s = ", s, Sys.time()))
print(paste("t.eval = ", t.eval, Sys.time()))
print(paste("id iter = ", id.iter, Sys.time()))

### Assign variables for model we will be fitting
if (model == 1){
  eq.RHS <- paste(do.call(paste0, expand.grid(c("match", "proph", "year1", "year2", "agecl1", "agecl2"), paste(".", 1:12, sep = ""))), collapse="+")
  eq <- paste("Surv(Tstart, Tstop, status) ~ ", eq.RHS,  "+ strata(trans)", sep = "")
  eq <- as.formula(eq)
} else if (model == 2){
  eq.RHS <- paste(do.call(paste0, expand.grid(c("agecl1", "agecl2"), paste(".", 1:12, sep = ""))), collapse="+")
  eq <- paste("Surv(Tstart, Tstop, status) ~ ", eq.RHS,  "+ strata(trans)", sep = "")
  eq <- as.formula(eq)
}

### Develop a model on entire dataset except individual of interest
cfull <- coxph(eq, data = subset(msebmt, id != id.iter), method = "breslow")

### Get location of individual in msebmt
pat.loc <- which(msebmt$id == id.iter)

### Create a miniture dataset, on which to generate predictions in (must be in mstate format and have a row for every transition)
pat.dat <- msebmt[rep(pat.loc[1], 12), 9:12]
pat.dat$trans <- 1:12
attr(pat.dat, "trans") <- tmat
pat.dat <- expand.covs(pat.dat, covs, longnames = FALSE)
pat.dat$strata <- pat.dat$trans

### Fit cause-specific hazards
msf.pat <- msfit(cfull, pat.dat, trans = tmat)

### Generate 5 year transition probabilities at time s
pt <- probtrans(msf.pat, predt = s)

### Write a function to extract the transition probabilities from state j into each state, after followup time f.time
extract.tp <- function(tp.object, j, t.eval){
  ### Create output object
  output.object <- return(subset(tp.object[[j]], time > t.eval) %>% slice(1) %>% select(-c(time)))
}

### Will generate risks out of every state j and store in tp.id
tp.id <- vector("list", 6)
for (j in 1:6){
  tp.id[[j]] <- data.frame("id" = id.iter, extract.tp(tp.object = pt, j, t.eval))
}

### Clean workspace
rm(list = setdiff(ls(), list("tp.id", "id.iter", "t.eval", "s", "model")))

### Save image
save.image(paste("data/calc_tp_leaveoo_model", model, "_s", s, "_t", t.eval, "_id", id.iter, ".RData", sep = ""))
print("IMAGE SAVED")