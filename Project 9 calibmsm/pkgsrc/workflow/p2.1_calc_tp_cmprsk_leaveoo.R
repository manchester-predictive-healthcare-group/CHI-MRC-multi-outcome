###
### This program will...
###

### Set working directory
rm(list=ls())
setwd("/mnt/bmh01-rds/mrc-multi-outcome/Project_5")
getwd()

### Load packages
source("code/z_load_packages.R")

### Assign j, s and t.eval from comand line
args <- commandArgs(trailingOnly = T)
j <- as.numeric(args[1])
s <- as.numeric(args[2])
t.eval <- as.numeric(args[3])
print(paste("j = ", j, Sys.time()))
print(paste("s = ", s, Sys.time()))
print(paste("t.eval = ", t.eval, Sys.time()))

### Load the prepped data
load("data/prep_ebmt.RData")

### Define tmat
if (j == 1){
  tmat <- transMat(x = list(c(2, 3, 5, 6), c(), c(), c(),
                            c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
} else if (j == 2){
  tmat <- transMat(x = list(c(), c(4, 5, 6), c(), c(),
                            c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
} else if (j == 3){
  tmat <- transMat(x = list(c(), c(), c(4, 5, 6), c(),
                            c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
} else if (j == 4){
  tmat <- transMat(x = list(c(), c(), c(), c(5, 6),
                            c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
}

### Also define the original transition matrix, so we can extract the appropriate predictors from msebmt
tmat.original <- transMat(x = list(c(2, 3, 5, 6), c(4,5,6), c(4,5,6), c(5, 6),
                                   c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
tmat.original.j <- tmat.original[j, ][!is.na(tmat.original[j, ])]

### Want to validate the competing risks model out of state j t time s, so remove all observations not in state j at time s
msebmt.j <- msebmt %>% subset(from == j & Tstart <= s & s < Tstop) 
ebmt.j <- ebmt %>% subset(id %in% msebmt.j$id)

### For mstate dataset, set Tstart to zero, and deduct s from Tstop and time variables
msebmt.j <- mutate(msebmt.j, Tstart = 0,  
                   time = Tstop - s,
                   Tstop = Tstop - s,
                   trans = rep(1:sum(!is.na(tmat)), nrow(ebmt.j)))


### Assign the predictor variables for transitions from tmat.original.j, to be for the first n transitions
msebmt.j[,paste(do.call(paste0, expand.grid(c("match", "proph", "year1", "year2", "agecl1", "agecl2"), paste(".", 1:sum(!is.na(tmat)), sep = ""))))] <-
  msebmt.j[,paste(do.call(paste0, expand.grid(c("match", "proph", "year1", "year2", "agecl1", "agecl2"), paste(".", tmat.original.j, sep = ""))))]

### Assign variables for model we will be fitting
eq.RHS <- paste(do.call(paste0, expand.grid(c("match", "proph", "year1", "year2", "agecl1", "agecl2"), paste(".", 1:sum(!is.na(tmat)), sep = ""))), collapse="+")
eq <- paste("Surv(Tstart, Tstop, status) ~ ", eq.RHS,  "+ strata(trans)", sep = "")
eq <- as.formula(eq)

### Loop through id.iter
for (id.iter in ebmt.j$id){
  
  print(paste("id.iter = ", id.iter, Sys.time()))
  
  ### Develop a model on entire dataset except individual of interest
  cfull <- coxph(eq, data = subset(msebmt.j, id != id.iter), method = "breslow")
  
  ### Get location of individual in msebmt.j
  pat.loc <- which(msebmt.j$id == id.iter)
  
  ### Create a miniture dataset, on which to generate predictions in (must be in mstate format and have a row for every transition)
  pat.dat <- msebmt.j[rep(pat.loc[1], sum(!is.na(tmat))), 9:12]
  pat.dat$trans <- 1:sum(!is.na(tmat))
  attr(pat.dat, "trans") <- tmat
  pat.dat <- expand.covs(pat.dat, covs, longnames = FALSE)
  pat.dat$strata <- pat.dat$trans
  
  ### Fit cause-specific hazards
  msf.pat <- msfit(cfull, pat.dat, trans = tmat)
  
  ### Generate 5 year transition probabilities for this patient from times s = 0
  pt <- probtrans(msf.pat, predt = 0)
  
  ### Write a function to extract the transition probabilities from state j into each state, after followup time f.time
  extract.tp <- function(tp.object, state.j, f.time){
    ### Create output object
    output.object <- return(subset(tp.object[[state.j]], time > f.time) %>% slice(1) %>% select(-c(time)))
  }
  
  ### Calculate required transition probabilities and store in output dataset
  tp.id <- data.frame("id" = id.iter, extract.tp(tp.object = pt, state.j = j, f.time = t.eval - s))
  
  ### Remove everything except id.iter and tp.id
  rm(list = setdiff(ls(), list("tp.id", "id.iter", "t.eval", "s", "j", "msebmt.j", "ebmt.j", "tmat", "eq", "covs")))
  
  ### Save image
  save.image(paste("data/calc_tp_cmprsk_leaveoo_j", j, "_s", s, "_t.eval", t.eval, "_iditer", id.iter, ".RData", sep = ""))
}

