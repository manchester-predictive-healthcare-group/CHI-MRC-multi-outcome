##################################################
### Code to prepare datasets provided with package
##################################################

### Load the ebmt dataset from mstate package and save it do /data
###
data("ebmt4")
ebmt <- ebmt4

### Define tmat
tmat <- transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
                          c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
tmat

### Create mstate format
msebmt <- msprep(data = ebmt, trans = tmat, time = c(NA, "rec", "ae","recae", "rel", "srv"),
                 status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"),
                 keep = c("match", "proph", "year", "agecl"))
msebmt[msebmt$id == 1, c(1:8, 10:12)]

### Define covariates for model
covs <- c("match", "proph", "year", "agecl")
msebmt <- expand.covs(msebmt, covs, longnames = FALSE)
msebmt[msebmt$id == 1, -c(9, 10, 12:48, 61:84)]

### Change time scale of model into years
msebmt[, c("Tstart", "Tstop", "time")] <- msebmt[, c("Tstart","Tstop", "time")]

### Create a variable which is maximum observed follow up time for all individuals, this is when they were either censored, relapsed or died
ebmt$dtcens <- pmin(ebmt$rel, ebmt$srv)
ebmt$dtcens.s <- 1 - pmax(ebmt$rel.s, ebmt$srv.s)

### Assign variables for model we will be fitting
eq.RHS <- paste(do.call(paste0, expand.grid(c("match", "proph", "year1", "year2", "agecl1", "agecl2"), paste(".", 1:12, sep = ""))), collapse="+")
eq <- paste("Surv(Tstart, Tstop, status) ~ ", eq.RHS,  "+ strata(trans)", sep = "")
eq <- as.formula(eq)

###################
### Create tps0 ###
###################

### Assign parameters
s <- 0
t.eval <- 1826

### Create dataframe to store predicted risks
tps0test <- vector("list", 6)
for (j in 1:6){
  tps0test[[j]] <-  data.frame(matrix(NA, ncol = 14, nrow = nrow(ebmt)))
  colnames(tps0test[[j]]) <- c("id", paste("pstate", 1:6, sep = ""), paste("se", 1:6, sep = ""), "j")
}

### Run through each patient id
for (id.iter in 1:nrow(ebmt)){

  print(paste("id.iter = ", id.iter, Sys.time()))

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
  extract.tp <- function(tp.object, state, f.time){
    ### Create output object
    output.object <- as.numeric(subset(tp.object[[state]], time > f.time) %>% slice(1) %>% select(-c(time)))
    return(output.object)
  }

  ### Calculate required transition probabilities and store in output dataset
  ### Will generate risks out of every state j and store in tp.id
  for (j in 1:6){
    tps0test[[j]][id.iter, ] <- c(id.iter, extract.tp(tp.object = pt, state = j, f.time = t.eval - s), j)
  }

}

### Combine into one dataset
tps0test <- do.call("rbind", tps0test)

### Use the p9 programs to get in the correct form....

### Use this data
usethis::use_data(ebmt, overwrite = TRUE)
usethis::use_data(msebmcal, overwrite = TRUE)
usethis::use_data(ebmcal, overwrite = TRUE)
usethis::use_data(tps0, overwrite = TRUE)
usethis::use_data(tps100, overwrite = TRUE)
