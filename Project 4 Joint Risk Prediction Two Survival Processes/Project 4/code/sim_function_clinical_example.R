####################################################################################################################################
####################################################################################################################################
### Function to take cohort and get data into a common data format, compatible with existing code, and meaning outcome variables ###
### can be easily changed
####################################################################################################################################
####################################################################################################################################
transform.data.cop <- function(data.in, variables.vec){
  ### Also need to create new variable called CVD  (time until either have occured)
  data.in$CVD_ev_t <- pmin(data.in$CHD_MI_ev_t, data.in$Stroke_TIA_ev_t, data.in$HF_ev_t)
  data.in$CVD_ev_c <- as.factor(pmax(as.numeric(data.in$CHD_MI_ev_c == 1), as.numeric(data.in$Stroke_TIA_ev_c == 1), as.numeric(data.in$HF_ev_c == 1)))
  
  ### Also need to create new variable called CVD_Diab  (time until both have occured)
  data.in$CVD_Diab_ev_t <- pmax(data.in$CVD_ev_t, data.in$Diab_t2_ev_t)
  data.in$CVD_Diab_ev_c <- as.factor(pmin(as.numeric(data.in$CVD_ev_c == 1), as.numeric(data.in$Diab_t2_ev_c == 1)))
  
  data.in.A <- data.in
  data.in.B <- data.in
  data.in.C <- data.in
  
  data.in.A$time <- data.in.A$CVD_ev_t
  data.in.A$status <- as.numeric(data.in.A$CVD_ev_c) - 1
  data.in.A$outcome <- 1
  data.in.A$outcome_char <- "A"
  data.in.A$id <- 1:nrow(data.in.A)
  
  data.in.B$time <- data.in.B$Diab_t2_ev_t
  data.in.B$status <- as.numeric(data.in.B$Diab_t2_ev_c) - 1
  data.in.B$outcome <- 2
  data.in.B$outcome_char <- "B"
  data.in.B$id <- 1:nrow(data.in.B)
  
  data.in.C$time <- data.in.C$CVD_Diab_ev_t
  data.in.C$status <- as.numeric(data.in.C$CVD_Diab_ev_c) - 1
  data.in.C$outcome <- 3
  data.in.C$outcome_char <- "C"
  data.in.C$id <- 1:nrow(data.in.C)
  
  data.in.A <- dplyr::select(data.in.A, outcome, outcome_char, id, time, status, 
                             all_of(variables.vec))
  
  data.in.B <- dplyr::select(data.in.B, outcome, outcome_char, id, time, status, 
                             all_of(variables.vec))
  
  data.in.C <- dplyr::select(data.in.C, outcome, outcome_char, id, time, status, 
                             all_of(variables.vec))
  
  data.in <- rbind(data.in.A, data.in.B, data.in.C)
  
  return(data.in)
}



###################################################################################################################
###################################################################################################################
### Functions to generate a risk for a given model after model has been fitted. These will be applied over the  ###
### validation dataset by row to gnerate predicted risks
###################################################################################################################
###################################################################################################################

##################################################################################################
### Create a function which takes two elements, and produces a joint risk using the copula     ###
### This will be used in combination with apply to generate risk scores for validation dataset ###
##################################################################################################
myfunc.jointrisk.cop <- function(u.in, copula.in){
  return(pCopula(u = c(u.in[1], u.in[2]), copula = copula.in))
}


##############################################################################################################
### Create a funciton to generate predicted risk for a row of data from validation dataset using msm model ###
### This will be used in combination with apply to generate risk scores for validation dataset             ###
##############################################################################################################
myfunc.jointrisk.msm <- function(row, msm.data.in, msm.model.in, msm.tmat.in, variables.vec){
  
  ## Start by extracting the corret row from msm.data.in, using the first entity of row (which is the id variable)
  pat <- msm.data.in[as.numeric(row[1]), ]
  
  # Turn into an msprep object
  pat <- msprep(pat, trans = msm.tmat.in, time = c(NA, "state2", "state3", "state4"),
                status = c(NA, "state2.s", "state3.s", "state4.s"), keep = variables.vec)
  # Make it four rows (one for each model)
  pat <- pat[rep(1, 4), variables.vec ]
  # Assign trans variable and associate transition matrix
  pat$trans <- 1:4
  attr(pat, "trans") <- msm.tmat.in
  # Expand covariates to allow different effects per transition, after creating covs object
  covs <- variables.vec
  pat <- expand.covs(pat, covs, longnames = FALSE)
  pat$strata <- pat$trans
  
  
  ## Apply the multistate model
  msf <- msfit(msm.model.in, pat, trans = msm.tmat.in)
  
  ## Generate risks for patients from time t = 0
  pt <- probtrans(msf, predt = 0)
  
  ## Calculate joint risks, joint survival, and marginal survival probabilities, at time t.eval
  
  # Joint risk is probability of being in state 4
  risk.joint.est <- pt[[1]]$pstate4[max(which(pt[[1]]$time < t.eval))]
  
  return(risk.joint.est)
}

##################################################################################################################
### Create a funciton to generate predicted risk for a row of data from validation dataset using frailty model ###
### This will be used in combination with apply to generate risk scores for validation dataset                 ###
##################################################################################################################
myfunc.jointrisk.frailty <- function(x.in, betas_A, betas_B, fit.in){
  
  if (fit.in[["frail.dist"]] == "gamma"){
    surv.marg.A.func <- function(frail.term){
      return(pweibull(t.eval, 
                      shape = fit.in[["bh.shape.A"]], 
                      scale = 1/(frail.term*exp(fit.in[["bh.int.A"]] + x.in %*% betas_A)),
                      lower.tail = FALSE))
    }
    surv.marg.B.func <- function(frail.term){
      return(pweibull(t.eval, 
                      shape = fit.in[["bh.shape.B"]], 
                      scale = 1/(frail.term*exp(fit.in[["bh.int.B"]] + x.in %*% betas_B)),
                      lower.tail = FALSE))
    }
    risk.marg.A.func <- function(frail.term){
      return(1-pweibull(t.eval, 
                        shape = fit.in[["bh.shape.A"]], 
                        scale = 1/(frail.term*exp(fit.in[["bh.int.A"]] + x.in %*% betas_A)),
                        lower.tail = FALSE))
    }
    risk.marg.B.func <- function(frail.term){
      return(1-pweibull(t.eval, 
                        shape = fit.in[["bh.shape.B"]], 
                        scale = 1/(frail.term*exp(fit.in[["bh.int.B"]] + x.in %*% betas_B)),
                        lower.tail = FALSE))
    }
  } else if (fit.in[["frail.dist"]] == "normal"){
    surv.marg.A.func <- function(frail.term){
      return(pweibull(t.eval, 
                      shape = fit.in[["bh.shape.A"]], 
                      scale = 1/(exp(frail.term + fit.in[["bh.int.A"]] + x.in %*% betas_A)),
                      lower.tail = FALSE))
    }
    surv.marg.B.func <- function(frail.term){
      return(pweibull(t.eval, 
                      shape = fit.in[["bh.shape.B"]], 
                      scale = 1/(exp(frail.term + fit.in[["bh.int.B"]] + x.in %*% betas_B)),
                      lower.tail = FALSE))
    }
    risk.marg.A.func <- function(frail.term){
      return(1-pweibull(t.eval, 
                        shape = fit.in[["bh.shape.A"]], 
                        scale = 1/(exp(frail.term + fit.in[["bh.int.A"]] + x.in %*% betas_A)),
                        lower.tail = FALSE))
    }
    risk.marg.B.func <- function(frail.term){
      return(1-pweibull(t.eval, 
                        shape = fit.in[["bh.shape.B"]], 
                        scale = 1/(exp(frail.term + fit.in[["bh.int.B"]] + x.in %*% betas_B)),
                        lower.tail = FALSE))
    }
  }
  
  ### Integrate these functions over the frailty distributions to get the marginal risks/surival probabilities,
  ### and Integrate the product of these functions over the frailty distributions to get the joint risks/survival probabilities
  
  ## First need to create functions which are the product of the entity we want to integrate over, and the density function of the 
  ## random effect (with variance estimated from the data)
  if (fit.in[["frail.dist"]] == "gamma"){
    surv.marg.est.A.for.int <- function(z.in){
      surv.marg.A.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    surv.marg.est.B.for.int <- function(z.in){
      surv.marg.B.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    risk.marg.est.A.for.int <- function(z.in){
      risk.marg.A.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    risk.marg.est.B.for.int <- function(z.in){
      risk.marg.B.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    surv.joint.est.for.int <- function(z.in){
      surv.marg.A.func(z.in)*surv.marg.B.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
    risk.joint.est.for.int <- function(z.in){
      risk.marg.A.func(z.in)*risk.marg.B.func(z.in)*dgamma(z.in, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    }
  } else if (fit.in[["frail.dist"]] == "normal"){
    surv.marg.est.A.for.int <- function(z.in){
      surv.marg.A.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    surv.marg.est.B.for.int <- function(z.in){
      surv.marg.B.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    risk.marg.est.A.for.int <- function(z.in){
      risk.marg.A.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    risk.marg.est.B.for.int <- function(z.in){
      risk.marg.B.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    surv.joint.est.for.int <- function(z.in){
      surv.marg.A.func(z.in)*surv.marg.B.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
    risk.joint.est.for.int <- function(z.in){
      risk.marg.A.func(z.in)*risk.marg.B.func(z.in)*dnorm(z.in, mean = 0, sd = fit.in[["frail.var.est"]])
    }
  }
  
  ## Now do the integration
  if (fit.in[["frail.dist"]] == "gamma"){
    ## Define upper limit for integration at the 99.999th percentile of the ditribution
    upper.lim <- qgamma(0.99999, shape = fit.in[["frail.var.est"]], scale = 1/fit.in[["frail.var.est"]])
    
    ## Now run the integration
    #surv.marg.est.A <- cubintegrate(f = surv.marg.est.A.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    #surv.marg.est.B <- cubintegrate(f = surv.marg.est.B.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    #risk.marg.est.A <- cubintegrate(f = risk.marg.est.A.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    #risk.marg.est.B <- cubintegrate(f = risk.marg.est.B.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    #surv.joint.est <- cubintegrate(f = surv.joint.est.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
    risk.joint.est <- cubintegrate(f = risk.joint.est.for.int, lower = 0, upper = upper.lim, method = "hcubature")$integral
  } else if (fit.in[["frail.dist"]] == "normal"){
    ## Define lower and upper limit for integration at the 0.01th 99.99th percentile of the ditribution
    lower.lim <- qnorm(0.00001, mean = 0, sd = fit.in[["frail.var.est"]])
    upper.lim <- qnorm(0.99999, mean = 0, sd = fit.in[["frail.var.est"]])
    
    ## Now run the imputation
    #surv.marg.est.A <- cubintegrate(f = surv.marg.est.A.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    #surv.marg.est.B <- cubintegrate(f = surv.marg.est.B.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    #risk.marg.est.A <- cubintegrate(f = risk.marg.est.A.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    #risk.marg.est.B <- cubintegrate(f = risk.marg.est.B.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    #surv.joint.est <- cubintegrate(f = surv.joint.est.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
    risk.joint.est <- cubintegrate(f = risk.joint.est.for.int, lower = lower.lim, upper = upper.lim, method = "pcubature")$integral
  }
  
  ## Create output object
  #     output.obj <- c("surv.marg.est.A" = surv.marg.est.A,
  #                     "surv.marg.est.B" = surv.marg.est.B,
  #                     "surv.joint.est" = surv.joint.est,
  #                     "risk.marg.est.A" = risk.marg.est.A,
  #                     "risk.marg.est.B" = risk.marg.est.B,
  #                     "risk.joint.est" = risk.joint.est)
  
  return(as.numeric(risk.joint.est))
}


########################################################################################################################################################
########################################################################################################################################################
### Function to create predicted observed risks for a given outcome, status indicator, predicted risk, number of knots and time point for evaluation ###
### (i.e. do the validation)
########################################################################################################################################################
########################################################################################################################################################
create.pred.obs <- function(time.in, status.in, risk.pred.in, n.knots.in, t.eval, model.name){
  
  #   time.in = data.valid$time_C
  #   status.in = data.valid$status_C
  #   risk.pred.in = predrisk.dual[["risk.joint.est"]]
  #   n.knots.in = 0
  #   t.eval = t.eval
  
  if (n.knots.in == 0){
    ## Create loglogsurv
    loglogsurv <- log(-log(1 - risk.pred.in))
    ## Fit validation model
    coxph.valid <- coxph(Surv(time.in, status.in) ~ loglogsurv)
    ## Create linear predictors
    obs.lp <- predict(coxph.valid, response = "lp")
    ## Baseline hazard
    basehaz.valid <- basehaz(coxph.valid, centered = TRUE) 
    ## baseline hazard at time t
    basehaz.valid.t <- basehaz.valid$hazard[max(which(basehaz.valid$time < t.eval))]
    ## Calculate predicted observed risks
    risk.obs <- 1 - exp(-basehaz.valid.t*exp(obs.lp))
    
    output <- arrange(data.frame("id" = 1:length(risk.obs), "risk.obs" = risk.obs, "risk.pred" = risk.pred.in, "model" = model.name), risk.pred.in)
    
    return(output)
    
  } else if (n.knots.in %in%c(3,4,5,6,7)){
    
    ## Create loglogsurv
    loglogsurv <- log(-log(1 - risk.pred.in))
    ## Calculate placement of knots
    Knots <- rcspline.eval(loglogsurv, nk = n.knots.in, knots.only = TRUE)
    ## Fit validation model
    coxph.valid <- coxph(Surv(time.in, status.in) ~ rcs(loglogsurv, Knots))
    ## Create linear predictors
    obs.lp <- predict(coxph.valid, response = "lp")
    ## Baseline hazard
    basehaz.valid <- basehaz(coxph.valid, centered = TRUE) 
    ## baseline hazard at time t
    basehaz.valid.t <- basehaz.valid$hazard[max(which(basehaz.valid$time < t.eval))]
    ## Calculate predicted observed risks
    risk.obs <- 1 - exp(-basehaz.valid.t*exp(obs.lp))
    
    output <- arrange(data.frame("id" = 1:length(risk.obs), "risk.obs" = risk.obs, "risk.pred" = risk.pred.in, "model" = model.name), risk.pred.in)
    
    return(output)
  }
}


#######################################################
### Function to create data for risk by decile plot ###
#######################################################
create.pred.obs.by.decile <- function(time.in, status.in, risk.pred.in, t.eval, model.name){
  
  ### Create cut offs for decile of predicted risk
  cutoffs <- quantile(risk.pred.in, probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  
  ### Create groups of individuals by decile of predicted risk
  temp.data <- data.frame("time" = time.in, "status" = status.in, "risk.pred" = risk.pred.in)
  temp.data.list <- vector("list", 10)
  for (cut in 1:10){
    if (cut == 1){
      temp.data.list[[cut]] <- temp.data[temp.data$risk.pred <= cutoffs[cut], ]
    } else if (cut != 1){
      temp.data.list[[cut]] <- temp.data[temp.data$risk.pred <= cutoffs[cut] & temp.data$risk.pred > cutoffs[(cut-1)], ]
    }
  }
  
  ### Calculate mean predicted risk within each group, and observed risk as kaplan meier estimate of survival within each group
  pred.risk.decile <- rep(NA, 10)
  obs.risk.decile <- rep(NA, 10)
  
  for (cut in 1:10){
    ### Mean predicted risk
    pred.risk.decile[cut] <- mean(temp.data.list[[cut]]$risk.pred)
    
    ### Observed risk
    ## Fit KM model
    km.model <- survfit(Surv(time, status) ~ 1, data = temp.data.list[[cut]])
    
    ## Get estimate of survival probability at first time after t.eval
    obs.risk.decile[cut] <- 1 - km.model$surv[which(km.model$time > t.eval)[1]]
    
  }
  
  return(data.frame("pred.risk.decile" = pred.risk.decile, "obs.risk.decile" = obs.risk.decile, "model" = rep(model.name, length(pred.risk.decile))))
}



#################################################################
#################################################################
### Function to fit copula model and generate predicted risks ###
#################################################################
#################################################################

calc.predrisk.copula <- function(data.devel, data.valid, t.eval, variables.vec, copula.in, rotate.in){
  
  #   ### Create validation and development datasets
  #   data.devel <- complete.data.1[1:2000, ]
  #   data.valid <- complete.data.1[2001:4000, ]
  #   str(data.devel)
  #   ### Define variables we want to adjust for
  #   variables.vec <- c("Age", "Smoking")
  #   
  #   ### Choose time point to evaluate calibration
  #   t.eval <- 3652.50
  #   
  #   ### Transform datasets ready for analysis
  #   data.devel <- transform.data.cop(data.devel, variables.vec = variables.vec)
  #   data.valid <- transform.data.cop(data.valid, variables.vec = variables.vec)
  #   
  #   data.devel = data.devel
  #   data.valid = data.valid
  #   t.eval = t.eval
  #   variables.vec = variables.vec
  #   copula.in = "clayton"
  #   rotate.in = 0
  
  ### First transform data.devel and data.valid approprately into wide format
  data.devel <- tidyr::pivot_wider(data.devel, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status))
  data.valid <- tidyr::pivot_wider(data.valid, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status))
  
  
  ### Create equations for fitting the copula model
  eqA.out <- "time_A ~ s(log(time_A), bs = 'mpi') + "
  eqB.out <- "time_B ~ s(log(time_B), bs = 'mpi') + "
  pred.formula <- paste(variables.vec, collapse = "+")
  eqA <- as.formula(paste(eqA.out, pred.formula, sep = " "))
  eqB <- as.formula(paste(eqB.out, pred.formula, sep = " "))
  
  if (copula.in == "clayton"){
    if (rotate.in == 0){
      gjrm.fit <- gjrm(list(eqA, eqB), data = data.devel, surv = TRUE, 
                       margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = "C0")
    } else if (rotate.in == 90){
      gjrm.fit <- gjrm(list(eqA, eqB), data = data.devel, surv = TRUE, 
                       margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = "C90")
    } else if (rotate.in == 180){
      gjrm.fit <- gjrm(list(eqA, eqB), data = data.devel, surv = TRUE, 
                       margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = "C180")
    } else if (rotate.in == 270){
      gjrm.fit <- gjrm(list(eqA, eqB), data = data.devel, surv = TRUE, 
                       margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = "C270")
    }
  } else if (copula.in == "gumbel"){
    if (rotate.in == 0){
      gjrm.fit <- gjrm(list(eqA, eqB), data = data.devel, surv = TRUE, 
                       margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = "G0")
    } else if (rotate.in == 90){
      gjrm.fit <- gjrm(list(eqA, eqB), data = data.devel, surv = TRUE, 
                       margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = "G90")
    } else if (rotate.in == 180){
      gjrm.fit <- gjrm(list(eqA, eqB), data = data.devel, surv = TRUE, 
                       margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = "G180")
    } else if (rotate.in == 270){
      gjrm.fit <- gjrm(list(eqA, eqB), data = data.devel, surv = TRUE, 
                       margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = "G270")
    }
  } else if (copula.in == "frank"){
    gjrm.fit <- gjrm(list(eqA, eqB), data = data.devel, surv = TRUE, 
                     margins = c("PH", "PH"), cens1 = status_A, cens2 = status_B, Model = "B", BivD = "F")
  }
  
  ### Fit cox models to extract the baseline hazards
  ## First define equations
  coxA.out <- "Surv(time_A, status_A) ~"
  coxB.out <- "Surv(time_B, status_B) ~"
  cox.pred <- paste(variables.vec, collapse = "+")
  coxformulaA <- as.formula(paste(coxA.out, cox.pred, sep = " "))
  coxformulaB <- as.formula(paste(coxB.out, cox.pred, sep = " "))
  
  ## Fit the coxph models
  coxph.cop.A <- coxph(coxformulaA, data = data.devel)
  coxph.cop.B <- coxph(coxformulaB, data = data.devel)
  
  ## Extract baseline hazards at time t.eval
  basehaz.A <- basehaz(coxph.cop.A, centered = FALSE)
  basehaz.B <- basehaz(coxph.cop.B, centered = FALSE)
  
  basehaz.A.t <- basehaz.A$hazard[max(which(basehaz.A$time < t.eval))]
  basehaz.B.t <- basehaz.B$hazard[max(which(basehaz.B$time < t.eval))]
  
  ## Extract predictors for data.valid
  model.matrix.formula <- as.formula(paste("~", paste(variables.vec, collapse = "+"), sep = " "))
  data.X <- model.matrix(model.matrix.formula, data.valid[, variables.vec])
  
  
  ## Extract coefficients from copula model
  betaA <- data.frame(summary(gjrm.fit)$tableP1)
  betaB <- data.frame(summary(gjrm.fit)$tableP2)
  eta.est <- summary(gjrm.fit)$theta
  
  ## Calc marginal survival probabilities and risks
  surv.marg.est.A <- exp(-basehaz.A.t*exp(data.X[, 2:ncol(data.X)] %*% betaA[2:nrow(betaA), "Estimate"]))
  surv.marg.est.B <- exp(-basehaz.B.t*exp(data.X[, 2:ncol(data.X)] %*% betaB[2:nrow(betaB), "Estimate"]))
  
  risk.marg.est.A <- 1 - surv.marg.est.A
  risk.marg.est.B <- 1 - surv.marg.est.B
  
  ## Create copula and rotate if appropriate
  cl.est <- archmCopula(family = tolower(copula.in), param = eta.est, dim = 2)
  if (rotate.in == 90){
    cl.est <- rotCopula(cl.est, flip = c(TRUE, FALSE))
  } else if (rotate.in == 180){
    cl.est <- rotCopula(cl.est, flip = c(TRUE, TRUE))
  } else if (rotate.in == 270){
    cl.est <- rotCopula(cl.est, flip = c(FALSE, TRUE))
  }
  
  
  ## Create a dataframe with marginal survival scores
  survs.data.frame <- cbind(surv.marg.est.A, surv.marg.est.B)
  
  ## And use these to estimate joint survival and risk, using the function defined earlier to prediction joint risk from a copula
  surv.joint.est <- apply(survs.data.frame, 1, myfunc.jointrisk.cop, copula.in = cl.est)
  risk.joint.est <- 1 - surv.marg.est.A - surv.marg.est.B + surv.joint.est
  
  #   return(list("surv.joint.est" = surv.joint.est, "risk.joint.est" = risk.joint.est, 
  #               "surv.marg.est.A" = surv.marg.est.A, "surv.marg.est.B", surv.marg.est.B))
  return(list("copula.fit" = gjrm.fit, "risk.joint.est" = as.numeric(risk.joint.est)))
}


#######################################################################
#######################################################################
### Function to fit dual outcome model and generate predicted risks ###
#######################################################################
#######################################################################

### Write a function to fit the marginal coxph models
calc.predrisk.dual <- function(data.devel, data.valid, t.eval, variables.vec){
  
  ### First transform data.devel and data.valid approprately into wide format
  data.devel <- tidyr::pivot_wider(data.devel, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status))
  data.valid <- tidyr::pivot_wider(data.valid, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status))
  
  ### Define equation for cox model
  coxC.out <- "Surv(time_C, status_C) ~"
  g <- paste(variables.vec, collapse = "+")
  coxformulaC <- as.formula(paste(coxC.out, g, sep = " "))
  
  ## Fit the coxph models
  coxph.C <- coxph(coxformulaC, data = data.devel)
  
  ## Extract baseline hazards
  basehaz.C <- basehaz(coxph.C, centered = TRUE)
  
  ## Extract baseline hazards at the time point of interest
  basehaz.C.t <- basehaz.C$hazard[max(which(basehaz.C$time < t.eval))]
  
  ## Calculate the linear predictor
  obs.lp.C <- predict(coxph.C, newdata = data.valid, type = "lp")
  
  ## Calculate risks
  risk.joint.est <- 1 - as.numeric(exp(-basehaz.C.t*exp(obs.lp.C)))
  
  ## Create output object
  output.obj <- list("risk.joint.est" = risk.joint.est)
  
  ## Return output object
  return(list("dual.fit" = coxph.C, "risk.joint.est" = as.numeric(risk.joint.est)))
  #return(output.obj)
  
}



##################################################################
##################################################################
### Function to fit product model and generate predicted risks ###
##################################################################
##################################################################

### Write a function to fit the marginal coxph models
calc.predrisk.product <- function(data.devel, data.valid, t.eval, variables.vec){
  
  ### First transform data.devel and data.valid approprately into wide format
  data.devel <- tidyr::pivot_wider(data.devel, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status))
  data.valid <- tidyr::pivot_wider(data.valid, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status))
  
  ### Define equation for cox models
  coxA.out <- "Surv(time_A, status_A) ~"
  coxB.out <- "Surv(time_B, status_B) ~"
  g <- paste(variables.vec, collapse = "+")
  coxformulaA <- as.formula(paste(coxA.out, g, sep = " "))
  coxformulaB <- as.formula(paste(coxB.out, g, sep = " "))
  
  ## Fit a coxph model for each outcome
  coxph.A <- coxph(coxformulaA, data = data.devel)
  coxph.B <- coxph(coxformulaB, data = data.devel)
  
  ## Extract baseline hazards
  basehaz.A <- basehaz(coxph.A, centered = TRUE)
  basehaz.B <- basehaz(coxph.B, centered = TRUE)
  
  ## Extract baseline hazards at the time point of interest
  basehaz.A.t <- basehaz.A$hazard[max(which(basehaz.A$time < t.eval))]
  basehaz.B.t <- basehaz.B$hazard[max(which(basehaz.B$time < t.eval))]
  
  ## Create linear preditor
  obs.lp.A <- predict(coxph.A, newdata = data.valid, type = "lp")
  obs.lp.B <- predict(coxph.B, newdata = data.valid, type = "lp")
  
  
  ## Calculate risks
  risk.marg.est.A <- 1 - as.numeric(exp(-basehaz.A.t*exp(obs.lp.A)))
  risk.marg.est.B <- 1 - as.numeric(exp(-basehaz.B.t*exp(obs.lp.B)))
  
  ## Calculate joint risk
  risk.joint.est <- risk.marg.est.A*risk.marg.est.B
  
  ## Return output object
  #return(risk.joint.est)
  return(list("risk.joint.est" = risk.joint.est))
  
}


##############################################################
##############################################################
### Function to fit msm model and generate predicted risks ###
##############################################################
##############################################################

### Write a function to fit the marginal coxph models
calc.predrisk.msm <- function(data.devel, data.valid, t.eval, variables.vec, msm.iter, msm.chunk){
  
  #   data.devel = data.devel
  #   data.valid = data.valid
  #   t.eval = t.eval
  #   variables.vec = variables.vec
  #   msm.iter = 1
  #   msm.chunk = 10
  
  ### First transform data.devel and data.valid approprately into wide format
  data.devel.wide <- tidyr::pivot_wider(data.devel, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status))
  data.valid.wide <- tidyr::pivot_wider(data.valid, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status))
  
  ### Write function that gets data into format for msm
  transform.data.msm <- function(data.in){
    ## Remove duplicate predictors variables and rename
    data.anal.pre <- data.in %>% 
      dplyr::select(id, time_A, time_B, status_A, status_B, all_of(variables.vec))
    
    ## Assign event times
    data.anal.pre$state1 <- rep(0, nrow(data.anal.pre))
    data.anal.pre$state2 <- data.anal.pre$time_A
    data.anal.pre$state3 <- data.anal.pre$time_B
    data.anal.pre$state4 <- pmax(data.anal.pre$time_A, data.anal.pre$time_B)
    
    ## Assign censoring indicator  
    ## First assign all 0 (censored)
    data.anal.pre$state2.s <- rep(0, nrow(data.anal.pre))
    data.anal.pre$state3.s <- rep(0, nrow(data.anal.pre))
    data.anal.pre$state4.s <- rep(0, nrow(data.anal.pre))
    
    ## Entry to state2 will be observed if time_A < time_B and status_A = 1 
    ## Entry to state3 will be observed if time_B < time_A and status_B = 1 
    ## Entry to state4 will be observed if status_A and status_B = 1
    data.anal.pre <- data.anal.pre %>% 
      mutate(state2.s = case_when(time_A < time_B & status_A == 1 ~ 1, TRUE ~ as.numeric(as.character(state2.s))),
             state3.s = case_when(time_B < time_A & status_B == 1 ~ 1, TRUE ~ as.numeric(as.character(state3.s))),
             state4.s = case_when(status_A == 1 & status_B == 1 ~ 1, TRUE ~ as.numeric(as.character(state4.s)))
      )
    
    ## Select only variables of interest
    data.anal.pre <- dplyr::select(data.anal.pre, id, state1, state2, state3, state4, state2.s, state3.s, state4.s, all_of(variables.vec))
    
    ## Note that we only observe the initial transition into A/B if it happens first 
    ## (i.e. the event time for entry into state 2 and state 3 should be the same, 
    ## with censoring on whichever event doesnt happen first)
    ## This will be dealt with using the msprep functionality to deal with this
    
    ## Change from tibble into data.frame
    data.anal.pre <- data.frame(data.anal.pre)
    
    return(data.anal.pre)
  }
  
  ### Get both development and validation dataset into appropriate format
  data.devel.msm <- transform.data.msm(data.devel.wide)
  data.valid.msm <- transform.data.msm(data.valid.wide)
  
  ## Create a transition matrix for the allowed transitions
  tmat <- transMat(x = list(c(2,3), c(4), c(4), c()), 
                   names = c("state1", "state2", "state3", "state4"))
  
  ## Now can prepr the data into format which can be analysed uing mstate
  data.anal <- msprep(data.devel.msm, trans = tmat, time = c(NA, "state2", "state3", "state4"),
                      status = c(NA, "state2.s", "state3.s", "state4.s"), keep = variables.vec)
  
  ## Note how the censoring transitions into state 2/3 are now done so at the correct time
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- variables.vec
  data.anal <- expand.covs(data.anal, covs, longnames = FALSE)
  
  ## For creating model formula, I want to extract a vector just with expanded variable names
  ## Do this by removing all the rows that wer ealso present on data.devel.msm (non expanded data frame)
  msmformula.pred <- paste(colnames(data.anal)[-(1:ncol(data.devel.msm))], collapse = "+")
  
  ## Create formula for msm
  msmformula.out <- "Surv(Tstart, Tstop, status) ~ "
  msmformula <- as.formula(paste(msmformula.out, msmformula.pred, "+ strata(trans)", sep = " "))
  
  ### Now to fit the cause-specific hazard models MSM and store for output
  msm.model <- coxph(msmformula, 
                     data = data.anal, method = "breslow")
  
  #   ### Also save the development data, the format of which will be used to create "newdata" arguments for when we generate risk scores
  #   msm.data <- data.anal.pre
  
  ### Also save transition matrix
  msm.tmat <- tmat
  
  ### Apply function defined previously to calculate risks in validation dataset
  risk.joint.est <- apply(data.valid.msm[((msm.iter-1)*(msm.chunk)+1):(msm.iter*msm.chunk), ], 1, myfunc.jointrisk.msm,
                          msm.data.in = data.valid.msm, 
                          msm.model.in = msm.model, 
                          msm.tmat.in = msm.tmat,
                          variables.vec = variables.vec)
  
  return(list("msm.fit" = msm.model, "risk.joint.est" = as.numeric(risk.joint.est)))
  
}


##################################################################
##################################################################
### Function to fit frailty model and generate predicted risks ###
##################################################################
##################################################################

calc.predrisk.frailty <- function(data.devel, data.valid, t.eval, variables.vec, frail.dist, baseline.dist, n.iter, seed.frail){
  
  stopifnot(frail.dist %in% c("gamma","normal"), local = TRUE)
  
  #   data.devel = data.devel
  #   data.valid = data.devel
  #   t.eval = t.eval
  #   variables.vec = variables.vec
  #   frail.dist = "normal"
  #   baseline.dist = "weibull"
  #   n.iter = 200
  
  ## Remove observations that are type C
  data.devel.removeC <- data.devel[data.devel$outcome_char %in% c("A", "B"), ]
  data.valid.removeC <- data.valid[data.valid$outcome_char %in% c("A", "B"), ]
  
  ## Get data in format for analysis with rstan
  data.devel.uncens.A <- filter(data.devel.removeC, status == 1 & outcome == 1)
  data.devel.uncens.B <- filter(data.devel.removeC, status == 1 & outcome == 2)
  data.devel.cens.A <- filter(data.devel.removeC, status == 0 & outcome == 1)
  data.devel.cens.B <- filter(data.devel.removeC, status == 0 & outcome == 2)
  
  X_cens_A.temp <- data.devel.cens.A %>% dplyr::select(all_of(variables.vec))
  X_uncens_A.temp <- data.devel.uncens.A %>% dplyr::select(all_of(variables.vec))
  X_cens_B.temp <- data.devel.cens.B %>% dplyr::select(all_of(variables.vec))
  X_uncens_B.temp <- data.devel.uncens.B %>% dplyr::select(all_of(variables.vec))
  
  X_cens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_A.temp)[, -1] %>% data.matrix()
  X_uncens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_A.temp)[, -1] %>% data.matrix()
  X_cens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_B.temp)[, -1] %>% data.matrix()
  X_uncens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_B.temp)[, -1] %>% data.matrix()
  
  stan_data <- list(N = length(unique(data.devel.removeC$id)),
                    N_cens_A = length(unique(data.devel.cens.A$id)),
                    N_uncens_A = length(unique(data.devel.uncens.A$id)),
                    N_cens_B = length(unique(data.devel.cens.B$id)),
                    N_uncens_B = length(unique(data.devel.uncens.B$id)),
                    #stacked_N = nrow(data.devel.cens),
                    
                    P = ncol(X_cens_A.temp),
                    
                    IDs_cens_A = data.devel.cens.A$id,
                    IDs_uncens_A = data.devel.uncens.A$id,
                    IDs_cens_B = data.devel.cens.B$id,
                    IDs_uncens_B = data.devel.uncens.B$id,
                    
                    X_cens_A = X_cens_A.temp,
                    X_uncens_A = X_uncens_A.temp,
                    X_cens_B = X_cens_B.temp,
                    X_uncens_B = X_uncens_B.temp,
                    
                    times_cens_A = data.devel.cens.A$time,
                    times_uncens_A = data.devel.uncens.A$time,
                    times_cens_B = data.devel.cens.B$time,
                    times_uncens_B = data.devel.uncens.B$time)
  
  ## Assign the appropriate stan model depending on what frailty distribution
  if (frail.dist == "normal" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_exp.stan")
    fit <- sampling(sm, data=stan_data, seed=seed.frail, chains=2, cores=1, iter=n.iter)
  } else if (frail.dist == "gamma" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_exp.stan")
    fit <- sampling(sm, data=stan_data, seed=seed.frail, chains=2, cores=1, iter=n.iter)
  } else if (frail.dist == "normal" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_weibull_inverted.stan")
    fit <- sampling(sm, data=stan_data, seed=seed.frail, chains=2, cores=1, iter=n.iter)
  } else if (frail.dist == "gamma" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_weibull_inverted.stan")
    fit <- sampling(sm, data=stan_data, seed=seed.frail, chains=2, cores=1, iter=n.iter)
  }
  
  ## Assign estimates to an object and output
  if (baseline.dist == "exp"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  } else if (baseline.dist == "weibull"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "shape_A", "shape_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  }
  
  ## Create an object with output that can be used in the myfunc.jointrisk.frailty to estimate joint risk
  frail.fit <- vector("list", 11)
  if (baseline.dist == "weibull"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- model.estimates$summary["shape_A","mean"]
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- model.estimates$summary["shape_B","mean"]
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  } else if (baseline.dist == "exp"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- 1
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- 1
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  }
  frail.fit[[9]] <- model.estimates$summary["frail_param","mean"]
  frail.fit[[10]] <- frail.dist
  frail.fit[[11]] <- baseline.dist
  
  names(frail.fit) <- c("betas_A", "betas_B", "model.summary", "NA", 
                        "bh.shape.A", "bh.int.A", "bh.shape.B", "bh.int.B",
                        "frail.var.est", 
                        "frail.dist", "bh.dist")
  
  ### Use the output from the frailty model to run numerical integration to calculate joint risk
  ### This will be done using previously defined function, myfunc.jointrisk.frailty
  
  ### Want to apply this to each observation in vlaidaiton dataset, so first get it into wide format, one row per individual,
  ### and just the predictor variables in model.matrix format (i.e. dummy variables for categorical vars)
  data.valid.for.pred <- tidyr::pivot_wider(data.valid.removeC, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status)) %>%
    dplyr::select(all_of(variables.vec))
  
  data.valid.for.pred <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), 
                                      data = data.valid.for.pred)[, -1]
  
  ### Now apply the function
  risk.joint.est <- apply(data.valid.for.pred, 1, myfunc.jointrisk.frailty, 
                          betas_A = frail.fit[["betas_A"]], betas_B = frail.fit[["betas_B"]], fit.in = frail.fit)
  
  return(list("rstan.model" = fit, "frail.fit" = frail.fit, "risk.joint.est" = as.numeric(risk.joint.est)))
  
}


#####################################################
#####################################################
### Same again but with parallelisation of chains ###
#####################################################
#####################################################
calc.predrisk.frailty.parallel.10 <- function(data.devel, data.valid, t.eval, variables.vec, frail.dist, baseline.dist, n.iter){
  
  stopifnot(frail.dist %in% c("gamma","normal"), local = TRUE)
  
  #   data.devel = data.devel
  #   data.valid = data.devel
  #   t.eval = t.eval
  #   variables.vec = variables.vec
  #   frail.dist = "normal"
  #   baseline.dist = "weibull"
  #   n.iter = 200
  
  ## Remove observations that are type C
  data.devel.removeC <- data.devel[data.devel$outcome_char %in% c("A", "B"), ]
  data.valid.removeC <- data.valid[data.valid$outcome_char %in% c("A", "B"), ]
  
  ## Get data in format for analysis with rstan
  data.devel.uncens.A <- filter(data.devel.removeC, status == 1 & outcome == 1)
  data.devel.uncens.B <- filter(data.devel.removeC, status == 1 & outcome == 2)
  data.devel.cens.A <- filter(data.devel.removeC, status == 0 & outcome == 1)
  data.devel.cens.B <- filter(data.devel.removeC, status == 0 & outcome == 2)
  
  X_cens_A.temp <- data.devel.cens.A %>% dplyr::select(all_of(variables.vec))
  X_uncens_A.temp <- data.devel.uncens.A %>% dplyr::select(all_of(variables.vec))
  X_cens_B.temp <- data.devel.cens.B %>% dplyr::select(all_of(variables.vec))
  X_uncens_B.temp <- data.devel.uncens.B %>% dplyr::select(all_of(variables.vec))
  
  X_cens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_A.temp)[, -1] %>% data.matrix()
  X_uncens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_A.temp)[, -1] %>% data.matrix()
  X_cens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_B.temp)[, -1] %>% data.matrix()
  X_uncens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_B.temp)[, -1] %>% data.matrix()
  
  stan_data <- list(N = length(unique(data.devel.removeC$id)),
                    N_cens_A = length(unique(data.devel.cens.A$id)),
                    N_uncens_A = length(unique(data.devel.uncens.A$id)),
                    N_cens_B = length(unique(data.devel.cens.B$id)),
                    N_uncens_B = length(unique(data.devel.uncens.B$id)),
                    #stacked_N = nrow(data.devel.cens),
                    
                    P = ncol(X_cens_A.temp),
                    
                    IDs_cens_A = data.devel.cens.A$id,
                    IDs_uncens_A = data.devel.uncens.A$id,
                    IDs_cens_B = data.devel.cens.B$id,
                    IDs_uncens_B = data.devel.uncens.B$id,
                    
                    X_cens_A = X_cens_A.temp,
                    X_uncens_A = X_uncens_A.temp,
                    X_cens_B = X_cens_B.temp,
                    X_uncens_B = X_uncens_B.temp,
                    
                    times_cens_A = data.devel.cens.A$time,
                    times_uncens_A = data.devel.uncens.A$time,
                    times_cens_B = data.devel.cens.B$time,
                    times_uncens_B = data.devel.uncens.B$time)
  
  ## Assign the appropriate stan model depending on what frailty distribution
  if (frail.dist == "normal" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_exp.stan")
    cl <- makeCluster(11)
    registerDoParallel(11)
    fit <- (foreach(input=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  } else if (frail.dist == "gamma" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_exp.stan")
    cl <- makeCluster(11)
    registerDoParallel(11)
    fit <- (foreach(input=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  } else if (frail.dist == "normal" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_weibull_inverted.stan")
    cl <- makeCluster(11)
    registerDoParallel(11)
    fit <- (foreach(input=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  } else if (frail.dist == "gamma" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_weibull_inverted.stan")
    cl <- makeCluster(11)
    registerDoParallel(11)
    fit <- (foreach(input=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  }
  
  ### Combine output from each chain into a single stan object
  fit <- sflist2stanfit(fit)
  
  ## Assign estimates to an object and output
  if (baseline.dist == "exp"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  } else if (baseline.dist == "weibull"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "shape_A", "shape_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  }
  
  ## Create an object with output that can be used in the myfunc.jointrisk.frailty to estimate joint risk
  frail.fit <- vector("list", 11)
  if (baseline.dist == "weibull"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- model.estimates$summary["shape_A","mean"]
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- model.estimates$summary["shape_B","mean"]
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  } else if (baseline.dist == "exp"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- 1
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- 1
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  }
  frail.fit[[9]] <- model.estimates$summary["frail_param","mean"]
  frail.fit[[10]] <- frail.dist
  frail.fit[[11]] <- baseline.dist
  
  names(frail.fit) <- c("betas_A", "betas_B", "model.summary", "NA", 
                        "bh.shape.A", "bh.int.A", "bh.shape.B", "bh.int.B",
                        "frail.var.est", 
                        "frail.dist", "bh.dist")
  
  ### Use the output from the frailty model to run numerical integration to calculate joint risk
  ### This will be done using previously defined function, myfunc.jointrisk.frailty
  
  ### Want to apply this to each observation in vlaidaiton dataset, so first get it into wide format, one row per individual,
  ### and just the predictor variables in model.matrix format (i.e. dummy variables for categorical vars)
  data.valid.for.pred <- tidyr::pivot_wider(data.valid.removeC, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status)) %>%
    dplyr::select(all_of(variables.vec))
  
  data.valid.for.pred <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), 
                                      data = data.valid.for.pred)[, -1]
  
  ### Now apply the function
  risk.joint.est <- apply(data.valid.for.pred, 1, myfunc.jointrisk.frailty, 
                          betas_A = frail.fit[["betas_A"]], betas_B = frail.fit[["betas_B"]], fit.in = frail.fit)
  
  return(list("rstan.model" = fit, "frail.fit" = frail.fit, "risk.joint.est" = as.numeric(risk.joint.est)))
  
}


#####################################################
#####################################################
### Same again but with parallelisation of chains ###
#####################################################
#####################################################
calc.predrisk.frailty.parallel.5 <- function(data.devel, data.valid, t.eval, variables.vec, frail.dist, baseline.dist, n.iter){
  
  stopifnot(frail.dist %in% c("gamma","normal"), local = TRUE)
  
  #   data.devel = data.devel
  #   data.valid = data.devel
  #   t.eval = t.eval
  #   variables.vec = variables.vec
  #   frail.dist = "normal"
  #   baseline.dist = "weibull"
  #   n.iter = 200
  
  ## Remove observations that are type C
  data.devel.removeC <- data.devel[data.devel$outcome_char %in% c("A", "B"), ]
  data.valid.removeC <- data.valid[data.valid$outcome_char %in% c("A", "B"), ]
  
  ## Get data in format for analysis with rstan
  data.devel.uncens.A <- filter(data.devel.removeC, status == 1 & outcome == 1)
  data.devel.uncens.B <- filter(data.devel.removeC, status == 1 & outcome == 2)
  data.devel.cens.A <- filter(data.devel.removeC, status == 0 & outcome == 1)
  data.devel.cens.B <- filter(data.devel.removeC, status == 0 & outcome == 2)
  
  X_cens_A.temp <- data.devel.cens.A %>% dplyr::select(all_of(variables.vec))
  X_uncens_A.temp <- data.devel.uncens.A %>% dplyr::select(all_of(variables.vec))
  X_cens_B.temp <- data.devel.cens.B %>% dplyr::select(all_of(variables.vec))
  X_uncens_B.temp <- data.devel.uncens.B %>% dplyr::select(all_of(variables.vec))
  
  X_cens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_A.temp)[, -1] %>% data.matrix()
  X_uncens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_A.temp)[, -1] %>% data.matrix()
  X_cens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_B.temp)[, -1] %>% data.matrix()
  X_uncens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_B.temp)[, -1] %>% data.matrix()
  
  stan_data <- list(N = length(unique(data.devel.removeC$id)),
                    N_cens_A = length(unique(data.devel.cens.A$id)),
                    N_uncens_A = length(unique(data.devel.uncens.A$id)),
                    N_cens_B = length(unique(data.devel.cens.B$id)),
                    N_uncens_B = length(unique(data.devel.uncens.B$id)),
                    #stacked_N = nrow(data.devel.cens),
                    
                    P = ncol(X_cens_A.temp),
                    
                    IDs_cens_A = data.devel.cens.A$id,
                    IDs_uncens_A = data.devel.uncens.A$id,
                    IDs_cens_B = data.devel.cens.B$id,
                    IDs_uncens_B = data.devel.uncens.B$id,
                    
                    X_cens_A = X_cens_A.temp,
                    X_uncens_A = X_uncens_A.temp,
                    X_cens_B = X_cens_B.temp,
                    X_uncens_B = X_uncens_B.temp,
                    
                    times_cens_A = data.devel.cens.A$time,
                    times_uncens_A = data.devel.uncens.A$time,
                    times_cens_B = data.devel.cens.B$time,
                    times_uncens_B = data.devel.uncens.B$time)
  
  ## Assign the appropriate stan model depending on what frailty distribution
  if (frail.dist == "normal" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_exp.stan")
    cl <- makeCluster(6)
    registerDoParallel(6)
    fit <- (foreach(input=c(1, 2, 3, 4, 5), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  } else if (frail.dist == "gamma" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_exp.stan")
    cl <- makeCluster(6)
    registerDoParallel(6)
    fit <- (foreach(input=c(1, 2, 3, 4, 5), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  } else if (frail.dist == "normal" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_weibull_inverted.stan")
    cl <- makeCluster(6)
    registerDoParallel(6)
    fit <- (foreach(input=c(1, 2, 3, 4, 5), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  } else if (frail.dist == "gamma" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_weibull_inverted.stan")
    cl <- makeCluster(6)
    registerDoParallel(6)
    fit <- (foreach(input=c(1, 2, 3, 4, 5), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  }
  
  ### Combine output from each chain into a single stan object
  fit <- sflist2stanfit(fit)
  
  ## Assign estimates to an object and output
  if (baseline.dist == "exp"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  } else if (baseline.dist == "weibull"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "shape_A", "shape_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  }
  
  ## Create an object with output that can be used in the myfunc.jointrisk.frailty to estimate joint risk
  frail.fit <- vector("list", 11)
  if (baseline.dist == "weibull"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- model.estimates$summary["shape_A","mean"]
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- model.estimates$summary["shape_B","mean"]
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  } else if (baseline.dist == "exp"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- 1
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- 1
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  }
  frail.fit[[9]] <- model.estimates$summary["frail_param","mean"]
  frail.fit[[10]] <- frail.dist
  frail.fit[[11]] <- baseline.dist
  
  names(frail.fit) <- c("betas_A", "betas_B", "model.summary", "NA", 
                        "bh.shape.A", "bh.int.A", "bh.shape.B", "bh.int.B",
                        "frail.var.est", 
                        "frail.dist", "bh.dist")
  
  ### Use the output from the frailty model to run numerical integration to calculate joint risk
  ### This will be done using previously defined function, myfunc.jointrisk.frailty
  
  ### Want to apply this to each observation in vlaidaiton dataset, so first get it into wide format, one row per individual,
  ### and just the predictor variables in model.matrix format (i.e. dummy variables for categorical vars)
  data.valid.for.pred <- tidyr::pivot_wider(data.valid.removeC, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status)) %>%
    dplyr::select(all_of(variables.vec))
  
  data.valid.for.pred <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), 
                                      data = data.valid.for.pred)[, -1]
  
  ### Now apply the function
  risk.joint.est <- apply(data.valid.for.pred, 1, myfunc.jointrisk.frailty, 
                          betas_A = frail.fit[["betas_A"]], betas_B = frail.fit[["betas_B"]], fit.in = frail.fit)
  
  return(list("rstan.model" = fit, "frail.fit" = frail.fit, "risk.joint.est" = as.numeric(risk.joint.est)))
  
}


#####################################################
#####################################################
### Same again but with parallelisation of chains, and allowing specification of initial values ###
#####################################################
#####################################################
calc.predrisk.frailty.parallel.5.init <- function(data.devel, data.valid, t.eval, variables.vec, frail.dist, baseline.dist, n.iter){
  
  stopifnot(frail.dist %in% c("gamma","normal"), local = TRUE)
  
  #   data.devel = data.devel
  #   data.valid = data.devel
  #   t.eval = t.eval
  #   variables.vec = variables.vec
  #   frail.dist = "normal"
  #   baseline.dist = "weibull"
  #   n.iter = 200
  
  ## Remove observations that are type C
  data.devel.removeC <- data.devel[data.devel$outcome_char %in% c("A", "B"), ]
  data.valid.removeC <- data.valid[data.valid$outcome_char %in% c("A", "B"), ]
  
  ## Get data in format for analysis with rstan
  data.devel.uncens.A <- filter(data.devel.removeC, status == 1 & outcome == 1)
  data.devel.uncens.B <- filter(data.devel.removeC, status == 1 & outcome == 2)
  data.devel.cens.A <- filter(data.devel.removeC, status == 0 & outcome == 1)
  data.devel.cens.B <- filter(data.devel.removeC, status == 0 & outcome == 2)
  
  X_cens_A.temp <- data.devel.cens.A %>% dplyr::select(all_of(variables.vec))
  X_uncens_A.temp <- data.devel.uncens.A %>% dplyr::select(all_of(variables.vec))
  X_cens_B.temp <- data.devel.cens.B %>% dplyr::select(all_of(variables.vec))
  X_uncens_B.temp <- data.devel.uncens.B %>% dplyr::select(all_of(variables.vec))
  
  X_cens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_A.temp)[, -1] %>% data.matrix()
  X_uncens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_A.temp)[, -1] %>% data.matrix()
  X_cens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_B.temp)[, -1] %>% data.matrix()
  X_uncens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_B.temp)[, -1] %>% data.matrix()
  
  stan_data <- list(N = length(unique(data.devel.removeC$id)),
                    N_cens_A = length(unique(data.devel.cens.A$id)),
                    N_uncens_A = length(unique(data.devel.uncens.A$id)),
                    N_cens_B = length(unique(data.devel.cens.B$id)),
                    N_uncens_B = length(unique(data.devel.uncens.B$id)),
                    #stacked_N = nrow(data.devel.cens),
                    
                    P = ncol(X_cens_A.temp),
                    
                    IDs_cens_A = data.devel.cens.A$id,
                    IDs_uncens_A = data.devel.uncens.A$id,
                    IDs_cens_B = data.devel.cens.B$id,
                    IDs_uncens_B = data.devel.uncens.B$id,
                    
                    X_cens_A = X_cens_A.temp,
                    X_uncens_A = X_uncens_A.temp,
                    X_cens_B = X_cens_B.temp,
                    X_uncens_B = X_uncens_B.temp,
                    
                    times_cens_A = data.devel.cens.A$time,
                    times_uncens_A = data.devel.uncens.A$time,
                    times_cens_B = data.devel.cens.B$time,
                    times_uncens_B = data.devel.uncens.B$time)
  
  ### Define initial values
  
  ### Create functions to generate initial values at random for each chain
  ## The only differences between specification of initial values is frail param, we have a shorter upper limit for normal frailty
  ## because this is applied after being exponentiated.
  
  ## Also the coefficient for age is initialised between -0.1 and 0.1, whereas all others are between -1.5 and 1.5
  set_inits_normal <- function(seed = 1){
    set.seed(seed)
    return(list(shape_A = runif(1, 0.5, 1.5), 
                intercept_A = runif(1, -16, -13), 
                betas_A = c(runif(1, -0.1, 0.1), runif((stan_data$P - 1), -1.5, 1.5)),
                shape_B = runif(1, 0.5, 1.5), 
                intercept_B = runif(1, -13, -9), 
                betas_B = c(runif(1, -0.1, 0.1), runif((stan_data$P - 1), -1.5, 1.5)),
                frail_param = runif(1, 0.5, 1.5))
           )
  }
  
  set_inits_gamma <- function(seed = 1){
    set.seed(seed)
    return(list(shape_A = runif(1, 0.5, 1.5), 
                intercept_A = runif(1, -16, -13), 
                betas_A = c(runif(1, -0.1, 0.1), runif((stan_data$P - 1), -1.5, 1.5)),
                shape_B = runif(1, 0.5, 1.5), 
                intercept_B = runif(1, -13, -9), 
                betas_B = c(runif(1, -0.1, 0.1), runif((stan_data$P - 1), -1.5, 1.5)),
                frail_param = runif(1, 0.5, 2))
    )
  }

  ### Create the initial values for each using above functions
  init.values.in.normal <- list(
    set_inits_normal(seed = 1),
    set_inits_normal(seed = 2),
    set_inits_normal(seed = 3),
    set_inits_normal(seed = 4),
    set_inits_normal(seed = 5)
    )
  
  init.values.in.gamma <- list(
    set_inits_gamma(seed = 1),
    set_inits_gamma(seed = 2),
    set_inits_gamma(seed = 3),
    set_inits_gamma(seed = 4),
    set_inits_gamma(seed = 5)
  )
  
  ## Assign the appropriate stan model depending on what frailty distribution
  if (frail.dist == "normal" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_exp.stan")
    cl <- makeCluster(6)
    registerDoParallel(6)
    fit <- (foreach(input=c(1, 2, 3, 4, 5), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  } else if (frail.dist == "gamma" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_exp.stan")
    cl <- makeCluster(6)
    registerDoParallel(6)
    fit <- (foreach(input=c(1, 2, 3, 4, 5), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter)
            })
    stopCluster(cl)
  } else if (frail.dist == "normal" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_weibull_inverted.stan")
    cl <- makeCluster(6)
    registerDoParallel(6)
    fit <- (foreach(input=c(1, 2, 3, 4, 5), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter, init = list(init.values.in.normal[[input]]))
            })
    stopCluster(cl)
  } else if (frail.dist == "gamma" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_weibull_inverted.stan")
    cl <- makeCluster(6)
    registerDoParallel(6)
    fit <- (foreach(input=c(1, 2, 3, 4, 5), .combine=list, .multicombine=TRUE, 
                    .packages=c("rstan"))
            %dopar%{fit <- sampling(sm, data=stan_data, seed=input, chains=1, cores=1, iter=n.iter, init = list(init.values.in.gamma[[input]]))
            })
    stopCluster(cl)
  }
  
  ### Combine output from each chain into a single stan object
  fit <- sflist2stanfit(fit)

  ## Assign estimates to an object and output
  if (baseline.dist == "exp"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  } else if (baseline.dist == "weibull"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "shape_A", "shape_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  }
  
  ## Create an object with output that can be used in the myfunc.jointrisk.frailty to estimate joint risk
  frail.fit <- vector("list", 11)
  if (baseline.dist == "weibull"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- model.estimates$summary["shape_A","mean"]
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- model.estimates$summary["shape_B","mean"]
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  } else if (baseline.dist == "exp"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- 1
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- 1
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  }
  frail.fit[[9]] <- model.estimates$summary["frail_param","mean"]
  frail.fit[[10]] <- frail.dist
  frail.fit[[11]] <- baseline.dist
  
  names(frail.fit) <- c("betas_A", "betas_B", "model.summary", "NA", 
                        "bh.shape.A", "bh.int.A", "bh.shape.B", "bh.int.B",
                        "frail.var.est", 
                        "frail.dist", "bh.dist")
  
  ### Use the output from the frailty model to run numerical integration to calculate joint risk
  ### This will be done using previously defined function, myfunc.jointrisk.frailty
  
  ### Want to apply this to each observation in vlaidaiton dataset, so first get it into wide format, one row per individual,
  ### and just the predictor variables in model.matrix format (i.e. dummy variables for categorical vars)
  data.valid.for.pred <- tidyr::pivot_wider(data.valid.removeC, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status)) %>%
    dplyr::select(all_of(variables.vec))
  
  data.valid.for.pred <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), 
                                      data = data.valid.for.pred)[, -1]
  
  ### Now apply the function
  risk.joint.est <- apply(data.valid.for.pred, 1, myfunc.jointrisk.frailty, 
                          betas_A = frail.fit[["betas_A"]], betas_B = frail.fit[["betas_B"]], fit.in = frail.fit)
  
  return(list("rstan.model" = fit, "frail.fit" = frail.fit, "risk.joint.est" = as.numeric(risk.joint.est)))
  
}

#####################################################
#####################################################
### Same again but with auto parallelisation of chains (not done through doParallel) ###
#####################################################
#####################################################
calc.predrisk.frailty.parallel.5.auto <- function(data.devel, data.valid, t.eval, variables.vec, frail.dist, baseline.dist, n.iter){
  
  stopifnot(frail.dist %in% c("gamma","normal"), local = TRUE)
  
  #   data.devel = data.devel
  #   data.valid = data.devel
  #   t.eval = t.eval
  #   variables.vec = variables.vec
  #   frail.dist = "normal"
  #   baseline.dist = "weibull"
  #   n.iter = 200
  
  ## Remove observations that are type C
  data.devel.removeC <- data.devel[data.devel$outcome_char %in% c("A", "B"), ]
  data.valid.removeC <- data.valid[data.valid$outcome_char %in% c("A", "B"), ]
  
  ## Get data in format for analysis with rstan
  data.devel.uncens.A <- filter(data.devel.removeC, status == 1 & outcome == 1)
  data.devel.uncens.B <- filter(data.devel.removeC, status == 1 & outcome == 2)
  data.devel.cens.A <- filter(data.devel.removeC, status == 0 & outcome == 1)
  data.devel.cens.B <- filter(data.devel.removeC, status == 0 & outcome == 2)
  
  X_cens_A.temp <- data.devel.cens.A %>% dplyr::select(all_of(variables.vec))
  X_uncens_A.temp <- data.devel.uncens.A %>% dplyr::select(all_of(variables.vec))
  X_cens_B.temp <- data.devel.cens.B %>% dplyr::select(all_of(variables.vec))
  X_uncens_B.temp <- data.devel.uncens.B %>% dplyr::select(all_of(variables.vec))
  
  X_cens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_A.temp)[, -1] %>% data.matrix()
  X_uncens_A.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_A.temp)[, -1] %>% data.matrix()
  X_cens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_cens_B.temp)[, -1] %>% data.matrix()
  X_uncens_B.temp <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), data = X_uncens_B.temp)[, -1] %>% data.matrix()
  
  stan_data <- list(N = length(unique(data.devel.removeC$id)),
                    N_cens_A = length(unique(data.devel.cens.A$id)),
                    N_uncens_A = length(unique(data.devel.uncens.A$id)),
                    N_cens_B = length(unique(data.devel.cens.B$id)),
                    N_uncens_B = length(unique(data.devel.uncens.B$id)),
                    #stacked_N = nrow(data.devel.cens),
                    
                    P = ncol(X_cens_A.temp),
                    
                    IDs_cens_A = data.devel.cens.A$id,
                    IDs_uncens_A = data.devel.uncens.A$id,
                    IDs_cens_B = data.devel.cens.B$id,
                    IDs_uncens_B = data.devel.uncens.B$id,
                    
                    X_cens_A = X_cens_A.temp,
                    X_uncens_A = X_uncens_A.temp,
                    X_cens_B = X_cens_B.temp,
                    X_uncens_B = X_uncens_B.temp,
                    
                    times_cens_A = data.devel.cens.A$time,
                    times_uncens_A = data.devel.uncens.A$time,
                    times_cens_B = data.devel.cens.B$time,
                    times_uncens_B = data.devel.uncens.B$time)
  
  ## Assign the appropriate stan model depending on what frailty distribution
  if (frail.dist == "normal" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_exp.stan")
    fit <- sampling(sm, data=stan_data, seed=1, chains=5, cores=5, iter=n.iter)
  } else if (frail.dist == "gamma" & baseline.dist == "exp"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_exp.stan")
    fit <- sampling(sm, data=stan_data, seed=1, chains=5, cores=5, iter=n.iter)
  } else if (frail.dist == "normal" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_normal_weibull_inverted.stan")
    fit <- sampling(sm, data=stan_data, seed=1, chains=5, cores=5, iter=n.iter)
  } else if (frail.dist == "gamma" & baseline.dist == "weibull"){
    sm <- rstan::stan_model("Project 4/code/ExpSurvModel_multivariate_vectorized_gamma_weibull_inverted.stan")
    fit <- sampling(sm, data=stan_data, seed=1, chains=5, cores=5, iter=n.iter)
  }
  
  ## Assign estimates to an object and output
  if (baseline.dist == "exp"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  } else if (baseline.dist == "weibull"){
    model.estimates <- summary(fit,  par = c("betas_A", "betas_B",
                                             "shape_A", "shape_B",
                                             "intercept_A", "intercept_B",
                                             "frail_param"))
  }
  
  ## Create an object with output that can be used in the myfunc.jointrisk.frailty to estimate joint risk
  frail.fit <- vector("list", 11)
  if (baseline.dist == "weibull"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- model.estimates$summary["shape_A","mean"]
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- model.estimates$summary["shape_B","mean"]
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  } else if (baseline.dist == "exp"){
    frail.fit[[1]] <- model.estimates$summary[grep("betas_A", rownames(model.estimates$summary)), "mean"]
    frail.fit[[2]] <- model.estimates$summary[grep("betas_B", rownames(model.estimates$summary)), "mean"]
    frail.fit[[3]] <- model.estimates
    frail.fit[[4]] <- NA
    frail.fit[[5]] <- 1
    frail.fit[[6]] <- model.estimates$summary["intercept_A","mean"]
    frail.fit[[7]] <- 1
    frail.fit[[8]] <- model.estimates$summary["intercept_B","mean"]
  }
  frail.fit[[9]] <- model.estimates$summary["frail_param","mean"]
  frail.fit[[10]] <- frail.dist
  frail.fit[[11]] <- baseline.dist
  
  names(frail.fit) <- c("betas_A", "betas_B", "model.summary", "NA", 
                        "bh.shape.A", "bh.int.A", "bh.shape.B", "bh.int.B",
                        "frail.var.est", 
                        "frail.dist", "bh.dist")
  
  ### Use the output from the frailty model to run numerical integration to calculate joint risk
  ### This will be done using previously defined function, myfunc.jointrisk.frailty
  
  ### Want to apply this to each observation in vlaidaiton dataset, so first get it into wide format, one row per individual,
  ### and just the predictor variables in model.matrix format (i.e. dummy variables for categorical vars)
  data.valid.for.pred <- tidyr::pivot_wider(data.valid.removeC, id_cols = c(id, all_of(variables.vec)), names_from = outcome_char, values_from = c(time, status)) %>%
    dplyr::select(all_of(variables.vec))
  
  data.valid.for.pred <- model.matrix(as.formula(paste("~ ", paste(variables.vec, collapse = "+"), sep = "")), 
                                      data = data.valid.for.pred)[, -1]
  
  ### Now apply the function
  risk.joint.est <- apply(data.valid.for.pred, 1, myfunc.jointrisk.frailty, 
                          betas_A = frail.fit[["betas_A"]], betas_B = frail.fit[["betas_B"]], fit.in = frail.fit)
  
  return(list("rstan.model" = fit, "frail.fit" = frail.fit, "risk.joint.est" = as.numeric(risk.joint.est)))
  
}