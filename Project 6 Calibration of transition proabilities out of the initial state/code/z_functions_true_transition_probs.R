### Storing these functions in a seperate workspace, so they don't clog up the other function file,
### as these are very long, and have lots of functions within them

###
###
### Function to calculate true transition probabilities for DGM1
###
###
calc.true.transition.probs.DGM1 <- function(u.eval, t.eval, x1.eval, x2.eval,
                                       shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                                       shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                                       shape15, scale15, #shape and scale for weibull baseline hazard for transition 1 -> 5
                                       shape24, scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                                       shape25, scale25, #shape and scale for weibull baseline hazard for transition 2 -> 5
                                       shape34, scale34, #shape and scale for weibull baseline hazard for transition 3 -> 4
                                       shape35, scale35, #shape and scale for weibull baseline hazard for transition 3 -> 5
                                       shape45, scale45, #shape and scale for weibull baseline hazard for transition 3 -> 5
                                       beta12.x1, beta12.x2, #covariate effects for transiion 12
                                       beta13.x1, beta13.x2, #covariate effects for transiion 13
                                       beta15.x1, beta15.x2, #covariate effects for transiion 15
                                       beta24.x1, beta24.x2, #covariate effects for transiion 24
                                       beta25.x1, beta25.x2, #covariate effects for transiion 25
                                       beta34.x1, beta34.x2, #covariate effects for transiion 34
                                       beta35.x1, beta35.x2, #covariate effects for transiion 35
                                       beta45.x1, beta45.x2 #covariate effects for transiion 45
                                       ){
  ###
  ### Start by defining the cause-specific hazards
  csh12 <- function(t.in, x1, x2){
    return(exp(beta12.x1*x1 + beta12.x2*x2)*(shape12/scale12)*((t.in)/scale12)^(shape12 - 1))
  }
  
  csh13 <- function(t.in, x1, x2){
    return(exp(beta13.x1*x1 + beta13.x2*x2)*(shape13/scale13)*((t.in)/scale13)^(shape13 - 1))
  }
  
  csh15 <- function(t.in, x1, x2){
    return(exp(beta15.x1*x1 + beta15.x2*x2)*(shape15/scale15)*((t.in)/scale15)^(shape15 - 1))
  }
  
  csh24 <- function(t.in, x1, x2){
    return(exp(beta24.x1*x1 + beta24.x2*x2)*(shape24/scale24)*((t.in)/scale24)^(shape24 - 1))
  }
  
  csh25 <- function(t.in, x1, x2){
    return(exp(beta25.x1*x1 + beta25.x2*x2)*(shape25/scale25)*((t.in)/scale25)^(shape25 - 1))
  }
  
  csh34 <- function(t.in, x1, x2){
    return(exp(beta34.x1*x1 + beta34.x2*x2)*(shape34/scale34)*((t.in)/scale34)^(shape34 - 1))
  }
  
  csh35 <- function(t.in, x1, x2){
    return(exp(beta35.x1*x1 + beta35.x2*x2)*(shape35/scale35)*((t.in)/scale35)^(shape35 - 1))
  }
  
  csh45 <- function(t.in, x1, x2){
    return(exp(beta45.x1*x1 + beta45.x2*x2)*(shape45/scale45)*((t.in)/scale45)^(shape45 - 1))
  }
  
  ###
  ### Define survival functions out of each state
  S1 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta12.x1*x1 + beta12.x2*x2)*(t.in/scale12)^shape12 +
      exp(beta13.x1*x1 + beta13.x2*x2)*(t.in/scale13)^shape13 +
      exp(beta15.x1*x1 + beta15.x2*x2)*(t.in/scale15)^shape15
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S2 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta24.x1*x1 + beta24.x2*x2)*(t.in/scale24)^shape24 +
      exp(beta25.x1*x1 + beta25.x2*x2)*(t.in/scale25)^shape25
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S3 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta34.x1*x1 + beta34.x2*x2)*(t.in/scale34)^shape34 +
      exp(beta35.x1*x1 + beta35.x2*x2)*(t.in/scale35)^shape35
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S4 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta45.x1*x1 + beta45.x2*x2)*(t.in/scale45)^shape45
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  
  ###
  ###
  ### Define transition probabilities
  ###
  ###
  
  ### In format P.ij.route (from state i, to state j, via route)
  
  ###
  ### Transitions out of state 4
  ###
  P.44 <- function(u.in, t.in, x1, x2){
    return(S4(t.in, x1, x2)/S4(u.in, x1, x2))
  }
  
  P.45 <- function(u.in, t.in, x1, x2){
    return(1 - S4(t.in, x1, x2)/S4(u.in, x1, x2))
  }
  
  
  ###
  ### Transitions out of state 3
  ###
  
  P.33 <- function(u.in, t.in, x1, x2){
    return(S3(t.in, x1, x2)/S3(u.in, x1, x2))
  }
  
  
  P.35.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh35(r.in, x1, x2)*S3(r.in, x1, x2)/S3(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.35.4 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh34(r.in, x1, x2)*(S3(r.in, x1, x2)/S3(u.in, x1, x2))*P.45(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.34 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh34(r.in, x1, x2)*(S3(r.in, x1, x2)/S3(u.in, x1, x2))*P.44(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  ###
  ### Transitions out of state 2
  ###
  P.22 <- function(u.in, t.in, x1, x2){
    return(S2(t.in, x1, x2)/S2(u.in, x1, x2))
  }
  
  P.25.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh25(r.in, x1, x2)*S2(r.in, x1, x2)/S2(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.25.4 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh24(r.in, x1, x2)*(S2(r.in, x1, x2)/S2(u.in, x1, x2))*P.45(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.24 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh24(r.in, x1, x2)*(S2(r.in, x1, x2)/S2(u.in, x1, x2))*P.44(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  ###
  ### Transitions out of state 1
  ###
  
  ###
  ### Staying in state 1
  P.11 <- function(u.in, t.in, x1, x2){
    return(S1(t.in, x1, x2)/S1(u.in, x1, x2))
  }
  
  ###
  ### Transitions to states 2 or 3 (and staying)
  P.12 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.22(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  P.13 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.33(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ###
  ### Transitions to state 4
  P.14.2 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.24(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  P.14.3 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.34(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Combine both for transition probability from 1 to 4
  P.14 <- function(u.in, t.in, x1, x2){
    return(P.14.2(u.in, t.in, x1, x2) + P.14.3(u.in, t.in, x1, x2))
  }
  
  ###
  ### Transitions to state 5
  
  ## Direct
  P.15.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh15(r.in, x1, x2)*S1(r.in, x1, x2)/S1(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 2 only
  P.15.2 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.25.direct(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 2 and 4
  P.15.24 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.25.4(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 3 only
  P.15.3 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.35.direct(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 3 and 4
  P.15.34 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.35.4(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Combine all for transition probability from 1 to 5
  P.15 <- function(u.in, t.in, x1, x2){
    return(P.15.direct(u.in, t.in, x1, x2) + 
             P.15.2(u.in, t.in, x1, x2) +
             P.15.24(u.in, t.in, x1, x2) + 
             P.15.3(u.in, t.in, x1, x2) +
             P.15.34(u.in, t.in, x1, x2))
  }
  
  ### Create output vector with transition probabilities
  output <- c("P.11" = P.11(u.eval, t.eval, x1.eval, x2.eval),
              "P.12" = P.12(u.eval, t.eval, x1.eval, x2.eval),
              "P.13" = P.13(u.eval, t.eval, x1.eval, x2.eval),
              "P.14" = P.14(u.eval, t.eval, x1.eval, x2.eval),
              "P.15" = P.15(u.eval, t.eval, x1.eval, x2.eval))
  return(output)
}

# 
# ### Storing these seperately so they don't clog up workspace, and make z_functions.R file difficult to navigate
# 
# ### Baseline hazards
# shape12 <- 1
# scale12 <- 1588.598
# 
# shape13 <- 1
# scale13 <- 0.5*1588.598
# 
# shape15 <- 1
# scale15 <- 5*1588.598
# 
# shape24 <- 1
# scale24 <- 1588.598
# 
# shape25 <- 1
# scale25 <- 5*1588.598
# 
# shape34 <- 1
# scale34 <- 0.5*1588.598
# 
# shape35 <- 1
# scale35 <- 5*1588.598
# 
# shape45 <- 1
# scale45 <- 5*1588.598
# 
# #qweibull(0.8, 1, 1588.598)
# 
# ## Covariate effects
# beta12.x1 <- 1
# beta12.x2 <- 1
# beta13.x1 <- 0.5
# beta13.x2 <- 0.5
# beta15.x1 <- 1
# beta15.x2 <- 0.5
# beta24.x1 <- 0.5
# beta24.x2 <- 1
# beta25.x1 <- 1
# beta25.x2 <- 1
# beta34.x1 <- 0.5
# beta34.x2 <- 0.5
# beta35.x1 <- 1
# beta35.x2 <- 0.5
# beta45.x1 <- 0.5
# beta45.x2 <- 1
# 
# ### Let's test this
# true.trans.probs <- calc.true.transition.probs.DGM1(u.eval = 0, t.eval = ceiling(7*365.25), x1.eval = 0, x2.eval = 0,
#                        shape12 = 1, scale12 = 1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 2
#                        shape13 = 1, scale13 = 0.5*1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 3
#                        shape15 = 1, scale15 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 5
#                        shape24 = 1, scale24 = 1588.598, #shape and scale for weibull baseline hazard for transition 2 -> 4
#                        shape25 = 1, scale25 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 2 -> 5
#                        shape34 = 1, scale34 = 0.5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 4
#                        shape35 = 1, scale35 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 5
#                        shape45 = 1, scale45 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 5
#                        beta12.x1 = 1, beta12.x2 = 1, #covariate effects for transiion 12
#                        beta13.x1 = 0.5, beta13.x2 = 0.5, #covariate effects for transiion 13
#                        beta15.x1 = 1, beta15.x2 = 0.5, #covariate effects for transiion 15
#                        beta24.x1 = 0.5, beta24.x2 = 1, #covariate effects for transiion 24
#                        beta25.x1 = 1, beta25.x2 = 1, #covariate effects for transiion 25
#                        beta34.x1 = 0.5, beta34.x2 = 0.5, #covariate effects for transiion 34
#                        beta35.x1 = 1, beta35.x2 = 0.5, #covariate effects for transiion 35
#                        beta45.x1 = 0.5, beta45.x2 = 1 #covariate effects for transiion 45
#                        )
# true.trans.probs
# 


###
###
### Function to calculate true transition probabilities for DGM2
###
###
calc.true.transition.probs.DGM2 <- function(u.eval, t.eval, x1.eval, x2.eval, out.num.states,
                                            shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                                            shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                                            shape16, scale16, #shape and scale for weibull baseline hazard for transition 1 -> 6
                                            shape24, scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                                            shape26, scale26, #shape and scale for weibull baseline hazard for transition 2 -> 6
                                            shape35, scale35, #shape and scale for weibull baseline hazard for transition 3 -> 5
                                            shape36, scale36, #shape and scale for weibull baseline hazard for transition 3 -> 6
                                            shape46, scale46, #shape and scale for weibull baseline hazard for transition 4 -> 6
                                            shape56, scale56, #shape and scale for weibull baseline hazard for transition 5 -> 6
                                            beta12.x1, beta12.x2, #covariate effects for transiion 1 -> 2
                                            beta13.x1, beta13.x2, #covariate effects for transiion 1 -> 3
                                            beta16.x1, beta16.x2, #covariate effects for transiion 1 -> 6
                                            beta24.x1, beta24.x2, #covariate effects for transiion 2 -> 4
                                            beta26.x1, beta26.x2, #covariate effects for transiion 2 -> 6
                                            beta35.x1, beta35.x2, #covariate effects for transiion 3 -> 5
                                            beta36.x1, beta36.x2, #covariate effects for transiion 3 -> 6
                                            beta46.x1, beta46.x2, #covariate effects for transiion 4 -> 6
                                            beta56.x1, beta56.x2 #covariate effects for transiion 5 -> 6
){
  ###
  ### Start by defining the cause-specific hazards
  csh12 <- function(t.in, x1, x2){
    return(exp(beta12.x1*x1 + beta12.x2*x2)*(shape12/scale12)*((t.in)/scale12)^(shape12 - 1))
  }
  
  csh13 <- function(t.in, x1, x2){
    return(exp(beta13.x1*x1 + beta13.x2*x2)*(shape13/scale13)*((t.in)/scale13)^(shape13 - 1))
  }
  
  csh16 <- function(t.in, x1, x2){
    return(exp(beta16.x1*x1 + beta16.x2*x2)*(shape16/scale16)*((t.in)/scale16)^(shape16 - 1))
  }
  
  csh24 <- function(t.in, x1, x2){
    return(exp(beta24.x1*x1 + beta24.x2*x2)*(shape24/scale24)*((t.in)/scale24)^(shape24 - 1))
  }
  
  csh26 <- function(t.in, x1, x2){
    return(exp(beta26.x1*x1 + beta26.x2*x2)*(shape26/scale26)*((t.in)/scale26)^(shape26 - 1))
  }
  
  csh35 <- function(t.in, x1, x2){
    return(exp(beta35.x1*x1 + beta35.x2*x2)*(shape35/scale35)*((t.in)/scale35)^(shape35 - 1))
  }
  
  csh36 <- function(t.in, x1, x2){
    return(exp(beta36.x1*x1 + beta36.x2*x2)*(shape36/scale36)*((t.in)/scale36)^(shape36 - 1))
  }
  
  csh46 <- function(t.in, x1, x2){
    return(exp(beta46.x1*x1 + beta46.x2*x2)*(shape46/scale46)*((t.in)/scale46)^(shape46 - 1))
  }
  
  csh56 <- function(t.in, x1, x2){
    return(exp(beta56.x1*x1 + beta56.x2*x2)*(shape56/scale56)*((t.in)/scale56)^(shape56 - 1))
  }
  
  ###
  ### Define survival functions out of each state
  S1 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta12.x1*x1 + beta12.x2*x2)*(t.in/scale12)^shape12 +
      exp(beta13.x1*x1 + beta13.x2*x2)*(t.in/scale13)^shape13 +
      exp(beta16.x1*x1 + beta16.x2*x2)*(t.in/scale16)^shape16
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S2 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta24.x1*x1 + beta24.x2*x2)*(t.in/scale24)^shape24 +
      exp(beta26.x1*x1 + beta26.x2*x2)*(t.in/scale26)^shape26
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S3 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta35.x1*x1 + beta35.x2*x2)*(t.in/scale35)^shape35 +
      exp(beta36.x1*x1 + beta36.x2*x2)*(t.in/scale36)^shape36
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S4 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta46.x1*x1 + beta46.x2*x2)*(t.in/scale46)^shape46
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S5 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta56.x1*x1 + beta56.x2*x2)*(t.in/scale56)^shape56
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  
  ###
  ###
  ### Define transition probabilities
  ###
  ###
  
  ### In format P.ij.route (from state i, to state j, via route)
  
  ###
  ### Transitions out of state 5
  ###
  P.55 <- function(u.in, t.in, x1, x2){
    return(S5(t.in, x1, x2)/S5(u.in, x1, x2))
  }
  
  P.56 <- function(u.in, t.in, x1, x2){
    return(1 - S5(t.in, x1, x2)/S5(u.in, x1, x2))
  }
  
  ###
  ### Transitions out of state 4
  ###
  P.44 <- function(u.in, t.in, x1, x2){
    return(S4(t.in, x1, x2)/S4(u.in, x1, x2))
  }
  
  P.46 <- function(u.in, t.in, x1, x2){
    return(1 - S4(t.in, x1, x2)/S4(u.in, x1, x2))
  }
  
  
  
  ###
  ### Transitions out of state 3
  ###
  
  P.33 <- function(u.in, t.in, x1, x2){
    return(S3(t.in, x1, x2)/S3(u.in, x1, x2))
  }
  
  
  P.36.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh36(r.in, x1, x2)*S3(r.in, x1, x2)/S3(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.36.5 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh35(r.in, x1, x2)*(S3(r.in, x1, x2)/S3(u.in, x1, x2))*P.56(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.35 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh35(r.in, x1, x2)*(S3(r.in, x1, x2)/S3(u.in, x1, x2))*P.55(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  ###
  ### Transitions out of state 2
  ###
  P.22 <- function(u.in, t.in, x1, x2){
    return(S2(t.in, x1, x2)/S2(u.in, x1, x2))
  }
  
  P.26.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh26(r.in, x1, x2)*S2(r.in, x1, x2)/S2(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.26.4 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh24(r.in, x1, x2)*(S2(r.in, x1, x2)/S2(u.in, x1, x2))*P.46(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.24 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh24(r.in, x1, x2)*(S2(r.in, x1, x2)/S2(u.in, x1, x2))*P.44(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  ###
  ### Transitions out of state 1
  ###
  
  ###
  ### Staying in state 1
  P.11 <- function(u.in, t.in, x1, x2){
    return(S1(t.in, x1, x2)/S1(u.in, x1, x2))
  }
  
  ###
  ### Transitions to states 2 or 3 (and staying)
  P.12 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.22(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  P.13 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.33(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ###
  ### Transition to state 4
  P.14 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.24(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ###
  ### Transition to state 5
  P.15 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.35(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }

  
  ###
  ### Transitions to state 6
  
  ## Direct
  P.16.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh16(r.in, x1, x2)*S1(r.in, x1, x2)/S1(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 2 only
  P.16.2 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.26.direct(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 2 and 4
  P.16.24 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.26.4(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 3 only
  P.16.3 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.36.direct(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 3 and 5
  P.16.35 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.36.5(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Combine all for transition probability from 1 to 5
  P.16 <- function(u.in, t.in, x1, x2){
    return(P.16.direct(u.in, t.in, x1, x2) + 
             P.16.2(u.in, t.in, x1, x2) +
             P.16.24(u.in, t.in, x1, x2) + 
             P.16.3(u.in, t.in, x1, x2) +
             P.16.35(u.in, t.in, x1, x2))
  }
  
  ### Create output vector with transition probabilities. Have two options, for whether we want the true probs of the underling states,
  ### or the probability of being in the states according to the 5 state model
  if (out.num.states == 5){
    output <- c("P.11" = P.11(u.eval, t.eval, x1.eval, x2.eval),
                "P.12" = P.12(u.eval, t.eval, x1.eval, x2.eval),
                "P.13" = P.13(u.eval, t.eval, x1.eval, x2.eval),
                "P.14" = P.14(u.eval, t.eval, x1.eval, x2.eval) + P.15(u.eval, t.eval, x1.eval, x2.eval),
                "P.15" = P.16(u.eval, t.eval, x1.eval, x2.eval))
  } else if (out.num.states == 6){
    output <- c("P.11" = P.11(u.eval, t.eval, x1.eval, x2.eval),
                "P.12" = P.12(u.eval, t.eval, x1.eval, x2.eval),
                "P.13" = P.13(u.eval, t.eval, x1.eval, x2.eval),
                "P.14" = P.14(u.eval, t.eval, x1.eval, x2.eval),
                "P.15" = P.15(u.eval, t.eval, x1.eval, x2.eval),
                "P.16" = P.16(u.eval, t.eval, x1.eval, x2.eval))
  }

  return(output)
}


# ### Storing these seperately so they don't clog up workspace, and make z_functions.R file difficult to navigate
# u.eval <- 0
# t.eval <- ceiling(7*365.25)
# x1.eval <- 0
# x2.eval <- 0
# ### Baseline hazards
# shape12 <- 1
# scale12 <- 1588.598
# 
# shape13 <- 1
# scale13 <- 0.5*1588.598
# 
# shape16 <- 1
# scale16 <- 5*1588.598
# 
# shape24 <- 1
# scale24 <- 1588.598
# 
# shape26 <- 1
# scale26 <- 5*1588.598
# 
# shape35 <- 1
# scale35 <- 0.5*1588.598
# 
# shape36 <- 1
# scale36 <- 5*1588.598
# 
# shape46 <- 1
# scale46 <- 5*1588.598
# 
# shape56 <- 1
# scale56 <- 5*1588.598
# 
# #qweibull(0.8, 1, 1588.598)
# 
# ## Covariate effects
# beta12.x1 <- 1
# beta12.x2 <- 1
# beta13.x1 <- 0.5
# beta13.x2 <- 0.5
# beta16.x1 <- 1
# beta16.x2 <- 0.5
# beta24.x1 <- 0.5
# beta24.x2 <- 1
# beta26.x1 <- 1
# beta26.x2 <- 1
# beta35.x1 <- 0.5
# beta35.x2 <- 0.5
# beta36.x1 <- 1
# beta36.x2 <- 0.5
# beta46.x1 <- 0.5
# beta46.x2 <- 1
# beta56.x1 <- 0.5
# beta56.x2 <- 1