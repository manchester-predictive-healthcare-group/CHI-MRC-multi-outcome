##########################################################################
### Program with functions to calculate 'true' transition proabilities ###
### Storing these functions in a seperate workspace, so they don't     ###
### clog up the other function file. This file follows process in      ###
### supplementary material part 1 of manuscript:                       ###
### "calibration plots for multistate models"                          ###
##########################################################################

###
### Function to calculate true transition probabilities for DGM1
###
calc.true.transition.probs.DGM1 <- function(u.eval, t.eval, x12.eval, x13.eval, x15.eval, x24.eval, x25.eval, x34.eval, x35.eval, x45.eval,
                                       shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                                       shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                                       shape15, scale15, #shape and scale for weibull baseline hazard for transition 1 -> 5
                                       shape24, scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                                       shape25, scale25, #shape and scale for weibull baseline hazard for transition 2 -> 5
                                       shape34, scale34, #shape and scale for weibull baseline hazard for transition 3 -> 4
                                       shape35, scale35, #shape and scale for weibull baseline hazard for transition 3 -> 5
                                       shape45, scale45, #shape and scale for weibull baseline hazard for transition 3 -> 5
                                       beta.x12,  #covariate effects for transiion 12
                                       beta.x13,  #covariate effects for transiion 13
                                       beta.x15,  #covariate effects for transiion 15
                                       beta.x24,  #covariate effects for transiion 24
                                       beta.x25,  #covariate effects for transiion 25
                                       beta.x34,  #covariate effects for transiion 34
                                       beta.x35,  #covariate effects for transiion 35
                                       beta.x45  #covariate effects for transiion 45
                                       ){
  ###
  ### Start by defining the cause-specific hazards
  csh12 <- function(t.in, x12){
    return(exp(beta.x12*x12)*(shape12/scale12)*((t.in)/scale12)^(shape12 - 1))
  }
  
  csh13 <- function(t.in, x13){
    return(exp(beta.x13*x13)*(shape13/scale13)*((t.in)/scale13)^(shape13 - 1))
  }
  
  csh15 <- function(t.in, x15){
    return(exp(beta.x15*x15)*(shape15/scale15)*((t.in)/scale15)^(shape15 - 1))
  }
  
  csh24 <- function(t.in, x24){
    return(exp(beta.x24*x24)*(shape24/scale24)*((t.in)/scale24)^(shape24 - 1))
  }
  
  csh25 <- function(t.in, x25){
    return(exp(beta.x25*x25)*(shape25/scale25)*((t.in)/scale25)^(shape25 - 1))
  }
  
  csh34 <- function(t.in, x34){
    return(exp(beta.x34*x34)*(shape34/scale34)*((t.in)/scale34)^(shape34 - 1))
  }
  
  csh35 <- function(t.in, x35){
    return(exp(beta.x35*x35)*(shape35/scale35)*((t.in)/scale35)^(shape35 - 1))
  }
  
  csh45 <- function(t.in, x45){
    return(exp(beta.x45*x45)*(shape45/scale45)*((t.in)/scale45)^(shape45 - 1))
  }
  
  ###
  ### Define survival functions out of each state
  S1 <- function(t.in, x12, x13, x15){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta.x12*x12)*(t.in/scale12)^shape12 +
      exp(beta.x13*x13)*(t.in/scale13)^shape13 +
      exp(beta.x15*x15)*(t.in/scale15)^shape15
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S2 <- function(t.in, x24, x25){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta.x24*x24)*(t.in/scale24)^shape24 +
      exp(beta.x25*x25)*(t.in/scale25)^shape25
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S3 <- function(t.in, x34, x35){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta.x34*x34)*(t.in/scale34)^shape34 +
      exp(beta.x35*x35)*(t.in/scale35)^shape35
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S4 <- function(t.in, x45){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta.x45*x45)*(t.in/scale45)^shape45
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
  P.44 <- function(u.in, t.in, x45){
    return(S4(t.in, x45)/S4(u.in, x45))
  }
  
  P.45 <- function(u.in, t.in, x45){
    return(1 - S4(t.in, x45)/S4(u.in, x45))
  }
  
  
  ###
  ### Transitions out of state 3
  ###
  
  P.33 <- function(u.in, t.in, x34, x35){
    return(S3(t.in, x34, x35)/S3(u.in, x34, x35))
  }
  
  
  P.35.direct <- function(u.in, t.in, x34, x35){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh35(r.in, x35)*S3(r.in, x34, x35)/S3(u.in, x34, x35))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.35.4 <- function(u.in, t.in, x34, x35, x45){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh34(r.in, x34)*(S3(r.in, x34, x35)/S3(u.in, x34, x35))*P.45(r.in, t.in, x45))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.34 <- function(u.in, t.in, x34, x35, x45){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh34(r.in, x34)*(S3(r.in, x34, x35)/S3(u.in, x34, x35))*P.44(r.in, t.in, x45))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  ###
  ### Transitions out of state 2
  ###
  P.22 <- function(u.in, t.in, x24, x25){
    return(S2(t.in, x24, x25)/S2(u.in, x24, x25))
  }
  
  P.25.direct <- function(u.in, t.in, x24, x25){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh25(r.in, x25)*S2(r.in, x24, x25)/S2(u.in, x24, x25))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.25.4 <- function(u.in, t.in, x24, x25, x45){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh24(r.in, x24)*(S2(r.in, x24, x25)/S2(u.in, x24, x25))*P.45(r.in, t.in, x45))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  P.24 <- function(u.in, t.in, x24, x25, x45){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh24(r.in, x24)*(S2(r.in, x24, x25)/S2(u.in, x24, x25))*P.44(r.in, t.in, x45))
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
  P.11 <- function(u.in, t.in, x12, x13, x15){
    return(S1(t.in, x12, x13, x15)/S1(u.in, x12, x13, x15))
  }
  
  ###
  ### Transitions to states 2 or 3 (and staying)
  P.12 <- function(u.in, t.in, x12, x13, x15, x24, x25){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x12)*(S1(r.in, x12, x13, x15)/S1(u.in, x12, x13, x15))*P.22(r.in, t.in, x24, x25))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  P.13 <- function(u.in, t.in, x12, x13, x15, x34, x35){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x13)*(S1(r.in, x12, x13, x15)/S1(u.in, x12, x13, x15))*P.33(r.in, t.in, x34, x35))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ###
  ### Transitions to state 4
  P.14.2 <- function(u.in, t.in, x12, x13, x15, x24, x25, x45){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x12)*(S1(r.in, x12, x13, x15)/S1(u.in, x12, x13, x15))*P.24(r.in, t.in, x24, x25, x45))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  P.14.3 <- function(u.in, t.in, x12, x13, x15, x34, x35, x45){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x13)*(S1(r.in, x12, x13, x15)/S1(u.in, x12, x13, x15))*P.34(r.in, t.in, x34, x35, x45))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Combine both for transition probability from 1 to 4
  P.14 <- function(u.in, t.in, x12, x13, x15, x24, x25, x34, x35, x45){
    return(P.14.2(u.in, t.in, x12, x13, x15, x24, x25, x45) + P.14.3(u.in, t.in, x12, x13, x15, x34, x35, x45))
  }
  
  ###
  ### Transitions to state 5
  
  ## Direct
  P.15.direct <- function(u.in, t.in, x12, x13, x15){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh15(r.in, x15)*S1(r.in, x12, x13, x15)/S1(u.in, x12, x13, x15))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 2 only
  P.15.2 <- function(u.in, t.in, x12, x13, x15, x24, x25){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x12)*(S1(r.in, x12, x13, x15)/S1(u.in, x12, x13, x15))*P.25.direct(r.in, t.in, x24, x25))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 2 and 4
  P.15.24 <- function(u.in, t.in, x12, x13, x15, x24, x25, x45){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x12)*(S1(r.in, x12, x13, x15)/S1(u.in, x12, x13, x15))*P.25.4(r.in, t.in, x24, x25, x45))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 3 only
  P.15.3 <- function(u.in, t.in, x12, x13, x15, x34, x35){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x13)*(S1(r.in, x12, x13, x15)/S1(u.in, x12, x13, x15))*P.35.direct(r.in, t.in, x34, x35))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 3 and 4
  P.15.34 <- function(u.in, t.in, x12, x13, x15, x34, x35, x45){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x13)*(S1(r.in, x12, x13, x15)/S1(u.in, x12, x13, x15))*P.35.4(r.in, t.in, x34, x35, x45))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Combine all for transition probability from 1 to 5
  P.15 <- function(u.in, t.in, x12, x13, x15, x24, x25, x34, x35, x45){
    return(P.15.direct(u.in, t.in, x12, x13, x15) + 
             P.15.2(u.in, t.in, x12, x13, x15, x24, x25) +
             P.15.24(u.in, t.in, x12, x13, x15, x24, x25, x45) + 
             P.15.3(u.in, t.in, x12, x13, x15, x34, x35) +
             P.15.34(u.in, t.in, x12, x13, x15, x34, x35, x45))
  }
  
  ### Create output vector with transition probabilities
  output <- c("P.11" = P.11(u.eval, t.eval, x12.eval, x13.eval, x15.eval),
              "P.12" = P.12(u.eval, t.eval, x12.eval, x13.eval, x15.eval, x24.eval, x25.eval),
              "P.13" = P.13(u.eval, t.eval, x12.eval, x13.eval, x15.eval, x34.eval, x35.eval),
              "P.14" = P.14(u.eval, t.eval, x12.eval, x13.eval, x15.eval, x24.eval, x25.eval, x34.eval, x35.eval, x45.eval),
              "P.15" = P.15(u.eval, t.eval, x12.eval, x13.eval, x15.eval, x24.eval, x25.eval, x34.eval, x35.eval, x45.eval))
  return(output)
}


