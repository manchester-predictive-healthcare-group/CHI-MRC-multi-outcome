### Deal with global variables note
utils::globalVariables(c("id", "from", "to", "Tstart", "Tstop", "status", "state.poly.fac", #variables from data inputted into calc_calib_blr
                         "state.poly", #variables from data inputted into calc_calib_mlr
                         "value", "pred", "obs", "obs.upper", "obs.lower", "line.group", "mapping")) #variables from data inputted into plot_calib_blr
