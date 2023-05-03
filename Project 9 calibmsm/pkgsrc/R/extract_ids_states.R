###
### Function to extract all individuals in state j, at time t.eval, from a dataset in mstate format
### Used in other functions which assess calibration using blr/mlr at specific time points (3.4 and 3.5)
###
extract_ids_states <- function(data.mstate, tmat, j, t.eval){

  ### Define maximum state number
  max.state <- max(data.mstate$to)

  ### Identify which states are absorbing states
  absorbing.states <- which(apply(tmat, 1, function(x) {sum(!is.na(x))}) == 0)

  ### For non-absorbing states, to be in state j at time t, you must have an observations from state j, where Tstart <= t.eval < Tstop
  if (!(j %in% absorbing.states)){
    ## Extract ids
    ids.state.j <- base::subset(data.mstate, from == j & Tstart <= t.eval & t.eval < Tstop) %>%
      dplyr::select(id) %>%
      dplyr::distinct(id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$id)
  } else if (j %in% absorbing.states){
    ### For absorbing state, just have to have moved into it
    ids.state.j <- base::subset(data.mstate, to == j & t.eval >= Tstop & status == 1) %>%
      dplyr::select(id) %>%
      dplyr::distinct(id)
    ## Put into numeric vector
    ids.state.j <- as.numeric(ids.state.j$id)
  }

  return(ids.state.j)
}
