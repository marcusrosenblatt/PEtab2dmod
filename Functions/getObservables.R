getObervables <- function(model){
  ## Load observables
  myobs <- read.csv(file = obs_file, sep = "\t") %>% as.data.frame()
  obsNames <- myobs$observableId %>% as.character()
  obsFormula <- myobs$observableFormula %>% as.character()
  names(obsFormula) <- obsNames
  observables <- obsFormula %>% as.eqnvec()
  return(observables)
}
