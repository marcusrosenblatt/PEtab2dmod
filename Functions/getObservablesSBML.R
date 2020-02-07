#' Import observables from PEtab. 
#' 
#' @description This function imports observables from the PEtab observable file.
#'  
#' @param model PEtab observable file as .tsv
#'   
#' @return Eqnvec of observables.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   

getObservablesSBML <- function(observables){
  ## Load observables
  myobs <- read.csv(file = observables, sep = "\t") %>% as.data.frame()
  obsNames <- myobs$observableId %>% as.character()
  obsFormula <- myobs$observableFormula %>% as.character()
  names(obsFormula) <- obsNames
  observables <- obsFormula %>% as.eqnvec()
  return(observables)
}
