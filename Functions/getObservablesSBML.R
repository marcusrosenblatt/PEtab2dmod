#' Import observables from PEtab. 
#' 
#' @description This function imports observables from the PEtab observable file.
#'  
#' @param observables PEtab observable file as .tsv
#'   
#' @return Eqnvec of observables.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   

getObservablesSBML <- function(observables){
  ## Load observables
  myobs <- read.csv(file = observables, sep = "\t") %>% as.data.frame()
  obsNames <- myobs$observableId %>% as.character()
  
  # # rename observables with _obs
  # obsNames <- paste0(obsNames,"_obs")
  
  obsFormula <- myobs$observableFormula %>% as.character()
  obsFormula[which(myobs$observableTransformation=="log")] <- paste0("log(", obsFormula[which(myobs$observableTransformation=="log")], ")")
  obsFormula[which(myobs$observableTransformation=="log10")] <- paste0("log10(", obsFormula[which(myobs$observableTransformation=="log10")], ")")
  names(obsFormula) <- obsNames
  observables <- obsFormula %>% as.eqnvec()
  return(observables)
}
