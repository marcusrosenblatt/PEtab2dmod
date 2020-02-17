#' Plot observables and states of an SBML/PEtab model
#' 
#' @description This function plots data and fits of an SBML/PEtab model after the import.
#' Note: Certain objects generated during the model import by importPEtabSBML have to be present in the global environment and are used as default variables if not specified differntly.
#'  
#' @param g observation function as obsfn
#' @param x prediction function as prdfn
#' @param p parameter function as parfn
#' @param mydata data as datalist
#' @param pars parameter as vector
#' @param times times as vector
#'   
#' @return NULL
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
plotPEtabSBML <- function(g1 = g,
                          x1 = x,
                          p1 = p0,
                          mydata1 = mydata,
                          pouter1 = pouter,
                          times1 = times){
  
  prediction <- (g1*x1*p1)(times1, pouter1)
  P1 <- plotCombined(prediction, mydata1) 
  print(P1)
  # plotPrediction(prediction)
  
  # return(modelname)
}