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
#' @param ... further arguments going to \link{plotCombined}
#'   
#' @return NULL
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
plotPEtabSBML <- function(g1 = g,
                          x1 = x,
                          p1 = p0,
                          mydata1 = mydata,
                          pouter1 = pouter,
                          times1 = times,
                          errfn = err,
                          ...){
  if(!is.null(errfn)){
    prediction <- as.data.frame((g1*x1*p1)(times1, pouter1), errfn = err)
    data <- as.data.frame(mydata1)
    P1 <- ggplot() + geom_line(data=prediction, aes(x=time, y=value, color=condition)) + 
      geom_ribbon(data=prediction, aes(x=time, ymin=value-sigma, ymax=value+sigma, fill=condition), color=NA, alpha=0.25) +
      geom_point(data=data, aes(x=time, y=value, color=condition)) +
      theme_dMod() + scale_color_dMod() + scale_fill_dMod() +
      facet_wrap(~name, scales="free")
  } else {
    prediction <- (g1*x1*p1)(times1, pouter1)
    P1 <- plotCombined(prediction, mydata1, ...) 
  }
  print(P1)
  # plotPrediction(prediction)
  
  # return(modelname)
}