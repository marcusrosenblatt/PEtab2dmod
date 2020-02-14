#' Import Parameters from PEtab 
#' 
#' @description This function imports fixed and fitted parameters from the PEtab parameter file as named vectors.
#'  
#' @param parameters PEtab parameter file as .tsv
#'   
#' @return constraints and pouter as list of named vectros.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
getParametersSBML <- function(parameters){
  mypars <- read.csv(file = parameters, sep = "\t")
  fixed <- mypars %>% filter(estimate == 0)
  constraints <- NULL
  if(nrow(fixed)>0){
    for(i in 1:length(fixed$parameterScale)) {
      parscale <- fixed$parameterScale[i]
      par <- fixed$parameterId[i] %>% as.character()
      value <- fixed$nominalValue[i]
      if(parscale == "lin") constraints <- c(constraints, value)
      if(parscale == "log10") constraints <- c(constraints, log10(value))
      if(parscale == "log") constraints <- c(constraints, log(value))
      else paste("This type of parameterScale is not supported.")
      names(constraints)[i] <- par
    } 
    parscales <- fixed$parameterScale %>% as.character()
    pars <- fixed$parameterId %>% as.character()
    names(parscales) <- pars
    attr(constraints,"parscale") <- parscales
  }
  estimated <- mypars %>% filter(estimate == 1)
  pouter <- NULL
  if(nrow(estimated)>0){
    for(i in 1:length(estimated$parameterScale)) {
      parscale <- estimated$parameterScale[i]
      par <- estimated$parameterId[i] %>% as.character()
      value <- estimated$nominalValue[i]
      if(parscale == "lin") pouter <- c(pouter, value)
      if(parscale == "log10") pouter <- c(pouter, log10(value))
      if(parscale == "log") pouter <- c(pouter, log(value))
      else paste("This type of parameterScale is not supported.")
      names(pouter)[i] <- par
    } 
    parscales <- estimated$parameterScale %>% as.character()
    pars <- estimated$parameterId %>% as.character()
    names(parscales) <- pars
    attr(pouter,"parscale") <- parscales
  }
  return(list(constraints=constraints,pouter=pouter))
}
