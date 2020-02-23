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
  parlower <- NULL
  parupper <- NULL
  if(nrow(estimated)>0){
    for(i in 1:length(estimated$parameterScale)) {
      parscale <- estimated$parameterScale[i]
      par <- estimated$parameterId[i] %>% as.character()
      value <- estimated$nominalValue[i]
      lowervalue <- estimated$lowerBound[i]
      uppervalue <- estimated$upperBound[i]
      if(parscale == "lin"){
        pouter <- c(pouter, value)
        parlower <- c(parlower, lowervalue)
        parupper <- c(parupper, uppervalue)
      } else if(parscale == "log10"){
        pouter <- c(pouter, log10(value))
        parlower <- c(parlower, log10(lowervalue))
        parupper <- c(parupper, log10(uppervalue))
      } else if(parscale == "log"){
        pouter <- c(pouter, log(value))
        parlower <- c(parlower, log(lowervalue))
        parupper <- c(parupper, log(uppervalue))
      } else paste("This type of parameterScale is not supported.")
      names(pouter)[i] <- par
      names(parlower)[i] <- par
      names(parupper)[i] <- par
    } 
    parscales <- estimated$parameterScale %>% as.character()
    pars <- estimated$parameterId %>% as.character()
    names(parscales) <- pars
    attr(pouter,"parscale") <- parscales
    attr(pouter,"lowerBound") <- parlower
    attr(pouter,"upperBound") <- parupper
  }
  return(list(constraints=constraints,pouter=pouter))
}
