#' Import Data from PEtab 
#' 
#' @description This function imports data from the PEtab data file as a data list and defines errors if an error model is required.
#'  
#' @param data PEtab data file as .tsv
#' @param observables observables as eqnvec
#'   
#' @return data as data list and errors (if required) as a named vector.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
getDataSBML <- function(data, observables){
  mydata <- read.csv(file = data, sep = "\t")
  errors <- NULL
  if(!mydata$noiseParameters %>% is.numeric) {
    # define errors
    obs <- observables %>% names()
    errors <- mydata$noiseParameters %>% levels
    which_err <- c(1:length(obs))
    if(length(errors) != length(obs)) errors <- rep(errors,length(obs))
    names(errors) <- obs[which_err]
    # set fixed sigmas to NA
    mydata$noiseParameters <- NA
  }
  #rename observables with _obs
  obs <- mydata$observableId %>% as.character()
  obs <- paste0(obs,"_obs")
  mydata$observableId <- obs
  # select necessary data columns
  data <- data.frame(name = mydata$observableId, time = mydata$time, 
                     value = mydata$measurement, sigma = mydata$noiseParameters,
                     condition = mydata$simulationConditionId) %>% as.datalist()
  return(list(data=data,errors=errors))
}
