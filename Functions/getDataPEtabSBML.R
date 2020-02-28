#' Import Data from PEtab 
#' 
#' @description This function imports data from the PEtab data file as a data list and defines errors if an error model is required.
#'  
#' @param data PEtab data file as .tsv
#' @param observables observables as eqnvec
#'   
#' @return data as data list and errors (if required) as eqnvec.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
getDataPEtabSBML <- function(data, observables){
  mydata <- read.csv(file = data, sep = "\t")
  myobs <- read.csv(file = observables, sep = "\t") %>% as.data.frame()
  obs <- myobs$observableId %>% as.character()
  errors <- NULL
  if(!mydata$noiseParameters %>% is.numeric) {
    # define errors
    errors <- myobs$noiseFormula %>% as.character()
    which_err <- c(1:length(obs))
    if(length(errors) != length(obs)) errors <- rep(errors,length(obs))
    names(errors) <- obs[which_err]
    errors <- as.eqnvec(errors)
    # set fixed sigmas to NA
    mydata$noiseParameters <- NA
  }
  
  # # rename observables with _obs
  # obs <- mydata$observableId %>% as.character() %>% paste0("_obs")
  # mydata$observableId <- obs
  # if(!is.null(errors)) names(errors) <- paste0(names(errors), "_obs")
  
  # select necessary data columns
  data <- data.frame(name = mydata$observableId, time = mydata$time, 
                     value = mydata$measurement, sigma = mydata$noiseParameters,
                     condition = mydata$simulationConditionId) 
  obs2log <- myobs$observableId[which(myobs$observableTransformation=="log")]
  data$value[which(data$name%in%obs2log)] <- log(data$value[which(data$name%in%obs2log)])
  obs2log10 <- myobs$observableId[which(myobs$observableTransformation=="log10")]
  data$value[which(data$name%in%obs2log10)] <- log10(data$value[which(data$name%in%obs2log10)])
  data <- data %>% as.datalist()
  
  return(list(data=data,errors=errors))
}