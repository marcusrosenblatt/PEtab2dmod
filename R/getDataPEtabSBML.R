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
#' @export
#' 
getDataPEtabSBML <- function(data, observables){
  mydata <- read.csv(file = data, sep = "\t")
  myobs <- read.csv(file = observables, sep = "\t") %>% as.data.frame()
  obs <- myobs$observableId %>% as.character()
  errors <- NULL
  
  if(!is.null(mydata$noiseParameters) & !is.null(myobs$noiseFormula)) cat(red("Warning: errors specified in data and observable file.\n"))
  if(!mydata$noiseParameters %>% is.numeric) {
    # define errors
    errors <- myobs$noiseFormula %>% as.character()
    which_err <- c(1:length(obs))
    if(length(errors) != length(obs)) errors <- rep(errors,length(obs))
    names(errors) <- obs[which_err]
    errors <- as.eqnvec(errors)
    # set fixed sigmas 
    if(is.null(mydata$noiseParameters)){
      mydata$noiseParameters <- errors %>% as.numeric()
    } else mydata$noiseParameters <- NA
  }
  
  # # rename observables with _obs
  # obs <- mydata$observableId %>% as.character() %>% paste0("_obs")
  # mydata$observableId <- obs
  # if(!is.null(errors)) names(errors) <- paste0(names(errors), "_obs")
  mydata$simulationConditionId <- as.character(mydata$simulationConditionId)
  if(!is.null(mydata$observableParameters)){
    for (observable in unique(mydata$observableId)){
      for(condition in unique(mydata$simulationConditionId)){
        sub <- subset(mydata, simulationConditionId==condition & observableId==observable)
        if(nrow(sub) > 0){
          if(length(unique(sub$observableParameters)) > 1){
            index <- which(mydata$simulationConditionId==condition & mydata$observableId==observable)
            mydata$simulationConditionId[index] <- paste0(mydata$simulationConditionId[index], "_",mydata$observableParameters[index])
          }
        }
      }
    }
  }
  
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