#' Import condition.grid from PEtab 
#' 
#' @description This function imports the experimental conditions from the PEtab condition file as a gondition.grid.
#'  
#' @param conditions PEtab condition file as .tsv
#'   
#' @return condition.grid as data frame.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
getConditionsSBML <- function(conditions,data){
  condition.grid_orig <- read.csv(file = conditions, sep = "\t")
  mydata <- read.csv(file = data, sep = "\t")
  
  # check which conditions are observed
  condis_obs <- mydata$simulationConditionId %>% unique()
  # check which observables exist
  observables <- mydata$observableId %>% unique()
  
  # replace "" by NA
  if(!is.null(mydata$observableParameters)){
    mydata$observableParameters <- mydata$observableParameters %>% as.character()
    mydata <- mydata %>% mutate(observableParameters = ifelse(observableParameters == "",NA,observableParameters))
  }
  if(!is.null(mydata$noiseParameters)){
  mydata$noiseParameters <- mydata$noiseParameters %>% as.character()
  mydata <- mydata %>% mutate(noiseParameters = ifelse(noiseParameters == "",NA,noiseParameters))
  }
  # generate columns for observableParameters
  if(!is.numeric(mydata$observableParameters) & !is.null(mydata$observableParameters)) 
    {
    condition.grid_obs <- data.frame(conditionId = condis_obs)
    for (obs in observables) 
      {
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs) 
        {
        if(condition %in% data_obs$simulationConditionId){
          row_pars <- NULL
          obs_par <- subset(data_obs, simulationConditionId == condition)$observableParameters %>% unique() %>% as.character()
          if(!is.na(obs_par)){
            # one or more observable parameters?
            if(str_detect(obs_par,";")){
              myobspars <- strsplit(obs_par,";")[[1]]
              for(i in 1:length(myobspars)) {
                row_pars <- c(row_pars, myobspars[i])
              }
            } else row_pars <- c(row_pars, obs_par)
          }
          if(!is.null(row_pars)) for (par in 1:length(row_pars)) {
            col_name <- paste0("observableParameter",par,"_",obs)
            condition.grid_obs[which(condition.grid_obs$conditionId==condition),col_name] <- row_pars[par]
          }
        }
      } 
      
    }
    mycondition.grid <- suppressWarnings(inner_join(condition.grid_orig,condition.grid_obs, by = "conditionId"))
    # avoid warning if not all conditions are observed
  }
  
  # generate columns for noiseParameters
  if(!is.numeric(mydata$noiseParameters) & !is.null(mydata$noiseParameters)) 
  {
    if(exists("mycondition.grid")) {condition.grid_orig <- mycondition.grid}
    condition.grid_noise <- data.frame(conditionId = condis_obs)
    for (obs in observables) 
    {
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs) 
      {
        if(condition %in% data_obs$simulationConditionId){
          row_pars <- NULL
          noise_par <- subset(data_obs, simulationConditionId == condition)$noiseParameters %>% unique() %>% as.character()
          if(!is.na(noise_par)){
            # one or more observable parameters?
            if(str_detect(noise_par,";")){
              myobspars <- strsplit(noise_par,";")[[1]]
              for(i in 1:length(myobspars)) {
                row_pars <- c(row_pars, myobspars[i])
              }
            } else row_pars <- c(row_pars, noise_par)
          }
          if(!is.null(row_pars)) for (par in 1:length(row_pars)) {
            col_name <- paste0("noiseParameter",par,"_",obs)
            condition.grid_noise[which(condition.grid_noise$conditionId==condition),col_name] <- row_pars[par]
          }
        }
      } 
      
    }
    mycondition.grid <- suppressWarnings(inner_join(condition.grid_orig,condition.grid_noise, by = "conditionId"))
    # avoid warning if not all conditions are observed
  }
  
  if(!exists("mycondition.grid")) mycondition.grid <- condition.grid_orig
  rownames(mycondition.grid) <- mycondition.grid$conditionId
  # mycondition.grid$conditionId <- NULL ## we need this column in cases with just one condition!
  
  # check if all conditions are observed
  if(nrow(mycondition.grid) < nrow(condition.grid_orig)) print("There exist non-observed conditions!")
  
  for(i in 1:nrow(mycondition.grid)){
    for(j in 1:ncol(mycondition.grid)){
      if(is.na(mycondition.grid[i,j])) mycondition.grid[i,j] <- "1"
    }
  }

  return(mycondition.grid)
}
