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
  
  # generate columns for observableParameters
  if(!is.numeric(mydata$observableParameters) & !Reduce("&",is.na(mydata$observableParameters))) 
    {
    condition.grid_obs <- data.frame(conditionId = condis_obs)
    for (obs in observables) 
      {
      col_pars <- NULL
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs) 
        {
        if(condition %in% data_obs$simulationConditionId){
          obs_par <- subset(data_obs, simulationConditionId == condition)$observableParameters %>% unique() %>% as.character()
          # one or more observable parameters?
          if(str_detect(obs_par,";")){
            myobspars <- strsplit(obs_par,";")[[1]]
            for(i in 1:length(myobspars)) {
              col_pars <- c(col_pars, myobspars[i])
            }
          } else col_pars <- c(col_pars, obs_par)
        }
      } 
      for (par in 1:length(col_pars)) {
        col_name <- paste0("observableParameter",par,"_",obs)
        condition.grid_obs[col_name] <- col_pars[par]
      }
    }
    mycondition.grid <- suppressWarnings(inner_join(condition.grid_orig,condition.grid_obs, by = "conditionId"))
    # avoid warning if not all conditions are observed
  }
  
  # generate columns for noiseParameters
  if(!is.numeric(mydata$noiseParameters)  & !Reduce("&",is.na(mydata$noiseParameters))) 
  {
    if(exists("mycondition.grid")) {condition.grid_orig <- mycondition.grid}
    condition.grid_noise <- data.frame(conditionId = condis_obs)
    for (obs in observables) 
    {
      col_pars <- NULL
      #print(obs)
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs) 
      {
        #print(condition)
        if(condition %in% data_obs$simulationConditionId){
          noise_par <- subset(data_obs, simulationConditionId == condition)$noiseParameters %>% unique() %>% as.character()
          
          # one or more noise parameters?
          if(str_detect(noise_par,";")){
            myobspars <- strsplit(noise_par,";")[[1]]
            for(i in 1:length(myobspars)) {
              col_pars <- c(col_pars, myobspars[i])
            }
          } else col_pars <- c(col_pars, noise_par)
        }
      } 
      for (par in 1:length(col_pars)) {
        col_name <- paste0("noiseParameter",par,"_",obs)
        condition.grid_noise[col_name] <- col_pars[par]
      }
    }
    mycondition.grid <- suppressWarnings(inner_join(condition.grid_orig,condition.grid_noise, by = "conditionId")) 
    # avoid warning if not all conditions are observed
  }
  
  if(!exists("mycondition.grid")) mycondition.grid <- condition.grid_orig
  rownames(mycondition.grid) <- mycondition.grid$conditionId
  mycondition.grid$conditionId <- NULL
  
  # check if all conditions are observed
  if(nrow(mycondition.grid) < nrow(condition.grid_orig)) print("There exist non-observed conditions!")

  return(mycondition.grid)
}
