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
      col_par <- NULL
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs) 
        {
        obs_par <- subset(data_obs, simulationConditionId == condition)$observableParameters %>% unique() %>% as.character()
        col_par <- c(col_par, obs_par)
      } 
      col_name <- paste0("observableParameter1_",obs)
      condition.grid_obs[col_name] <- col_par
    }
    condition.grid <- suppressWarnings(inner_join(condition.grid_orig,condition.grid_obs, by = "conditionId"))
    # avoid warning if not all conditions are observed
  }
  
  # generate columns for noiseParameters
  if(!mydata$noiseParameters %>% is.numeric  & !Reduce("&",is.na(mydata$noiseParameters))) 
  {
    condition.grid_noise <- data.frame(conditionId = condis_obs)
    for (obs in observables) 
    {
      col_par <- NULL
      data_obs <- subset(mydata, observableId == obs)
      for (condition in condis_obs) 
      {
        noise_par <- subset(data_obs, simulationConditionId == condition)$noiseParameters %>% unique() %>% as.character()
        col_par <- c(col_par, noise_par)
      } 
      col_name <- paste0("noiseParameter1_",obs)
      condition.grid_noise[col_name] <- col_par
    }
    condition.grid <- suppressWarnings(inner_join(condition.grid_orig,condition.grid_noise, by = "conditionId")) 
    # avoid warning if not all conditions are observed
  }
  
  rownames(condition.grid) <- condition.grid$conditionId
  condition.grid$conditionId <- NULL
  
  # check if all conditions are observed
  if(nrow(condition.grid) < nrow(condition.grid_orig)) print("There exist non-observed conditions!")

  return(condition.grid)
}
