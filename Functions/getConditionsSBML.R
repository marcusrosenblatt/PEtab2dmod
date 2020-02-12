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
getConditionsSBML <- function(conditions){
  condition.grid <- read.csv(file = conditions, sep = "\t")
  rownames(condition.grid) <- condition.grid$conditionId
  condition.grid$conditionId <- NULL
  return(condition.grid)
}
