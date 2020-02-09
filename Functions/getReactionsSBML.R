#' Import reactions from SBML. 
#' 
#' @description This function imports reactions from SBML. Reactions are written to an eqnlist object. 
#' Assignment rules for input functions are substituted in the rate column. Time is introduced as an additionel state t.
#'  
#' @param model SBML file as .xml
#'   
#' @return Eqnlist of reactions.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
getReactionsSBML <- function(model){
  m = readSBML(model)$getModel()
  
  reactions <- NULL
  events <- NULL
  
  N_reactions <- m$getNumReactions()
  for (reaction in 0:(N_reactions-1)){
    Reactantstring <- ""
    Productstring <- ""
    eq <- m$getModel()$getReaction(reaction)
    Reactantnr <- eq$getNumReactants(reaction)
    if(Reactantnr > 0)Reactantstring <- paste0( eq$getReactant(0)$getStoichiometry(), "*", eq$getReactant(0)$getSpecies())
    if(Reactantnr > 1) for (s in 1:(Reactantnr-1)) {
      Reactantstring <- paste0(Reactantstring, " + ",
                               paste0(eq$getReactant(s)$getStoichiometry(), "*", eq$getReactant(s)$getSpecies()))
    }
    Productnr <- eq$getNumProducts(reaction)
    if(Productnr > 0) Productstring <- paste0( eq$getProduct(0)$getStoichiometry(), "*", eq$getProduct(0)$getSpecies())
    if(Productnr > 1) for (s in 1:(Productnr-1)) {
      Productstring <- paste0(Productstring, " + ",
                              paste0(eq$getProduct(s)$getStoichiometry(), "*", eq$getProduct(s)$getSpecies()))
    }
    rate <- sub("pow", "", sub(", ", "**", eq$getKineticLaw()$getFormula())) # to be double checked
    #rate <- replaceOperation("pow", "**", eq$getKineticLaw()$getFormula())
    reactions <- reactions %>% addReaction(Reactantstring, Productstring, rate) 
  }
  
  N_rules <- m$getNumRules()
  if (N_rules > 0){
    reactions <- reactions %>% addReaction("", "t", "1")
    reactions <- reactions %>% addReaction("", "time", "1")
    for (rule in 0:(N_rules-1)){
      reactions$rates <- replaceSymbols(m$getRule(rule)$getVariable(),
                                        paste0("(",m$getRule(rule)$getFormula(), ")"), reactions$rates)  # substitute m$getRule(0)$getVariable() by m$getRule(0)$getFormula()
    }
  }
  
  
  for(fun in c("piecewise")){
    for(reaction in reactions$rates){
      if(str_detect(reaction, fun)){
        split <- str_split(reaction, fun)[[1]][2]
        pos <- which(strsplit(split, "")[[1]]==")")[2]
        event <- paste0(fun, substr(split, 1, pos))
        events <- c(events, event)
      }
    }
  }
  events <- unique(events)
  if(!is.null(events)) for(i in 1:length(events)){
    replace <- gsub("\\(", "\\\\\\(", events[i])
    #replace <- gsub("\\)", "\\\\\\)", replace)
    reactions$rates <- gsub(replace, paste0("event", i), reactions$rates)
  }
  
  return(list(reactions=reactions, events=events))
}
