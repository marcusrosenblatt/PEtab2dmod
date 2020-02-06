#' Import reactions from SBML. 
#' 
#' @description This function imports reactions from SBML.
#'  
#' @param model SBML file as .xml
#'   
#' @return Eqnlist of reactions.
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
getReactionsSBML <- function(model){
  m = readSBML(model)$getModel()
  N_reactions <- m$getNumReactions()
  reactions <- NULL
  for (r in 0:(N_reactions-1)){
    Reactantstring <- ""
    Productstring <- ""
    eq <- m$getModel()$getReaction(r)
    Reactantnr <- eq$getNumReactants(r)
    if(Reactantnr > 0)Reactantstring <- paste0( eq$getReactant(0)$getStoichiometry(), "*", eq$getReactant(0)$getSpecies())
    if(Reactantnr > 1) for (s in 1:(Reactantnr-1)) {
      Reactantstring <- paste0(Reactantstring, " + ",
                               paste0(eq$getReactant(s)$getStoichiometry(), "*", eq$getReactant(s)$getSpecies()))
    }
    Productnr <- eq$getNumProducts(r)
    if(Productnr > 0) Productstring <- paste0( eq$getProduct(0)$getStoichiometry(), "*", eq$getProduct(0)$getSpecies())
    if(Productnr > 1) for (s in 1:(Productnr-1)) {
      Productstring <- paste0(Productstring, " + ",
                              paste0(eq$getProduct(s)$getStoichiometry(), "*", eq$getProduct(s)$getSpecies()))
    }
    rate <- sub("pow", "", sub(", ", "**", eq$getKineticLaw()$getFormula()))
    reactions <- reactions %>% addReaction(Reactantstring, Productstring, rate) 
  }
  return(reactions)
}
