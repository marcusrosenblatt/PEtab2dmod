#' Calculate analytical steady states. 
#' 
#' @description This function follows the method published in [1]. The determined steady-state solution is tailored to parameter estimation. Please note that kinetic parameters might be fixed for solution of steady-state equations. Note that additional parameters might be introduced to ensure positivity of the solution.
#' @description The function calls a python script via rPython. Usage problems might occur when different python versions are used. The script was written and tested for python 2.7.12, sympy 0.7.6 and numpy 1.8.2.
#' @description Recently, users went into problems with RJSONIO when rPython was used. Unless a sound solution is available, please try to reinstall RJSONIO in these cases.
#' 
#' 
#' @param model Either name of the csv-file or the eqnlist of the model. If NULL, specify smatrix, states and rates by hand.
#' @param file Name of the file to which the steady-state equations are saved.
#'   Read this file with \code{\link{readRDS}}.
#' @param smatrix Numeric matrix, stiochiometry matrix of the system 
#' @param states Character vector, state vector of the system
#' @param rates Character vector, flux vector of the system
#' @param forcings Character vector with the names of the forcings
#' @param givenCQs Character vector with conserved quantities. Use the format c("A + pA = totA", "B + pB = totB"). If NULL, conserved quantities are automatically calculated.
#' @param neglect Character vector with names of states and parameters that must not be used for solving the steady-state equations
#' @param sparsifyLevel numeric, Upper bound for length of linear combinations used for simplifying the stoichiometric matrix
#' @param outputFormat Define the output format. By default "R" generating dMod 
#'   compatible output. To obtain an output appropriate for d2d [2] "M" must be 
#'   selected.
#'   
#' @return Character vector of steady-state equations.
#'   
#' @references [1]
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4863410/}
#' @references [2]
#' \url{https://github.com/Data2Dynamics/d2d}
#' 
#' @author Marcus Rosenblatt, \email{marcus.rosenblatt@@fdm.uni-freiburg.de}
#'   
#' @export
#' @importFrom utils write.table
#' @example inst/examples/steadystates.R
getReactions <- function(model){
  m = readSBML(model)$getModel()
  N_reactions <- m$getNumReactions()
  reactions <- NULL
  for (r in 0:(N_reactions-1)){
    # addReaction(
    eq <- m$getModel()$getReaction(r)
    Reactantnr <- eq$getNumReactants(r)
    Reactantstring <- paste0( eq$getReactant(0)$getStoichiometry(), "*", eq$getReactant(0)$getSpecies())
    if(Reactantnr > 1) for (s in 1:(Reactantnr-1)) {
      Reactantstring <- paste0(Reactantstring, " + ",
                               paste0(eq$getReactant(s)$getStoichiometry(), "*", eq$getReactant(s)$getSpecies()))
    }
    Productnr <- eq$getNumProducts(r)
    Productstring <- paste0( eq$getProduct(0)$getStoichiometry(), "*", eq$getProduct(0)$getSpecies())
    if(Productnr > 1) for (s in 1:(Productnr-1)) {
      Productstring <- paste0(Productstring, " + ",
                              paste0(eq$getProduct(s)$getStoichiometry(), "*", eq$getProduct(s)$getSpecies()))
    }
    rate <- sub("pow", "", sub(", ", "**", eq$getKineticLaw()$getFormula()))
    reactions <- reactions %>% addReaction(Reactantstring, Productstring, rate) 
  }
  return(reactions)
}
