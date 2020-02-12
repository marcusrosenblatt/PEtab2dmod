#' Import initials from SBML. 
#' 
#' @description This function imports initial values or equations describing the same from SBML and writes them in a named vector. 
#'  
#' @param model SBML file as .xml
#'   
#' @return Named vector of initials.
#'   
#' @author Marcus Rosenblatt, Svenja Kemmer and Frank Bergmann
#'   
getInitialsSBML <- function(model){
  
  model = readSBML(model)$getModel()
  initials <- NULL
  species <- NULL
  
  # now we can go through all species
  for ( i in 0:(model$getNumSpecies()-1) ){
    current <- model$getSpecies(i)
    
    # now the species can have several cases that determine
    # their initial value 
    
    # it could be that the species is fully determined by an assignment rule
    # (that apply at all times), so we have to check rules first
    rule <- model$getRule(current$getId())
    if (!is.null(rule))
    {
      # ok there is a rule for this species so lets figure out what its type
      # is as that determines whether it applies at t0
      rule_type <- rule$getTypeCode()
      type_name <- libSBML::SBMLTypeCode_toString(rule_type, 'core')
      if (type_name == "AssignmentRule")
      {
        # the initial value is determined by the formula
        math <- rule$getMath()
        if (!is.null(math))
        {
          formula <- libSBML::formulaToL3String(math)
          # print(paste('Species: ', current$getId(), ' is determined at all times by formula: ', formula))
          initials <- c(initials,formula)
          species <- c(species,current$getId())
          
          # no need to look at other values so continue for another one
          next
        }
      }
      
      if (type_name == "RateRule")
      {
        math <- rule$getMath()
        if (!is.null(math))
        {
          formula <- libSBML::formulaToL3String(math)
          # print(paste('Species: ', current$getId(), ' has an ode rule with formula: ', formula))
          initials <- c(initials,formula)
          species <- c(species,current$getId())
          
          # even though there is an ODE attached to the species, its initial value is needed
        }
      }
    }
    
    
    # it could have an initial assignment
    ia <- model$getInitialAssignment(current$getId())
    if (!is.null(ia))
    {
      math <- ia$getMath()
      if (!is.null(math))
      {
        formula <- libSBML::formulaToL3String(math)
        # print(paste("Species: ", current$getId(), " has an initial assignment with formula: ", formula))
        initials <- c(initials,formula)
        species <- c(species,current$getId())
        
        # as soon as you have that formula, no initial concentration / amount applies
        # so we don't have to look at anything else for this species
        next
      }
    }
    
    
    # it could have an initial amount
    if (current$isSetInitialAmount())
    {
      # print (paste("Species: ", current$getId(), "has initial amount: ", current$getInitialAmount()))
      initials <- c(initials,current$getInitialAmount())
      species <- c(species,current$getId())
    }
    
    # it could have an initial concentration
    if (current$isSetInitialConcentration())
    {
      # print (paste("Species: ", current$getId(), "has initial concentration: ", current$getInitialConcentration()))
      initials <- c(initials,current$getInitialConcentration())
      species <- c(species,current$getId())
    }
  }
  
  names(initials) <- species
  
  
  ## extract compartments

  # check if compartments exist
  if(model$getNumCompartments()>0)
    {

    #initialize vectors
    comp_name <- NULL
    comp_size <- NULL
    
    for ( i in 0:(model$getNumCompartments()-1) )
      {
      # get compartment name and size
      which <- model$getCompartment(i)$getId()
      size <- model$getCompartment(i)$getSize()
      
      comp_name <- c(comp_name,which)
      comp_size <- c(comp_size,size)
      
    }
    compartments <- comp_size
    names(compartments) <- comp_name
  }
  
  return(list(initials = initials, compartments = compartments))
}


