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
  
  # Initialization
  reactions <- NULL
  events <- NULL
  compartments <- NULL
  
  # import compartments
  N_species <- m$getNumSpecies()
  compartments <- do.call(c, lapply(0:(N_species-1), function(i){ m$getSpecies(i)$getCompartment()}))
  names_compartments <- do.call(c, lapply(0:(N_species-1), function(i){ m$getSpecies(i)$getId() }))
  names(compartments) <- names_compartments
  
  if(unique(compartments)[1]=="default") compartments <- NULL
  
  # import reactions and adjust by means of compartments
  N_reactions <- m$getNumReactions()
  for (reaction in 0:(N_reactions-1)){
    Reactantstring <- ""
    Productstring <- ""
    eq <- m$getModel()$getReaction(reaction)
    Reactantnr <- eq$getNumReactants(reaction)
    if(Reactantnr > 0) Reactantstring <- paste0( eq$getReactant(0)$getStoichiometry(), "*", eq$getReactant(0)$getSpecies())
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
    if(FALSE){
      if(Reactantnr > 0) ReactantCompartment <- compartments$compartment[which(compartments$name==eq$getReactant(0)$getSpecies())]
        else ReactantCompartment <- "1"
      if(Productnr > 0) ProductCompartment <- compartments$compartment[which(compartments$name==eq$getProduct(0)$getSpecies())]
        else ProductCompartment <- "1"
      if(ReactantCompartment==ProductCompartment) 
        reactions <- reactions %>% addReaction(Reactantstring, Productstring, paste0("(",rate, ")/", ReactantCompartment))
      else if(ReactantCompartment=="1"){
        reactions <- reactions %>% addReaction(Reactantstring, Productstring, paste0("(",rate, ")/", ProductCompartment))
        } else if (ProductCompartment=="1"){
          reactions <- reactions %>% addReaction(Reactantstring, Productstring, paste0("(",rate, ")/", ReactantCompartment))
          } else reactions <- reactions %>% 
                  addReaction(Reactantstring, "", paste0("(",rate, ")/", ReactantCompartment)) %>% 
                  addReaction("", Productstring, paste0("(",rate, ")/", ProductCompartment))
    } else {
      if(!is.null(compartments)){
        pos <- which(strsplit(rate, "")[[1]]=="*")[1]
        rate <- substr(rate,pos+1,length(strsplit(rate, "")[[1]]=="*"))
      } 
      reactions <- reactions %>% addReaction(Reactantstring, Productstring, rate)
    }
  }
  # import inputs
  N_rules <- m$getNumRules()
  if (N_rules > 0){
    for (rule in 0:(N_rules-1)){
      reactions$rates <- replaceSymbols(m$getRule(rule)$getVariable(),
                                        paste0("(",m$getRule(rule)$getFormula(), ")"), reactions$rates)  # substitute m$getRule(0)$getVariable() by m$getRule(0)$getFormula()
    }
  }
  
  # replace function based inputs by events (done in reactions)
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
    reactions <- reactions %>% addReaction("", paste0("event", i), "0")
  }
  
  TransformEvents <- function(events){
    if(!is.null(events)){
      do.call(rbind, lapply(1:length(events), function(i){
        myevent <- events[i]
        if(str_detect(myevent, "piecewise") & (str_detect(myevent, "leq") | str_detect(myevent, "lt"))){
          timepoint <- strsplit(strsplit(myevent, ",")[[1]][3], ")")[[1]][1]
          first <- strsplit(strsplit(myevent, "\\(")[[1]][2], ",")[[1]][1]
          second <- strsplit(strsplit(myevent, ",")[[1]][4], ")")[[1]][1]
          if(!is.na(suppressWarnings(as.numeric(timepoint)))) timepoint <- as.numeric(timepoint) # avoid warning if variable is not numeric
          if(!is.na(suppressWarnings(as.numeric(first)))) first <- as.numeric(first)
          if(!is.na(suppressWarnings(as.numeric(second)))) second <- as.numeric(second)
          return(data.frame(var=paste0("event",i), time=c(0,timepoint), value=c(first, second), method="replace"))
        } else {print("Warning. Event not yet supported"); return(myevent)}
      }))
    } else return(NULL)
  }
  events <- TransformEvents(events)
  
  reactions$rates <- replaceSymbols(c("t", "TIME", "T"), "time", reactions$rates)
  
  mydata <- as.data.frame(reactions)
  reactions <- as.eqnlist(mydata, compartments)
  
  return(list(reactions=reactions, events=events))
}


