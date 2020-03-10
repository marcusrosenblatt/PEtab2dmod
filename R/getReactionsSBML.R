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
#' @export
#' 
getReactionsSBML <- function(model, conditions){
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
  
  # if(unique(compartments)[1]=="default") compartments <- NULL
  
  # import reactions and adjust by means of compartments
  N_reactions <- m$getNumReactions()
  for (reaction in 0:(N_reactions-1)){
    Reactantstring <- ""
    Productstring <- ""
    eq <- m$getReaction(reaction)
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
    formula <- eq$getKineticLaw()$getFormula()
    if(str_detect(formula, "Function"))
      rate <- formula # to be double checked  # works for Borghans now
    else 
      rate <- gsub("pow", "", gsub(", ", "**", formula))
    #rate <- replaceOperation("pow", "**", eq$getKineticLaw()$getFormula())
    if(!is.null(compartments)){
      if(Reduce("|", str_detect(rate, unique(compartments)))){
        pos <- which(strsplit(rate, "")[[1]]=="*")[1]
        rate <- substr(rate,pos+1,length(strsplit(rate, "")[[1]]=="*"))
      }
    }
    reactions <- reactions %>% addReaction(Reactantstring, Productstring, rate)
  }
  reactions$rates <- gsub(" ","",reactions$rates)
  # import functions
  N_fundefs <- m$getNumFunctionDefinitions()
  if (N_fundefs > 0){
    for (fun in 0:(N_fundefs-1)){
      mymath <- m$getFunctionDefinition(fun)$getMath()
      # print(fun)
      string <- function_def_to_string(m$getFunctionDefinition(fun)) %>% gsub(" ","",.)
      # print(string)
      first <- strsplit(string, "=")[[1]][1]
      second <- strsplit(string, "=")[[1]][2]
      # print(first)
      # print(second)
      first <- gsub("\\(", "\\\\\\(", first)
      second <- gsub("\\(", "\\\\\\(", second)
      reactions$rates <- gsub(first,
                              second, reactions$rates)  # substitute m$getRule(0)$getVariable() by m$getRule(0)$getFormula()
      #print(formulaToL3String(mymath$getChild(mymath$getNumChildren()-1)))
    }
  }
  
  # import inputs
  N_rules <- m$getNumRules()
  if (N_rules > 0){
    for (rule in 0:(N_rules-1)){
      # substitute m$getRule(0)$getVariable() by m$getRule(0)$getFormula()
      reactions$rates <- replaceSymbols(m$getRule(rule)$getVariable(),
                                        paste0("(",m$getRule(rule)$getFormula(), ")"), reactions$rates)  
    }
  }
  
  reactions$rates <- gsub(" ","",reactions$rates)
  
  # replace function based inputs by events (done in reactions)
  for(fun in c("piecewise")){
    for(reaction in reactions$rates){
      if(str_detect(reaction, fun)){
        split <- str_split(reaction, fun)[[1]][2]
        count_bracket <- 0
        done <- F
        for(z in 1:nchar(split)){
          if(substr(split, z, z)=="(") count_bracket <- count_bracket+1
          if(substr(split, z, z)==")") count_bracket <- count_bracket-1
          if(count_bracket==0 & !done) {done <- T; pos <- z}
        }
        #pos <- which(strsplit(split, "")[[1]]==")")[2]
        event <- paste0(fun, substr(split, 1, pos))
        events <- c(events, event)
      }
    }
  }
  events <- unique(events)
  if(!is.null(events)) for(i in 1:length(events)){
    replace <- gsub("\\(", "\\\\\\(", events[i])
    replace <- gsub("\\*", "\\\\\\*", replace)
    replace <- gsub("\\+", "\\\\\\+", replace)
    #replace <- gsub("\\)", "\\\\\\)", replace)
    reactions$rates <- gsub(replace, paste0("event", i), reactions$rates)
    reactions <- reactions %>% addReaction("", paste0("event", i), "0")
  }
  
  # replace mathematical expressions 
  reactions$rates <- replaceSymbols(c("t", "TIME", "T"), "time", reactions$rates)
  
  TransformEvents <- function(events){
    if(!is.null(events)){
      do.call(rbind, lapply(1:length(events), function(i){
        myevent <- events[i]
        if(str_detect(myevent, "piecewise") & (str_detect(myevent, "leq") | str_detect(myevent, "lt"))){
          expr1 <- strsplit(myevent, ",")[[1]][2]
          expr1 <- gsub(paste0(strsplit(expr1, "\\(")[[1]][1],"\\("), "", expr1)
          expr2 <- strsplit(strsplit(myevent, ",")[[1]][3], ")")[[1]][1]
          if(expr1=="time") timepoint <- expr2 else 
            if(str_detect(expr1, "time-")) timepoint <- gsub("time-", "", expr1) else cat("Warning: Event not yet supported.")
          first <- strsplit(strsplit(myevent, "\\(")[[1]][2], ",")[[1]][1]
          second <- strsplit(strsplit(myevent, ",")[[1]][4], ")")[[1]][1]
          if(!is.na(suppressWarnings(as.numeric(timepoint)))) timepoint <- as.numeric(timepoint) # avoid warning if variable is not numeric
          if(!is.na(suppressWarnings(as.numeric(first)))) first <- as.numeric(first)
          if(!is.na(suppressWarnings(as.numeric(second)))) second <- as.numeric(second)
          return(data.frame(var=paste0("event",i), time=c(0,timepoint), value=c(first, second), method="replace"))
        } else {cat("Warning: Event not yet supported"); return(myevent)}
      }))
    } else return(NULL)
  }
  events <- TransformEvents(events)
  
  
  ## check for preequilibration conditions and handle them via events
  preeqEvents <- NULL
  myconditions <- read.csv(file = conditions, sep = "\t")
  myCons <- myconditions$conditionId
  mypreeqCons <- NULL
  attrib <- NULL
  for (con in myCons){if(paste0("preeq_", con)%in%myCons) mypreeqCons <- c(mypreeqCons, con)}
  if(!is.null(mypreeqCons)){
    for (con in mypreeqCons){
      mycongrid <- filter(myconditions, conditionId==con | conditionId==paste0("preeq_", con))
      if(ncol(mycongrid)>1){
        for(i in 2:ncol(mycongrid)){
          preeqEvents <- addEvent(preeqEvents, var=names(mycongrid)[i], time=0, value=mycongrid[[which(mycongrid$conditionId==con),i]], method="replace")
          attrib <- c(attrib, mycongrid[[which(mycongrid$conditionId==paste0("preeq_",con)),i]])
        }
      }
    }
  }
  mystates <- reactions$states
  reactions_orig <- reactions
  attr(preeqEvents, "initials") <- attrib
  if(!is.null(preeqEvents)) for(i in 1:nrow(preeqEvents)){
    events <- rbind(events, preeqEvents[i,])
    reactions <- reactions %>% addReaction("", preeqEvents[[i,"var"]], "0")
  }
  

  # for(i in 1:length(reactions$rates)){
  #   reaction <- reactions$rates[i]
  #   if(str_detect(reaction, "pow")){
  #     reaction_new <- gsub("pow", "", reaction)
  #     # reaction_new <- gsub(", ", "**", reaction_new)
  #     reaction_new <- gsub(",", "**", reaction_new)
  #     reactions$rates[i] <- reaction_new
  #   } else reactions$rates[i] <- reaction
  # }
  
  mydata <- as.data.frame(reactions)
  reactions <- as.eqnlist(mydata, compartments)
  
  return(list(reactions=reactions, events=events, reactions_orig=reactions_orig, preeqEvents=preeqEvents, mystates=mystates))
}


## function from Frank
function_def_to_string <- function(fun)
{
  if (is.null(fun)) return;
  id <- fun$getId()
  
  math <- fun$getMath()
  if (is.null(math)) return;
  
  # the function will be of the form lambda(a, b, formula)
  # so the last child of the AST_Node is the actual math everything in front of it 
  # the arguments
  #
  # So to generate the desired output function_definition_id(args) = math
  # we use: 
  
  
  num_children <- math$getNumChildren()
  
  result <- paste(id, '(', sep="")
  for (i in 1:(num_children-1) ) # this leaves out the last one
  {
    result <- paste(result, libSBML::formulaToL3String(math$getChild(i-1)), sep="")
    
    if (i < (num_children-1))
      result <- paste(result, ', ', sep="")
    
  }
  
  result <- paste(result, ') = ', libSBML::formulaToL3String(math$getChild(num_children - 1)), sep="")
  
  return (result)
  
}

