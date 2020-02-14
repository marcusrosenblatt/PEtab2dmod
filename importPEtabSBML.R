#' Import an SBML model and corresponding PEtab objects 
#' 
#' @description This function imports an SBML model and corresponding PEtab objects, e.g. from the Benchmark collection.
#' Note: Objects such as model equations, parameters or data are automatically assigned to standard variables (reactions, observables, g, x, p0) and written to your current working directory (via <<-).
#' You can obtain your common variable names by the additional assign arguments.
#'  
#' @param parameters PEtab parameter file as .tsv
#' @param parameters PEtab parameter file as .tsv
#'   
#' @return NULL
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
importPEtabSBML <- function(modelname = "Boehm_JProteomeRes2014",
                       path2BC = "BenchmarkModels/",
                       SBML_file = NULL,
                       observable_file = NULL,
                       condition_file = NULL,
                       data_file = NULL,
                       parameter_file = NULL,
                       assign_reactions = NULL,
                       assign_observables = NULL,
                       assign_g = NULL,
                       assign_data = NULL,
                       assign_odemodel = NULL){
  starttime <- Sys.time()
  if(is.null(SBML_file))       SBML_file <- paste0(path2BC, modelname, "/model_", modelname, ".xml")
  if(is.null(observable_file)) observable_file <- paste0(path2BC, modelname, "/observables_", modelname, ".tsv")
  if(is.null(condition_file))  condition_file <- paste0(path2BC, modelname, "/experimentalCondition_", modelname, ".tsv")
  if(is.null(data_file))       data_file <- paste0(path2BC, modelname, "/measurementData_", modelname, ".tsv")
  if(is.null(parameter_file))  parameter_file <- paste0(path2BC, modelname, "/parameters_", modelname, ".tsv")
  mywd <- getwd()
  if(!file.exists(SBML_file)){cat(paste0("The file ",mywd,SBML_file, " does not exist. Please check spelling or provide the file name via the SBML_file argument.")); return(NULL)}
  if(!file.exists(observable_file)){cat(paste0("The file ",mywd,observable_file, " does not exist. Please check spelling or provide the file name via the observable_file argument.")); return(NULL)}
  if(!file.exists(condition_file)){cat(paste0("The file ",mywd,condition_file, " does not exist. Please check spelling or provide the file name via the condition_file argument.")); return(NULL)}
  if(!file.exists(data_file)){cat(paste0("The file ",mywd,data_file, " does not exist. Please check spelling or provide the file name via the data_file argument.")); return(NULL)}
  if(!file.exists(parameter_file)){cat(paste0("The file ",mywd,parameter_file, " does not exist. Please check spelling or provide the file name via the parameter_file argument.")); return(NULL)}
  cat("Reading SBML file ...\n")
  mylist <- getReactionsSBML(SBML_file)
  myreactions <- mylist$reactions
  myevents <- mylist$events
  if(is.null(assign_reactions)){reactions <<- myreactions} else {cat("Manual assignment not yet provided.")}
  cat("Reading observables ...\n")
  myobservables <<- getObservablesSBML(observable_file)
  if(is.null(assign_observables)){observables <<- myobservables} else {cat("Manual assignment not yet provided.")}
  cat("Compiling observable function ...\n")
  myg <- Y(myobservables, myreactions, compile=TRUE, modelname=paste0("g_",modelname))
  if(is.null(assign_g)){g <<- myg} else {cat("Manual assignment not yet provided.")}
  cat("Reading data file ...\n")
  mydataSBML <- getDataSBML(data_file, myobservables)
  mydata <- mydataSBML$data
  if(is.null(assign_data)){data <<- mydata} else {cat("Manual assignment not yet provided.")}
  cat("Compiling ODE model ...\n")
  myodemodel <- odemodel(reactions, forcings = NULL, events = myevents, fixed=NULL, modelname = paste0("odemodel_", modelname), jacobian = "inz.lsodes", compile = TRUE)
  if(is.null(assign_odemodel)){myodemodel <<- myodemodel} else {cat("Manual assignment not yet provided.")}
  cat("Reading parameters and initials ...\n")
  myparameters <- getParametersSBML(parameter_file)
  myconstraints <- myparameters$constraints
  myfit_values <- myparameters$pouter
  myinitialsSBML <- getInitialsSBML(SBML_file)
  mycompartments <- myinitialsSBML$compartments
  myinitials <- myinitialsSBML$initials
  cat("Check and compile error model ...\n")
  myerrors <- mydataSBML$errors
  if(!is.null(errors)){
    err <- Y(myerrors, f = c(as.eqnvec(myreactions), myobservables), states = names(myobservables), attach.input = FALSE, compile = F, modelname = "errfn")
  }
  cat("Generate parameter transformations ...\n")
  myinnerpars <- unique(c(getParameters(myodemodel), getSymbols(myobservables), getSymbols(myerrors)))
  names(myinnerpars) <- myinnerpars
  trafo <- as.eqnvec(myinnerpars, names = myinnerpars)
  trafo <- replaceSymbols(names(initials), initials, trafo)
  trafo <- replaceSymbols(names(compartments), compartments, trafo)
  trafo <- replaceSymbols(names(constraints), constraints, trafo)
  endtime <- Sys.time()
  cat(paste0(modelname, " imported in ",as.character(format(as.numeric(endtime-starttime), digits=3)), " seconds.\n"))
  return(modelname)
}