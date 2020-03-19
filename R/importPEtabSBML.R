#' Import an SBML model and corresponding PEtab objects 
#' 
#' @description This function imports an SBML model and corresponding PEtab files, e.g. from the Benchmark collection.
#'  
#' @param modelname name of folder containing all PEtab files of the model to be imported. NULL if file paths are defined separately (see below).
#' @param path2model path to model folder
#' @param TestCases TRUE to load feature test cases
#' @param path2TestCases path to feature test case folder
#' @param compile if FALSE, g, ODEmodel and err are loaded from .RData (if present) and compilation time is saved
#' @param SBML_file SBML model as .xml
#' @param observable_file PEtab observable file as .tsv
#' @param condition_file PEtab condition file as .tsv
#' @param data_file PEtab data file as .tsv
#' @param parameter_file PEtab parameter file as .tsv
#'
#' @details Objects such as model equations, parameters or data are automatically assigned to the following standard variables and written to your current working directory (via <<-):
#' reactions, observables, errors, g, x, p0, err, obj, mydata, ODEmodel, condition.grid, trafoL, pouter, times.
#' Compiled objects (g, ODEmodel and err) are saved in .RData.
#' 
#' @return name of imported model
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#' 
#' @export
#' 
importPEtabSBML <- function(modelname = "Boehm_JProteomeRes2014",
                            path2model = "BenchmarkModels/",
                            testCases = FALSE,
                            path2TestCases = "PEtabTests/",
                            compile = TRUE,
                            SBML_file = NULL,
                            observable_file = NULL,
                            condition_file = NULL,
                            data_file = NULL,
                            parameter_file = NULL
                            )
  {
  
  ## Define path to SBML and PEtab files --------------------
  
  starttime <- Sys.time()
  if(testCases == FALSE){
    if(is.null(SBML_file))       SBML_file <- paste0(path2model, modelname, "/model_", modelname, ".xml")
    if(is.null(observable_file)) observable_file <- paste0(path2model, modelname, "/observables_", modelname, ".tsv")
    if(is.null(condition_file))  condition_file <- paste0(path2model, modelname, "/experimentalCondition_", modelname, ".tsv")
    if(is.null(data_file))       data_file <- paste0(path2model, modelname, "/measurementData_", modelname, ".tsv")
    if(is.null(parameter_file))  parameter_file <- paste0(path2model, modelname, "/parameters_", modelname, ".tsv")
  } else{
    SBML_file <- paste0(path2TestCases, modelname, "/_model.xml")
    observable_file <- paste0(path2TestCases, modelname, "/_observables.tsv")
    condition_file <- paste0(path2TestCases, modelname, "/_conditions.tsv")
    data_file <- paste0(path2TestCases, modelname, "/_measurements.tsv")
    parameter_file <- paste0(path2TestCases, modelname, "/_parameters.tsv")
  }
  mywd <- getwd()
  if(!file.exists(SBML_file)){cat(paste0("The file ",mywd,SBML_file, " does not exist. Please check spelling or provide the file name via the SBML_file argument.")); return(NULL)}
  if(!file.exists(observable_file)){cat(paste0("The file ",mywd,observable_file, " does not exist. Please check spelling or provide the file name via the observable_file argument.")); return(NULL)}
  if(!file.exists(condition_file)){cat(paste0("The file ",mywd,condition_file, " does not exist. Please check spelling or provide the file name via the condition_file argument.")); return(NULL)}
  if(!file.exists(data_file)){cat(paste0("The file ",mywd,data_file, " does not exist. Please check spelling or provide the file name via the data_file argument.")); return(NULL)}
  if(!file.exists(parameter_file)){cat(paste0("The file ",mywd,parameter_file, " does not exist. Please check spelling or provide the file name via the parameter_file argument.")); return(NULL)}
  
  if(is.null(modelname)) modelname <- "mymodel"
  ## Load shared objects --------------------
  
  dir.create(paste0(mywd,"/CompiledObjects/"), showWarnings = FALSE)
  setwd(paste0(mywd,"/CompiledObjects/"))
  files_loaded <- FALSE
  if(compile == FALSE & file.exists(paste0(modelname,".RData"))){
    load(paste0(modelname,".RData"))
    files_loaded <- TRUE
  } 
  setwd(mywd)
  
  ## Model Definition - Equations --------------------
  
  cat("Reading SBML file ...\n")
  mylist <- getReactionsSBML(SBML_file, condition_file)
  myreactions <- mylist$reactions
  myreactions_orig <- mylist$reactions_orig
  myevents <- mylist$events
  mypreeqEvents <- mylist$preeqEvents
  mystates <- mylist$mystates
  reactions <<- myreactions
  
  
  ## Model Definition - Observables --------------------
  
  cat("Reading observables ...\n")
  myobservables <- getObservablesSBML(observable_file)
  observables <<- myobservables
  
  
  cat("Compiling observable function ...\n")
  if(!files_loaded) {
    setwd(paste0(mywd,"/CompiledObjects/"))
    myg <- Y(myobservables, myreactions, compile=TRUE, modelname=paste0("g_",modelname))
    setwd(mywd)
  }
  g <<- myg
  
  
  ## Get Data ------------
  
  cat("Reading data file ...\n")
  mydataSBML <- getDataPEtabSBML(data_file, observable_file)
  mydata <- mydataSBML$data
  mydata <<- mydata
  
  
  ## Model Generation ---------------------
  
  cat("Compiling ODE model ...\n")
  
  if(!files_loaded) {
    setwd(paste0(mywd,"/CompiledObjects/"))
    myodemodel <- odemodel(myreactions, forcings = NULL, events = myevents, fixed=NULL, modelname = paste0("odemodel_", modelname), jacobian = "inz.lsodes", compile = TRUE)
    setwd(mywd)
  }
  ODEmodel <<- myodemodel
  
  
  ## Check and define error model ------------ 
  
  cat("Check and compile error model ...\n")
  myerrors <- mydataSBML$errors
  errors <<- myerrors
  
  myerr <- NULL
  if(!files_loaded) {
    if(!is_empty(getSymbols(myerrors))){
      setwd(paste0(mywd,"/CompiledObjects/"))
      myerr <- Y(myerrors, f = c(as.eqnvec(myreactions), myobservables), states = names(myobservables), attach.input = FALSE, compile = TRUE, modelname = paste0("errfn_", modelname))
      setwd(mywd)
    }
  }
  err <<- myerr
  
  ## Define constraints, initials, parameters and compartments --------------
  
  cat("Reading parameters and initials ...\n")
  myparameters <- getParametersSBML(parameter_file, SBML_file)
  myconstraints <- myparameters$constraints
  SBMLfixedpars <- myparameters$SBMLfixedpars
  myfit_values <- myparameters$pouter
  myinitialsSBML <- getInitialsSBML(SBML_file, condition_file)
  mycompartments <- myinitialsSBML$compartments
  myinitials <- myinitialsSBML$initials
  
  ## Parameter transformations -----------
  
  # Generate condition.grid
  grid <- getConditionsSBML(conditions = condition_file, data = data_file) 
  mypreeqCons <- grid$preeqCons
  mycondition.grid <- grid$condition_grid
  
  if(!is.null(SBMLfixedpars)){
    for (i in 1:length(SBMLfixedpars)) {
      if(!names(SBMLfixedpars)[i] %in% names(mycondition.grid))  mycondition.grid[names(SBMLfixedpars)[i]] <- SBMLfixedpars[i]
    } 
  }
  condi_pars <- names(mycondition.grid)[!names(mycondition.grid) %in% c("conditionName","conditionId")]
  condition.grid <<- mycondition.grid
  
  cat("Generate parameter transformations ...\n")
  myinnerpars <- unique(c(getParameters(myodemodel), getParameters(myg), getSymbols(myerrors)))
  names(myinnerpars) <- myinnerpars
  trafo <- as.eqnvec(myinnerpars, names = myinnerpars)
  trafo <- replaceSymbols(names(mycompartments), mycompartments, trafo)
  # only overwrite intial if it's not defined in condition.grid
  for (i in 1:length(myinitials)) {
    if(!names(myinitials)[i] %in% names(mycondition.grid)){
      trafo <- replaceSymbols(names(myinitials)[i], myinitials[i], trafo)
    }
  }
  trafo <- replaceSymbols(names(myconstraints), myconstraints, trafo)
  
  # branch trafo for different conditions
  mytrafoL <- branch(trafo, table=mycondition.grid)
  # set preequilibration event initials to corresponding values
  if(!is.null(mypreeqEvents)){
    mypreeqEvents2replace <- filter(mypreeqEvents, !var%in%mystates)
    mytrafoL <- repar("x~y", mytrafoL , x = unique(mypreeqEvents2replace$var), y = attr(mypreeqEvents2replace, "initials"))
  }
  # set remaining event initial to 0
  mytrafoL <- repar("x~0", mytrafoL , x = setdiff(unique(myevents$var), unique(mypreeqEvents$var)))  
  
  # condition-specific assignment of parameters from condition grid
  if(length(condi_pars) > 0){
    for (j in 1:length(names(mytrafoL))) {
      for (i in 1:length(condi_pars)) {
        mytrafoL[[j]] <- repar(x~y, mytrafoL[[j]], x=condi_pars[i], y=mycondition.grid[j,condi_pars[i]])
      }
    }
  }
  
  
  # transform parameters according to scale defined in the parameter PEtab file
  parscales <- attr(myfit_values,"parscale")
  mynames <- names(parscales)
  for(i in 1:length(parscales)){
    par <- parscales[i]
    par[par=="lin"] <- ""
    par[par=="log10"] <- "10**"
    par[par=="log"] <- "exp"
    parameter <- mynames[i]
    mytrafoL <- repar(paste0("x~",par,"(x)"), mytrafoL, x = parameter)
  }
  trafoL <<- mytrafoL
  
  
  ## Specify prediction functions ------
  
  # Get numeric steady state for preequilibration conditions
  if(!is.null(mypreeqCons)){
    myf <- as.eqnvec(myreactions_orig)[mystates]
    cq <- conservedQuantities(myreactions_orig$smatrix)
    if(!is.null(cq)){
      for(i in 1:nrow(cq)){
        myf[getSymbols(cq)[1]] <- paste0(as.character(conservedQuantities(myreactions_orig$smatrix)[1,]),"-1")
      }
    }
    setwd(paste0(mywd,"/CompiledObjects/"))
    pSS <- P(myf, condition = "c0", method = "implicit", compile = TRUE, modelname = paste0("preeq_", modelname))
    setwd(mywd)
  } else pSS <- NULL
  
  cat("Generate prediction function ...\n")
  tolerances <- 1e-7
  myp0 <- myx <- NULL
  for (C in names(mytrafoL)) {
    if(C%in%mypreeqCons){
      myp0 <- myp0 + pSS*P(mytrafoL[[C]], condition = C)
    } else {
      myp0 <- myp0 + P(mytrafoL[[C]], condition = C)
    }
    myx <- myx + Xs(myodemodel, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
                    optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                    condition = C)
  }
  
  p0 <<- myp0
  x <<- myx
  
  ## Generate objective function and initial parameter set -------
  
  myouterpars <- getParameters(myp0)
  mypouter <- structure(rep(NA,length(myouterpars)), names = myouterpars)
  
  common <- intersect(names(mypouter),names(myfit_values))
  mypouter[common] <- myfit_values[common]
  attr(mypouter, "parscales") <- parscales[which(names(myfit_values)%in%names(mypouter))]
  attr(mypouter, "lowerBound") <- attr(myfit_values,"lowerBound")[which(names(myfit_values)%in%names(mypouter))]
  attr(mypouter, "upperBound") <- attr(myfit_values,"upperBound")[which(names(myfit_values)%in%names(mypouter))]
  pouter <<- mypouter
  
  
  ## Define objective function -------
  
  cat("Generate objective function ...\n")
  if(!is.null(myerrors)){
    myobj <- normL2(mydata, myg*myx*myp0, errmodel = myerr) #+ constraintL2(prior, sigma=16)
  } else myobj <- normL2(mydata, myg*myx*myp0)
  obj <<- myobj
  
  mytimes <- seq(0,max(do.call(c, lapply(1:length(mydata), function(i) max(mydata[[i]]$time)))), len=501)
  times <<- mytimes
  
  if(!files_loaded){
    setwd(paste0(mywd,"/CompiledObjects/"))
    save(list = c("myg","myodemodel","myerr"),file = paste0(modelname,".RData"))
    setwd(mywd)
  }
  
  model_name <<- modelname
  
  endtime <- Sys.time()
  mytimediff <- as.numeric(difftime(endtime, starttime, unit="secs"))
  if(mytimediff > 3600) cat(paste0(modelname, " imported in ",as.character(format(as.numeric(difftime(endtime, starttime, unit="hours")), digits=3)), " hours.\n")) else
    if(mytimediff > 60) cat(paste0(modelname, " imported in ",as.character(format(as.numeric(difftime(endtime, starttime, unit="mins")), digits=3)), " minutes.\n")) else
      cat(paste0(modelname, " imported in ",as.character(format(as.numeric(difftime(endtime, starttime, unit="secs")), digits=3)), " seconds.\n"))
  return(modelname)
}