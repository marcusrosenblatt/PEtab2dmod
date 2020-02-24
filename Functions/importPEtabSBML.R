#' Import an SBML model and corresponding PEtab objects 
#' 
#' @description This function imports an SBML model and corresponding PEtab objects, e.g. from the Benchmark collection.
#' Note: Objects such as model equations, parameters or data are automatically assigned to standard variables (reactions, observables, g, x, p0) and written to your current working directory (via <<-).
#' You can obtain your common variable names by the additional assign arguments.
#'  
#' @param parameters PEtab parameter file as .tsv
#' @param parameters PEtab parameter file as .tsv
#'   
#' @return name of imported model
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
importPEtabSBML <- function(modelname = "Boehm_JProteomeRes2014",
                            path2BC = "BenchmarkModels/",
                            compile = FALSE,
                            SBML_file = NULL,
                            observable_file = NULL,
                            condition_file = NULL,
                            data_file = NULL,
                            parameter_file = NULL,
                            assign_reactions = NULL,
                            assign_observables = NULL,
                            assign_g = NULL,
                            assign_data = NULL,
                            assign_ODEmodel = NULL,
                            assign_condition.grid = NULL,
                            assign_trafo = NULL,
                            assign_p = NULL,
                            assign_x = NULL,
                            assign_pouter = NULL,
                            assign_obj = NULL,
                            assign_times = NULL,
                            assign_err = NULL,
                            assign_errors = NULL){
  
  ## Define path to SBML and PEtab files --------------------
  
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
  mylist <- getReactionsSBML(SBML_file)
  myreactions <- mylist$reactions
  myevents <- mylist$events
  if(is.null(assign_reactions)){reactions <<- myreactions} else {cat("Manual assignment not yet provided.")}
  
  
  ## Model Definition - Observables --------------------
  
  cat("Reading observables ...\n")
  myobservables <- getObservablesSBML(observable_file)
  if(is.null(assign_observables)){observables <<- myobservables} else {cat("Manual assignment not yet provided.")}
  

  cat("Compiling observable function ...\n")
  if(!files_loaded) {
    setwd(paste0(mywd,"/CompiledObjects/"))
    myg <- Y(myobservables, myreactions, compile=TRUE, modelname=paste0("g_",modelname))
    setwd(mywd)
  }
  if(is.null(assign_g)){g <<- myg} else {cat("Manual assignment not yet provided.")}
  
  
  ## Get Data ------------
  
  cat("Reading data file ...\n")
  mydataSBML <- getDataPEtabSBML(data_file, observable_file)
  mydata <- mydataSBML$data
  if(is.null(assign_data)){mydata <<- mydata} else {cat("Manual assignment not yet provided.")}
  
  
  ## Model Generation ---------------------
  
  cat("Compiling ODE model ...\n")
  
  if(!files_loaded) {
    setwd(paste0(mywd,"/CompiledObjects/"))
    myodemodel <- odemodel(myreactions, forcings = NULL, events = myevents, fixed=NULL, modelname = paste0("odemodel_", modelname), jacobian = "inz.lsodes", compile = TRUE)
    setwd(mywd)
  }
  if(is.null(assign_ODEmodel)){ODEmodel <<- myodemodel} else {cat("Manual assignment not yet provided.")}

  
  ## Define constraints, initials, parameters and compartments --------------
  
  cat("Reading parameters and initials ...\n")
  myparameters <- getParametersSBML(parameter_file)
  myconstraints <- myparameters$constraints
  myfit_values <- myparameters$pouter
  myinitialsSBML <- getInitialsSBML(SBML_file)
  mycompartments <- myinitialsSBML$compartments
  myinitials <- myinitialsSBML$initials
  
  
  ## Check and define error model ------------ 
  
  cat("Check and compile error model ...\n")
  myerrors <- mydataSBML$errors
  if(is.null(assign_errors)){errors <<- myerrors} else {cat("Manual assignment not yet provided.")}
  
  myerr <- NULL
  if(!files_loaded) {
    if(!is.null(myerrors)){
      setwd(paste0(mywd,"/CompiledObjects/"))
      myerr <- Y(myerrors, f = c(as.eqnvec(myreactions), myobservables), states = names(myobservables), attach.input = FALSE, compile = TRUE, modelname = paste0("errfn_", modelname))
      setwd(mywd)
    }
    if(is.null(assign_err)){err <<- myerr} else {cat("Manual assignment not yet provided.")}
  }
  
  
  ## Parameter transformations -----------
  
  cat("Generate parameter transformations ...\n")
  myinnerpars <- unique(c(getParameters(myodemodel), getSymbols(myobservables), setdiff(getSymbols(myerrors),names(observables))))
  names(myinnerpars) <- myinnerpars
  trafo <- as.eqnvec(myinnerpars, names = myinnerpars)
  trafo <- replaceSymbols(names(myinitials), myinitials, trafo)
  trafo <- replaceSymbols(names(mycompartments), mycompartments, trafo)
  trafo <- replaceSymbols(names(myconstraints), myconstraints, trafo)
  
  # Generate condition.grid
  mycondition.grid <- getConditionsSBML(conditions = condition_file, data = data_file) 
  condi_pars <- names(mycondition.grid)[!names(mycondition.grid) %in% c("conditionName")]
  if(is.null(assign_condition.grid)){condition.grid <<- mycondition.grid} else {cat("Manual assignment not yet provided.")}
  
  # branch trafo for different conditions
  # set event initial to 0
  mytrafoL <- branch(trafo, table=mycondition.grid)
  mytrafoL <- repar("x~0", mytrafoL , x = unique(myevents$var))  
  
  # condition-specific assignment of parameters from condition grid
  if(length(condi_pars) > 0){
    for (j in 1:length(names(mytrafoL))) {
      for (i in 1:length(condi_pars)) {
        mytrafoL[[j]] <- repar(x~y, mytrafoL[[j]], x=condi_pars[i], y=mycondition.grid[j,i+1])
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
  if(is.null(assign_trafo)){trafoL <<- mytrafoL} else {cat("Manual assignment not yet provided.")}
  
  
  ## Specify prediction functions ------
  
  cat("Generate prediction function ...\n")
  tolerances <- 1e-7
  myp0 <- myx <- NULL
  for (C in names(trafoL)) {
    myp0 <- myp0 + P(trafoL[[C]], condition = C)
    myx <- myx + Xs(myodemodel, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
                optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                condition = C)
  }
  if(is.null(assign_p)){p0 <<- myp0} else {cat("Manual assignment not yet provided.")}
  if(is.null(assign_x)){x <<- myx} else {cat("Manual assignment not yet provided.")}
  
  ## Generate objective function and initial parameter set -------
  
  myouterpars <- getParameters(myp0)
  mypouter <- structure(rep(NA,length(myouterpars)), names = myouterpars)
  
  common <- intersect(names(mypouter),names(myfit_values))
  mypouter[common] <- myfit_values[common]
  attr(mypouter, "parscales") <- parscales[which(names(myfit_values)%in%names(mypouter))]
  attr(mypouter, "lowerBound") <- attr(myfit_values,"lowerBound")[which(names(myfit_values)%in%names(mypouter))]
  attr(mypouter, "upperBound") <- attr(myfit_values,"upperBound")[which(names(myfit_values)%in%names(mypouter))]
  if(is.null(assign_pouter)){pouter <<- mypouter} else {cat("Manual assignment not yet provided.")}
  
  
  ## Define objective function -------
  
  cat("Generate objective function ...\n")
  if(!is.null(myerrors)){
    myobj <- normL2(mydata, myg*myx*myp0, errmodel = myerr) #+ constraintL2(prior, sigma=16)
  } else myobj <- normL2(mydata, myg*myx*myp0)
  if(is.null(assign_obj)){obj <<- myobj} else {cat("Manual assignment not yet provided.")}
  
  mytimes <- seq(0,max(do.call(c, lapply(1:length(mydata), function(i) max(mydata[[i]]$time)))), len=501)
  if(is.null(assign_times)){times <<- mytimes} else {cat("Manual assignment not yet provided.")}
  
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