## Loading Libraries --------------------

library(PEtab2dmod)

## Import model--------------------
#1 Boehm_JProteomeRes2014 
#2 Fujita_SciSignal2010       - takes long
#3 Zheng_PNAS2012             - takes long
#4 Borghans_BiophysChem1997
#5 Elowitz_Nature2000
#6 Sneyd_PNAS2002
#7 Lucarelli_CellSystems2018  - takes long
#8 Crauste_CellSystems2017
#9 Schwen_PONE2014
#10 Raia_CancerResearch2011
#11 Fiedler_BMC2016

importPEtabSBML(modelname = "Boehm_JProteomeRes2014", compile = TRUE) 
plotPEtabSBML(name%in%names(observables))

testPEtabSBML(models = c("Elowitz_Nature2000"))

testPEtabSBML(tests = T)

modelname = "0002"
path2TestCases = "PEtabTests/"
SBML_file <- paste0(path2TestCases, modelname, "/_model.xml")
observable_file <- paste0(path2TestCases, modelname, "/_observables.tsv")
condition_file <- paste0(path2TestCases, modelname, "/_conditions.tsv")
data_file <- paste0(path2TestCases, modelname, "/_measurements.tsv")
parameter_file <- paste0(path2TestCases, modelname, "/_parameters.tsv")
mywd <- getwd()
files_loaded <- FALSE


modelname = "Boehm_JProteomeRes2014"
path2BC = "BenchmarkModels/"
SBML_file <- paste0(path2BC, modelname, "/model_", modelname, ".xml")
observable_file <- paste0(path2BC, modelname, "/observables_", modelname, ".tsv")
condition_file <- paste0(path2BC, modelname, "/experimentalCondition_", modelname, ".tsv")
data_file <- paste0(path2BC, modelname, "/measurementData_", modelname, ".tsv")
parameter_file <- paste0(path2BC, modelname, "/parameters_", modelname, ".tsv")
mywd <- getwd()
files_loaded <- FALSE


cat("Reading SBML file ...\n")
mylist <- getReactionsSBML(SBML_file)
myreactions <- mylist$reactions
myevents <- mylist$events

cat("Reading observables ...\n")
myobservables <- getObservablesSBML(observable_file)

cat("Compiling observable function ...\n")
if(!files_loaded) {
  setwd(paste0(mywd,"/CompiledObjects/"))
  myg <- Y(myobservables, myreactions, compile=TRUE, modelname=paste0("g_",modelname))
  setwd(mywd)
}

cat("Reading data file ...\n")
mydataSBML <- getDataPEtabSBML(data_file, observable_file)
mydata <- mydataSBML$data

cat("Compiling ODE model ...\n")

if(!files_loaded) {
  setwd(paste0(mywd,"/CompiledObjects/"))
  myodemodel <- odemodel(myreactions, forcings = NULL, events = myevents, fixed=NULL, modelname = paste0("odemodel_", modelname), jacobian = "inz.lsodes", compile = TRUE)
  setwd(mywd)
}

cat("Check and compile error model ...\n")
myerrors <- mydataSBML$errors
myerr <- NULL
if(!files_loaded) {
  if(!is_empty(getSymbols(myerrors))){
    setwd(paste0(mywd,"/CompiledObjects/"))
    myerr <- Y(myerrors, f = c(as.eqnvec(myreactions), myobservables), states = names(myobservables), attach.input = FALSE, compile = TRUE, modelname = paste0("errfn_", modelname))
    setwd(mywd)
  }
}
