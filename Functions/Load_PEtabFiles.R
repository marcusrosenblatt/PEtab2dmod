library(libSBML)
library(dplyr)
library(tidyr)
library(dMod)

#### 2. Set model path and load files ####
modelpath1 <- 'BenchmarkModels/Boehm_JProteomeRes2014/'
modelpath2 <- 'BenchmarkModels/Fujita_SciSignal2010/'
modelpath3 <- 'BenchmarkModels/Zheng_PNAS2012/'

mymodel <- paste0(modelpath1,"model_Boehm_JProteomeRes2014.xml")
mymodel2 <- paste0(modelpath2,"model_Fujita_SciSignal2010.xml")
mymodel3 <- paste0(modelpath3,"model_Zheng_PNAS2012_original.xml")

myobservables1 <- paste0(modelpath1,"observables_Boehm_JProteomeRes2014.tsv")
myobservables2 <- paste0(modelpath2,"observables_Fujita_SciSignal2010.tsv")
myobservables3 <- paste0(modelpath3,"observables_Zheng_PNAS2012.tsv")

myconditions1 <- paste0(modelpath1,"experimentalCondition_Boehm_JProteomeRes2014.tsv")
myconditions2 <- paste0(modelpath2,"experimentalCondition_Fujita_SciSignal2010.tsv")
myconditions3 <- paste0(modelpath3,"experimentalCondition_Zheng_PNAS2012.tsv")

mydata1 <- paste0(modelpath1,"measurementData_Boehm_JProteomeRes2014.tsv")
mydata2 <- paste0(modelpath2,"measurementData_Fujita_SciSignal2010.tsv")
mydata3 <- paste0(modelpath3,"measurementData_Zheng_PNAS2012.tsv")

myparameters1 <- paste0(modelpath1,"parameters_Boehm_JProteomeRes2014.tsv")
myparameters2 <- paste0(modelpath2,"parameters_Fujita_SciSignal2010.tsv")
myparameters3 <- paste0(modelpath3,"parameters_Zheng_PNAS2012.tsv")
mypars <- read.csv(file = myparameters1, sep = "\t") 
  


### Model import

model <- paste0(modelpath,"model_Boehm_JProteomeRes2014.xml") #args <- commandArgs(trailingOnly = TRUE) #
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

## Define parameters, initial parameter trafo and steady-state transformations --------------
constraints <- resolveRecurrence(
  c(cyt = "1e9",
                                   t_inc = "14",
                                   #scale = 0.01,
                                   N0 = "S",
                                   Iacc = "0",
                                   E = "0",
                                   #I = "1",
                                   Q = "0",
                                   R = "0",
                                   D = "0"))

## compartments
m$getCompartment(0)$getId()
m$getCompartment(1)$getId()
m$getCompartment(0)$getSize()
m$getCompartment(1)$getSize()


