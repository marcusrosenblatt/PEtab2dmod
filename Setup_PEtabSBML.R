## Loading Libraries --------------------

library(libSBML)
library(tidyverse)
library(dMod)
library(stringr)
source("Functions/getReactionsSBML.R")
source("Functions/getObservablesSBML.R")
source("Functions/getParametersSBML.R")
source("Functions/getConditionsSBML.R")
source("Functions/getDataPEtabSBML.R")
source("Functions/getInitialsSBML.R")
source("Functions/importPEtabSBML.R")
source("Functions/plotPEtabSBML.R")

## Import model--------------------
#1 Boehm_JProteomeRes2014 
#2 Fujita_SciSignal2010 
#3 Zheng_PNAS2012
#4 Borghans_BiophysChem1997
#5 Elowitz_Nature2000
#6 Sneyd_PNAS2002
#7 Lucarelli_CellSystems2018
#8 Crauste_CellSystems2017
#9 Schwen_PONE2014

importPEtabSBML(modelname = "Lucarelli_CellSystems2018", compile = TRUE) 
plotPEtabSBML(name%in%names(observables))
