## Loading Libraries --------------------

library(libSBML)
library(tidyverse)
library(dMod)
library(stringr)
library(crayon)
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
try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf){
  y <- try({setTimeLimit(cpu, elapsed); expr}, silent = TRUE) 
  if(inherits(y, "try-error")) NULL else y 
}

testPEtabSBML <- function(timelimit = 100, 
                          models=c("Boehm_JProteomeRes2014", "Fujita_SciSignal2010", "Zheng_PNAS2012")){
  for(model in models){
    fgh <- try_with_time_limit(
      {test <- try(importPEtabSBML(model, compile=T), silent = T)
      if(inherits(test, "try-error")) "import error" else test}
      , timelimit)
    if(fgh=="import error") cat(blue("Import error in", model)) else
      if(!is.null(fgh)){
        pdf(file = paste0("Test/plot_",model,".pdf"))
        plotPEtabSBML()
        dev.off()
        cat(green("Import and plot test for ",fgh, " successful!\n"))
      } else cat(blue("Time limit for",model, "exceeded."))
  }
}

testPEtabSBML()
