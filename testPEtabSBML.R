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
source("Functions/fitModelPEtabSBML.R")

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


testPEtabSBML <- function(models=c("Boehm_JProteomeRes2014",
                                   "Fujita_SciSignal2010",
                                   "Borghans_BiophysChem1997"#,
                                   # "Elowitz_Nature2000",
                                   # "Sneyd_PNAS2002",
                                   # "Crauste_CellSystems2017",
                                   # "Schwen_PONE2014",
                                   # "Raia_CancerResearch2011"#,
                                   #"Zheng_PNAS2012"#,
                                   #"Bachmann_MSB2011"
                                   #"Lucarelli_CellSystems2018",
), testFit = TRUE, timelimit = 1000){
  cat(green("Start test function...\n"))
  teststarttime <- Sys.time()
  output <- NULL 
  for(model in models){
    importtest <- F
    plottest <- F
    bestfit <- NA
    cat(blue(paste0("Testing ", model, "\n")))
    fgh <- try_with_time_limit(
      {test <- try(importPEtabSBML(model, compile=T), silent = T)
      if(inherits(test, "try-error")) "import error" else test}
      , timelimit)
    if(fgh=="import error") cat(yellow("Import error or time limit exceeded for", model, "\n")) else {
        testobj <- obj(pouter)
        if(is.numeric(testobj$value)) cat(green("Calculation of objective function successful.\n")) else cat(red("Warning: obj(pouter) is not numeric.\n"))
        if(testFit){
          fitstarttime <- Sys.time()
          myframe <- fitModelPEtabSBML()
          fitendtime <- Sys.time()
          if(is.parframe(myframe) & !is.null(myframe))
            if(is.numeric(obj(myframe[1,])$value)){
              cat(green("Fit test successful.\n"))
              bestfit <- obj(myframe[1,])$value
            } else cat(red("Warning: obj(myframe) is not numeric.\n"))
          else cat(red("Warning: Fit test not successful..\n"))
          mytimediff <- as.numeric(difftime(fitendtime, fitstarttime, unit="secs"))
          if(mytimediff > 3600) cat(green(paste0("Fitting done in ",as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit="hours")), digits=3)), " hours.\n"))) else
            if(mytimediff > 60) cat(green(paste0("Fitting done in ",as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit="mins")), digits=3)), " minutes.\n"))) else
              cat(green(paste0("Fitting done in ",as.character(format(as.numeric(difftime(fitendtime, fitstarttime, unit="secs")), digits=3)), " seconds.\n")))
          
        }
        pdf(file = paste0("Test/",model,"_plotAll.pdf"))
        plotPEtabSBML()
        dev.off()
        pdf(file = paste0("Test/",model,"_plotTargetsObserved.pdf"))
        plotPEtabSBML(name%in%names(observables))
        dev.off()
        pdf(file = paste0("Test/",model,"_plotConditionsObserved.pdf"))
        plotPEtabSBML(condition%in%names(mydata))
        dev.off()
        importtest <- T
        plottest <- T
        cat(green("Import and plot test for ",fgh, " successful!\n\n\n"))
    }
    output <- rbind(output, data.frame(modelname = model, import = importtest, plot=plottest, objective_value = testobj$value, bestfit=bestfit))
  }
  testendtime <- Sys.time()
  mytimediff <- as.numeric(difftime(testendtime, teststarttime, unit="secs"))
  if(mytimediff > 3600) cat(green(paste0("Test done in ",as.character(format(as.numeric(difftime(testendtime, teststarttime, unit="hours")), digits=3)), " hours.\n"))) else
    if(mytimediff > 60) cat(green(paste0("Test done in ",as.character(format(as.numeric(difftime(testendtime, teststarttime, unit="mins")), digits=3)), " minutes.\n"))) else
      cat(green(paste0("Test done in ",as.character(format(as.numeric(difftime(testendtime, teststarttime, unit="secs")), digits=3)), " seconds.\n")))
  return(output)
}

testPEtabSBML()
