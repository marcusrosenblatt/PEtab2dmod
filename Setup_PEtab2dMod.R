## Loading Libraries --------------------

library(libSBML)
library(tidyverse)
library(dMod)
library(stringr)
source("Functions/getReactionsSBML.R")
source("Functions/getObservablesSBML.R")
source("Functions/getParametersSBML.R")
source("Functions/getConditionsSBML.R")
source("Functions/getDataSBML.R")
source("Functions/getInitialsSBML.R")

## Load SBML and PEtab files --------------------

modelpath1 <- 'BenchmarkModels/Boehm_JProteomeRes2014/'
modelpath2 <- 'BenchmarkModels/Fujita_SciSignal2010/'
modelpath3 <- 'BenchmarkModels/Zheng_PNAS2012/'

model_Boehm <- list(model = paste0(modelpath1,"model_Boehm_JProteomeRes2014.xml"), 
                    observables = paste0(modelpath1,"observables_Boehm_JProteomeRes2014.tsv"),
                    conditions = paste0(modelpath1,"experimentalCondition_Boehm_JProteomeRes2014.tsv"),
                    data = paste0(modelpath1,"measurementData_Boehm_JProteomeRes2014.tsv"),
                    parameters = paste0(modelpath1,"parameters_Boehm_JProteomeRes2014.tsv"))

model_Fujita <- list(model = paste0(modelpath2,"model_Fujita_SciSignal2010.xml"), 
                     observables = paste0(modelpath2,"observables_Fujita_SciSignal2010.tsv"),
                     conditions = paste0(modelpath2,"experimentalCondition_Fujita_SciSignal2010.tsv"),
                     data = paste0(modelpath2,"measurementData_Fujita_SciSignal2010.tsv"),
                     parameters = paste0(modelpath2,"parameters_Fujita_SciSignal2010.tsv"))

model_Zheng <- list(model = paste0(modelpath3,"model_Zheng_PNAS2012_original.xml"), 
                     observables = paste0(modelpath3,"observables_Zheng_PNAS2012.tsv"),
                     conditions = paste0(modelpath3,"experimentalCondition_Zheng_PNAS2012.tsv"),
                     data = paste0(modelpath3,"measurementData_Zheng_PNAS2012.tsv"),
                     parameters = paste0(modelpath3,"parameters_Zheng_PNAS2012.tsv"))


## Model Definition - Equations --------------------

mymodel <- model_Fujita
model_name <- "test"
reactions <- getReactionsSBML(mymodel$model)$reactions
events <- getReactionsSBML(mymodel$model)$events

## Model Definition - Observables --------------------

observables <- getObservablesSBML(mymodel$observables)

# Observation function
g <- Y(observables, reactions, compile=TRUE, modelname=paste0("g_",model_name))


## Get Data ------------
data <- getDataSBML(mymodel$data,observables)$data


## Model Generation ---------------------

modeltest <- odemodel(reactions, forcings = NULL,
                     events = events,
                     fixed=NULL, modelname = paste0("odemodel_", model_name),
                     jacobian = "inz.lsodes", compile = TRUE)


## Define constraints, initials, parameters and compartments --------------

constraints <- getParametersSBML(mymodel$parameters)$constraints
initials <- getInitialsSBML(mymodel$model)$initials
compartments <- getInitialsSBML(mymodel$model)$compartments
fit_values <- getParametersSBML(mymodel$parameters)$pouter

## Define error model if applicable ------------

# check whether errors are estimated
errors <- getDataSBML(mymodel$data,observables)$errors
if(!is.null(errors))
{
  # gernerate error function
  err <- Y(errors, f = c(as.eqnvec(reactions), observables), states = names(observables), 
           attach.input = FALSE, compile = F, modelname = "errfn")
}

## Define inner parameters ----------

innerpars <- unique(c(getParameters(modeltest), getSymbols(observables), getSymbols(errors)))
names(innerpars) <- innerpars

## Parameter transformation -----------

trafo <- as.eqnvec(innerpars, names = innerpars)#replaceSymbols(names(observables), observables, innerpars)
trafo <- replaceSymbols(names(initials), initials, trafo)
trafo <- replaceSymbols(names(compartments), compartments, trafo)
#trafo <- replaceSymbols(c("t","time"), 0, trafo)
trafo <- replaceSymbols(names(constraints), constraints, trafo)

# Generate condition.grid
condition.grid <- getConditionsSBML(mymodel$conditions, mymodel$data)
# remove NAs
vec <- NULL
for (i in 1:nrow(condition.grid)) {if(Reduce("&",!is.na(condition.grid[i,]))) vec <- c(vec, i)}
condition.grid <- condition.grid[vec,]

parameters <- names(condition.grid)[!names(condition.grid) %in% c("conditionName")]

# branch trafo for different conditions
trafoL <- branch(trafo, table=condition.grid) %>%
  insert("x~0", x = unique(events$var)) %>%
  insert("x~10**(x)", x = .currentSymbols) 
  
for (j in 1:length(names(trafoL))) {
  for (i in 1:length(parameters)) {
    trafoL[[j]] <- repar(x~y, trafoL[[j]], x=parameters[i], y=condition.grid[j,i+1])
  }
}

## Specify prediction functions ------
tolerances <- 1e-7
p0 <- x <- NULL
for (C in names(trafoL)) {
  p0 <- p0 + P(trafoL[[C]], condition = C)
  x <- x + Xs(modeltest, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
              optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
              condition = C)
}

# p <- P(trafoL, modelname = "p")
# compile(p,g,err, output = "g_err_p", cores = detectFreeCores())

## Generate objective function and initial parameter set -------
outerpars <- getParameters(p0)
pouter <- structure(rep(NA,length(outerpars)), names = outerpars)

common <- intersect(names(pouter),names(fit_values))
pouter[common] <- fit_values[common]

## Define objective function -------

if(!is.null(errors))
  {
  obj <- normL2(data, g*x*p0, err) #+ constraintL2(prior, sigma=16)
} else obj <- normL2(data, g*x*p0)

times <- 0:max(data[[1]]$time)

# times <- 0:300
# 
# 
# 
# ## Perform some fits and plot result of best fit
# out <- mstrust(objfun=obj, center=prior, studyname="example", rinit = 0.1, rmax = 10,
#                fits = 8, cores = 2, samplefun = "rnorm", resultPath = ".",
#                stats = FALSE, narrowing = NULL, iterlim=200, sd = 3)
# 
# fitlist <- as.parframe(out)
# 
# bestfit <- unlist(fitlist[1,-c(1:4)])

prediction <- (g*x*p0)(times, pouter)
plotCombined(prediction, data) 



