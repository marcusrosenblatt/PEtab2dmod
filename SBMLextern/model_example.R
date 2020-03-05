## Loading Libraries --------------------

library(deSolve)
library(trust)
library(parallel)
library(ggplot2)
library(ggthemes)
library(cOde)
library(dMod)
library(blotIt)
library(cowplot)
library(magrittr)
library(tidyr)

## Model Definition - Equations --------------------

model_name <- "corona_27_01"
flist <- NULL %>% addReaction("S", "E", "R0*gam*S*I/N0") %>%
  addReaction("E", "I", "1/t_inc*E") %>%
  addReaction("", "Iacc", "1/t_inc*E") %>%
  addReaction("I", "Q", "gam*I/(1+ratio)") %>%
  addReaction("Q", "R", "Q*k0") %>%
  addReaction("I", "D", "gam*ratio*I/(1+ratio)") 


## Model Definition - Observables --------------------

observables <- eqnvec(#Iobs = "Iacc*scale/(1+scale)+Ioffset")
                      Iobs = "Iacc*scale+Ioffset",
                      Robs = "R*scale+Roffset",
                      Dobs = "D")

## Model Generation ---------------------

modelCorona <- odemodel(flist, forcings = NULL,
                     events = NULL,
                     fixed=NULL, modelname = paste0("odemodel_", model_name),
                     jacobian = "inz.lsodes", compile = TRUE)


## Define parameters, initial parameter trafo and steady-state transformations --------------
constraints <- resolveRecurrence(c(S = "1e9",
                                   t_inc = "14",
                                   #scale = 0.01,
                                   N0 = "S",
                                   Iacc = "0",
                                   E = "0",
                                   #I = "1",
                                   Q = "0",
                                   R = "0",
                                   D = "0"))

innerpars <- unique(c(getParameters(modelCorona), getSymbols(observables)))
names(innerpars) <- innerpars

trafo <- replaceSymbols(names(observables), observables, innerpars)
trafo <- replaceSymbols(names(constraints), constraints, trafo)

## Generate condition.grid
condition.grid <- data.frame(MGE = c(0),
                       row.names = c("seir"))

## Get Data ------------
my <- as.datalist(
  rbind(
  data.frame(name="Iobs", time=c(3,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33, 34, 35, 36),
             value=c(45,62,121,198,278,440,573,834,1294,1982,2757,4530, 5989, 7801, 9692, 11791, 14380, 17205, 20438, 24324, 28018),
             sigma=1, lloq=-1, condition="seir"),
  data.frame(name="Dobs", time=c(3,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34, 35, 36), 
             value=c(0,0,0,0,4,6,17,26,41,56,80,106,132,170,213,259, 304, 361, 425, 490, 563),
             sigma=1, lloq=-1, condition="seir"),
  data.frame(name="Robs", time=c(15,16,17,18,19,20,23,24,25,26,27,28,29,30,31,32,33,34, 35, 36), 
             value=c(7,12,15,19,25,25,34, 38, 49, 51, 60, 103, 124, 171, 243, 328, 475, 632, 892, 1153),
             sigma=1, lloq=-1, condition="seir")
  ))
data[[1]]$sigma <- data[[1]]$value**0.5 + 1

# ## Parameter transformation -----------
trafoL <- branch(trafo, table=condition.grid) %>%
#   insert(x~y, x=c("MGE"), y=MGE) %>%
#   insert(x~y, x=c("scale_HBx"), y=c("scale_HBx_TC"), conditionMatch="infected") %>%
#   insert(x~1, x=c("scale_HBx")) %>%
insert("x~exp(x)", x = .currentSymbols)

## Specify prediction functions ------
tolerances <- 1e-7
p0 <- x <- NULL
for (C in names(trafoL)) {
  p0 <- p0 + P(trafoL[[C]], condition = C)
  x <- x + Xs(modelCorona, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 5000),
              optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
              condition = C)
}

## Generate objective function and initial parameter set -------
outerpars <- attr(p0, "parameters")
pouter <- structure(rnorm(length(outerpars)), names = outerpars)
prior <- rep(0, length(outerpars)); names(prior) <- outerpars
g <- Y(observables, flist, compile=TRUE, modelname=paste0("g_",model_name))

obj <- normL2(data, g*x*p0) + constraintL2(prior, sigma=16)

sigma <- rep(3, length(outerpars)); names(sigma) <- outerpars
times <- 0:100

# Error model
which_err <- c(1:length(observables))
errors <- paste("sigma_rel", names(observables)[which_err], sep = "_")
names(errors) <- names(observables[which_err])

err <- Y(errors, f = c(as.eqnvec(reactions), observables), states = names(observables), 
         attach.input = FALSE, compile = F, modelname = "errfn")


## Perform some fits and plot result of best fit
out <- mstrust(objfun=obj, center=prior, studyname="example", rinit = 0.1, rmax = 10,
               fits = 8, cores = 2, samplefun = "rnorm", resultPath = ".",
               stats = FALSE, narrowing = NULL, iterlim=200, sd = 3)

fitlist <- as.parframe(out)

bestfit <- unlist(fitlist[1,-c(1:4)])

prediction <- (g*x*p0)(seq(0,100,by=0.1), bestfit)
plotCombined(prediction, data) 



