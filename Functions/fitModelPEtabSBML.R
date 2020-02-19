fitModelPEtabSBML <- function(objfun=obj){
  prior <- structure(rep(0,length(pouter)))
  names(prior) <- names(pouter)
  mstrust(objfun=objfun, center=prior, studyname="corona", rinit = 0.1, rmax = 10,
          fits = 8, cores = 4, samplefun = "rnorm", resultPath = ".",
          parlower = -5, parupper=3,
          stats = FALSE, narrowing = NULL, iterlim=400, sd = 3)
}