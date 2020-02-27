#' Fit a model imported via importPEtabSBML 
#' 
#' @description A wrapper function to use \link{mstrust} with imported PEtabSBML models. Some reasonable standard arguments for mstrust are used. Results of mstrust are written to Results folder.
#'  
#' @param objfun Objective function to be minimized as created by \link{importPEtabSBML}.
#' @param nrfits numeric, Number of fits to be performed
#' @param nrcores numeric, Number of cores to be used
#' @param useBounds  boolean, if TRUE, parameter bounds are taken as provided in PEtab format, if FALSE no parameter bounds are applied
#'   
#' @return parframe with the parameter estimated of the multi-start optimization
#'   
#' @author Marcus Rosenblatt and Svenja Kemmer
#'   
fitModelPEtabSBML <- function(objfun=obj, nrfits=4, nrcores=4, useBounds=TRUE){
  prior <- structure(rep(0,length(pouter)))
  names(prior) <- names(pouter)
  mywd <- getwd()
  dir.create(paste0(mywd,"/Test/mstrust/"), showWarnings = FALSE)
  if(useBounds) out <- mstrust(objfun=objfun, center=msParframe(prior, n = nrfits+1, seed=47)[-1,], studyname=model_name, rinit = 0.1, rmax = 10,
          fits = nrfits, cores = nrcores, samplefun = "rnorm", resultPath = "Test/mstrust/",
          parlower = attr(pouter, "lowerBound"), parupper=attr(pouter, "upperBound"),
          stats = FALSE, narrowing = NULL, iterlim=400, sd = 3)
  else out <- mstrust(objfun=objfun, center=msParframe(prior, n = nrfits, seed=47), studyname=model_name, rinit = 0.1, rmax = 10,
               fits = nrfits, cores = nrcores, samplefun = "rnorm", resultPath = "Test/mstrust/",
               stats = FALSE, narrowing = NULL, iterlim=400, sd = 3)
    if(any(lapply(out, function(fgh) fgh$converged)==TRUE)) return(as.parframe(out)) else {cat("No fit converged."); return(NULL)}
}