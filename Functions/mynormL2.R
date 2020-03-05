mynormL2 <- function(data, x, errmodel = NULL, times = NULL, attr.name = "data", outputLL = FALSE) {
  
  timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time)))))
  if (!is.null(times)) timesD <- sort(union(times, timesD))
  
  x.conditions <- names(attr(x, "mappings"))
  data.conditions <- names(data)
  if (!all(data.conditions %in% x.conditions)) 
    stop("The prediction function does not provide predictions for all conditions in the data.")
  e.conditions <- names(attr(errmodel, "mappings"))
  
  controls <- list(times = timesD, attr.name = attr.name, conditions = intersect(x.conditions, data.conditions))
  
  # might be necessary to "store" errmodel in the objective function (-> runbg)
  force(errmodel)  
  
  myfn <- function(..., fixed = NULL, deriv=TRUE, conditions = controls$conditions, env = NULL) {
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    # Generate output template
    pars_out <- colnames(getDerivs(as.parvec(pouter)))
    template <- objlist(
      value = 0,
      gradient = structure(rep(0, length(pars_out)), names = pars_out),
      hessian = matrix(0, nrow = length(pars_out), ncol = length(pars_out), dimnames = list(pars_out, pars_out))
    )
    
    # Import from controls
    timesD <- controls$times
    attr.name <- controls$attr.name
    
    # Create new environment if necessary
    if (is.null(env)) env <- new.env()
    
    prediction <- x(times = timesD, pars = pouter, fixed = fixed, deriv = deriv, conditions = conditions)
    
    # Apply res() and wrss() to compute residuals and the weighted residual sum of squares
    out.data <- lapply(conditions, function(cn) {
      err <- NULL
      if ((!is.null(errmodel) & is.null(e.conditions)) | (!is.null(e.conditions) && (cn %in% e.conditions))) {
        err <- errmodel(out = prediction[[cn]], pars = getParameters(prediction[[cn]]), conditions = cn)
        mywrss <- nll(res(data[[cn]], prediction[[cn]], err[[cn]]))
      } else {
        if(outputLL) 
          mywrss <- nll(res(data[[cn]], prediction[[cn]]))
        else mywrss <- wrss(res(data[[cn]], prediction[[cn]]))  
      }
      available <- intersect(pars_out, names(mywrss$gradient))
      result <- template
      result$value <- mywrss$value
      if (deriv) {
        result$gradient[available] <- mywrss$gradient[available]
        result$hessian[available, available] <- mywrss$hessian[available, available]  
      } else {
        result$gradient <- result$hessian <- NULL
      }
      
      return(result)
    })
    out.data <- Reduce("+", out.data)
    
    # Combine contributions and attach attributes
    out <- out.data
    attr(out, controls$attr.name) <- out.data$value
    
    if (!is.null(env)) {
      assign("prediction", prediction, envir = env)
    }
    
    attr(out, "env") <- env
    return(out)
    
    
  }
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- data.conditions
  attr(myfn, "parameters") <- attr(x, "parameters")
  attr(myfn, "modelname") <- modelname(x, errmodel)
  return(myfn)
  
}
