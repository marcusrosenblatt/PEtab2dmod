exportPEtab <- function(modelname="test",
                        myODEmodel=ODEmodel,
                        myobservables=observables,
                        mydata=mydata,
                        myparameters=pouter,
                        myconditions=condition.grid){
  dir.create("Export/", showWarnings = FALSE)
  exportwd <- paste0("Export/",modelname)
  dir.create(exportwd, showWarnings = FALSE)
  
  cat("Writing observable file...\n")
  myobs <- as.data.frame(myobservables)
  myobs$myobservables <- as.character(myobs$myobservables)
  formula <- NULL
  trafo <- NULL
  for(i in 1:length(myobs$myobservables)){
    obs <- myobs$myobservables[i]
    if(str_detect(obs, "log10\\(")){
      formula <- c(formula, gsub(" ", "", substr(obs,7,length(strsplit(obs, "")[[1]])-1)))
      trafo <- c(trafo, "log10")
    } else if (str_detect(obs, "log\\(")){
      formula <- c(formula, gsub(" ", "", substr(obs,5,length(strsplit(obs, "")[[1]])-1)))
      trafo <- c(trafo, "log")
    } else {
      formula <- c(formula, gsub(" ", "", obs))
      trafo <- c(trafo, "lin")
    }
  }
  out <- data.frame(observableId = rownames(myobs), observableFormula = formula, observableTransformation = trafo)
  write_tsv(out, path = paste0(exportwd,"/observables_",modelname,".tsv"))
  
  
  cat("Writing condition file...\n")
  out <- myconditions
  if(!"conditionId"%in%colnames(myconditions)){
    out <- cbind(conditionId=paste0("model1_data",1:nrow(myconditions)), out)
  }
  if(!"conditionName"%in%colnames(myconditions)){
    out <- cbind(conditionName=rownames(myconditions), out)
  }
  write_tsv(out, path = paste0(exportwd,"/experimentalCondition_",modelname,".tsv"))
  
  
  cat("Writing parameter file...\n")
  if(!is.null(attr(pouter, "parscales"))){
    out <- data.frame(parameterId = names(pouter), parameterScale=attr(pouter, "parscales"))
  } else out <- data.frame(parameterId = names(pouter), parameterScale="log")
  if(!is.null(attr(pouter, "lowerBound"))){
    out <- cbind(out, lowerBound = attr(pouter, "lowerBound"))  
  } else out <- cbind(out, data.frame(parameterId = names(pouter), lowerBound="-23"))
  if(!is.null(attr(pouter, "upperBound"))){
    out <- cbind(out, upperBound = attr(pouter, "upperBound"))  
  } else out <- cbind(out, data.frame(parameterId = names(pouter), upperBound="23"))
  out <- cbind(out, estimate=1)
  write_tsv(out, path = paste0(exportwd,"/parameters_",modelname,".tsv"))
  
  cat(green(paste0("PEtab files written to Export/", modelname, "/\n")))
}
