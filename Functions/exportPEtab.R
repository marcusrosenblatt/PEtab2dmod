exportPEtab <- function(modelname=modelname, ODEmodel=ODEmodel, observables=observables, mydata=mydata){
  dir.create(paste0("Export/",modelname))
  write.table(orgs, file = "orgs_updated.tsv", row.names=FALSE, sep="\t")
  cat(green(paste0("PEtab files written to Export/", modelname)))
  return(NULL)
}