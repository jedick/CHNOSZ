# CHNOSZ/palply.R

palply <- function(varlist, X, FUN, ...) {
  # A wrapper function to run parLapply if length(X) >= thermo()$opt$paramin
  # and package 'parallel' is available, otherwise run lapply
  if(length(X) >= get("thermo", CHNOSZ)$opt$paramin) {
    # Use option mc.cores to choose an appropriate cluster size.
    # and set max at 2 for now (per CRAN policies)
    nCores <- min(getOption("mc.cores"), get("thermo", CHNOSZ)$opt$maxcores)
    # Don't load methods package, to save startup time - ?makeCluster
    cl <- parallel::makeCluster(nCores, methods=FALSE)
    # Export the variables and notify the user
    if(! "" %in% varlist) {
      parallel::clusterExport(cl, varlist)
      message(paste("palply:", caller.name(4), "running", length(X), "calculations on",
        nCores, "cores with variable(s)", paste(varlist, collapse=", ")))
    } else {
      message(paste("palply:", caller.name(4), "running", length(X), "calculations on",
        nCores, "cores"))
    }
    # Run the calculations
    out <- parallel::parLapply(cl, X, FUN, ...)
    parallel::stopCluster(cl)
  } else out <- lapply(X, FUN, ...)
  return(out)
}
