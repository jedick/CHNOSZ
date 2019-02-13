# CHNOSZ/zzz.R
# this has the .onAttach function for package startup message and data initialization

# the CHNOSZ environment is made here in open code 20190214
# https://stackoverflow.com/questions/41954302/where-to-create-package-environment-variables
CHNOSZ <- new.env()

.onAttach <- function(libname,pkgname) {
  # version figuring adapted from package mgcv
  pkghelp <- library(help=CHNOSZ)$info[[1]]
  # things are different for older versions of R
  if(length(pkghelp)==1) pkghelp <- library(help=CHNOSZ)$info[[2]]
  version <- pkghelp[pmatch("Version:", pkghelp)]
  um <- strsplit(version, " ")[[1]]
  version <- um[nchar(um)>0][2]
  date <- pkghelp[pmatch("Date:", pkghelp)]
  um <- strsplit(date, " ")[[1]]
  date <- um[nchar(um)>0][2]
  # identify the program and version
  packageStartupMessage(paste("CHNOSZ version ", version, " (", date, ")", sep=""))
  # initialize the 'thermo' data object
  reset()
}
