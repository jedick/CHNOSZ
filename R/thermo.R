# CHNOSZ/data/thermo.R
# create or restore, and access the 'thermo' data object

# 20190213: move from data/thermo.R to R/thermo.R
# --> invocation changes from data(thermo) to reset()
# --> invocation changes from data(OBIGT) to OBIGT()

reset <- function() {
  # create thermo list
  thermodir <- system.file("extdata/thermo/", package="CHNOSZ")
  thermo <- list(
    # as.is: keep character values as character and not factor
    opt = as.list(read.csv(file.path(thermodir, "opt.csv"), as.is=TRUE)),
    element = read.csv(file.path(thermodir, "element.csv"), as.is=1:3),
    OBIGT = NULL,
    refs = NULL,
    buffer = read.csv(file.path(thermodir, "buffer.csv"), as.is=1:3),
    protein = read.csv(file.path(thermodir, "protein.csv"), as.is=1:4),
    groups = read.csv(file.path(thermodir, "groups.csv"), row.names=1, check.names=FALSE),
    stoich = read.csv(file.path(thermodir, "stoich.csv.xz"), as.is=TRUE),
    basis = NULL,
    species = NULL,
    opar = NULL
  )
  # store stoich as matrix (with non-unique row names), not data frame
  formula <- thermo$stoich[, 1]
  thermo$stoich <- as.matrix(thermo$stoich[, 2:ncol(thermo$stoich)])
  rownames(thermo$stoich) <- formula
  # give a summary of what we are doing
  if(!"thermo" %in% ls(CHNOSZ)) message("reset: creating \"thermo\" object")
  else message("reset: resetting \"thermo\" object")
  # place thermo in CHNOSZ environment
  assign("thermo", thermo, CHNOSZ)
  # run OBIGT() to add the thermodynamic data
  OBIGT()
}

# load default thermodynamic data (OBIGT) in thermo
OBIGT <- function() {
  # we only work if thermo is already in the CHNOSZ environment
  if(!"thermo" %in% ls(CHNOSZ)) stop("The CHNOSZ environment doesn't have a \"thermo\" object. Try running reset()")
  # create OBIGT data frame
  sources_aq <- paste0(c("H2O", "inorganic", "organic", "biotic"), "_aq")
  sources_cr <- paste0(c("inorganic", "organic", "Berman"), "_cr")
  sources_liq <- paste0(c("organic"), "_liq")
  sources_gas <- paste0(c("inorganic", "organic"), "_gas")
  sources <- c(sources_aq, sources_cr, sources_gas, sources_liq)
  OBIGTdir <- system.file("extdata/OBIGT/", package="CHNOSZ")
  # need explicit "/" for Windows
  sourcefiles <- paste0(OBIGTdir, "/", c(sources_aq, sources_cr, sources_gas, sources_liq), ".csv")
  datalist <- lapply(sourcefiles, read.csv, as.is=TRUE)
  OBIGT <- do.call(rbind, datalist)
  # also read references file
  refs <- read.csv(file.path(OBIGTdir, "refs.csv"), as.is=TRUE)
  # get thermo from CHNOSZ environment
  thermo <- get("thermo", CHNOSZ)
  # set OBIGT and refs
  thermo$OBIGT <- OBIGT
  thermo$refs <- refs
  # place modified thermo in CHNOSZ environment
  assign("thermo", thermo, CHNOSZ)
  # give a summary of some of the data
  message(paste("OBIGT: loading default database with",
    nrow(thermo$OBIGT[thermo$OBIGT$state=="aq",]),
    "aqueous,", nrow(thermo$OBIGT), "total species"))
  # warn if there are duplicated species
  idup <- duplicated(paste(thermo$OBIGT$name, thermo$OBIGT$state))
  if(any(idup)) warning("OBIGT: duplicated species: ", 
    paste(thermo$OBIGT$name[idup], "(", thermo$OBIGT$state[idup], ")", sep="", collapse=" "))
}

# a function to access or modify the thermo object 20190214
thermo <- function(...) {
  args <- list(...)
  # get the object
  thermo <- get("thermo", CHNOSZ)
  if(length(args) > 0) {
    # assign into the object
    slots <- names(args)
    for(i in 1:length(slots)) {
      # parse the name of the slot
      names <- strsplit(slots[i], "$", fixed=TRUE)[[1]]
      if(length(names) == 1) thermo[[names]] <- args[[i]]
      if(length(names) == 2) thermo[[names[1]]][[names[2]]] <- args[[i]]
    }
    assign("thermo", thermo, CHNOSZ)
  } else {
    # return the object
    thermo
  }
}
