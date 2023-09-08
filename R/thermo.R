# CHNOSZ/data/thermo.R
# Create or restore, and access the 'thermo' data object

# 20190213: Move from data/thermo.R to R/thermo.R
# --> invocation changes from data(thermo) to reset()
# --> invocation changes from data(OBIGT) to OBIGT()

reset <- function() {
  # Create thermo list
  thermodir <- system.file("extdata/thermo/", package = "CHNOSZ")
  thermo <- list(
    # Use as.is = TRUE to keep character values as character and not factor
    opt = as.list(read.csv(file.path(thermodir, "opt.csv"), as.is = TRUE)),
    element = read.csv(file.path(thermodir, "element.csv"), as.is = 1:3),
    OBIGT = NULL,
    refs = NULL,
    Berman = NULL,
    buffer = read.csv(file.path(thermodir, "buffer.csv"), as.is = 1:3),
    protein = read.csv(file.path(thermodir, "protein.csv"), as.is = 1:4),
    groups = read.csv(file.path(thermodir, "groups.csv"), row.names = 1, check.names = FALSE),
    stoich = read.csv(file.path(thermodir, "stoich.csv.xz"), as.is = TRUE),
    Bdot_acirc = read.csv(file.path(thermodir, "Bdot_acirc.csv"), as.is = TRUE),
    basis = NULL,
    species = NULL,
    opar = NULL
  )

  # Store stoich as matrix (with non-unique row names), not data frame
  formula <- thermo$stoich[, 1]
  thermo$stoich <- as.matrix(thermo$stoich[, 2:ncol(thermo$stoich)])
  rownames(thermo$stoich) <- formula
  
  # Make a named numeric vector for Bdot_acirc 20230309
  Bdot_acirc <- thermo$Bdot_acirc[, "acirc"]
  names(Bdot_acirc) <- thermo$Bdot_acirc[, "species"]
  thermo$Bdot_acirc <- Bdot_acirc

  # Get parameters in Berman equations from data files 20220203
  path <- system.file("extdata/Berman/", package = "CHNOSZ")
  files <- dir(path, "\\.csv$")
  # Put files in reverse chronological order (youngest first)
  files <- rev(files[order(sapply(strsplit(files, "_"), "[", 2))])
  # Read the parameters from each file
  Berman <- list()
  for(i in 1:length(files)) Berman[[i]] <- read.csv(file.path(path, files[i]), as.is = TRUE)
  # Assemble the parameters in a single data frame
  thermo$Berman <- do.call(rbind, Berman)

  # Message about what we are doing
  if(!"thermo" %in% ls(CHNOSZ)) packageStartupMessage("reset: creating \"thermo\" object")
  else packageStartupMessage("reset: resetting \"thermo\" object")
  # Place thermo in CHNOSZ environment
  assign("thermo", thermo, CHNOSZ)
  # Run OBIGT() to add the thermodynamic data
  OBIGT()
}

# Load default thermodynamic data (OBIGT) in thermo
OBIGT <- function(no.organics = FALSE) {
  # We only work if thermo is already in the CHNOSZ environment
  if(!"thermo" %in% ls(CHNOSZ)) stop("The CHNOSZ environment doesn't have a \"thermo\" object. Try running reset()")
  # Identify OBIGT data files
  sources_aq <- paste0(c("H2O", "inorganic", "organic"), "_aq")
  sources_cr <- paste0(c("Berman", "inorganic", "organic"), "_cr")
  sources_liq <- paste0(c("organic"), "_liq")
  sources_gas <- paste0(c("inorganic", "organic"), "_gas")
  sources <- c(sources_aq, sources_cr, sources_gas, sources_liq)
  OBIGTdir <- system.file("extdata/OBIGT/", package = "CHNOSZ")
  # Read data from each file
  datalist <- lapply(sources, function(source) {
    # Need explicit "/" for Windows
    file <- paste0(OBIGTdir, "/", source, ".csv")
    dat <- read.csv(file, as.is = TRUE)
    # Include CH4 with no.organics = TRUE 20220411
    if(no.organics & grepl("^organic", source)) dat <- subset(dat, formula == "CH4")
    dat
  })
  # Create OBIGT data frame
  OBIGT <- do.call(rbind, datalist)
  # Also read references file
  refs <- read.csv(file.path(OBIGTdir, "refs.csv"), as.is = TRUE)
  # Get thermo from CHNOSZ environment
  thermo <- get("thermo", CHNOSZ)
  # Set OBIGT and refs
  thermo$OBIGT <- OBIGT
  thermo$refs <- refs
  # Place modified thermo in CHNOSZ environment
  assign("thermo", thermo, CHNOSZ)
  # Message with brief summary of the data
  packageStartupMessage(paste("OBIGT: loading", ifelse(no.organics, "inorganic", "default"), "database with",
    nrow(thermo$OBIGT[thermo$OBIGT$state == "aq",]),
    "aqueous,", nrow(thermo$OBIGT), "total species"))
  # Warn if there are duplicated species
  idup <- duplicated(paste(thermo$OBIGT$name, thermo$OBIGT$state))
  if(any(idup)) warning("OBIGT: duplicated species: ", 
    paste(thermo$OBIGT$name[idup], "(", thermo$OBIGT$state[idup], ")", sep = "", collapse = " "))
}

# A function to access or modify the thermo object 20190214
# Revised for argument handling more like par() 20230310
thermo <- function (...) {

  # Get the arguments
  args <- list(...)

  # This part is taken from graphics::par()
  # To handle only names in c(), e.g. thermo(c("basis", "species"))
  if (all(unlist(lapply(args, is.character))))
      args <- as.list(unlist(args))
#  # To handle an argument analogous to old.par in the example for ?par
#  # ... but it breaks the "Adding an element" example in ?thermo
#  if (length(args) == 1) {
#      if (is.list(args[[1L]]) || is.null(args[[1L]])) 
#          args <- args[[1L]]
#  }

  # Get the name of the arguments
  argnames <- names(args)
  # Use "" for the name of each unnamed argument
  if(is.null(argnames)) argnames <- character(length(args))

  # Get the current 'thermo' object
  value <- original <- thermo <- get("thermo", CHNOSZ)
  # Loop over arguments
  if(length(args) > 0) {

    # Start with an empty return value with the right length
    value <- vector("list", length(args))
    for(i in 1:length(argnames)) {
      if(argnames[i] == "") {
        # For an unnnamed argument, retrieve the parameter from thermo
        # Parse the argument value to get the slots
        slots <- strsplit(args[[i]], "$", fixed = TRUE)[[1]]
        names <- args[[i]]
      } else {
        # For a named argument, assign the parameter in thermo
        # Parse the name of the argument to get the slots
        slots <- strsplit(argnames[i], "$", fixed = TRUE)[[1]]
        names <- argnames[i]
        # Perform the assignment in the local 'thermo' object
        if(length(slots) == 1) thermo[[slots[1]]] <- args[[i]]
        if(length(slots) == 2) thermo[[slots[1]]][[slots[2]]] <- args[[i]]
      }
      # Get the (original) parameter value
      if(length(slots) == 1) orig <- original[[slots[1]]]
      if(length(slots) == 2) orig <- original[[slots[1]]][[slots[2]]]
      # Put the parameter into the output value
      if(!is.null(orig)) value[[i]] <- orig
      names(value)[i] <- names
    }
    # Finally perform the assignment to 'thermo' in the CHNOSZ environment
    assign("thermo", thermo, CHNOSZ)

  }

  # Don't encapsulate a single unassigned parameter in a list
  if(is.null(names(args)) & length(value) == 1) value <- value[[1]]
  # Return the value
  value

}
