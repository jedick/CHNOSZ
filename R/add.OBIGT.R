# CHNOSZ/add.OBIGT.R
# add or change entries in the thermodynamic database

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("info.R")
#source("util.data.R")

mod.OBIGT <- function(..., zap = FALSE) {
  # Add or modify species in thermo()$OBIGT
  thermo <- get("thermo", CHNOSZ)
  # The names and values are in the arguments
  # This works for providing arguments via do.call
  args <- list(...)
  # This is needed if we are called with a list as the actual argument
  if(is.list(args[[1]])) args <- args[[1]]
  if(length(args) < 2) stop("please supply at least a species name and a property to update")
  if(is.null(names(args))) stop("all arguments after the first should be named")
  if(any(tail(nchar(names(args)), -1)==0)) stop("all arguments after the first should be named")
  # If the first argument is numeric, it's the species index
  if(is.numeric(args[[1]][1])) {
    ispecies <- args[[1]]
    args <- args[-1]
    speciesname <- info(ispecies, check.it = FALSE)$name
  } else {
    # If the name of the first argument is missing, assume it's the species name
    if(names(args)[1]=="") names(args)[1] <- "name"
    speciesname <- args$name
    # Search for this species, use check.protein = FALSE to avoid infinite loop when adding proteins
    # and suppressMessages to not show messages about matches of this name to other states
    if("state" %in% names(args)) ispecies <- suppressMessages(mapply(info.character, 
      species = args$name, state = args$state, check.protein = FALSE, SIMPLIFY = TRUE, USE.NAMES = FALSE))
    else ispecies <- suppressMessages(mapply(info.character, 
      species = args$name, check.protein = FALSE, SIMPLIFY = TRUE, USE.NAMES = FALSE))
  }
  # The column names of thermo()$OBIGT, split at the "."
  cnames <- c(do.call(rbind, strsplit(colnames(thermo$OBIGT), ".", fixed=TRUE)), colnames(thermo$OBIGT))
  # The columns we are updating
  icol <- match(names(args), cnames)
  if(any(is.na(icol))) stop(paste("properties not in thermo$OBIGT:", paste(names(args)[is.na(icol)], collapse=" ")) )
  # The column numbers for properties that matched after the split
  icol[icol > 44] <- icol[icol > 44] - 44
  icol[icol > 22] <- icol[icol > 22] - 22
  # Which species are new and which are old
  inew <- which(is.na(ispecies))
  iold <- which(!is.na(ispecies))
  # The arguments as data frame
  args <- data.frame(args, stringsAsFactors=FALSE)
  if(length(inew) > 0) {
    # The right number of blank rows of thermo()$OBIGT
    newrows <- thermo$OBIGT[1:length(inew), ]
    # If we don't know something it's NA
    newrows[] <- NA
    # Put in a default state
    newrows$state <- thermo$opt$state
    # The formula defaults to the name
    newrows$formula <- args$name[inew]
    # The units should also be set 20190530
    newrows$E_units <- thermo$opt$E.units
    # Fill in the columns
    newrows[, icol] <- args[inew, ]
    # Now check the formulas
    e <- tryCatch(makeup(newrows$formula), error=function(e) e)
    if(inherits(e, "error")) {
      warning("please supply a valid chemical formula as the species name or in the 'formula' argument")
      # Transmit the error from makeup
      stop(e)
    }
    # Assign to thermo()$OBIGT
    thermo$OBIGT <- rbind(thermo$OBIGT, newrows)
    rownames(thermo$OBIGT) <- NULL
    assign("thermo", thermo, CHNOSZ)
    # Update ispecies
    ntotal <- nrow(thermo$OBIGT)
    ispecies[inew] <- (ntotal-length(inew)+1):ntotal
    # Inform user
    message(paste("mod.OBIGT: added ", newrows$name, "(", newrows$state, ")", " with energy units of ", newrows$E_units, sep="", collapse="\n"))
  }
  if(length(iold) > 0) {
    # Loop over species
    for(i in 1:length(iold)) {
      # The old values and the state
      oldprop <- thermo$OBIGT[ispecies[iold[i]], icol]
      state <- thermo$OBIGT$state[ispecies[iold[i]]]
      # Zap (clear) all preexisting values except for state 20220324
      if(zap) {
        thermo$OBIGT[ispecies[iold[i]], ] <- NA
        thermo$OBIGT$state[ispecies[iold[i]]] <- state
      }
      # Tell user if they're the same, otherwise update the data entry
      if(isTRUE(all.equal(oldprop, args[iold[i], ], check.attributes = FALSE))) 
        message("mod.OBIGT: no change for ", speciesname[iold[i]], "(", state, ")")
      else {
        thermo$OBIGT[ispecies[iold[i]], icol] <- args[iold[i], ]
        assign("thermo", thermo, CHNOSZ)
        message("mod.OBIGT: updated ", speciesname[iold[i]], "(", state, ")")
      }
    }
  }
  return(ispecies)
}

add.OBIGT <- function(file, species=NULL, force=TRUE) {
  # Add/replace entries in thermo$OBIGT from values saved in a file
  # Only replace if force==TRUE
  thermo <- get("thermo", CHNOSZ)
  to1 <- thermo$OBIGT
  id1 <- paste(to1$name,to1$state)
  # We match system files with the file suffixes (.csv) removed
  sysfiles <- dir(system.file("extdata/OBIGT/", package="CHNOSZ"))
  sysnosuffix <- sapply(strsplit(sysfiles, "\\."), "[", 1)
  isys <- match(file, sysnosuffix)
  if(!is.na(isys)) file <- system.file(paste0("extdata/OBIGT/", sysfiles[isys]), package="CHNOSZ")
#  else {
#    # We also match single system files with the state suffix removed
#    # (e.g. "DEW" for "DEW_aq", but not "organic" because we have "organic_aq", "organic_cr", etc.)
#    sysnostate <- sapply(strsplit(sysnosuffix, "_"), "[", 1)
#    isys <- which(file==sysnostate)
#    if(length(isys)==1) file <- system.file(paste0("extdata/OBIGT/", sysfiles[isys]), package="CHNOSZ")
#  }
  # Read data from the file
  to2 <- read.csv(file, as.is=TRUE)
  Etxt <- paste(unique(to2$E_units), collapse = " and ")
  # Load only selected species if requested
  if(!is.null(species)) {
    idat <- match(species, to2$name)
    ina <- is.na(idat)
    if(!any(ina)) to2 <- to2[idat, ]
    else stop(paste("file", file, "doesn't have", paste(species[ina], collapse=", ")))
  }
  id2 <- paste(to2$name,to2$state)
  # Check if the data is compatible with thermo$OBIGT
  tr <- tryCatch(rbind(to1, to2), error = identity)
  if(inherits(tr, "error")) stop(paste(file, "is not compatible with thermo$OBIGT data frame."))
  # Match the new species to existing ones
  does.exist <- id2 %in% id1
  ispecies.exist <- na.omit(match(id2, id1))
  nexist <- sum(does.exist)
  # Keep track of the species we've added
  inew <- numeric()
  if(force) {
    # Replace existing entries
    if(nexist > 0) {
      to1[ispecies.exist, ] <- to2[does.exist, ]
      to2 <- to2[!does.exist, ]
      inew <- c(inew, ispecies.exist)
    }
  } else {
    # Ignore any new entries that already exist
    to2 <- to2[!does.exist, ]
    nexist <- 0
  }
  # Add new entries
  if(nrow(to2) > 0) {
    to1 <- rbind(to1, to2)
    inew <- c(inew, (length(id1)+1):nrow(to1))
  }
  # Commit the change
  thermo$OBIGT <- to1
  rownames(thermo$OBIGT) <- 1:nrow(thermo$OBIGT)
  assign("thermo", thermo, CHNOSZ)
  # Give the user a message
  message("add.OBIGT: read ", length(does.exist), " rows; made ", 
    nexist, " replacements, ", nrow(to2), " additions [energy units: ", Etxt, "]")
  #message("add.OBIGT: use OBIGT() or reset() to restore default database")
  return(invisible(inew))
}
