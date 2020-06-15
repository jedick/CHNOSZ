# CHNOSZ/add.obigt.R
# add or change entries in the thermodynamic database

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("info.R")

today <- function() {
  # write today's date in the format used in SUPCRT data files
  # e.g. 13.May.12 for 2012-05-13
  t <- date()
  tt <- unlist(strsplit(t, " "))
  # for single-digit days there is an extra space
  tt <- tt[!tt==""]
  tday <- tt[3]
  tmonth <- tt[2]
  tyear <- substr(tt[5], start=3, stop=4)
  return(paste(tday, tmonth, tyear, sep="."))
}

mod.obigt <- function(...) {
  # add or modify species in thermo$obigt
  thermo <- get("thermo", CHNOSZ)
  # the names and values are in the arguments
  # this works for providing arguments via do.call
  args <- list(...)
  # this is needed if we are called with a list as the actual argument
  if(is.list(args[[1]])) args <- args[[1]]
  if(length(args) < 2) stop("please supply at least a species name and a property to update")
  if(is.null(names(args))) stop("all arguments after the first should be named")
  if(any(tail(nchar(names(args)), -1)==0)) stop("all arguments after the first should be named")
  # if the first argument is numeric, it's the species index
  if(is.numeric(args[[1]][1])) {
    ispecies <- args[[1]]
  } else {
    # if the name of the first argument is missing, assume it's the species name
    if(names(args)[1]=="") names(args)[1] <- "name"
    # search for this species, use check.protein=FALSE to avoid infinite loop when adding proteins
    # and suppressMessages to not show messages about matches of this name to other states
    if("state" %in% names(args)) ispecies <- suppressMessages(mapply(info.character, 
      species=args$name, state=args$state, check.protein=FALSE, SIMPLIFY=TRUE, USE.NAMES=FALSE))
    else ispecies <- suppressMessages(mapply(info.character, 
      species=args$name, check.protein=FALSE, SIMPLIFY=TRUE, USE.NAMES=FALSE))
  }
  # the column names of thermo$obigt, split at the "."
  cnames <- c(do.call(rbind, strsplit(colnames(thermo$obigt), ".", fixed=TRUE)), colnames(thermo$obigt))
  # the columns we are updating
  icol <- match(names(args), cnames)
  if(any(is.na(icol))) stop(paste("properties not in thermo$obigt:", paste(names(args)[is.na(icol)], collapse=" ")) )
  # the column numbers for properties that matched after the split
  icol[icol > 42] <- icol[icol > 42] - 42
  icol[icol > 21] <- icol[icol > 21] - 21
  # which species are new and which are old
  inew <- which(is.na(ispecies))
  iold <- which(!is.na(ispecies))
  # the arguments as data frame
  args <- data.frame(args, stringsAsFactors=FALSE)
  if(length(inew) > 0) {
    # the right number of blank rows of thermo$obigt
    newrows <- thermo$obigt[1:length(inew), ]
    # if we don't know something it's NA
    newrows[] <- NA
    # put in a default state
    newrows$state <- thermo$opt$state
    # the formula defaults to the name
    newrows$formula <- args$name[inew]
    # the units should also be set 20190530
    newrows$E_units <- thermo$opt$E.units
    # fill in the columns
    newrows[, icol] <- args[inew, ]
    # now check the formulas
    e <- tryCatch(makeup(newrows$formula), error=function(e) e)
    if(inherits(e, "error")) {
      warning("please supply a valid chemical formula as the species name or in the 'formula' argument")
      # transmit the error from makeup
      stop(e)
    }
    # for aqueous species, supply a value for Z if it is missing, otherwise NA triggers AkDi model 20190224
    isaq <- newrows$state == "aq"
    if(any(isaq)) {
      mnrf <- makeup(newrows$formula)
      if(nrow(newrows)==1) mnrf <- list(mnrf)
      Z <- sapply(mnrf, "[", "Z")
      Z[is.na(Z)] <- 0
      newrows$z.T[isaq] <- Z[isaq]
    }
    # assign to thermo$obigt
    thermo$obigt <- rbind(thermo$obigt, newrows)
    rownames(thermo$obigt) <- NULL
    assign("thermo", thermo, CHNOSZ)
    # update ispecies
    ntotal <- nrow(thermo$obigt)
    ispecies[inew] <- (ntotal-length(inew)+1):ntotal
    # inform user
    message(paste("mod.obigt: added ", newrows$name, "(", newrows$state, ")", " with energy units of ", newrows$E_units, sep="", collapse="\n"))
  }
  if(length(iold) > 0) {
    # loop over species
    for(i in 1:length(iold)) {
      # the old values and the state
      oldprop <- thermo$obigt[ispecies[iold[i]], icol]
      state <- thermo$obigt$state[ispecies[iold[i]]]
      # tell user if they're the same, otherwise update the data entry
      if(isTRUE(all.equal(oldprop, args[iold[i], ], check.attributes=FALSE))) 
        message("mod.obigt: no change for ", args$name[iold[i]], "(", state, ")")
      else {
        thermo$obigt[ispecies[iold[i]], icol] <- args[iold[i], ]
        assign("thermo", thermo, CHNOSZ)
        message("mod.obigt: updated ", args$name[iold[i]], "(", state, ")")
      }
    }
  }
  return(ispecies)
}

add.obigt <- function(file, species=NULL, force=TRUE) {
  # add/replace entries in thermo$obigt from values saved in a file
  # only replace if force==TRUE
  thermo <- get("thermo", CHNOSZ)
  to1 <- thermo$obigt
  id1 <- paste(to1$name,to1$state)
  # we match system files with the file suffixes removed (e.g. "CHNOSZ_aq" or "DEW_aq")
  sysfiles <- dir(system.file("extdata/OBIGT/", package="CHNOSZ"))
  sysnosuffix <- sapply(strsplit(sysfiles, "\\."), "[", 1)
  isys <- match(file, sysnosuffix)
  if(!is.na(isys)) file <- system.file(paste0("extdata/OBIGT/", sysfiles[isys]), package="CHNOSZ")
  else {
    # we also match single system files with the state suffix removed
    # (e.g. "DEW" for "DEW_aq", but not "organic" because we have "organic_aq", "organic_cr", etc.)
    sysnostate <- sapply(strsplit(sysnosuffix, "_"), "[", 1)
    isys <- which(file==sysnostate)
    if(length(isys)==1) file <- system.file(paste0("extdata/OBIGT/", sysfiles[isys]), package="CHNOSZ")
  }
  # read data from the file
  to2 <- read.csv(file, as.is=TRUE)
  # add E_units column if it's missing 20190529
  if(!"E_units" %in% colnames(to2)) to2 <- data.frame(to2[, 1:7], E_units = "cal", to2[, 8:20], stringsAsFactors = FALSE)
  Etxt <- paste(unique(to2$E_units), collapse = " and ")
  # load only selected species if requested
  if(!is.null(species)) {
    idat <- match(species, to2$name)
    ina <- is.na(idat)
    if(!any(ina)) to2 <- to2[idat, ]
    else stop(paste("file", file, "doesn't have", paste(species[ina], collapse=", ")))
  }
  id2 <- paste(to2$name,to2$state)
  # check if the data is compatible with thermo$obigt
  tr <- tryCatch(rbind(to1, to2), error = identity)
  if(inherits(tr, "error")) stop(paste(file, "is not compatible with thermo$obigt data table."))
  # match the new species to existing ones
  does.exist <- id2 %in% id1
  ispecies.exist <- na.omit(match(id2, id1))
  nexist <- sum(does.exist)
  # keep track of the species we've added
  inew <- numeric()
  if(force) {
    # replace existing entries
    if(nexist > 0) {
      to1[ispecies.exist, ] <- to2[does.exist, ]
      to2 <- to2[!does.exist, ]
      inew <- c(inew, ispecies.exist)
    }
  } else {
    # ignore any new entries that already exist
    to2 <- to2[!does.exist, ]
    nexist <- 0
  }
  # add new entries
  if(nrow(to2) > 0) {
    to1 <- rbind(to1, to2)
    inew <- c(inew, (length(id1)+1):nrow(to1))
  }
  # commit the change
  thermo$obigt <- to1
  rownames(thermo$obigt) <- 1:nrow(thermo$obigt)
  assign("thermo", thermo, CHNOSZ)
  # give the user a message
  message("add.obigt: read ", length(does.exist), " rows; made ", 
    nexist, " replacements, ", nrow(to2), " additions [energy units: ", Etxt, "]")
  #message("add.obigt: use obigt() or reset() to restore default database")
  return(invisible(inew))
}
