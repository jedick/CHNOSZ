# CHNOSZ/info.R
# Search database for species names or formulas
# and retrieve thermodynamic properties of species
# 20061024 Extracted from species.R jmd
# 20120507 Code rewrite and split into info.[character,approx,numeric];
#   These functions expect arguments of length 1; 
#   info() handles longer arguments

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.data.R")

info <- function(species = NULL, state = NULL, check.it = TRUE) {

  ## Return information for one or more species in thermo()$OBIGT
  thermo <- get("thermo", CHNOSZ)
  # That should give us the data, not the thermo() function 20190928
  if(is.function(thermo)) stop("CHNOSZ package data is not available; use reset() or library(CHNOSZ) to load it")
  ## If no species are requested, summarize the available data  20101129
  if(is.null(species)) {
    message("info: 'species' is NULL; summarizing information about thermodynamic data...")
    message(paste("thermo()$OBIGT has", nrow(thermo$OBIGT[thermo$OBIGT$state == "aq", ]), "aqueous,",
      nrow(thermo$OBIGT), "total species"))
    message(paste("number of literature sources: ", nrow(thermo$refs), ", elements: ",
      nrow(thermo$element), ", buffers: ", length(unique(thermo$buffer$name)), sep = ""))
    message(paste("number of proteins in thermo()$protein is", nrow(thermo$protein), "from",
      length(unique(thermo$protein$organism)), "organisms"))
    return()
  }

  ## Run info.numeric or info.character depending on the input type
  if(is.numeric(species)) {

    out <- lapply(species, info.numeric, check.it)
    # If we have different states the column names could be different
    if(length(unique(unlist(lapply(out, names)))) > ncol(thermo$OBIGT)) {
      # make them the same as thermo$OBIGT
      out <- lapply(out, function(row) {
        colnames(row) <- colnames(thermo$OBIGT); return(row)
      }) 
    }
    # Turn the list into a data frame
    out <- do.call(rbind, out)
    # Ensure that the rownames are numeric values (not names possibly inherited from retrieve()) 20190224
    if(!is.null(attr(species, "names"))) row.names(out) <- species

  } else {

    # State and species should be same length
    if(!is.null(state)) {
      lmax <- max(length(species), length(state))
      state <- rep(state, length.out = lmax)
      species <- rep(species, length.out = lmax)
    }
    # Loop over the species
    out <- sapply(seq_along(species), function(i) {
      # First look for exact match
      ispecies <- info.character(species[i], state[i])
      # If no exact match and it's not a protein, show approximate matches (side effect of info.approx)
      if(identical(ispecies, NA) & !grepl("_", species[i])) ispecies.notused <- info.approx(species[i], state[i])
      # Do not accept multiple matches
      if(length(ispecies) > 1) ispecies <- NA
      return(ispecies)
    })

  }

  ## All done!
  return(out)
}

### Unexported functions ###

info.text <- function(ispecies, withsource = TRUE) {
  # a textual description of species name, formula, source, e.g.
  # CO2 [CO2(aq)] (SSW01, SHS89, 11.Oct.07)
  this <- get("thermo", CHNOSZ)$OBIGT[ispecies, ]
  out <- paste(this$name, " [", this$formula, "(", this$state, ")", "]", sep = "")
  if(!withsource) return(out)
  sourcetext <- this$ref1
  ref2 <- this$ref2
  if(!is.na(ref2)) sourcetext <- paste(sourcetext, ref2, sep = ", ")
  date <- this$date
  if(!is.na(date)) sourcetext <- paste(sourcetext, date, sep = ", ")
  out <- paste(out, " (", sourcetext, ")", sep = "")
  return(out)
}

info.character <- function(species, state = NULL, check.protein = TRUE) {
  # Return the rownumbers of thermo()$OBIGT having an exact match of 'species' to
  # thermo()$OBIGT$[species|abbrv|formula] or NA otherwise
  # A match to thermo()$OBIGT$state is also required if 'state' is not NULL
  # (first occurence of a match to species is returned otherwise)

  thermo <- get("thermo", CHNOSZ)
  # Find matches for species name, abbreviation or formula
  matches.species <- thermo$OBIGT$name == species | thermo$OBIGT$abbrv == species | thermo$OBIGT$formula == species
  # Since thermo()$OBIGT$abbrv contains NAs, convert NA results to FALSE
  matches.species[is.na(matches.species)] <- FALSE
  # Turn it in to no match if it's a protein in the wrong state
  ip <- pinfo(species)
  if(any(matches.species) & !is.na(ip) & !is.null(state)) {
    matches.state <- matches.species & grepl(state, thermo$OBIGT$state)
    if(!any(matches.state)) matches.species <- FALSE
  }
  # No match, not available
  if(!any(matches.species)) {
    # Unless it's a protein
    if(check.protein) {
      # Did we find a protein? add its properties to OBIGT
      if(!is.na(ip)) {
        # Here we use a default state from thermo()$opt$state
        if(is.null(state)) state <- thermo$opt$state
        # Add up protein properties
        eos <- protein.OBIGT(ip, state = state)
        # The real assignment work 
        nrows <- suppressMessages(mod.OBIGT(eos))
        thermo <- get("thermo", CHNOSZ)
        matches.species <- rep(FALSE, nrows)
        matches.species[nrows] <- TRUE
      } else return(NA)
    } else return(NA)
  }
  # Do we demand a particular state
  if(!is.null(state)) {
    # Special treatment for H2O: aq retrieves the liq
    if(species %in% c("H2O", "water") & state == "aq") state <- "liq"
    # The matches for both species and state
    matches.state <- matches.species & state == thermo$OBIGT$state
    if(!any(matches.state)) {
      # The requested state is not available for this species
      available.states <- thermo$OBIGT$state[matches.species]
      if(length(available.states) == 1) a.s.verb <- "is" else a.s.verb <- "are"
      a.s.text <- paste("'", available.states, "'", sep = "", collapse = " ")
      message("info.character: requested state '", state, "' for ", species, 
        " but only ", a.s.text, " ", a.s.verb, " available")
      # Warn about looking for aqueous methane (changed to CH4) 20200707
      if(identical(species, "methane") & identical(state, "aq")) {
        warning("'methane' is not an aqueous species; use 'CH4' instead\nTo revert to the old behavior, run mod.OBIGT(info('CH4'), name = 'methane')")
      }
      return(NA)
    }
    matches.species <- matches.state
  }
  # All of the species that match
  ispecies.out <- ispecies <- which(matches.species)
  # Processing for more than one match
  if(length(ispecies) > 1) {
    # If a single name matches, use that one (useful for distinguishing pseudo-H4SiO4 and H4SiO4) 20171020
    matches.name <- matches.species & thermo$OBIGT$name == species
    if(sum(matches.name) == 1) ispecies.out <- which(matches.name)
    else ispecies.out <- ispecies[1]  # otherwise, return only the first species that matches
    # Let user know if there is more than one state for this species
    mystate <- thermo$OBIGT$state[ispecies.out]
    ispecies.other <- ispecies[!ispecies %in% ispecies.out]
    otherstates <- thermo$OBIGT$state[ispecies.other]
    # For non-aqueous species, list other substances (isomers) in the same state 20210321
    if(mystate != "aq" & sum(otherstates == mystate) > 0) {
      otherstates[otherstates == mystate] <- thermo$OBIGT$name[ispecies.other[otherstates == mystate]]
    }
    transtext <- othertext <- ""
    # We count, but don't show the states for phase transitions (cr2, cr3, etc)
    istrans <- otherstates %in% c("cr2", "cr3", "cr4", "cr5", "cr6", "cr7", "cr8", "cr9")
    if(mystate == "cr") {
      # If we are "cr" we show the number of phase transitions
      ntrans <- sum(istrans)
      if(ntrans == 1) transtext <- paste(" with", ntrans, "phase transition")
      else if(ntrans > 1) transtext <- paste(" with", ntrans, "phase transitions")
    }
    myname <- NULL
    if(mystate != "aq") {
      # If it's not already in the species name, append the substance name 20210323
      myname <- thermo$OBIGT$name[ispecies.out]
      if(species == myname) myname <- NULL
    }
    otherstates <- unique(otherstates[!istrans])
    if(length(otherstates) == 1) othertext <- paste0("; also available in ", otherstates)
    if(length(otherstates) > 1) othertext <- paste0("; also available in ", paste(otherstates, collapse = ", "))
    if(transtext != "" | othertext != "") {
      starttext <- paste0("info.character: found ", species, "(", mystate, ")")
      if(!is.null(myname)) starttext <- paste0(starttext, " [", myname, "]")
      message(starttext, transtext, othertext)
    }
  }
  return(ispecies.out)
}

info.numeric <- function(ispecies, check.it = TRUE) {
  # From a numeric species index in 'ispecies' return the 
  #   thermodynamic properties and equations-of-state parameters
  thermo <- get("thermo", CHNOSZ)
  # if we're called with NA, return an empty row
  if(is.na(ispecies)) {
    this <- thermo$OBIGT[1,]
    this[] <- NA
    return(this)
  }
  this <- thermo$OBIGT[ispecies,]
  # Species indices must be in range
  ispeciesmax <- nrow(thermo$OBIGT)
  if(ispecies > ispeciesmax | ispecies < 1) 
    stop(paste("species index", ispecies, "not found in thermo()$OBIGT\n"))

  # Remove scaling factors on EOS parameters depending on state
  # Use new OBIGT2eos function here
  this <- OBIGT2eos(this, this$state)

  # Stop with an informative message if species don't have a model 20220929
  namodel <- is.na(this$model)
  if(namodel) stop(paste("Species has NA model:", info.text(ispecies, FALSE)), call. = FALSE)

  if(tolower(this$model) == "berman") { # this is Berman
    # Get G, H, S, and V for minerals with Berman parameters 20220203
    Bermandat <- Berman()
    Bermandat <- Bermandat[Bermandat$name == this$name, ]
    this[, c("G", "H", "S", "V")] <- Bermandat[, c("GfPrTr", "HfPrTr", "SPrTr", "VPrTr")] * c(1, 1, 1, 10)
    isBerman <- TRUE
  } else isBerman <- FALSE

  # Identify any missing GHS values
  naGHS <- is.na(this[10:12])
  # A missing one of G, H or S can cause problems for subcrt calculations at high T
  if(sum(naGHS) == 1) {
    # calculate a single missing one of G, H, or S from the others
    GHS <- as.numeric(GHS(as.character(this$formula), G = this[, 10], H = this[, 11], S = this[, 12], E_units = this$E_units))
    message("info.numeric: ", colnames(this)[10:12][naGHS], " of ",
      this$name, "(", this$state, ") is NA; set to ", round(GHS[naGHS],2), " ", this$E_units, " mol-1")
    this[, which(naGHS) + 9] <- GHS[naGHS]
  } 

  # Perform consistency checks for GHS and EOS parameters if check.it = TRUE
  # Don't do it for the AD species 20190219
  if(check.it & this$model != "AD") {
    # Check GHS if they are all present
    if(sum(naGHS) == 0) calcG <- check.GHS(this, return.difference = FALSE)
    # Check heat capacities in database against EOS parameters
    calcCp <- check.EOS(this, this$model, "Cp", return.difference = FALSE)
    # Fill in NA heat capacity
    if(!is.na(calcCp) & is.na(this$Cp)) {
      message("info.numeric: Cp of ", this$name, "(", this$state, ") is NA; set by EOS parameters to ", round(calcCp, 2), " ", this$E_units, " K-1 mol-1")
      this$Cp <- as.numeric(calcCp)
    } else if(isBerman) {
      # Calculate Cp for Berman minerals 20220208
      calcCp <- Berman(this$name)$Cp
      this$Cp <- calcCp
    }
    # Check volumes in database - only for HKF model (aq species)
    if(this$model %in% c("HKF", "DEW")) {
      calcV <- check.EOS(this, this$model, "V", return.difference = FALSE)
      # Fill in NA volume
      if(!is.na(calcV) & is.na(this$V)) {
        message("info.numeric: V of ", this$name, "(", this$state, ") is NA; set by EOS parameters to ", round(calcV, 2), " cm3 mol-1")
        this$V <- as.numeric(calcV)
      }
    }
  } # Done checking
  # All done!
  return(this)
}

info.approx <- function(species, state = NULL) {
  # Returns species indices that have an approximate match of 'species'
  # to thermo$OBIGT$[name|abbrv|formula], possibly restricted to a given state
  thermo <- get("thermo", CHNOSZ)
  if(!is.null(state)) this <- thermo$OBIGT[thermo$OBIGT$state == state, ]
  else this <- thermo$OBIGT
  # Only look for fairly close matches
  max.distance <- 0.1
  approx.name <- agrep(species, this$name, max.distance)
  approx.abbrv <- agrep(species, this$abbrv, max.distance)
  approx.formula <- agrep(species, this$formula, max.distance)
  approx.species <- unique(c(approx.name, approx.abbrv, approx.formula))
  if(!is.na(approx.species[1])) {
    # Show the names of the species
    if(length(approx.species) == 1) {
      message("info.approx: '", species, "' is similar to ", info.text(approx.species))
    } else {
      napprox.max <- 100
      exttext <- ":"
      if(length(approx.species) > napprox.max) exttext <- paste(" (showing first ", napprox.max, ")", sep = "")
      message("info.approx: '", species, "' is ambiguous; has approximate matches to ", 
        length(approx.species), " species", exttext)
      printout <- capture.output(print(thermo$OBIGT$name[head(approx.species, napprox.max)]))
      message(paste(printout, collapse = "\n"))
    }
    return(approx.species)
  }
  # If we got here there were no approximate matches
  # 20190127 look for the species in optional data files 
  for(opt in c("SLOP98", "SUPCRT92", "AD")) {
    optdat <- read.csv(system.file(paste0("extdata/OBIGT/", opt, ".csv"), package = "CHNOSZ"), as.is = TRUE)
    if(species %in% optdat$name) {
      message('info.approx: ', species, ' is in an optional database; use add.OBIGT("', opt, '", "', species, '") to load it')
      return(NA)
    }
  }
  message("info.approx: '", species, "' has no approximate matches")
  return(NA)
}
