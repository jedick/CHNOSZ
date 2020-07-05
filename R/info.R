# CHNOSZ/info.R
# search database for species names or formulas
# and retrieve thermodynamic properties of species
# 20061024 extracted from species.R jmd
# 20120507 code rewrite and split into info.[character,approx,numeric];
#   these functions expect arguments of length 1; 
#   info() handles longer arguments

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.data.R")

info <- function(species=NULL, state=NULL, check.it=TRUE) {
  ## return information for one or more species in thermo()$OBIGT
  thermo <- get("thermo", CHNOSZ)
  # that should give us the data, not the thermo() function 20190928
  if(is.function(thermo)) stop("CHNOSZ package data is not available; use reset() or library(CHNOSZ) to load it")
  ## if no species are requested, summarize the available data  20101129
  if(is.null(species)) {
    message("info: 'species' is NULL; summarizing information about thermodynamic data...")
    message(paste("thermo()$OBIGT has", nrow(thermo$OBIGT[thermo$OBIGT$state=="aq", ]), "aqueous,",
      nrow(thermo$OBIGT), "total species"))
    message(paste("number of literature sources: ", nrow(thermo$refs), ", elements: ",
      nrow(thermo$element), ", buffers: ", length(unique(thermo$buffers$name)), sep=""))
    message(paste("number of proteins in thermo()$protein is", nrow(thermo$protein), "from",
      length(unique(thermo$protein$organism)), "organisms"))
    return()
  }
  ## run info.numeric or info.character depending on the input type
  if(is.numeric(species)) {
    out <- lapply(species, info.numeric, check.it)
    # if we have different states the column names could be different
    if(length(unique(unlist(lapply(out, names)))) > ncol(thermo$OBIGT)) {
      # make them the same as thermo$OBIGT
      out <- lapply(out, function(row) {
        colnames(row) <- colnames(thermo$OBIGT); return(row)
      }) 
    }
    # turn the list into a data frame
    out <- do.call(rbind, out)
    # ensure that the rownames are numeric values (not names possibly inherited from retrieve()) 20190224
    if(!is.null(attr(species, "names"))) row.names(out) <- species
  } else {
    # state and species should be same length
    if(!is.null(state)) {
      lmax <- max(length(species), length(state))
      state <- rep(state, length.out=lmax)
      species <- rep(species, length.out=lmax)
    }
    # loop over the species
    out <- sapply(seq_along(species), function(i) {
      # first look for exact match
      ispecies <- info.character(species[i], state[i])
      # if no exact match and it's not a protein, show approximate matches (side effect of info.approx)
      if(identical(ispecies, NA) & !grepl("_", species[i])) ispecies.notused <- info.approx(species[i], state[i])
      # do not accept multiple matches
      if(length(ispecies) > 1) ispecies <- NA
      return(ispecies)
    })
  }
  ## all done!
  return(out)
}

### unexported functions ###

info.text <- function(ispecies) {
  # a textual description of species name, formula, source, e.g.
  # CO2 [CO2(aq)] (SSW01, SHS89, 11.Oct.07)
  this <- get("thermo", CHNOSZ)$OBIGT[ispecies, ]
  sourcetext <- this$ref1
  ref2 <- this$ref2
  if(!is.na(ref2)) sourcetext <- paste(sourcetext, ref2, sep=", ")
  date <- this$date
  if(!is.na(date)) sourcetext <- paste(sourcetext, date, sep=", ")
  out <- paste(this$name, " [", this$formula, "(", this$state, ")", "] (", sourcetext, ")", sep="")
  return(out)
}

info.character <- function(species, state=NULL, check.protein=TRUE) {
  # returns the rownumbers of thermo()$OBIGT having an exact match of 'species' to
  # thermo()$OBIGT$[species|abbrv|formula] or NA otherwise
  # a match to thermo()$OBIGT$state is also required if 'state' is not NULL
  # (first occurence of a match to species is returned otherwise)
  thermo <- get("thermo", CHNOSZ)
  # find matches for species name, abbreviation or formula
  matches.species <- thermo$OBIGT$name==species | thermo$OBIGT$abbrv==species | thermo$OBIGT$formula==species
  # since thermo()$OBIGT$abbrv contains NAs, convert NA results to FALSE
  matches.species[is.na(matches.species)] <- FALSE
  # turn it in to no match if it's a protein in the wrong state
  ip <- pinfo(species)
  if(any(matches.species) & !is.na(ip) & !is.null(state)) {
    matches.state <- matches.species & grepl(state, thermo$OBIGT$state)
    if(!any(matches.state)) matches.species <- FALSE
  }
  # no match, not available
  if(!any(matches.species)) {
    # unless it's a protein
    if(check.protein) {
      # did we find a protein? add its properties to OBIGT
      if(!is.na(ip)) {
        # here we use a default state from thermo()$opt$state
        if(is.null(state)) state <- thermo$opt$state
        # add up protein properties
        eos <- protein.OBIGT(ip, state=state)
        # the real assignment work 
        nrows <- suppressMessages(mod.OBIGT(eos))
        thermo <- get("thermo", CHNOSZ)
        matches.species <- rep(FALSE, nrows)
        matches.species[nrows] <- TRUE
      } else return(NA)
    } else return(NA)
  }
  # do we demand a particular state
  if(!is.null(state)) {
    # special treatment for H2O: aq retrieves the liq
    if(species %in% c("H2O", "water") & state=="aq") state <- "liq"
    # the matches for both species and state
    matches.state <- matches.species & state == thermo$OBIGT$state
    if(!any(matches.state)) {
      # the requested state is not available for this species
      available.states <- thermo$OBIGT$state[matches.species]
      if(length(available.states)==1) a.s.verb <- "is" else a.s.verb <- "are"
      a.s.text <- paste("'", available.states, "'", sep="", collapse=" ")
      message("info.character: requested state '", state, "' for ", species, 
        " but only ", a.s.text, " ", a.s.verb, " available")
      return(NA)
    }
    matches.species <- matches.state
  }
  # all of the species that match
  ispecies.out <- ispecies <- which(matches.species)
  # processing for more than one match
  if(length(ispecies) > 1) {
    # if a single name matches, use that one (useful for distinguishing pseudo-H4SiO4 and H4SiO4) 20171020
    matches.name <- matches.species & thermo$OBIGT$name==species
    if(sum(matches.name)==1) ispecies.out <- which(matches.name)
    else ispecies.out <- ispecies[1]  # otherwise, return only the first species that matches
    # let user know if there is more than one state for this species
    mystate <- thermo$OBIGT$state[ispecies.out]
    ispecies.other <- ispecies[!ispecies %in% ispecies.out]
    otherstates <- thermo$OBIGT$state[ispecies.other]
    # for minerals (cr), use the word "phase"; otherwise, use "state" 20190209
    word <- "state"
    # substitute the mineral name for "cr" 20190121
    if(mystate == "cr" | sum(otherstates=="cr") > 1) {
      word <- "phase"
      otherstates[otherstates=="cr"] <- thermo$OBIGT$name[ispecies.other[otherstates=="cr"]]
    }
    transtext <- othertext <- ""
    # we count, but don't show the states for phase transitions (cr2, cr3, etc)
    istrans <- otherstates %in% c("cr2", "cr3", "cr4", "cr5", "cr6", "cr7", "cr8", "cr9")
    if(mystate=="cr") {
      # if we are "cr" we show the number of phase transitions
      ntrans <- sum(istrans)
      if(ntrans == 1) transtext <- paste(" with", ntrans, "phase transition")
      else if(ntrans > 1) transtext <- paste(" with", ntrans, "phase transitions")
      # if it's not already in the species name, substitute the mineral name for "cr" 20190121
      if(species != thermo$OBIGT$name[ispecies.out]) mystate <- thermo$OBIGT$name[ispecies.out]
    }
    otherstates <- otherstates[!istrans]
    if(length(otherstates) == 1) othertext <- paste0("; other available ", word, " is ", otherstates)
    if(length(otherstates) > 1) othertext <- paste0("; other available ", word, "s are ", paste(otherstates, collapse=", "))
    if(transtext != "" | othertext != "") {
      starttext <- paste0("info.character: found ", species, "(", mystate, ")")
      message(starttext, transtext, othertext)
    }
  }
  return(ispecies.out)
}

info.numeric <- function(ispecies, check.it=TRUE) {
  # from a numeric species index in 'ispecies' return the 
  # thermodynamic properties and equations-of-state parameters
  thermo <- get("thermo", CHNOSZ)
  # if we're called with NA, return an empty row
  if(is.na(ispecies)) {
    this <- thermo$OBIGT[1,]
    this[] <- NA
    return(this)
  }
  this <- thermo$OBIGT[ispecies,]
  # species indices must be in range
  ispeciesmax <- nrow(thermo$OBIGT)
  if(ispecies > ispeciesmax | ispecies < 1) 
    stop(paste("species index", ispecies, "not found in thermo()$OBIGT\n"))
  # remove scaling factors on EOS parameters depending on state
  # use new OBIGT2eos function here
  this <- OBIGT2eos(this, this$state)
  # identify any missing GHS values
  naGHS <- is.na(this[9:11])
  # a missing one of G, H or S can cause problems for subcrt calculations at high T
  if(sum(naGHS)==1) {
    # calculate a single missing one of G, H, or S from the others
    GHS <- as.numeric(GHS(as.character(this$formula), G=this[,9], H=this[,10], S=this[,11], E_units=this$E_units))
    message("info.numeric: ", colnames(this)[9:11][naGHS], " of ",
      this$name, "(", this$state, ") is NA; set to ", round(GHS[naGHS],2), " ", this$E_units, " mol-1")
    this[, which(naGHS)+8] <- GHS[naGHS]
  } 
  # now perform consistency checks for GHS and EOS parameters if check.it=TRUE
  # don't do it for the AkDi species 20190219
  if(check.it & !"xi" %in% colnames(this)) {
    # check GHS if they are all present
    if(sum(naGHS)==0) calcG <- checkGHS(this)
    # check tabulated heat capacities against EOS parameters
    calcCp <- checkEOS(this, this$state, "Cp")
    # fill in NA heat capacity
    if(!is.na(calcCp) & is.na(this$Cp)) {
      message("info.numeric: Cp of ", this$name, "(", this$state, ") is NA; set by EOS parameters to ", round(calcCp, 2), " ", this$E_units, " K-1 mol-1")
      this$Cp <- as.numeric(calcCp)
    }
    # check tabulated volumes - only for aq (HKF equation)
    if(identical(this$state, "aq")) {
      calcV <- checkEOS(this, this$state, "V")
      # fill in NA volume
      if(!is.na(calcV) & is.na(this$V)) {
        message("info.numeric: V of ", this$name, "(", this$state, ") is NA; set by EOS parameters to ", round(calcV, 2), " cm3 mol-1")
        this$V <- as.numeric(calcV)
      }
    }
  } # done checking
  # all done!
  return(this)
}

info.approx <- function(species, state=NULL) {
  # returns species indices that have an approximate match of 'species'
  # to thermo$OBIGT$[name|abbrv|formula], 
  # possibly restricted to a given state
  thermo <- get("thermo", CHNOSZ)
  if(!is.null(state)) this <- thermo$OBIGT[thermo$OBIGT$state==state, ]
  else this <- thermo$OBIGT
  # only look for fairly close matches
  max.distance <- 0.1
  approx.name <- agrep(species, this$name, max.distance)
  approx.abbrv <- agrep(species, this$abbrv, max.distance)
  approx.formula <- agrep(species, this$formula, max.distance)
  approx.species <- unique(c(approx.name, approx.abbrv, approx.formula))
  if(!is.na(approx.species[1])) {
    # show the names of the species
    if(length(approx.species)==1) {
      message("info.approx: '", species, "' is similar to ", info.text(approx.species))
    } else {
      napprox.max <- 100
      exttext <- ":"
      if(length(approx.species) > napprox.max) exttext <- paste(" (showing first ", napprox.max, ")", sep="")
      message("info.approx: '", species, "' is ambiguous; has approximate matches to ", 
        length(approx.species), " species", exttext)
      printout <- capture.output(print(thermo$OBIGT$name[head(approx.species, napprox.max)]))
      message(paste(printout, collapse="\n"))
    }
    return(approx.species)
  }
  # if we got here there were no approximate matches
  # 20190127 look for the species in optional data files 
  for(opt in c("SLOP98", "SUPCRT92", "OldAA", "AkDi")) {
    optdat <- read.csv(system.file(paste0("extdata/OBIGT/", opt, ".csv"), package="CHNOSZ"), as.is=TRUE)
    if(species %in% optdat$name) {
      message('info.approx: ', species, ' is in an optional database; use add.OBIGT("', opt, '", "', species, '") to load it')
      return(NA)
    }
  }
  message("info.approx: '", species, "' has no approximate matches")
  return(NA)
}
