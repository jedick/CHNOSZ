# CHNOSZ/basis.R
# set up the basis species of a thermodynamic system

basis <- function(species=NULL, state=NULL, logact=NULL, delete=FALSE) {
  thermo <- get("thermo")
  oldbasis <- thermo$basis
  ## delete the basis species if requested
  if(delete | identical(species, "")) {
    thermo$basis <- NULL
    thermo$species <- NULL
    assign("thermo", thermo, "CHNOSZ")
    return(invisible(oldbasis))
  }
  ## return the basis definition if requested
  if(is.null(species)) return(oldbasis)
  ## from now on we need something to work with
  if(length(species)==0) stop("species argument is empty")
  # is the species one of the preset keywords?
  if(species[1] %in% preset.basis()) return(preset.basis(species[1]))
  # the species names/formulas have to be unique
  if(!length(unique(species))==length(species)) stop("species names are not unique")
  ## processing 'state' and 'logact' arguments
  # they should be same length as species
  if(!is.null(state)) state <- rep(state, length.out=length(species))
  if(!is.null(logact)) logact <- rep(logact, length.out=length(species))
  # results should be identical for
  # basis(c('H2O','CO2','H2'), rep('aq',3), c(0,-3,-3))
  # basis(c('H2O','CO2','H2'), c(0,-3,-3), rep('aq',3))
  # first of all, do we have a third argument?
  if(!is.null(logact)) {
    # does the 3rd argument look like states?
    if(is.character(logact[1])) {
      # swap the arguments into their correct places
      tmp <- logact
      logact <- state
      state <- tmp
    }
  } else {
    # if the second argument is numeric, treat it like logacts
    if(is.numeric(state[1])) {
      logact <- state
      state <- NULL
    }
  }
  ## processing 'species' argument
  # pH transformation
  if("pH" %in% species) {
    logact[species=="pH"] <- -logact[species=="pH"]
    if(!is.null(logact)) species[species=="pH"] <- "H+"
  }
  # Eh and pe transformations
  if("pe" %in% species) {
    logact[species=="pe"] <- -logact[species=="pe"]
    if(!is.null(logact)) species[species=="pe"] <- "e-"
  }
  if("Eh" %in% species) {
    # 20090209 should be careful with this conversion as it's only for 25 deg C
    # to be sure, just don't call species("Eh")
    if(!is.null(logact)) logact[species=="Eh"] <- -convert(logact[species=="Eh"],"pe")
    species[species=="Eh"] <- "e-"
  }
  ## if all species are in the existing basis definition, 
  ## *and* at least one of state or logact is not NULL
  ## modify the states and/or logacts of the existing basis species
  if(all(species %in% rownames(oldbasis)) | all(species %in% oldbasis$ispecies)) 
    if(!is.null(state) | !is.null(logact))
      return(mod.basis(species, state, logact))
  ## we're on to making a new basis definition
  # use default logacts if they aren't present
  if(is.null(logact)) logact <- rep(0, length(species))
  # if species argument is numeric, it's species indices
  if(is.numeric(species[1])) {
    ispecies <- species
    ina <- ispecies > nrow(thermo$obigt)
  } else {
    # get species indices using states from the argument, or default states
    if(!is.null(state)) ispecies <- suppressMessages(info(species, state, check.it=FALSE))
    else ispecies <- suppressMessages(info(species, check.it=FALSE))
    # check if we got all the species
    ina <- is.na(ispecies)
    # info() returns a list if any of the species had multiple approximate matches
    # we don't accept any of those
    if(is.list(ispecies)) ina <- ina | sapply(ispecies,length) > 1
  }
  if(any(ina)) stop(paste("species not available:", paste(species[ina], "(", state[ina], ")", sep="", collapse=" ")))
  # load new basis species
  return(put.basis(ispecies, logact))
}

### unexported functions ###

# to add the basis to thermo$obigt
put.basis <- function(ispecies, logact = rep(NA, length(ispecies))) {
  thermo <- get("thermo")
  state <- thermo$obigt$state[ispecies]
  # make the basis matrix, revised 20120114
  # get the elemental makeup of each species,
  # counting zero for any element that only appears in other species in the set
  comp <- makeup(ispecies, count.zero=TRUE)
  # turn the list into a matrix
  comp <- sapply(comp, c)
  # transpose to get put basis species on the rows
  comp <- t(comp)
  # note, makeup(count.zero=TRUE) above gave elements (colnames) sorted alphabetically
  # rownames identify the species
  rownames(comp) <- as.character(thermo$obigt$formula[ispecies])
  # FIXME: the electron doesn't look like a chemical formula
  # this is needed for affinity() to understand a 'pe' or 'Eh' variable
  if("(Z-1)" %in% rownames(comp)) rownames(comp)[rownames(comp)=="(Z-1)"] <- "e-"
  # now check it for validity of basis species
  # the first test: matrix is square
  if( nrow(comp) > ncol(comp) ) {
    if("Z" %in% colnames(comp)) stop("the number of basis species is greater than the number of elements and charge")
    else stop("the number of basis species is greater than the number of elements")
  }
  if( nrow(comp) < ncol(comp) ) {
    if("Z" %in% colnames(comp)) stop("the number of basis species is less than the number of elements and charge")
    else stop("the number of basis species is less than the number of elements")
  }
  # the second test: matrix is invertible
  if(class(try(solve(comp), silent=TRUE))=='try-error') 
    stop("singular stoichiometric matrix")
  # store the basis definition in thermo$basis, including
  # both numeric and character data, so we need to use a data frame
  comp <- cbind(as.data.frame(comp), ispecies, logact, state, stringsAsFactors=FALSE)
  # ready to assign to the global thermo object
  thermo$basis <- comp
  assign("thermo", thermo, "CHNOSZ")
  # remove the species since there's no guarantee the
  # new basis includes all their elements
  species(delete=TRUE)
  return(thermo$basis)
}

# modify the states or logact values in the existing basis definition
mod.basis <- function(species, state=NULL, logact=NULL) {
  thermo <- get("thermo")
  # the basis must be defined
  if(is.null(thermo$basis)) stop("basis is not defined")
  # loop over each species to modify
  for(i in 1:length(species)) {
    # which basis species are we looking at?
    if(is.numeric(species)) {
      ib <- match(species[i], thermo$basis$ispecies)
      if(is.na(ib)) stop(paste(species[i],"is not a species index of one of the basis species"))
    } else {
      ib <- match(species[i], rownames(thermo$basis))
      if(is.na(ib)) stop(paste(species[i],"is not a formula of one of the basis species"))
    }
    # first modify the state
    if(!is.null(state)) {
      if(state[i] %in% thermo$buffer$name) {
        # this is the name of a buffer
        ibuff <- which(as.character(thermo$buffers$name)==state[i])
        # check that each species in the buffer is compositionally
        # compatible with the basis definition
        for(k in 1:length(ibuff)) {
          ispecies <- suppressMessages(info(as.character(thermo$buffers$species)[ibuff[k]],
            as.character(thermo$buffers$state)[ibuff[k]], check.it=FALSE))
          bufmakeup <- makeup(ispecies)
          inbasis <- names(bufmakeup) %in% colnames(basis()) 
          if(FALSE %in% inbasis) {
            stop(paste("the elements '",c2s(names(bufmakeup)[!inbasis]),
              "' of species '",thermo$buffers$species[ibuff[k]],"' in buffer '",state[i],
              "' are not in the basis\n",sep=""))
          }
        }
        thermo$basis$logact[ib] <- state[i]
      } else {
        # first, look for a species with the same _name_ in the requested state
        myname <- thermo$obigt$name[thermo$basis$ispecies[ib]]
        ispecies <- suppressMessages(info(myname, state[i], check.it=FALSE))
        if(is.na(ispecies) | is.list(ispecies)) {
          # if that failed, look for a species with the same _formula_ in the requested state
          myformula <- thermo$obigt$formula[thermo$basis$ispecies[ib]]
          ispecies <- suppressMessages(info(myformula, state[i], check.it=FALSE))
          if(is.na(ispecies) | is.list(ispecies)) {
            # if that failed, we're out of luck
            if(myname==myformula) nametxt <- myname else nametxt <- paste(myname, "or", myformula)
            stop(paste0("state or buffer '", state[i], "' not found for ", nametxt, "\n"))
          }
        }
        thermo$basis$ispecies[ib] <- ispecies
        thermo$basis$state[ib] <- state[i]
      }
    } 
    # then modify the logact
    if(!is.null(logact)) {
      # allow this to be non-numeric in case we're called by swap.basis() while a buffer is active  20181109
      if(can.be.numeric(logact[i])) thermo$basis$logact[ib] <- as.numeric(logact[i])
      else thermo$basis$logact[ib] <- logact[i]
    }
    # assign the result to the CHNOSZ environment
    assign("thermo", thermo, "CHNOSZ")
  }
  return(thermo$basis)
} 

# to load a preset basis definition by keyword
preset.basis <- function(key=NULL) {
  # the available keywords
  basis.key <- c("CHNOS", "CHNOS+", "CHNOSe", "CHNOPS+", "CHNOPSe", "MgCHNOPS+", "MgCHNOPSe", "FeCHNOS", "FeCHNOS+", "QEC4", "QEC", "QEC+")
  # just list the keywords if none is specified
  if(is.null(key)) return(basis.key)
  # delete any previous basis definition
  basis("")
  # match the keyword to the available ones
  ibase <- match(key, basis.key)
  if(is.na(ibase)) stop(paste(key, "is not a keyword for preset basis species"))
  if(ibase==1) species <- c("CO2", "H2O", "NH3", "H2S", "oxygen")
  else if(ibase==2) species <- c("CO2", "H2O", "NH3", "H2S", "oxygen", "H+")
  else if(ibase==3) species <- c("CO2", "H2O", "NH3", "H2S", "e-", "H+")
  else if(ibase==4) species <- c("CO2", "H2O", "NH3", "H3PO4", "H2S", "oxygen", "H+")
  else if(ibase==5) species <- c("CO2", "H2O", "NH3", "H3PO4", "H2S", "e-", "H+")
  else if(ibase==6) species <- c("Mg+2", "CO2", "H2O", "NH3", "H3PO4", "H2S", "oxygen", "H+")
  else if(ibase==7) species <- c("Mg+2", "CO2", "H2O", "NH3", "H3PO4", "H2S", "e-", "H+")
  else if(ibase==8) species <- c("Fe2O3", "CO2", "H2O", "NH3", "H2S", "oxygen")
  else if(ibase==9) species <- c("Fe2O3", "CO2", "H2O", "NH3", "H2S", "oxygen", "H+")
  else if(ibase %in% c(10, 11)) species <- c("glutamine", "glutamic acid", "cysteine", "H2O", "oxygen")
  else if(ibase==12) species <- c("glutamine", "glutamic acid", "cysteine", "H2O", "oxygen", "H+")
  # get the preset logact
  logact <- preset.logact(species)
  # for QEC4, we use logact = -4 for the amino acids
  if(key=="QEC4") logact[1:3] <- -4
  # load the species and return the result
  return(basis(species, logact))
}

# logarithms of activities for preset basis definitions
preset.logact <- function(species) {
  bases <- c("H2O", "CO2", "NH3", "H2S", "oxygen", "H+", "e-", "Fe2O3",
             "glutamine", "glutamic acid", "cysteine")
  # values for QEC amino acids from Dick, 2017 (http://doi.org/10.1101/097667)
  logact <- c(0, -3, -4, -7, -80, -7, -7, 0,
              -3.2, -4.5, -3.6)
  ibase <- match(species, bases)
  logact <- logact[ibase]
  # any unmatched species gets a logarithm of activity of -3
  logact[is.na(logact)] <- -3
  return(logact)
}

