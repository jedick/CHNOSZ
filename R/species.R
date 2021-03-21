# CHNOSZ/species.R
# define species of interest 

# Function to create or add to the species data frame
species <- function(species=NULL, state=NULL, delete=FALSE, add=FALSE, index.return=FALSE) {
  # 20080925 default quiet=TRUE 20101003 default quiet=FALSE
  # 20120128 remove 'quiet' argument (messages can be hidden with suppressMessages())
  # 20120523 return thermo()$species instead of rownumbers therein, and remove message showing thermo()$species
  # 20200714 add 'add' argument (by default, previous species are deleted)
  thermo <- get("thermo", CHNOSZ)
  ## argument processing
  # we can't deal with NA species
  if(identical(species, NA)) {
    cn <- caller.name()
    if(length(cn) > 0) ctext <- paste("(calling function was ", cn, ")", sep="") else ctext <- ""
    stop(paste("'species' is NA", ctext))
  }
  # delete the entire species definition or only selected species
  if(delete) {
    # delete the entire definition if requested
    if(is.null(species)) {
      thermo$species <- NULL
      assign("thermo", thermo, CHNOSZ)
      return(thermo$species)
    }
    # from here we're trying to delete already defined species
    if(is.null(thermo$species)) stop("nonexistent species definition")
    # match species to delete by name and species number
    isp <- rep(NA, length(species))
    ispname <- match(species, thermo$species$name)
    ispnum <- match(species, 1:nrow(thermo$species))
    isp[!is.na(ispname)] <- ispname[!is.na(ispname)]
    isp[!is.na(ispnum)] <- ispnum[!is.na(ispnum)]
    # filter out non-matching species
    ina <- is.na(isp)
    if(any(ina)) warning(paste("species:",
      paste(species[ina], collapse=" "), "not present, so can not be deleted"))
    isp <- isp[!ina]
    # go on to delete this/these species
    if(length(isp) > 0) {
      thermo$species <- thermo$species[-isp,]
      if(nrow(thermo$species)==0) thermo$species <- NULL
      else rownames(thermo$species) <- 1:nrow(thermo$species)
      assign("thermo", thermo, CHNOSZ)
    }
    return(thermo$species)
  }
  # if no species or states are given, just return the species list
  if(is.null(species) & is.null(state)) return(thermo$species)
#  # delete the previous species definition unless add = TRUE 20200714
#  if(!add) {
#    thermo$species <- NULL
#    assign("thermo", thermo, CHNOSZ)
#  }
  # if no species are given use all of them if available
  if(is.null(species) & !is.null(thermo$species)) species <- 1:nrow(thermo$species)
  # parse state argument
  state <- state.args(state)
  # make species and state arguments the same length
  if(length(species) > length(state) & !is.null(state)) state <- rep(state,length.out=length(species)) else 
  if(length(state) > length(species) & !is.null(species)) species <- rep(species,length.out=length(state))
  # rename state as logact if it's numeric, or get default values of logact
  if(is.numeric(state[1])) {
    logact <- state
    state <- NULL
  } else logact <- NULL
  # if they don't look like states (aq,gas,cr) or activities (numeric), 
  # use them as a suffix for species name (i.e., a protein_organism)
  allstates <- unique(thermo$OBIGT$state)
  if( sum(state %in% allstates) < length(state) & !can.be.numeric(state[1]) & !can.be.numeric(species[1]) ) {
    species <- paste(species, state, sep="_")
    state <- rep(thermo$opt$state, length.out=length(state))
  }
  # parse species argument
  iOBIGT <- NULL
  if(is.character(species[1])) {
    # look for named species in species definition
    ispecies <- match(species, thermo$species$name)
    # if all species names match, and logact is given, re-call the function with the species indices
    if(!any(is.na(ispecies)) & !is.null(logact)) return(species(ispecies, state=logact, index.return=index.return))
    # look for species in thermo()$OBIGT
    iOBIGT <- suppressMessages(info(species, state))
    # since that could have updated thermo()$OBIGT (with proteins), re-read thermo
    thermo <- get("thermo", CHNOSZ)
    # check if we got all the species
    ina <- is.na(iOBIGT)
    if(any(ina)) stop(paste("species not available:", paste(species[ina], collapse=" ")))
  } else {
    # if species is numeric and a small number it refers to the species definition,
    # else to species index in thermo()$OBIGT
    nspecies <- nrow(thermo$species)
    if(is.null(thermo$species)) nspecies <- 0
    if(max(species) > nspecies) iOBIGT <- species
  }
  ## done with argument processing ... now to do work
  # create or add to species definition
  if(!is.null(iOBIGT)) {
    if(is.null(thermo$basis)) stop("basis species are not defined")
    # the coefficients in reactions to form the species from basis species
    # wrap values in unname in case they have names from retrieve(), otherwise makeup() doesn't work as intended 20190225
    f <- (species.basis(unname(iOBIGT)))
    # the states and species names
    state <- as.character(thermo$OBIGT$state[iOBIGT])
    name <- as.character(thermo$OBIGT$name[iOBIGT])
    # get default values of logact
    if(is.null(logact)) {
      logact <- rep(0, length(species))
      logact[state=="aq"] <- -3
    }
    # create the new species
    newspecies <- data.frame(f, ispecies=iOBIGT, logact=logact, state=state, name=name, stringsAsFactors=FALSE)
    # "H2PO4-" looks better than "H2PO4."
    colnames(newspecies)[1:nrow(thermo$basis)] <- rownames(thermo$basis)
    # initialize or add to species data frame
    if(is.null(thermo$species) | !add) {
      thermo$species <- newspecies
      ispecies <- 1:nrow(thermo$species)
    } else {
      # don't add species that already exist
      idup <- newspecies$ispecies %in% thermo$species$ispecies
      thermo$species <- rbind(thermo$species, newspecies[!idup, ])
      ispecies <- match(newspecies$ispecies, thermo$species$ispecies)
    }
    rownames(thermo$species) <- seq(nrow(thermo$species))
  } else {
    # update activities or states of existing species
    # first get the rownumbers in thermo()$species
    if(is.numeric(species[1])) {
      ispecies <- species
      # if state and logact are both NULL we don't do anything but return the selected species
      if(is.null(state) & is.null(logact)) {
        if(index.return) return(ispecies)
        else return(thermo$species[ispecies, ])
      }
    } else ispecies <- match(species, thermo$species$name)
    # replace activities?
    if(!is.null(logact)) {
      thermo$species$logact[ispecies] <- logact
    } else {
      # change states, checking for availability of the desired state
      for(i in 1:length(ispecies)) {
        myform <- thermo$OBIGT$formula[thermo$species$ispecies[ispecies[i]]]
        #iOBIGT <- which(thermo$OBIGT$name==thermo$species$name[ispecies[k]] | thermo$OBIGT$formula==myform)
        # 20080925 don't match formula -- two proteins might have the
        # same formula (e.g. YLR367W and YJL190C)
        #iOBIGT <- which(thermo$OBIGT$name==thermo$species$name[ispecies[k]])
        # 20091112 do match formula if it's not a protein -- be able to 
        # change "carbon dioxide(g)" to "CO2(aq)"
        if(length(grep("_",thermo$species$name[ispecies[i]])) > 0)  
          iOBIGT <- which(thermo$OBIGT$name==thermo$species$name[ispecies[i]])
        else {
          iOBIGT <- which(thermo$OBIGT$name==thermo$species$name[ispecies[i]] & thermo$OBIGT$state==state[i])
          if(length(iOBIGT)==0)
            iOBIGT <- which(thermo$OBIGT$name==thermo$species$name[ispecies[i]] | thermo$OBIGT$formula==myform)
        }
        if(!state[i] %in% thermo$OBIGT$state[iOBIGT]) 
          warning(paste("can't update state of species", ispecies[i], "to", state[i], "\n"), call.=FALSE)
        else {
          ii <- match(state[i], thermo$OBIGT$state[iOBIGT])
          thermo$species$state[ispecies[i]] <- state[i]
          thermo$species$name[ispecies[i]] <- thermo$OBIGT$name[iOBIGT[ii]]
          thermo$species$ispecies[ispecies[i]] <- as.numeric(rownames(thermo$OBIGT)[iOBIGT[ii]])
        }
      }
    }
  }
  assign("thermo", thermo, CHNOSZ)
  # return the new species definition or the index(es) of affected species
  if(index.return) return(ispecies)
  else return(thermo$species)
}

### unexported functions ###

# to retrieve the coefficients of reactions to form the species from the basis species
species.basis <- function(species=get("thermo", CHNOSZ)$species$ispecies, mkp = NULL) {
  # current basis elements
  bmat <- basis.elements()
  tbmat <- t(bmat)
  # what are the elements?
  belem <- rownames(tbmat)
  # get the species makeup into a matrix
  # If 'species' is NULL, get it from the 'mkp' argument (for use by subcrt()) 20210321
  if(!is.null(species)) mkp <- as.matrix(sapply(makeup(species, count.zero=TRUE), c))
  # the positions of the species elements in the basis elements
  ielem <- match(rownames(mkp), belem)
  # the elements of the species must be contained by the basis species
  if(any(is.na(ielem))) stop(paste("element(s) not in the basis:", 
    paste(rownames(mkp)[is.na(ielem)], collapse=" ")))
  # the positions of the basis elements in the species elements
  jelem <- match(belem, rownames(mkp))
  # keep track of which ones are NA's; 
  # index them as one here but turn them to zero later
  ina <- is.na(jelem)
  jelem[ina] <- 1
  # now put the species matrix into the same order as the basis
  mkp <- mkp[jelem, , drop=FALSE]
  # fill zeros for any basis element not in the species
  mkp[ina, ] <- 0
  # solve for the basis coefficients and transpose
  nbasis <- t(apply(mkp, 2, function(x) solve(tbmat, x)))
  # very small numbers are probably a floating point artifact
  # can cause problems in situations where zeros are needed
  # (manifests as issue in longex("phosphate"), where which.balance()
  #  identifies H2O as conserved component)
  # 20140201 set digits (to R default) becuase getOption("digits") is changed in knitr
  out <- zapsmall(nbasis, digits=7)
  # add names of species and basis species
  colnames(out) <- colnames(tbmat)
  # add names of species only if it was a character argument
  if(all(is.character(species))) rownames(out) <- species
  return(out)
} 
