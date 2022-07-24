# CHNOSZ/mosaic.R
# calculate affinities with changing basis species
# 20141220 jmd initial version
# 20190129 complete rewrite to use any number of groups of changing basis species
#   and improve speed by pre-calculating subcrt values (sout)
# 20190505 bug fix: adjust affinities of species formation reactions for mole fractions of basis species

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("basis.R")
#source("util.character.R")
#source("util.args.R")

# Function to calculate affinities with mosaic of basis species
mosaic <- function(bases, blend = TRUE, stable = list(), loga_aq = NULL, ...) {

  # Argument recall 20190120
  # if the first argument is the result from a previous mosaic() calculation,
  # just update the remaining arguments
  if(is.list(bases)) {
    if(identical(bases[1], list(fun = "mosaic"))) {
      aargs <- bases$args
      # we can only update arguments given in ...
      ddd <- list(...)
      if(length(ddd) > 0) {
        for(i in 1:length(ddd)) {
          if(names(ddd)[i] %in% names(aargs)) aargs[[names(ddd)[i]]] <- ddd[[i]]
          else aargs <- c(aargs, ddd[i])
        }
      }
      return(do.call(mosaic, aargs))
    }
  }

  # Backward compatibility 20190131:
  # bases can be a vector instead of a list
  if(!is.list(bases)) {
    bases <- list(bases)
    otherargs <- list(...)
    allargs <- c(list(bases = bases, blend = blend, stable = stable), otherargs)
    out <- do.call(mosaic, allargs)
    # Replace A.bases (affinity calculations for all groups of basis species) with backwards-compatible A.bases
    out$A.bases <- out$A.bases[[1]]
    return(out)
  }

  if(length(stable) == 2) {
    # Use only predominant basis species for mosaic stacking 20220723
    stable2.orig <- stable2 <- stable[[2]]
    # Why c(1, ?
    # The first basis species should always be included (because it has to be swapped out for the others)
    istable2 <- sort(unique(c(1, as.numeric(stable2))))
    for(i in seq_along(istable2)) stable2[stable2.orig == istable2[i]] <- i
    bases2 <- bases[[2]][istable2]
    bases[[2]] <- bases2
    stable[[2]] <- stable2
  }

  # Save starting basis and species definition
  basis0 <- get("thermo", CHNOSZ)$basis
  species0 <- get("thermo", CHNOSZ)$species
  # Get species indices of requested basis species
  ispecies <- lapply(bases, info)
  if(any(is.na(unlist(ispecies)))) stop("one or more of the requested basis species is unavailable")
  # Identify starting basis species
  ispecies0 <- sapply(ispecies, "[", 1)
  ibasis0 <- match(ispecies0, basis0$ispecies)
  # Quit if starting basis species are not present
  ina <- is.na(ibasis0)
  if(any(ina)) {
    names0 <- unlist(lapply(bases, "[", 1))
    stop("the starting basis does not have ", paste(names0[ina], collapse = " and "))
  }

  ddd <- list(...)
  if("sout" %in% names(ddd)) {
    ddd_has_sout <- TRUE
    # Get sout from ... (from solubility()) 20210322
    sout <- ddd$sout
  } else {
    ddd_has_sout <- FALSE
    # Run subcrt() calculations for all basis species and formed species 20190131
    # this avoids repeating the calculations in different calls to affinity()
    # add all the basis species here - the formed species are already present
    lapply(bases, species, add = TRUE)
    sout <- affinity(..., return.sout = TRUE)
  }

  # Calculate affinities of the basis species themselves
  A.bases <- list()
  for(i in 1:length(bases)) {
    message("mosaic: calculating affinities of basis species group ", i, ": ", paste(bases[[i]], collapse=" "))
    mysp <- species(bases[[i]])
    # 20191111 include only aq species in total activity
    iaq <- mysp$state == "aq"
    # Use as.numeric in case a buffer is active 20201014
    if(any(iaq)) species(which(iaq), as.numeric(basis0$logact[ibasis0[i]]))
    if(ddd_has_sout) A.bases[[i]] <- suppressMessages(affinity(...))
    else A.bases[[i]] <- suppressMessages(affinity(..., sout = sout))
  }

  # Get all combinations of basis species
  newbases <- as.matrix(expand.grid(ispecies))
  allbases <- matrix(basis0$ispecies, nrow = 1)[rep(1, nrow(newbases)), , drop = FALSE]
  allbases[, ibasis0] <- newbases

  # Calculate affinities of species for all combinations of basis species
  aff.species <- list()
  message("mosaic: calculating affinities of species for all ", nrow(allbases), " combinations of the basis species")
  # Run backwards so that we put the starting basis species back at the end
  for(i in nrow(allbases):1) {
    # Get default loga from starting basis species
    thislogact <- basis0$logact
    # Use logact = 0 for solids 20191111
    states <- sout$species$state[match(allbases[i, ], sout$species$ispecies)]
    icr <- grepl("cr", states)
    thislogact[icr] <- 0
    # Use loga_aq for log(activity) of mosaiced aqueous basis species 20220722
    if(!is.null(loga_aq)) {
      if(length(loga_aq) != length(ibasis0)) stop("'loga_aq' should have same length as 'bases'")
      iaq <- grep("aq", states)
      # Loop over sets of mosaiced basis species
      for(j in 1:length(ibasis0)) {
        if(is.na(loga_aq[j])) next
        iaqbasis0 <- intersect(iaq, ibasis0[j])
        thislogact[iaqbasis0] <- loga_aq[j]
      }
    }
    put.basis(allbases[i, ], thislogact)
    # We have to define the species using the current basis
    species(species0$ispecies, species0$logact)
    if(ddd_has_sout) aff.species[[i]] <- suppressMessages(affinity(...))
    else aff.species[[i]] <- suppressMessages(affinity(..., sout = sout))
  }

  # Calculate equilibrium mole fractions for each group of basis species
  group.fraction <- list()
  blend <- rep(blend, length(A.bases))
  E.bases <- list()
  for(i in 1:length(A.bases)) {
    if(blend[i] & is.null(stable[i][[1]])) {
      # This isn't needed (and doesn't work) if all the affinities are NA 20180925
      if(any(!sapply(A.bases[[1]]$values, is.na))) {
        # 20190504: when equilibrating the changing basis species, use a total activity equal to the activity from the basis definition
        # 20191111 use equilibrate(loga.balance = ) instead of setting activities in species definition
        e <- equilibrate(A.bases[[i]], loga.balance = as.numeric(basis0$logact[ibasis0[i]]))
        # Exponentiate to get activities then divide by total activity
        a.equil <- lapply(e$loga.equil, function(x) 10^x)
        a.tot <- Reduce("+", a.equil)
        group.fraction[[i]] <- lapply(a.equil, function(x) x / a.tot)
        # Include the equilibrium activities in the output of this function 20190504
        E.bases[[i]] <- e
      } else {
        group.fraction[[i]] <- A.bases[[i]]$values
      }
    } else {
      # For blend = FALSE, we just look at whether
      # a basis species predominates within its group
      if(is.null(stable[i][[1]])) {
        d <- diagram(A.bases[[i]], plot.it = FALSE, limit.water = FALSE)
        predom <- d$predominant
      } else {
        # Get the stable species from the argument 20200715
        predom <- stable[i][[1]]
      }
      group.fraction[[i]] <- list()
      for(j in 1:length(bases[[i]])) {
        # If a basis species predominates, it has a mole fraction of 1, or 0 otherwise
        yesno <- predom
        yesno[yesno != j] <- 0
        yesno[yesno == j] <- 1
        group.fraction[[i]][[j]] <- yesno
      }
    }
  }

  # Make an indexing matrix for all combinations of basis species
  ind.mat <- list()
  for(i in 1:length(ispecies)) ind.mat[[i]] <- 1:length(ispecies[[i]])
  ind.mat <- as.matrix(expand.grid(ind.mat))

  # Loop over combinations of basis species
  for(icomb in 1:nrow(ind.mat)) {
    # Loop over groups of changing basis species
    for(igroup in 1:ncol(ind.mat)) {
      # Get mole fractions for this particular basis species
      basisx <- group.fraction[[igroup]][[ind.mat[icomb, igroup]]]
      # Loop over species
      for(jspecies in 1:length(aff.species[[icomb]]$values)) {
        # Get coefficient of this basis species in the formation reaction for this species
        nbasis <- aff.species[[icomb]]$species[jspecies, ibasis0[igroup]]
        # Adjust affinity of species for mole fractions (i.e. lower activity) of basis species 20190505
        aff.adjust <- nbasis * log10(basisx)
        # Avoid infinite values (from log10(0))
        isfin <- is.finite(aff.adjust)
        aff.species[[icomb]]$values[[jspecies]][isfin] <- aff.species[[icomb]]$values[[jspecies]][isfin] + aff.adjust[isfin]
      }
      # Multiply fractions of basis species from each group to get overall fraction
      if(igroup==1) groupx <- basisx
      else groupx <- groupx * basisx
    }
    # Multiply affinities by the mole fractions of basis species
    aff.species[[icomb]]$values <- lapply(aff.species[[icomb]]$values, function(values) values * groupx)
  }
  
  # Get total affinities for the species
  A.species <- aff.species[[1]]
  for(i in 1:length(A.species$values)) {
    # Extract the affinity contributions from each basis species
    A.values <- lapply(lapply(aff.species, "[[", "values"), "[[", i)
    # Sum them to get total affinities for this species
    A.species$values[[i]] <- Reduce("+", A.values)
  }

  # For argument recall, include all arguments in output 20190120
  allargs <- c(list(bases = bases, blend = blend), list(...))
  # Return the affinities for the species and basis species
  return(list(fun = "mosaic", args = allargs, A.species = A.species, A.bases = A.bases, E.bases = E.bases))
}
