# CHNOSZ/mosaic.R
# calculate affinities with changing basis species
# 20141220 jmd initial version
# 20190129 complete rewrite to use any number of groups of changing basis species
#   and improve speed by pre-calculating subcrt values (sout)
# 20190505 bug fix: adjust affinities of species formation reactions for mole fractions of basis species

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("basis.R")
#source("util.character.R")
#source("util.args.R")

# function to calculate affinities with mosaic of basis species
mosaic <- function(bases, bases2 = NULL, blend = TRUE, predominant = list(), ...) {

  # argument recall 20190120
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

  # backward compatibility 20190131:
  # bases can be a vector instead of a list
  # bases2 can be present
  if(!is.list(bases)) {
    bases <- list(bases)
    hasbases2 <- FALSE
    if(!is.null(bases2)) {
      bases <- c(bases, list(bases2))
      hasbases2 <- TRUE
    }
    otherargs <- list(...)
    allargs <- c(list(bases = bases, blend = blend, predominant = predominant), otherargs)
    out <- do.call(mosaic, allargs)
    # replace A.bases (affinity calculations for all groups of basis species) with backwards-compatbile A.bases and A.bases2
    if(hasbases2) A.bases2 <- out$A.bases[[2]]
    A.bases <- out$A.bases[[1]]
    out$A.bases <- A.bases
    if(hasbases2) out <- c(out, list(A.bases2 = A.bases2))
    return(out)
  }

  # save starting basis and species definition
  basis0 <- get("thermo", CHNOSZ)$basis
  species0 <- get("thermo", CHNOSZ)$species
  # get species indices of requested basis species
  ispecies <- lapply(bases, info)
  if(any(is.na(unlist(ispecies)))) stop("one or more of the requested basis species is unavailable")
  # identify starting basis species
  ispecies0 <- sapply(ispecies, "[", 1)
  ibasis0 <- match(ispecies0, basis0$ispecies)
  # quit if starting basis species are not present
  ina <- is.na(ibasis0)
  if(any(ina)) {
    names0 <- unlist(lapply(bases, "[", 1))
    stop("the starting basis does not have ", paste(names0[ina], collapse = " and "))
  }

  # run subcrt() calculations for all basis species and formed species 20190131
  # this avoids repeating the calculations in different calls to affinity()
  # add all the basis species here - the formed species are already present
  lapply(bases, species, add = TRUE)
  sout <- affinity(..., return.sout = TRUE)

  # calculate affinities of the basis species themselves
  A.bases <- list()
  for(i in 1:length(bases)) {
    message("mosaic: calculating affinities of basis species group ", i, ": ", paste(bases[[i]], collapse=" "))
    mysp <- species(bases[[i]])
    # 20191111 include only aq species in total activity
    iaq <- mysp$state == "aq"
    if(any(iaq)) species(which(iaq), basis0$logact[ibasis0[i]])
    A.bases[[i]] <- suppressMessages(affinity(..., sout = sout))
  }

  # get all combinations of basis species
  newbases <- as.matrix(expand.grid(ispecies))
  allbases <- matrix(basis0$ispecies, nrow = 1)[rep(1, nrow(newbases)), , drop = FALSE]
  allbases[, ibasis0] <- newbases

  # calculate affinities of species for all combinations of basis species
  aff.species <- list()
  message("mosaic: calculating affinities of species for all ", nrow(allbases), " combinations of the basis species")
  # run backwards so that we put the starting basis species back at the end
  for(i in nrow(allbases):1) {
    # use logact = 0 for solids 20191111
    thislogact <- basis0$logact
    states <- sout$species$state[match(allbases[i, ], sout$species$ispecies)]
    icr <- grepl("cr", states)
    thislogact[icr] <- 0
    put.basis(allbases[i, ], thislogact)
    # we have to define the species using the current basis
    species(species0$ispecies, species0$logact)
    aff.species[[i]] <- suppressMessages(affinity(..., sout = sout))
  }

  # calculate equilibrium mole fractions for each group of basis species
  group.fraction <- list()
  blend <- rep(blend, length(A.bases))
  E.bases <- list()
  for(i in 1:length(A.bases)) {
    if(blend[i] & is.null(predominant[i][[1]])) {
      # this isn't needed (and doesn't work) if all the affinities are NA 20180925
      if(any(!sapply(A.bases[[1]]$values, is.na))) {
        # 20190504: when equilibrating the changing basis species, use a total activity equal to the activity from the basis definition
        # 20191111 use equilibrate(loga.balance = ) instead of setting activities in species definition
        e <- equilibrate(A.bases[[i]], loga.balance = basis0$logact[ibasis0[i]])
        # exponentiate to get activities then divide by total activity
        a.equil <- lapply(e$loga.equil, function(x) 10^x)
        a.tot <- Reduce("+", a.equil)
        group.fraction[[i]] <- lapply(a.equil, function(x) x / a.tot)
        # include the equilibrium activities in the output of this function 20190504
        E.bases[[i]] <- e
      } else {
        group.fraction[[i]] <- A.bases[[i]]$values
      }
    } else {
      # for blend = FALSE, we just look at whether
      # a basis species predominates within its group
      if(is.null(predominant[i][[1]])) {
        d <- diagram(A.bases[[i]], plot.it = FALSE, limit.water = FALSE)
        predom <- d$predominant
      } else {
        # get the predominances from the argument 20200715
        predom <- predominant[i][[1]]
      }
      group.fraction[[i]] <- list()
      for(j in 1:length(bases[[i]])) {
        # if a basis species predominates, it has a mole fraction of 1, or 0 otherwise
        yesno <- predom
        yesno[yesno != j] <- 0
        yesno[yesno == j] <- 1
        group.fraction[[i]][[j]] <- yesno
      }
    }
  }

  # make an indexing matrix for all combinations of basis species
  ind.mat <- list()
  for(i in 1:length(ispecies)) ind.mat[[i]] <- 1:length(ispecies[[i]])
  ind.mat <- as.matrix(expand.grid(ind.mat))

  # loop over combinations of basis species
  for(icomb in 1:nrow(ind.mat)) {
    # loop over groups of changing basis species
    for(igroup in 1:ncol(ind.mat)) {
      # get mole fractions for this particular basis species
      basisx <- group.fraction[[igroup]][[ind.mat[icomb, igroup]]]
      # loop over species
      for(jspecies in 1:length(aff.species[[icomb]]$values)) {
        # get coefficient of this basis species in the formation reaction for this species
        nbasis <- aff.species[[icomb]]$species[jspecies, ibasis0[igroup]]
        # adjust affinity of species for mole fractions (i.e. lower activity) of basis species 20190505
        aff.adjust <- nbasis * log10(basisx)
        # avoid infinite values (from log10(0))
        isfin <- is.finite(aff.adjust)
        aff.species[[icomb]]$values[[jspecies]][isfin] <- aff.species[[icomb]]$values[[jspecies]][isfin] + aff.adjust[isfin]
      }
      # multiply fractions of basis species from each group to get overall fraction
      if(igroup==1) groupx <- basisx
      else groupx <- groupx * basisx
    }
    # multiply affinities by the mole fractions of basis species
    aff.species[[icomb]]$values <- lapply(aff.species[[icomb]]$values, function(values) values * groupx)
  }
  
  # get total affinities for the species
  A.species <- aff.species[[1]]
  for(i in 1:length(A.species$values)) {
    # extract the affinity contributions from each basis species
    A.values <- lapply(lapply(aff.species, "[[", "values"), "[[", i)
    # sum them to get total affinities for this species
    A.species$values[[i]] <- Reduce("+", A.values)
  }

  # for argument recall, include all arguments in output 20190120
  allargs <- c(list(bases = bases, blend = blend), list(...))
  # return the affinities for the species and basis species
  return(list(fun = "mosaic", args = allargs, A.species = A.species, A.bases = A.bases, E.bases = E.bases))
}
