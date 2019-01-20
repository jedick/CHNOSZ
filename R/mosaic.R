# CHNOSZ/mosaic.R
# calculate affinities with changing basis species
# 20141220 jmd

# function to calculate affinities with mosaic of basis species
mosaic <- function(bases, bases2=NULL, blend=FALSE, ...) {

  # argument recall 20190120
  # if the first argument is the result from a previous mosaic() calculation,
  # just update the remaining arguments
  if(is.list(bases)) {
    if(identical(bases[1], list(fun="mosaic"))) {
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

  if(is.null(bases2)) {
    # the arguments for affinity()
    myargs <- list(...)
  } else {
    # the arguments for affinity (first set of basis species; outer loop)
    myargs1 <- list(...)
    # the arguments for mosaic() (second set of basis species; inner loop)
    myargs <- list(bases=bases2, blend=blend, ...)
  }

  # are the swapped basis species on the plot?
  # (the first one should be present in the starting basis set)
  iswap <- match(bases[1], names(myargs))
  # the log activity of the starting basis species
  logact.swap <- basis()$logact[ibasis(bases[1])]

  # a list where we'll keep the affinity calculations
  affs <- list()
  for(i in seq_along(bases)) {
    message(paste("mosaic: current basis species is", bases[i], sep=" "))
    # set up argument list: name of swapped-in basis species
    if(!is.na(iswap)) names(myargs)[iswap] <- bases[i]
    # calculate affinities
    if(is.null(bases2)) {
      affs[[i]] <- do.call(affinity, myargs)
    } else {
      mcall <- do.call(mosaic, myargs)
      affs[[i]] <- mcall$A.species
      A.bases2 <- mcall$A.bases
    }
    # change the basis species; restore the original at the end of the loop
    if(can.be.numeric(logact.swap)) logact.swap <- as.numeric(logact.swap)
    if(i < length(bases)) {
      swap.basis(bases[i], bases[i+1]) 
      # TODO: basis() requires the formula to identify the basis species,
      # would be nicer to just use the ibasis here
      bformula <- rownames(basis())[ibasis(bases[i+1])]
      basis(bformula, logact.swap)
    } else {
      swap.basis(bases[i], bases[1])
      bformula <- rownames(basis())[ibasis(bases[1])]
      basis(bformula, logact.swap)
    }
  }

  # calculate affinities of formation of basis species
  message(paste("mosaic: combining diagrams for", paste(bases, collapse=" "), sep=" "))
  ispecies <- species()$ispecies
  species.logact <- species()$logact
  species(delete=TRUE)
  species(bases)
  if(is.null(bases2)) A.bases <- do.call(affinity, myargs)
  else A.bases <- do.call(affinity, myargs1)
  # restore original species with original activities
  species(delete=TRUE)
  species(ispecies, species.logact)

  # affinities calculated using the first basis species
  A.species <- affs[[1]]
  if(blend) {
    # calculate affinities using relative abundances of basis species
    # this isn't needed (and doesn't work) if all the affinities are NA 20180925
    if(any(!sapply(A.species$values, is.na))) {
      e <- equilibrate(A.bases)
      # what is the total activity of the basis species?
      a.tot <- Reduce("+", lapply(e$loga.equil, function(x) 10^x))
      for(j in seq_along(affs)) {
        for(i in seq_along(A.species$values)) {
          # start with zero affinity
          if(j==1) A.species$values[[i]][] <- 0
          # add affinity scaled by __relative__ abundance of this basis species
          A.species$values[[i]] <- A.species$values[[i]] + affs[[j]]$values[[i]] * 10^e$loga.equil[[j]]/a.tot
        }
      }
    }
  } else {
    # use affinities from the single predominant basis species
    d <- diagram(A.bases, plot.it=FALSE)
    # merge affinities using the second, third, ... basis species
    for(j in tail(seq_along(affs), -1)) {
      is.predominant <- d$predominant==j
      # diagram() produces NA beyond water limits on Eh-pH diagrams (but we can't use NA for indexing, below)
      is.predominant[is.na(is.predominant)] <- FALSE
      for(i in seq_along(A.species$values)) {
        A.species$values[[i]][is.predominant] <- affs[[j]]$values[[i]][is.predominant]
      }
    }
  }

  # for argument recall, include all arguments in output 20190120
  allargs <- c(list(bases=bases, bases2=bases2, blend=blend), list(...))
  # return the affinities for the species and basis species
  if(is.null(bases2)) return(list(fun="mosaic", args=allargs, A.species=A.species, A.bases=A.bases))
  else return(list(fun="mosaic", args=allargs, A.species=A.species, A.bases=A.bases, A.bases2=A.bases2))
}
