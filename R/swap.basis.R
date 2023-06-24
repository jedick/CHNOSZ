# CHNOSZ/swap.basis.R
# Functions related to swapping basis species
# Extracted from basis() 20120114 jmd

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("basis.R")

# Return the current basis elements
basis.elements <- function(basis = thermo()$basis) {
  if(is.null(basis)) stop("basis species are not defined")
  return(as.matrix(basis[, 1:nrow(basis), drop=FALSE]))
}

# Calculate chemical potentials of elements from logarithms of activity of basis species
element.mu <- function(basis = thermo()$basis, T = 25) {
  # Matrix part of the basis definition
  basis.mat <- basis.elements(basis)
  # The standard Gibbs energies of the basis species
  # Don't take it from thermo()$OBIGT, even at 25 degC, because G for H2O is NA there
  # the sapply(..., "[", 1) is needed to get the first value, in case subcrt appends a polymorph column (i.e. for S(cr))  20171105
  TK <- convert(T, "K")
  G <- unlist(sapply(subcrt(basis$ispecies, T = TK, property="G", convert = FALSE)$out, "[", 1))
  # Chemical potentials of the basis species
  species.mu <- G - convert(basis$logact, "G", T = TK)
  # Chemical potentials of the elements
  element.mu <- solve(basis.mat, species.mu)
  # Give them useful names
  names(element.mu) <- colnames(basis.mat)
  return(element.mu)
}

# Calculate logarithms of activity of basis species from chemical potentials of elements
basis.logact <- function(emu, basis = thermo()$basis, T = 25) {
  # Matrix part of the basis definition
  basis.mat <- basis.elements(basis)
  # Elements in emu can't be less than the number in the basis
  if(length(emu) < ncol(basis.mat)) stop("number of elements in 'emu' is less than those in basis")
  # Sort names of emu in order of those in basis.mat
  ielem <- match(names(emu), colnames(basis.mat))
  # Check that elements of basis.mat and emu are identical
  if(any(is.na(ielem))) stop(paste("element(s)", paste(names(emu)[is.na(ielem)], collapse = " "), "not found in basis"))
  # The standard Gibbs energies of the basis species
  # Don't take it from thermo()$OBIGT, even at 25 degC, because G for H2O is NA there
  # The sapply(..., "[", 1) is needed to get the first value, in case subcrt appends a polymorph column (i.e. for S(cr))  20171105
  TK <- convert(T, "K")
  G <- unlist(sapply(subcrt(basis$ispecies, T = TK, property = "G", convert = FALSE)$out, "[", 1))
  # The chemical potentials of the basis species in equilibrium
  # with the chemical potentials of the elements
  basis.mu <- colSums((t(basis.mat)*emu)) - G
  # Convert chemical potentials to logarithms of activity
  basis.logact <- -convert(basis.mu, "logK", T = TK)
  # Give them useful names
  names(basis.logact) <- rownames(basis.mat)
  return(basis.logact)
}

ibasis <- function(species) {
  # Get the index of a basis species from a species index, name or formula
  basis <- basis()
  if(is.numeric(species)) ib <- match(species, basis$ispecies)
  else {
    # Character: first look for formula of basis species
    ib <- match(species, rownames(basis))
    # If that doesn't work, look for name of basis species
    if(is.na(ib)) ib <- match(species, get("thermo", CHNOSZ)$OBIGT$name[basis$ispecies])
  }
  return(ib)
}

# Swap in one basis species for another
swap.basis <- function(species, species2, T = 25) {
  # Before we do anything, remember the old basis and species definitions
  oldbasis <- get("thermo", CHNOSZ)$basis
  ts <- species()
  if(is.null(oldbasis)) 
    stop("swapping basis species requires an existing basis definition")
  # Both arguments must have length 1
  if(missing(species) | missing(species2))
    stop("two species must be identified")
  if(length(species) > 1 | length(species2) > 2)
    stop("can only swap one species for one species")
  # Replace pH with H+ and pe and Eh with e- 20170205
  if(species == "pH") species <- "H+"
  if(species %in% c("Eh", "pe")) species <- "e-"
  if(species2 == "pH") species2 <- "H+"
  if(species2 %in% c("Eh", "pe")) species2 <- "e-"
  # Arguments are good, now find the basis species to swap out
  ib <- ibasis(species)
  if(is.na(ib)) stop(paste("basis species '",species,"' is not defined",sep = ""))
  # Find species2 in the thermodynamic database
  if(is.numeric(species2)) ispecies2 <- species2
  else ispecies2 <- suppressMessages(info(species2))
  if(is.na(ispecies2) | is.list(ispecies2))
    stop(paste("a species matching '",species2,"' is not available in thermo()$OBIGT",sep = ""))
  # Try to load the new basis species
  ispecies <- oldbasis$ispecies
  ispecies[ib] <- ispecies2
  newbasis <- put.basis(ispecies)$basis
  # Now deal with the activities
  if(!all(can.be.numeric(oldbasis$logact))) {
    # If there are any buffers, just set the old activities
    bl <- oldbasis$logact
  } else {
    # No buffers, so we can recalculate activities to maintain the chemical potentials of the elements
    # What were the original chemical potentials of the elements?
    emu <- element.mu(oldbasis, T = T)
    # The corresponding logarithms of activities of the new basis species
    bl <- basis.logact(emu, newbasis, T = T)
  }
  # Update the basis with these logacts
  mb <- mod.basis(ispecies, state = newbasis$state, logact = bl)
  # Delete, then restore species if they were defined
  species(delete = TRUE)
  if(!is.null(ts)) {
    suppressMessages(species(ts$ispecies))
    suppressMessages(species(1:nrow(get("thermo", CHNOSZ)$species), ts$logact))
  }
  # All done, return the new basis definition
  return(mb)
}
