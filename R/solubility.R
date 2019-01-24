# solubility.R: vectorized solubility calculations without uniroot
# 20181031 jmd
# 20181106 work on the output from affinity(); no "equilibrate()" needed!
# 20190117 add find.IS and test for dissociation reaction

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("equilibrate.R")
#source("util.misc.R")
#source("species.R")

solubility <- function(aout, dissociation=NULL, find.IS=FALSE, in.terms.of=NULL) {
  ## concept: the logarithms of activities of species at equilibrium are equal to
  ## Astar, the affinities calculated for unit activities of species

  ## is aout the output from mosaic() instead of affinity()?
  aout.save <- aout
  thisfun <- aout$fun
  if(thisfun=="mosaic") aout <- aout$A.species

  ## does the system involve a dissociation reaction?
  if(is.null(dissociation)) {
    # assume FALSE unless determined otherwise
    dissociation <- FALSE
    # if the reaction to form the first species involves the second basis species, we consider it to be a dissociation reaction
    if(aout$species[1, 2] != 0) {
      # 20190123 (corundum calculation): if there are only H2O, H+, and e-
      # besides the first basis species, it's not a dissociation reaction
      nbasis <- nrow(aout$basis)
      nH2O <- sum(rownames(aout$basis) %in% c("H2O", "H+", "e-", "O2", "H2"))
      if(nbasis > (nH2O + 1)) {
        # if we got here, and the second basis species is H2O, H+, e-, O2 (or others?), we don't have enough information, so stop
        if(colnames(aout$species)[2] %in% c("H2O", "H+", "e-", "O2", "H2")) {
          stop("Unsure whether the first formation reaction is a dissociation reaction.\nSet the 'dissociation' argument to TRUE or FALSE, or redefine the basis to put a product ion second.")
        }
        dissociation <- TRUE
      }
    }
    message("solubility: test for dissociation reaction returns ", dissociation)
  } else message("solubility: argument for dissociation reaction is ", dissociation)

  # get starting ionic strength (probably zero, but could be anything set by user)
  IS <- aout$IS
  if(is.null(IS)) IS <- 0

  ## to output loga.balance we need the balancing coefficients
  bout <- balance(aout)
  n.balance <- bout$n.balance
  balance <- bout$balance
  # get logarithm of total activity of the conserved basis species
  logabfun <- function(loga.equil, n.balance) {
    # exponentiate, multiply by n.balance, sum, logarithm
    a.equil <- mapply("^", 10, loga.equil, SIMPLIFY = FALSE)
    a.balance <- mapply("*", a.equil, n.balance, SIMPLIFY = FALSE)
    a.balance <- Reduce("+", a.balance)
    log10(a.balance)
  }
  
  # for find.IS=TRUE, iterate to converge on ionic strength
  niter <- 1
  while(TRUE) {
    ## the values in aout can be calculated for other than
    ## unit activities of species, so we have to take away the activites
    Astar <- function(i) aout$values[[i]] + aout$species$logact[i]
    loga.equil <- lapply(1:length(aout$values), Astar)

    ## for a dissociation on a *per reaction* (not system) basis,
    ## apply the divisor here and skip the if(dissociation){} part below
    ## (can be used to reproduce Fig. 4 of Manning et al., 2013)
    if(is.numeric(dissociation)) loga.equil <- lapply(loga.equil, "/", dissociation)

    loga.balance <- logabfun(loga.equil, bout$n.balance)

    # recalculate things for a dissociation reaction (like CaCO3 = Ca+2 + CO3+2)
    if(isTRUE(dissociation)) {
      # the multiplicity becomes the exponent in the reaction quotient
      loga.split <- loga.balance / 2
      # the contribution to affinity
      Asplit <- lapply(n.balance, "*", loga.split)
      # adjust the affinity and get new equilibrium activities
      aout$values <- mapply("-", aout$values, Asplit, SIMPLIFY=FALSE)
      loga.equil <- lapply(1:length(aout$values), Astar)
      # check that the new loga.balance == loga.split
      # TODO: does this work for non-1:1 species?
      loga.balance <- logabfun(loga.equil, n.balance)
      stopifnot(all.equal(loga.balance, loga.split))
    }

    # get the old IS
    IS.old <- rep(IS, length.out = length(aout$values[[1]]))
    # calculate the ionic strength for unpaired ions: IS = sum(molality * Z^2) / 2
    sum.mZ2 <- 0
    ion.names <- character()
    # is there an ion in the basis species?
    if(dissociation) {
      basis.ion <- rownames(aout$basis)[2]
      Z.basis.ion <- makeup(basis.ion)["Z"]
      if(!is.na(Z.basis.ion)) {
        sum.mZ2 <- sum.mZ2 + 10^loga.balance * Z.basis.ion^2
        ion.names <- c(ion.names, basis.ion)
      }
    }
    # add ions present in the species of interest
    for(i in 1:length(loga.equil)) {
      species.ion <- aout$species$name[i]
      Z.species.ion <- makeup(species.ion)["Z"]
      if(!is.na(Z.species.ion)) {
        sum.mZ2 <- sum.mZ2 + 10^loga.equil[[i]] * Z.species.ion^2
        ion.names <- c(ion.names, species.ion)
      }
    }
    IS <- sum.mZ2 / 2
    # report ions used in the ionic strength calculation
    if(find.IS & niter==1) message("solubility: ionic strength calculated for ", paste(ion.names, collapse=" "))
    # report the current ionic strength
    if(find.IS) message("solubility: (iteration ", niter, ") ionic strength range is ", paste(round(range(IS), 4), collapse=" "))
    # stop iterating if we reached the tolerance (or find.IS=FALSE)
    if(!find.IS | all(IS - IS.old < 1e-4)) break
    # on the first iteration, expand argument values for affinity() or mosaic()
    if(niter==1) {
      if(thisfun=="affinity") for(i in 1:length(aout$vals)) {
        aout.save$args[[i]] <- aout$vals[[i]]
      }
      else if(thisfun=="mosaic") {
        for(i in 1:length(aout$vals)) {
          argname <- names(aout$args)[i]
          aout.save$args[[argname]] <- aout$vals[[i]]
        }
      }
    }
    # recalculate the affinity using the new IS
    aout <- suppressMessages(do.call(thisfun, list(aout.save, IS = IS)))
    if(thisfun=="mosaic") aout <- aout$A.species
    niter <- niter + 1
  }

  # do we want the solubility expressed in terms of
  # something other than the first basis species? 20190123
  if(!is.null(in.terms.of)) {
    # write the reaction between the basis species and the new species
    sbasis <- species.basis(in.terms.of)
    # divide the activity of the conserved basis species by the coefficient in the formation reaction
    ibalance <- which.balance(aout$species)
    coeff <- sbasis[, ibalance][1]
    loga.balance <- loga.balance - log10(coeff)
  }

  # make the output
  # (we don't deal with normalized formulas yet, so for now m.balance==n.balance)
  c(aout, list(balance=bout$balance, m.balance=bout$n.balance, n.balance=bout$n.balance,
    loga.balance=loga.balance, Astar=loga.equil, loga.equil=loga.equil))
}
