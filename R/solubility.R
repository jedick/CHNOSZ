# solubility.R: vectorized solubility calculations without uniroot
# 20181031 jmd
# 20181106 work on the output from affinity(); no "equilibrate()" needed!
# 20190117 add find.IS and test for dissociation reaction

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("equilibrate.R")
#source("util.misc.R")

solubility <- function(aout, dissociation=NULL, find.IS=FALSE) {
  ## concept: the logarithms of activities of species at equilibrium are equal to
  ## Astar, the affinities calculated for unit activities of species

  ## does the system involve a dissociation reaction?
  if(is.null(dissociation)) {
    # assume FALSE unless determined otherwise
    dissociation <- FALSE
    # react the first basis species to form the first species
    sres <- suppressMessages(subcrt(c(rownames(aout$basis)[1], aout$species$name[1]), c(-1, 1)))
    # note that the reaction is auto-balanced using all the basis species
    # if the reaction involves the second basis species, we consider it to be a dissociation reaction
    if(rownames(aout$basis)[2] %in% sres$reaction$formula) {
      # however, if the second basis species is H2O, H+, e-, O2 (or others?), we don't have enough information, so stop
      if(rownames(aout$basis)[2] %in% c("H2O", "H+", "e-", "O2")) {
        print(sres$reaction)
        stop("unsure whether this is a dissociation reaction. if it is, redefine the basis to put a product ion second")
      }
      dissociation <- TRUE
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
  # get logarithm of total activity of the balancing basis species
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
    # calculate the ionic strength (assuming no ion pairing)
    # i.e. for Sr+2 and SO4-2
    # TODO: determine the correct values for the actual species in the system
    mol.cation <- mol.anion <- 10^loga.balance
    Z <- 2
    IS <- (mol.cation * Z^2 + mol.anion * Z^2) / 2
    # report the current ionic strength
    if(find.IS) message("solubility: (iteration ", niter, ") ionic strength range is ", paste(round(range(IS), 4), collapse=" "))
    # stop iterating if we reached the tolerance (or find.IS=FALSE)
    if(!find.IS | all(IS - IS.old < 1e-4)) break
    # recalculate the affinity using the new IS
    aout <- suppressMessages(do.call(aout$fun, list(aout, IS = IS)))
    niter <- niter + 1
  }

  # make the output
  # (we don't deal with normalized formulas yet, so for now m.balance==n.balance)
  c(aout, list(balance=bout$balance, m.balance=bout$n.balance, n.balance=bout$n.balance,
    loga.balance=loga.balance, Astar=loga.equil, loga.equil=loga.equil))
}
