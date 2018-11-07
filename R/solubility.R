# solubility.R: vectorized solubility calculations without uniroot
# 20181031 jmd
# 20181106 work on the output from affinity(); no "equilibrate()" needed!

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("equilibrate.R")
#source("util.misc.R")

solubility <- function(aout, balance=NULL, split=FALSE) {
  ## concept: the logarithms of activities of species at equilibrium are equal to
  ## Astar, the affinities calculated for unit activities of species
  
  ## however, the values in aout can be calculated for other than
  ## unit activities of species, so we have to take away the activites
  Astar <- function(i) aout$values[[i]] + aout$species$logact[i]
  loga.equil <- lapply(1:length(aout$values), Astar)

  ## for a dissociation (split) on a *per reaction* (not system) basis,
  ## apply the divisor here and skip the if(split){} part below
  ## (can be used to reproduce Fig. 4 of Manning et al., 2013)
  if(is.numeric(split)) loga.equil <- lapply(loga.equil, "/", split)

  ## to output loga.balance we need the balancing coefficients
  bout <- balance(aout, balance)
  n.balance <- bout$n.balance
  balance <- bout$balance
  # get logarithm of total activity of the balancing basis species
  logabfun <- function(loga.equil, n.balance) {
    # exponentiate, multiply by n.balance, sum, logarithm
    a.equil <- mapply("^", 10, loga.equil, SIMPLIFY = FALSE)
    a.balance <- mapply("*", a.equil, n.balance, SIMPLIFY=FALSE)
    a.balance <- Reduce("+", a.balance)
    log10(a.balance)
  }
  loga.balance <- logabfun(loga.equil, bout$n.balance)

  # recalculate things for a 1:1 split species (like CaCO3 = Ca+2 + CO3+2)
  if(isTRUE(split)) {
    # the multiplicity becomes the exponent in the reaction quotient
    loga.split <- loga.balance / 2
    # the contribution to affinity
    Asplit <- lapply(n.balance, "*", loga.split)
    # adjust the affinity and get new equilibrium activities
    aout$values <- mapply("-", aout$values, Asplit, SIMPLIFY=FALSE)
    loga.equil <- lapply(1:length(aout$values), Astar)
    # check that the new loga.balance == loga.split; this might not work for non-1:1 species
    loga.balance <- logabfun(loga.equil, n.balance)
    stopifnot(all.equal(loga.balance, loga.split))
  }

  # make the output
  # (we don't deal with normalized formulas yet, so for now m.balance==n.balance)
  c(aout, list(balance=bout$balance, m.balance=bout$n.balance, n.balance=bout$n.balance,
    loga.balance=loga.balance, Astar=loga.equil, loga.equil=loga.equil))
}
