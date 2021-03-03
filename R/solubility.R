# solubility.R: vectorized solubility calculations without uniroot
# 20181031 jmd
# 20181106 work on the output from affinity(); no "equilibrate()" needed!
# 20190117 add find.IS and test for dissociation reaction
# 20190730 add codeanal; test all species for dissociation reaction

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("equilibrate.R")
#source("util.misc.R")
#source("species.R")
#source("util.args.R")
#source("util.character.R")

solubility <- function(aout, dissociation = NULL, find.IS = FALSE, in.terms.of = NULL, codeanal = FALSE) {
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
    # we can only test for dissociation with more than one basis species 20190730
    nbasis <- nrow(aout$basis)
if(codeanal) print(paste("number of basis species:", nbasis))
    if(nbasis > 1) {
if(codeanal) print("-- testing for dissociation reactions:")
      # if the reaction to form the first species involves the second basis species, we consider it to be a dissociation reaction
      #if(aout$species[1, 2] != 0) {
      # change this to test for all species 20190730
      if(all(aout$species[2] != 0)) {
if(codeanal) print("   reactions to form all species involve the second basis species")
        # 20190123 (corundum calculation): if there are only H2O, H+, and e-
        # besides the first basis species, it's not a dissociation reaction
        nH2O <- sum(rownames(aout$basis) %in% c("H2O", "H+", "e-", "O2", "H2"))
        if(nbasis > (nH2O + 1)) {
          # if we got here, and the second basis species is H2O, H+, e-, O2 (or others?), we don't have enough information, so stop
          if(rownames(aout$basis)[2] %in% c("H2O", "H+", "e-", "O2", "H2")) {
            stop("Unsure whether the first formation reaction is a dissociation reaction.\nSet the 'dissociation' argument to TRUE or FALSE, or redefine the basis to put a product ion second.")
          }
          dissociation <- TRUE
        } else if(codeanal) print("   the basis consists of only the conserved species and H2O-derived species")
      } else if(codeanal) print("   at least one formation reaction doesn't involve the second basis species")
    } else if(codeanal) print("-- no dissociation reactions are possible")
    message("solubility: test for dissociation reaction returns ", dissociation)
  } else message("solubility: from argument, dissociation ratio is ", dissociation)

  # get starting ionic strength (probably zero, but could be anything set by user)
  IS <- aout$IS
  if(is.null(IS)) IS <- 0

  ## to output loga.balance we need the balancing coefficients
  bout <- balance(aout)
  n.balance <- bout$n.balance
  balance <- bout$balance
if(codeanal) print(paste0("balancing coefficients [", balance, "]: ", paste(n.balance, collapse = " ")))
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
    if(is.numeric(dissociation)) {
if(codeanal) for(ii in 1:length(loga.equil)) print(paste0("loga.equil0 [", aout$species$name[ii], "]: ", round(loga.equil[[ii]], 3)))
if(codeanal) print(paste0("applying ", dissociation, "-fold dissociation correction (no interaction)"))
      loga.equil <- lapply(loga.equil, "/", dissociation)
    }


    # get the total activity of the balancing basis species
    loga.balance <- logabfun(loga.equil, n.balance)
if(codeanal & !isTRUE(dissociation)) print(paste0("loga.balance [", balance, "]: ", round(loga.balance, 3)))

    # recalculate things for a dissociation reaction (like CaCO3 = Ca+2 + CO3+2)
    if(isTRUE(dissociation)) {
      ndissoc <- 2
if(codeanal) print(paste0("loga.balance0 [", balance, "]: ", round(loga.balance, 3)))
if(codeanal) for(ii in 1:length(loga.equil)) print(paste0("loga.equil0 [", aout$species$name[ii], "]: ", round(loga.equil[[ii]], 3)))
if(codeanal) print(paste0("applying ", ndissoc, "-fold dissociation correction (with interaction)"))
      # the number of dissociated products is the exponent in the activity product
      loga.dissoc <- loga.balance / ndissoc
      # the contribution to affinity
      Adissoc <- lapply(n.balance, "*", loga.dissoc)
      # adjust the affinity and get new equilibrium activities
      aout$values <- mapply("-", aout$values, Adissoc, SIMPLIFY = FALSE)
      loga.equil <- lapply(1:length(aout$values), Astar)
      loga.balance <- logabfun(loga.equil, n.balance)
if(codeanal) print(paste0("loga.balance [", balance, "]: ", round(loga.balance, 3)))
      # check that the new loga.balance == loga.dissoc
      # TODO: does this work for non-1:1 species?
      stopifnot(all.equal(loga.balance, loga.dissoc))
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
    # instead of using aout$species$name, use info() to get formulas 20190309
    species.formulas <- suppressMessages(info(aout$species$ispecies)$formula)
    for(i in 1:length(loga.equil)) {
      species.ion <- species.formulas[i]
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
    # (e.g. convert pH = c(0, 14) to pH = seq(0, 14, 256) so that it has the same length as the IS values)
    # we don't do this if aout$vals is NA 20190731
    if(niter==1 & !all(is.na(aout$vals))) {
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
if(codeanal) print(paste0("loga.balance [in terms of ", in.terms.of, "]: ", round(loga.balance, 3)))
  }

  # make the output (we don't deal with normalized formulas yet, so for now m.balance==n.balance)
  # indicate the function used to make this output 20190525
  aout$fun <- "solubility"
  # add names to loga.equil 20190731
  names(loga.equil) <- aout$species$name
  c(aout, list(balance=bout$balance, m.balance=bout$n.balance, n.balance=bout$n.balance, in.terms.of=in.terms.of,
    loga.balance=loga.balance, Astar=loga.equil, loga.equil=loga.equil))
}

# Calculate solubilities of multiple minerals 20210303
# a_cr: affinities for minerals (all bearing the same metal)
# i_aq: aqueous species that can be produced by dissolution of the minerals
# FIXME: what to do about 'dissociation' argument?
solubilities <- function(a_cr, i_aq, in.terms.of = NULL, dissociation = NULL) {
  # If a_cr is the output from mosaic(), just get the species' affinities
  is.mosaic <- FALSE
  m_cr <- a_cr
  if(identical(a_cr$fun, "mosaic")) {
    a_cr <- a_cr$A.species
    is.mosaic <- TRUE
  }
  # Find all stable minerals across diagram
  d_cr <- diagram(a_cr, plot.it = FALSE)
  d_cr.stable <- d_cr$species$ispecies[unique(as.vector(d_cr$predominant))]
  # Use basis species in a_cr as a template for the solubility calculations
  ispecies <- a_cr$basis$ispecies
  logact <- a_cr$basis$logact

  # Make a list to store the calculated solubilities for each mineral
  slist <- list()
  # Loop over stable minerals
  for(i in seq_along(d_cr.stable)) {
    # Define basis species with the mineral first (so it will be dissolved)
    ispecies[1] <- d_cr.stable[i]
    basis(ispecies, logact)
    # Add aqueous species (no need to define activities here - they will be calculated)
    species(i_aq)
    # Calculate affinities of formation reactions of this mineral at same conditions as a_cr (argument recall)
    if(is.mosaic) a <- mosaic(m_cr)$A.species else a <- affinity(a_cr)
    # Calculate solubility of this mineral
    s <- solubility(a, in.terms.of = in.terms.of, dissociation = dissociation)
    # Store the solubilities in the list
    slist[[i]] <- s$loga.balance
  }

  # The overall solubility is the *minimum* among all the minerals
  smin <- do.call(pmin, slist)
  # Put this into the last-computed 'solubility' object
  s$loga.balance <- smin
  # Change the function name stored in the object so diagram() plots loga.balance automatically
  s$fun <- "solubilities"
  # Return the object
  s
}

