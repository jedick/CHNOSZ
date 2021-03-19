# solubility.R: vectorized solubility calculations without uniroot
# 20181031 first version jmd
# 20181106 work on the output from affinity(); no "equilibrate()" needed!
# 20190117 add find.IS and test for dissociation reaction
# 20210319 use list of aqueous species as main argument (with back-compatibility for affinity output) and handle multiple minerals

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("equilibrate.R")
#source("util.misc.R")
#source("species.R")
#source("util.args.R")
#source("util.character.R")

# Function to calculate solubilities of multiple minerals 20210303
# species() should be used first to load the minerals (all bearing the same metal)
# 'iaq' lists aqueous species that can be produced by dissolution of the minerals
# '...' contains arguments for affinity() or mosaic() (i.e. plotting variables)
solubility <- function(iaq, ..., in.terms.of = NULL, dissociate = FALSE, find.IS = FALSE) {

  # If iaq is the output of affinity(), use old method 20210318
  if(is.list(iaq)) return(solubility_calc(aout = iaq, in.terms.of = in.terms.of, dissociate = dissociate, find.IS = find.IS))
  # Check whether to use affinity() or mosaic()
  ddd <- list(...)
  if(identical(names(ddd)[1], "bases")) is.mosaic <- TRUE else is.mosaic <- FALSE

  # Save current thermodynamic system settings
  thermo <- get("thermo", CHNOSZ)
  # Use current basis species as a template for the solubility calculations
  ispecies <- basis()$ispecies
  logact <- basis()$logact
  # The current formed species are the minerals to be dissolved
  mineral <- species()
  if(is.null(mineral)) stop("please load minerals or gases with species()")

  # Make a list to store the calculated solubilities for each mineral
  slist <- list()
  # Loop over minerals
  for(i in seq_along(mineral$ispecies)) {
    # Print message
    message(paste("solubility: calculating for", mineral$name[i]))
    # Define basis species with the mineral first (so it will be dissolved)
    ispecies[1] <- mineral$ispecies[i]
    logact[1] <- mineral$logact[i]
    # Use numeric values first and put in buffer names second (needed for demo/gold.R)
    loga.numeric <- suppressWarnings(as.numeric(logact))
    basis(ispecies, loga.numeric)
    is.na <- is.na(loga.numeric)
    if(any(is.na)) basis(rownames(basis())[is.na], logact[is.na])
    # Add aqueous species (no need to define activities here - they will be calculated by solubility_calc)
    species(iaq)
    if(is.mosaic) a <- suppressMessages(mosaic(...)$A.species) else a <- suppressMessages(affinity(...))
    # Calculate solubility of this mineral
    scalc <- solubility_calc(a, in.terms.of = in.terms.of, dissociate = dissociate, find.IS = find.IS)
    # Store the solubilities in the list
    slist[[i]] <- scalc$loga.balance
  }
  
  # Restore the original thermodynamic system settings
  assign("thermo", thermo, CHNOSZ)

  if(length(mineral$ispecies) == 1) {
    # For one mineral, return the results of the solubility calculation
    scalc
  } else {
    # For multiple minerals, the overall solubility is the *minimum* among all the minerals
    smin <- do.call(pmin, slist)
    # Put this into the last-computed 'solubility' object
    scalc$loga.balance <- smin
    scalc$loga.equil <- slist
    scalc$species <- mineral
    # Change the function name stored in the object so diagram() plots loga.balance automatically
    scalc$fun <- "solubilities"
    # Return the object
    scalc
  }
}


# The "nuts and bolts" of solubility calculations
# Moved from solubility() to solubility_calc() 20210318
solubility_calc <- function(aout, dissociate = NULL, find.IS = FALSE, in.terms.of = NULL) {
  ## concept: the logarithms of activities of species at equilibrium are equal to
  ## Astar, the affinities calculated for unit activities of species

  ## is aout the output from mosaic() instead of affinity()?
  aout.save <- aout
  thisfun <- aout$fun
  if(thisfun=="mosaic") aout <- aout$A.species

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

    ## for a dissociation on a *per reaction* (not system) basis, apply the divisor here
    ## (can be used to reproduce Fig. 4 of Manning et al., 2013)
    if(is.numeric(dissociate)) {
      loga.equil <- lapply(loga.equil, "/", dissociate)
    }

    # get the total activity of the balancing basis species
    loga.balance <- logabfun(loga.equil, n.balance)

    # recalculate things for a dissociation reaction (like CaCO3 = Ca+2 + CO3+2)
    if(isTRUE(dissociate)) {
      ndissoc <- 2
      # the number of dissociated products is the exponent in the activity product
      loga.dissoc <- loga.balance / ndissoc
      # the contribution to affinity
      Adissoc <- lapply(n.balance, "*", loga.dissoc)
      # adjust the affinity and get new equilibrium activities
      aout$values <- mapply("-", aout$values, Adissoc, SIMPLIFY = FALSE)
      loga.equil <- lapply(1:length(aout$values), Astar)
      loga.balance <- logabfun(loga.equil, n.balance)
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
    if(dissociate) {
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
  }

  # make the output (we don't deal with normalized formulas yet, so for now m.balance==n.balance)
  # indicate the function used to make this output 20190525
  aout$fun <- "solubility"
  # add names to loga.equil 20190731
  names(loga.equil) <- aout$species$name
  c(aout, list(balance=bout$balance, m.balance=bout$n.balance, n.balance=bout$n.balance, in.terms.of=in.terms.of,
    loga.balance=loga.balance, Astar=loga.equil, loga.equil=loga.equil))
}

