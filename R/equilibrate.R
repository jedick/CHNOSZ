# CHNOSZ/equilibrate.R
# functions to calculation logarithm of activity
# of species in (metastable) equilibrium

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.misc.R")
#source("util.character.R")

equilibrate <- function(aout, balance=NULL, loga.balance=NULL, 
  ispecies = !grepl("cr", aout$species$state), normalize=FALSE, as.residue=FALSE,
  method=c("boltzmann", "reaction"), tol=.Machine$double.eps^0.25) {
  ### calculate equilibrium activities of species from the affinities 
  ### of their formation reactions from basis species at given activities
  ### split from diagram() 20120925 jmd
  ## if aout is the output from mosaic(), combine the equilibrium activities of basis species
  ## and formed species into an object that can be plotted with diagram() 20190505
  if(aout$fun == "mosaic") {
    # calculate equilibrium activities of species
    if(missing(ispecies)) ispecies <- 1:length(aout$A.species$values)
    if(missing(method)) eqc <- equilibrate(aout$A.species, balance = balance, loga.balance = loga.balance,
      ispecies = ispecies, normalize = normalize, as.residue = as.residue, tol = tol)
    else eqc <- equilibrate(aout$A.species, balance = balance, loga.balance = loga.balance,
      ispecies = ispecies, normalize = normalize, as.residue = as.residue, method = method, tol = tol)
    # make combined object for basis species and species:
    # put together the species matrix and logarithms of equilibrium activity
    eqc$species <- rbind(aout$E.bases[[1]]$species, eqc$species)
    eqc$loga.equil <- c(aout$E.bases[[1]]$loga.equil, eqc$loga.equil)
    # we also need to combine 'values' (values of affinity) because diagram() uses this to get the number of species
    eqc$values <- c(aout$E.bases[[1]]$values, eqc$values)
    return(eqc)
  }
  ## number of possible species
  nspecies <- length(aout$values)
  ## get the balancing coefficients
  bout <- balance(aout, balance)
  n.balance.orig <- n.balance <- bout$n.balance
  balance <- bout$balance
  ## if solids (cr) species are present, find them on a predominance diagram 20191111
  iscr <- grepl("cr", aout$species$state)
  ncr <- sum(iscr)
  if(ncr > 0) dout <- diagram(aout, balance = balance, normalize = normalize, as.residue = as.residue, plot.it = FALSE, limit.water = FALSE)
  if(ncr == nspecies) {
    ## we get here if there are only solids 20200714
    m.balance <- NULL
    Astar <- NULL
    loga.equil <- aout$values
    for(i in 1:length(loga.equil)) loga.equil[[i]][] <- NA
  } else {
    ## we get here if there are any aqueous species 20200714
    ## take selected species in 'ispecies'
    if(length(ispecies)==0) stop("the length of ispecies is zero")
    if(is.logical(ispecies)) ispecies <- which(ispecies)
    # take out species that have NA affinities
    ina <- sapply(aout$values, function(x) all(is.na(x)))
    ispecies <- ispecies[!ina[ispecies]]
    if(length(ispecies)==0) stop("all species have NA affinities")
    if(!identical(ispecies, 1:nspecies)) {
      message(paste("equilibrate: using", length(ispecies), "of", nspecies, "species"))
      aout$species <- aout$species[ispecies, ]
      aout$values <- aout$values[ispecies]
      n.balance <- n.balance[ispecies]
    }
    ## number of species that are left
    nspecies <- length(aout$values)
    ## say what the balancing coefficients are
    if(length(n.balance) < 100) message(paste("equilibrate: n.balance is", c2s(n.balance)))
    ## logarithm of total activity of the balancing basis species
    if(is.null(loga.balance)) {
      # sum up the activities, then take absolute value
      # in case n.balance is negative
      sumact <- abs(sum(10^aout$species$logact * n.balance))
      loga.balance <- log10(sumact)
    }
    # make loga.balance the same length as the values of affinity
    loga.balance <- unlist(loga.balance)
    nvalues <- length(unlist(aout$values[[1]]))
    if(length(loga.balance) == 1) {
      # we have a constant loga.balance
      message(paste0("equilibrate: loga.balance is ", loga.balance))
      loga.balance <- rep(loga.balance, nvalues)
    } else {
      # we are using a variable loga.balance (supplied by the user)
      if(!identical(length(loga.balance), nvalues)) stop("length of loga.balance (", length(loga.balance), ") doesn't match the affinity values (", nvalues, ")")
      message(paste0("equilibrate: loga.balance has same length as affinity values (", length(loga.balance), ")"))
    }
    ## normalize -- normalize the molar formula by the balance coefficients
    m.balance <- n.balance
    isprotein <- grepl("_", as.character(aout$species$name))
    if(any(normalize) | as.residue) {
      if(any(n.balance < 0)) stop("one or more negative balancing coefficients prohibit using normalized molar formulas")
      n.balance[normalize|as.residue] <- 1
      if(as.residue) message(paste("equilibrate: using 'as.residue' for molar formulas"))
      else message(paste("equilibrate: using 'normalize' for molar formulas"))
      # set the formula divisor (m.balance) to 1 for species whose formulas are *not* normalized
      m.balance[!(normalize|as.residue)] <- 1
    } else m.balance[] <- 1
    ## Astar: the affinities/2.303RT of formation reactions with
    ## formed species in their standard-state activities
    Astar <- lapply(1:nspecies, function(i) { 
      # 'starve' the affinity of the activity of the species,
      # and normalize the value by the nor-molar ratio
      (aout$values[[i]] + aout$species$logact[i]) / m.balance[i] 
    })
    ## chose a method and compute the equilibrium activities of species
    if(missing(method)) {
      if(all(n.balance==1)) method <- method[1]
      else method <- method[2]
    }
    message(paste("equilibrate: using", method[1], "method"))
    if(method[1]=="boltzmann") loga.equil <- equil.boltzmann(Astar, n.balance, loga.balance)
    else if(method[1]=="reaction") loga.equil <- equil.reaction(Astar, n.balance, loga.balance, tol)
    ## if we normalized the formulas, get back to activities to species
    if(any(normalize) & !as.residue) {
      loga.equil <- lapply(1:nspecies, function(i) {
        loga.equil[[i]] - log10(m.balance[i])
      })
    }
  }
  ## process cr species 20191111
  if(ncr > 0) {
    # cr species were excluded from equilibrium calculation, so get values back to original lengths
    norig <- length(dout$values)
    n.balance <- n.balance.orig
    imatch <- match(1:norig, ispecies)
    m.balance <- m.balance[imatch]
    Astar <- Astar[imatch]
    loga.equil1 <- loga.equil[[1]]
    loga.equil <- loga.equil[imatch]
    # replace NULL loga.equil with input values (cr only)
    ina <- which(is.na(imatch))
    for(i in ina) {
      loga.equil[[i]] <- loga.equil1
      loga.equil[[i]][] <- dout$species$logact[[i]]
    }
    aout$species <- dout$species
    aout$values <- dout$values
    # find the grid points where any cr species is predominant
    icr <- which(grepl("cr", dout$species$state))
    iscr <- lapply(icr, function(x) dout$predominant == x)
    iscr <- Reduce("|", iscr)
    # at those grid points, make the aqueous species' activities practically zero
    for(i in 1:norig) {
      if(i %in% icr) next
      loga.equil[[i]][iscr] <- -999
    }
    # at the other grid points, make the cr species' activities practically zero
    for(i in icr) {
      ispredom <- dout$predominant == i
      loga.equil[[i]][!ispredom] <- -999
    }
  }
  ## put together the output
  out <- c(aout, list(balance=balance, m.balance=m.balance, n.balance=n.balance,
    loga.balance=loga.balance, Astar=Astar, loga.equil=loga.equil))
  # done!
  return(out)
}

equil.boltzmann <- function(Astar, n.balance, loga.balance) {
  # 20090217 new "abundance" function
  # return equilibrium logarithms of activity of species
  # works using Boltzmann distribution
  # A/At = e^(Astar/n.balance) / sum(e^(Astar/n.balance))
  # where A is activity of the ith residue and
  # At is total activity of residues
  # advantages over equil.reaction
  # 2) loops over species only - much faster
  # 3) no root finding - those games might fail at times
  # disadvantage:
  # 1) only works for per-residue reactions
  # 2) can create NaN logacts if the Astars are huge/small
  if(any(n.balance!=1)) stop("won't run equil.boltzmann for balance <> 1")
  # initialize output object
  A <- Astar
  # remember the dimensions of elements of Astar (could be vector,matrix)
  Astardim <- dim(Astar[[1]])
  Anames <- names(Astar)
  # first loop: make vectors
  A <- palply("", 1:length(A), function(i) as.vector(A[[i]]))
  loga.balance <- as.vector(loga.balance)
  # second loop: get the exponentiated Astars (numerators)
  # need to convert /2.303RT to /RT
  #A[[i]] <- exp(log(10)*Astar[[i]]/n.balance[i])/n.balance[i]
  A <- palply("", 1:length(A), function(i) exp(log(10)*Astar[[i]]/n.balance[i]))
  # third loop: accumulate the denominator
  # initialize variable to hold the sum
  At <- A[[1]]
  At[] <- 0
  for(i in 1:length(A)) At <- At + A[[i]]*n.balance[i]
  # fourth loop: calculate log abundances
  A <- palply("", 1:length(A), function(i) loga.balance + log10(A[[i]]/At))
  # fifth loop: replace dimensions
  for(i in 1:length(A)) dim(A[[i]]) <- Astardim
  # add names and we're done!
  names(A) <- Anames
  return(A)
}

equil.reaction <- function(Astar, n.balance, loga.balance, tol=.Machine$double.eps^0.25) {
  # to turn the affinities/RT (A) of formation reactions into 
  # logactivities of species (logact(things)) at metastable equilibrium
  # 20090217 extracted from diagram and renamed to abundance.old
  # 20080217 idea / 20120128 cleaned-up strategy
  # for any reaction stuff = thing,
  #   A = logK - logQ 
  #     = logK - logact(thing) + logact(stuff)
  # given Astar = A + logact(thing),
  # given Abar = A / n.balance,
  #   logact(thing) = Astar - Abar * n.balance  [2]
  # where n.balance is the number of the balanced quanitity
  # (conserved component) in each species
  # equilibrium values of logact(thing) satifsy:
  # 1) Abar is equal for all species
  # 2) log10( sum of (10^logact(thing) * n.balance) ) = loga.balance  [1]
  #
  # because of the logarithms, we can't solve the equations directly
  # instead, use uniroot() to compute Abar satisfying [1]

  # we can't run on one species
  if(length(Astar)==1) stop("at least two species needed for reaction-based equilibration")
  # remember the dimensions and names
  Adim <- dim(Astar[[1]])
  Anames <- names(Astar)
  # make a matrix out of the list of Astar
  Astar <- list2array(lapply(Astar, c))
  if(length(loga.balance) != nrow(Astar)) stop("length of loga.balance must be equal to the number of conditions for affinity()")
  # that produces the same result (other than colnames) and is much faster than
  #Astar <- sapply(Astar, c)  
  # also, latter has NULL nrow for length(Astar[[x]])==1
  # some function definitions
  # to calculate log of activity of balanced quantity from logact(thing) of all species [1]
  logafun <- function(logact) log10(sum(10^logact * n.balance))
  # to calculate logact(thing) from Abar for the ith condition [2]
  logactfun <- function(Abar, i) Astar[i,] - Abar * n.balance
  # to calculate difference between logafun and loga.balance for the ith condition
  logadiff <- function(Abar, i) loga.balance[i] - logafun(logactfun(Abar, i))
  # to calculate a range of Abar that gives negative and positive values of logadiff for the ith condition
  Abarrange <- function(i) {
    # starting guess of Abar (min/max) from range of Astar / n.balance
    Abar.range <- range(Astar[i, ] / n.balance)
    # diff(Abar.range) can't be 0 (dlogadiff.dAbar becomes NaN)
    if(diff(Abar.range)==0) {
      Abar.range[1] <- Abar.range[1] - 0.1
      Abar.range[2] <- Abar.range[2] + 0.1
    }
    # the range of logadiff
    logadiff.min <- logadiff(Abar.range[1], i)
    logadiff.max <- logadiff(Abar.range[2], i)
    # we're out of luck if they're both infinite
    if(is.infinite(logadiff.min) & is.infinite(logadiff.max))
      stop("FIXME: there are no initial guesses for Abar that give
        finite values of the differences in logarithm of activity
        of the conserved component")
    # if one of them is infinite we might have a chance
    if(is.infinite(logadiff.min)) {
      # decrease the Abar range by increasing the minimum
      Abar.range[1] <- Abar.range[1] + 0.99 * diff(Abar.range)
      logadiff.min <- logadiff(Abar.range[1], i)
      if(is.infinite(logadiff.min)) stop("FIXME: the second initial guess for Abar.min failed")
    }
    if(is.infinite(logadiff.max)) {
      # decrease the Abar range by decreasing the maximum
      Abar.range[2] <- Abar.range[2] - 0.99 * diff(Abar.range)
      logadiff.max <- logadiff(Abar.range[2], i)
      if(is.infinite(logadiff.max)) stop("FIXME: the second initial guess for Abar.max failed")
    }
    iter <- 0
    while(logadiff.min > 0 | logadiff.max < 0) {
      # the change of logadiff with Abar
      # it's a weighted mean of the n.balance
      dlogadiff.dAbar <- (logadiff.max - logadiff.min) / diff(Abar.range)
      # change Abar to center logadiff (min/max) on zero
      logadiff.mean <- mean(c(logadiff.min, logadiff.max))
      Abar.range <- Abar.range - logadiff.mean / dlogadiff.dAbar
      # one iteration is enough for the examples in the package
      # but there might be a case where the range of logadiff doesn't cross zero
      # (e.g. for the carboxylic acid example previously in revisit.Rd)
      logadiff.min <- logadiff(Abar.range[1], i)
      logadiff.max <- logadiff(Abar.range[2], i)
      iter <- 1
      if(iter > 5) {
        stop("FIXME: we seem to be stuck! This function (Abarrange() in
          equil.reaction()) can't find a range of Abar such that the differences
          in logarithm of activity of the conserved component cross zero")
      }
    }
    return(Abar.range)
  }
  # to calculate an equilibrium Abar for the ith condition
  Abarfun <- function(i) {
    # get limits of Abar where logadiff brackets zero
    Abar.range <- Abarrange(i)
    # now for the real thing: uniroot!
    Abar <- uniroot(logadiff, interval=Abar.range, i=i, tol=tol)$root
    return(Abar)
  }
  # calculate the logact(thing) for each condition
  logact <- palply("", 1:nrow(Astar), function(i) {
    # get the equilibrium Abar for each condition
    Abar <- Abarfun(i)
    return(logactfun(Abar, i))
  }) 
  # restore the dimensions and names
  if(length(Adim)==1) logact <- list2array(logact)
  else logact <- sapply(logact, c)
  logact <- lapply(1:nrow(logact), function(i) {
    thisla <- list(logact[i,])[[1]]
    dim(thisla) <- Adim
    return(thisla)
  })
  names(logact) <- Anames
  # all done!
  return(logact)
}

# a function to calculate the total moles of the elements in the output from equilibrate 20190505
moles <- function(eout) {
  # exponentiate loga.equil to get activities
  act <- lapply(eout$loga.equil, function(x) 10^x)
  # initialize list for moles of basis species
  nbasis <- rep(list(act[[1]] * 0), nrow(eout$basis))
  # loop over species
  for(i in 1:nrow(eout$species)) {
    # loop over basis species
    for(j in 1:nrow(eout$basis)) {
      # the coefficient of this basis species in the formation reaction of this species
      n <- eout$species[i, j]
      # accumulate the number of moles of basis species
      nbasis[[j]] <- nbasis[[j]] + act[[i]] * n
    }
  }
  # initialize list for moles of elements (same as number of basis species)
  nelem <- rep(list(act[[1]] * 0), nrow(eout$basis))
  # loop over basis species
  for(i in 1:nrow(eout$basis)) {
    # loop over elements
    for(j in 1:nrow(eout$basis)) {
      # the coefficient of this element in the formula of this basis species
      n <- eout$basis[i, j]
      # accumulate the number of moles of elements
      nelem[[j]] <- nelem[[j]] + nbasis[[i]] * n
    }
  }
  # add element names
  names(nelem) <- colnames(eout$basis)[1:nrow(eout$basis)]
  nelem
}

### unexported functions ###

# return a list containing the balancing coefficients (n.balance) and a textual description (balance)
balance <- function(aout, balance=NULL) {
  ## generate n.balance from user-given or automatically identified basis species
  ## extracted from equilibrate() 20120929
  # 'balance' can be:
  #   NULL                    - autoselect using which.balance
  #   name of basis species   - balanced on this basis species
  #   "length"                   - balanced on sequence length of proteins 
  #                             (default if balance is missing and all species are proteins)
  #   1                       - balanced on one mole of species (formula units)
  #   numeric vector          - user-defined n.balance
  #   "volume"                - standard-state volume listed in thermo()$OBIGT
  # the index of the basis species that might be balanced
  ibalance <- numeric()
  # deal with proteins
  isprotein <- grepl("_", as.character(aout$species$name))
  if(is.null(balance) & all(isprotein)) balance <- "length"
  # try to automatically find a balance
  if(is.null(balance)) {
    ibalance <- which.balance(aout$species)
    # no shared basis species and balance not specified by user - an error
    if(length(ibalance) == 0) stop("no basis species is present in all formation reactions")
  } 
  # change "1" to 1 (numeric) 20170206
  if(identical(balance, "1")) balance <- 1
  if(is.numeric(balance[1])) {
    # a numeric vector
    n.balance <- rep(balance, length.out=length(aout$values))
    msgtxt <- paste0("balance: on supplied numeric argument (", paste(balance, collapse = ","), ")")
    if(identical(balance, 1)) msgtxt <- paste(msgtxt, "[1 means balance on formula units]")
    message(msgtxt)
  } else {
    # "length" for balancing on protein length
    if(identical(balance, "length")) {
      if(!all(isprotein)) stop("'length' was the requested balance, but some species are not proteins")
      n.balance <- protein.length(aout$species$name)
      message("balance: on protein length")
    } else if(identical(balance, "volume")) {
      n.balance <- info(aout$species$ispecies, check.it=FALSE)$V
      message("balance: on volume")
    } else {
      # is the balance the name of a basis species?
      if(length(ibalance)==0) {
        ibalance <- match(balance, rownames(aout$basis))
        if(is.na(ibalance)) stop("basis species (", balance, ") not available to balance reactions")
      }
      # the name of the basis species (need this if we got ibalance which which.balance, above)
      balance <- colnames(aout$species)[ibalance[1]]
      message(paste("balance: on moles of", balance, "in formation reactions"))
      # the balancing coefficients
      n.balance <- aout$species[, ibalance[1]]
      # we check if that all formation reactions contain this basis species
      if(any(n.balance==0)) stop("some species have no ", balance, " in the formation reaction")
    }
  }
  return(list(n.balance=n.balance, balance=balance))
}
