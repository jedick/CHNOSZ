# CHNOSZ/affinity.R
# Calculate affinities of formation reactions

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.affinity.R")
#source("util.units.R")
#source("util.character.R")
#source("util.list.R")
#source("subcrt.R")
#source("buffer.R")
#source("util.args.R")
#source("util.data.R")
#source("species.R")
#source("info.R")
#source("hkf.R")
#source("cgl.R")

affinity <- function(..., property = NULL, sout = NULL, exceed.Ttr = FALSE, exceed.rhomin = FALSE,
  return.buffer = FALSE, return.sout = FALSE, balance = "PBB", iprotein = NULL, loga.protein = 0, transect = NULL) {
  # ...: variables over which to calculate
  # property: what type of energy
  #   (G.basis, G.species, logact.basis, logK, logQ, A)
  # return.buffer: return buffered activities
  # balance: balance protein buffers on PBB
  # exceed.Ttr: extrapolate Gibbs energies
  #   of minerals beyond their T-limits?
  # sout: provide a previously calculated output from subcrt
  # iprotein: build these proteins from residues (speed optimization)

  # History: 20061027 jmd version 1
  # this is where energy.args() used to sit
  # this is where energy() used to sit

  # Argument recall 20190117
  # If the first argument is the result from a previous affinity() calculation,
  # just update the remaining arguments
  args.orig <- list(...)
  # We can only do anything with at least one argument
  if(length(args.orig) > 0) {
    if(identical(args.orig[[1]][1], list(fun = "affinity"))) {
      aargs <- args.orig[[1]]$args
      # We can only update arguments given after the first argument
      if(length(args.orig) > 1) {
        for(i in 2:length(args.orig)) {
          if(names(args.orig)[i] %in% names(aargs)) aargs[[names(args.orig)[i]]] <- args.orig[[i]]
          else aargs <- c(aargs, args.orig[i])
        }
      }
      return(do.call(affinity, aargs))
    }
  }

  # The argument list
  args <- energy.args(args.orig, transect = transect)
  args <- c(args, list(sout = sout, exceed.Ttr = exceed.Ttr, exceed.rhomin = exceed.rhomin))

  # The species we're given
  thermo <- get("thermo", CHNOSZ)
  mybasis <- thermo$basis
  myspecies <- thermo$species

  if(!is.null(property)) {
    # The user just wants an energy property
    buffer <- FALSE
    args$what <- property
    out <- do.call("energy", args)
    a <- out$a
    sout <- out$sout
  } else {

    # Affinity calculations
    property <- args$what

    # Protein stuff
    # Note that affinities of the residues are modified by ionization calculations in energy(), not here
    if(!is.null(iprotein)) {
      # Check all proteins are available
      if(any(is.na(iprotein))) stop("`iprotein` has some NA values")
      if(!all(iprotein %in% 1:nrow(thermo$protein))) stop("some value(s) of `iprotein` are not rownumbers of thermo()$protein")
      # Add protein residues to the species list
      resnames <- c("H2O", aminoacids(3))
      # Residue activities set to zero; account for protein activities later
      resprot <- paste(resnames,"RESIDUE", sep = "_")
      species(resprot, 0)
      # Re-read thermo because the preceding command changed the species
      thermo <- get("thermo", CHNOSZ)
      ires <- match(resprot, thermo$species$name)
    }

    # Buffer stuff
    buffer <- FALSE
    hasbuffer <- !can.be.numeric(mybasis$logact)
    # However, variables specified in the arguments take precedence over the buffer 20171013
    isnotvar <- !rownames(mybasis) %in% args$vars
    ibufbasis <- which(hasbuffer & isnotvar)
    if(!is.null(mybasis) & length(ibufbasis) > 0) {
      buffer <- TRUE
      message('affinity: loading buffer species')
      if(!is.null(thermo$species)) is.species <- 1:nrow(thermo$species) else is.species <- numeric()
      # Load the species in the buffer and get their species number(s)
      is.buffer <- buffer(logK = NULL)
      # Re-read thermo because the preceding command changed the species
      thermo <- get("thermo", CHNOSZ)
      buffers <- names(is.buffer)
      # Find species that are only in the buffer, not in the starting species list
      is.only.buffer <- setdiff(unlist(is.buffer), is.species)
    }

    # Here we call 'energy'
    aa <- do.call("energy", args)
    a <- aa$a
    sout <- aa$sout

    if(return.sout) return(sout)

    # More buffer stuff
    if(buffer) {
      args$what <- "logact.basis"
      args$sout <- sout
      logact.basis.new <- logact.basis <- do.call("energy", args)$a
      ibasis.new <- numeric()
      for(k in 1:length(buffers)) {
        ibasis <- which(as.character(mybasis$logact) == buffers[k])
        # Calculate the logKs from the affinities
        logK <- a
        for(i in 1:length(logK)) {
          logK[[i]] <- logK[[i]] + thermo$species$logact[i]
          for(j in 1:length(logact.basis.new)) {
            logK[[i]] <- logK[[i]] - logact.basis.new[[j]] * thermo$species[i, j]
          }
        }
        lbn <- buffer(logK = logK, ibasis = ibasis, logact.basis = logact.basis.new,
          is.buffer = is.buffer[[k]], balance = balance)
        for(j in 1:length(logact.basis.new)) if(j %in% ibasis) logact.basis.new[[j]] <- lbn[[2]][[j]]
        # Calculation of the buffered activities' effect on chemical affinities
        is.only.buffer.new <- is.only.buffer[is.only.buffer %in% is.buffer[[k]]]
        for(i in 1:length(a)) {
          if(i %in% is.only.buffer.new) next
          for(j in 1:nrow(mybasis)) {
            # Let's only do this for the basis species specified by the user
            # even if others could be buffered
            if(!j %in% ibufbasis) next
            if(!j %in% ibasis) next
            aa <- a[[i]]
            a[[i]] <- aa + (logact.basis.new[[j]] - logact.basis[[j]]) * thermo$species[i, j]
            #if(!identical(a[[i]], aa)) print(paste(i, j))
          }
        }
        if(k == length(buffers) & return.buffer) {
          logact.basis.new <- lbn[[2]]
          ibasis.new <- c(ibasis.new, lbn[[1]])
        } else ibasis.new <- c(ibasis.new, ibasis)
      }
      species(is.only.buffer, delete = TRUE)
      if(length(is.only.buffer) > 0) a <- a[-is.only.buffer]
      # To return the activities of buffered basis species
      tb <- logact.basis.new[unique(ibasis.new)]
      if(!is.null(ncol(tb[[1]]))) {
        nd <- sum(dim(tb[[1]]) > 1)
        # TODO: apply names for more than two dimensions
        if(nd < 3) {
          for(i in 1:length(tb)) {
            #tb[[i]] <- as.data.frame(tb[[i]])
            if(nd > 0) rownames(tb[[i]]) <- 
              seq(args$lims[[1]][1], args$lims[[1]][2], length.out = args$lims[[1]][3])
            if(nd > 1) colnames(tb[[i]]) <- 
              seq(args$lims[[2]][1], args$lims[[2]][2], length.out = args$lims[[2]][3])
          }
        }
      }
    }

    # More iprotein stuff
    if(!is.null(iprotein)) {
      # Fast protein calculations 20090331
      # Function to calculate affinity of formation reactions from those of residues
      loga.protein <- rep(loga.protein, length.out = length(iprotein))
      protein.fun <- function(ip) {
        tpext <- as.numeric(thermo$protein[iprotein[ip], 5:25])
        return(Reduce("+", pprod(a[ires], tpext)) - loga.protein[ip])
      }
      # Use another level of indexing to let the function report on its progress
      jprotein <- 1:length(iprotein)
      protein.affinity <- palply("", jprotein, protein.fun)
      ## Update the species list
      # We use negative values for ispecies to denote that
      # they index thermo$protein and not thermo$species
      ispecies <- -iprotein
      # The current species list, containing the residues
      resspecies <- thermo$species
      # Now we can delete the residues from the species list
      species(ires, delete = TRUE)
      # State and protein names
      state <- resspecies$state[1]
      name <- paste(thermo$protein$protein[iprotein], thermo$protein$organism[iprotein], sep = "_")
      # The numbers of basis species in formation reactions of the proteins
      protbasis <- t(t((resspecies[ires, 1:nrow(mybasis)])) %*% t((thermo$protein[iprotein, 5:25])))
      # Put them together
      protspecies <- cbind(protbasis, data.frame(ispecies = ispecies, logact = loga.protein, state = state, name = name))
      myspecies <- rbind(myspecies, protspecies)
      rownames(myspecies) <- 1:nrow(myspecies)
      ## Update the affinity values
      names(protein.affinity) <- ispecies
      a <- c(a, protein.affinity)
      a <- a[-ires]
    }

  }

  # Put together return values
  # Constant T and P, in Kelvin and bar
  T <- args$T
  P <- args$P
  # The variable names and values
  vars <- args$vars
  vals <- args$vals
  # If T or P is variable it's not constant;
  # and set it to user's units
  iT <- match("T", vars)
  if(!is.na(iT)) {
    T <- numeric()
    vals[[iT]] <- outvert(vals[[iT]], "K")
  }
  iP <- match("P", vars)
  if(!is.na(iP) > 0) {
    P <- numeric()
    vals[[iP]] <- outvert(vals[[iP]], "bar")
  }
  # Get Eh
  iEh <- match("Eh", names(args.orig))
  if(!is.na(iEh)) {
    vars[[iEh]] <- "Eh"
    # We have to reconstruct the values used because
    # they got converted to log_a(e-) at an unknown temperature
    Eharg <- args.orig[[iEh]]
    if(length(Eharg) > 3) Ehvals <- Eharg
    else if(length(Eharg) == 3) Ehvals <- seq(Eharg[1], Eharg[2], length.out = Eharg[3])
    else if(length(Eharg) == 2) Ehvals <- seq(Eharg[1], Eharg[2], length.out = 256)
    vals[[iEh]] <- Ehvals
  }
  # Get pe and pH
  ipe <- match("pe", names(args.orig))
  if(!is.na(ipe)) {
    ie <- match("e-", names(args$lims))
    vars[[ie]] <- "pe"
    vals[[ie]] <- -args$vals[[ie]]
  }
  ipH <- match("pH", names(args.orig))
  if(!is.na(ipH)) {
    iH <- match("H+", names(args$lims))
    vars[[iH]] <- "pH"
    vals[[iH]] <- -args$vals[[iH]]
  }
  # Use the variable names for the vals list 20190415
  names(vals) <- vars

  # Content of return value depends on buffer request
  if(return.buffer) return(c(tb, list(vars = vars, vals = vals)))
  # For argument recall, include all arguments (except sout) in output 20190117
  allargs <- c(args.orig, list(property = property, exceed.Ttr = exceed.Ttr, exceed.rhomin = exceed.rhomin,
    return.buffer = return.buffer, balance = balance, iprotein = iprotein, loga.protein = loga.protein))
  # Add IS value only if it given as an argument 20171101
  # (even if its value is 0, the presence of IS will trigger diagram() to use "m" instead of "a" in axis labels)
  iIS <- match("IS", names(args.orig))
  if(!is.na(iIS)) a <- list(fun = "affinity", args = allargs, sout = sout, property = property,
                            basis = mybasis, species = myspecies, T = T, P = P, IS = args$IS, vars = vars, vals = vals, values = a)
  else a <- list(fun = "affinity", args = allargs, sout = sout, property = property,
                 basis = mybasis, species = myspecies, T = T, P = P, vars = vars, vals = vals, values = a)
  if(buffer) a <- c(a, list(buffer = tb))
  return(a)
}
