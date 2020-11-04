# CHNOSZ/affinity.R
# calculate affinities of formation reactions

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.affinity.R")
#source("util.units.R")
#source("util.character.R")
#source("util.list.R")
#source("subcrt.R")
#source("buffer.R")
#source("util.args.R")
#source("util.data.R")
#source("species.R")

affinity <- function(..., property=NULL, sout=NULL, exceed.Ttr=FALSE, exceed.rhomin=FALSE,
  return.buffer=FALSE, return.sout=FALSE, balance="PBB", iprotein=NULL, loga.protein=-3, transect = NULL) {
  # ...: variables over which to calculate
  # property: what type of energy
  #   (G.basis,G.species,logact.basis,logK,logQ,A)
  # return.buffer: return buffered activities
  # balance: balance protein buffers on PBB
  # exceed.Ttr: extrapolate Gibbs energies
  #   of minerals beyond their T-limits?
  # sout: provide a previously calculated output from subcrt
  # iprotein: build these proteins from residues (speed optimization)

  # history: 20061027 jmd version 1
  # this is where energy.args() used to sit
  # this is where energy() used to sit

  # argument recall 20190117
  # if the first argument is the result from a previous affinity() calculation,
  # just update the remaining arguments
  args.orig <- list(...)
  # we can only do anything with at least one argument
  if(length(args.orig) > 0) {
    if(identical(args.orig[[1]][1], list(fun="affinity"))) {
      aargs <- args.orig[[1]]$args
      # we can only update arguments given after the first argument
      if(length(args.orig) > 1) {
        for(i in 2:length(args.orig)) {
          if(names(args.orig)[i] %in% names(aargs)) aargs[[names(args.orig)[i]]] <- args.orig[[i]]
          else aargs <- c(aargs, args.orig[i])
        }
      }
      return(do.call(affinity, aargs))
    }
  }

  # the argument list
  args <- energy.args(args.orig, transect = transect)
  args <- c(args, list(sout=sout, exceed.Ttr=exceed.Ttr, exceed.rhomin=exceed.rhomin))

  # the species we're given
  thermo <- get("thermo", CHNOSZ)
  mybasis <- thermo$basis
  myspecies <- thermo$species

  if(!is.null(property)) {
    # the user just wants an energy property
    buffer <- FALSE
    args$what <- property
    out <- do.call("energy",args)
    a <- out$a
    sout <- out$sout
  } else {

    # affinity calculations
    property <- args$what

    # iprotein stuff
    # note that this affinities of the residues are subject
    # to ionization calculations in energy() so no explicit accounting
    # is needed here
    if(!is.null(iprotein)) {
      # check all proteins are available
      if(any(is.na(iprotein))) stop("`iprotein` has some NA values")
      if(!all(iprotein %in% 1:nrow(thermo$protein))) stop("some value(s) of `iprotein` are not rownumbers of thermo()$protein")
      # add protein residues to the species list
      resnames <- c("H2O",aminoacids(3))
      # residue activities set to zero;
      # account for protein activities later
      resprot <- paste(resnames,"RESIDUE",sep="_")
      species(resprot, 0)
      thermo <- get("thermo", CHNOSZ)
      ires <- match(resprot, thermo$species$name)
    }

    # buffer stuff
    buffer <- FALSE
    hasbuffer <- !can.be.numeric(mybasis$logact)
    # however, variables specified in the arguments take precedence over the buffer 20171013
    isnotvar <- !rownames(mybasis) %in% args$vars
    ibufbasis <- which(hasbuffer & isnotvar)
    if(!is.null(mybasis) & length(ibufbasis) > 0) {
      buffer <- TRUE
      message('affinity: loading buffer species')
      if(!is.null(thermo$species)) is.species <- 1:nrow(thermo$species) else is.species <- numeric()
      is.buffer <- buffer(logK=NULL)
      thermo <- get("thermo", CHNOSZ)
      is.buff <- numeric()
      for(i in 1:length(is.buffer)) is.buff <- c(is.buff,as.numeric(is.buffer[[i]]))
      is.only.buffer <- is.buff[!is.buff %in% is.species]
      buffers <- names(is.buffer)
      # reorder the buffers according to thermo$buffer
      buffers <- buffers[order(match(buffers,thermo$buffer$name))]
    }

    # here we call 'energy'
    aa <- do.call("energy",args)
    a <- aa$a
    sout <- aa$sout

    if(return.sout) return(sout)

    # more buffer stuff
    if(buffer) {
      args$what <- "logact.basis"
      args$sout <- sout
      logact.basis.new <- logact.basis <- do.call("energy",args)$a
      ibasis.new <- numeric()
      for(k in 1:length(buffers)) {
        ibasis <- which(as.character(mybasis$logact)==buffers[k])
        # calculate the logKs from the affinities
        logK <- a
        for(i in 1:length(logK)) {
          logK[[i]] <- logK[[i]] + thermo$species$logact[i]
          for(j in 1:length(logact.basis.new)) {
            logK[[i]] <- logK[[i]] - logact.basis.new[[j]] * thermo$species[i,j]
          }
        }
        lbn <- buffer(logK=logK,ibasis=ibasis,logact.basis=logact.basis.new,
          is.buffer=as.numeric(is.buffer[[which(names(is.buffer)==buffers[k])]]),balance=balance)
        for(j in 1:length(logact.basis.new)) if(j %in% ibasis) logact.basis.new[[j]] <- lbn[[2]][[j]]
        # calculation of the buffered activities' effect on chemical affinities
        is.only.buffer.new <- is.only.buffer[is.only.buffer %in% is.buffer[[k]]]
        for(i in 1:length(a)) {
          if(i %in% is.only.buffer.new) next
          for(j in 1:nrow(mybasis)) {
            # let's only do this for the basis species specified by the user
            # even if others could be buffered
            if(!j %in% ibufbasis) next
            if(!j %in% ibasis) next
            aa <- a[[i]]
            a[[i]] <- aa + (logact.basis.new[[j]] - logact.basis[[j]]) * thermo$species[i,j]
            #if(!identical(a[[i]],aa)) print(paste(i,j))
          }
        }
        if(k==length(buffers) & return.buffer) {
          logact.basis.new <- lbn[[2]]
          ibasis.new <- c(ibasis.new,lbn[[1]])
        } else ibasis.new <- c(ibasis.new,ibasis)
      }
      species(is.only.buffer,delete=TRUE)
      if(length(is.only.buffer) > 0) a <- a[-is.only.buffer]
      # to return the activities of buffered basis species
      tb <- logact.basis.new[unique(ibasis.new)]
      if(!is.null(ncol(tb[[1]]))) {
        nd <- sum(dim(tb[[1]]) > 1)
        # TODO: apply names for more than two dimensions
        if(nd < 3) {
          for(i in 1:length(tb)) {
            #tb[[i]] <- as.data.frame(tb[[i]])
            if(nd > 0) rownames(tb[[i]]) <- 
              seq(args$lims[[1]][1],args$lims[[1]][2],length.out=args$lims[[1]][3])
            if(nd > 1) colnames(tb[[i]]) <- 
              seq(args$lims[[2]][1],args$lims[[2]][2],length.out=args$lims[[2]][3])
          }
        }
      }
    }

    # more iprotein stuff
    if(!is.null(iprotein)) {
      # 20090331 fast protein calculations
      # function to calculate affinity of formation reactions
      # from those of residues
      loga.protein <- rep(loga.protein,length.out=length(iprotein))
      protein.fun <- function(ip) {
        tpext <- as.numeric(thermo$protein[iprotein[ip],5:25])
        return(Reduce("+", pprod(a[ires],tpext)) - loga.protein[ip])
      }
      # use another level of indexing to let the function
      # report on its progress
      jprotein <- 1:length(iprotein)
      protein.affinity <- palply("", jprotein, protein.fun)
      ## update the species list
      # we use negative values for ispecies to denote that
      # they index thermo$protein and not thermo$species
      ispecies <- -iprotein
      # the current species list, containing the residues
      resspecies <- thermo$species
      # now we can delete the residues from the species list
      species(ires,delete=TRUE)
      # state and protein names
      state <- resspecies$state[1]
      name <- paste(thermo$protein$protein[iprotein],thermo$protein$organism[iprotein],sep="_")
      # the numbers of basis species in formation reactions of the proteins
      protbasis <- t(t((resspecies[ires,1:nrow(mybasis)])) %*% t((thermo$protein[iprotein,5:25])))
      # put them together
      protspecies <- cbind(protbasis,data.frame(ispecies=ispecies,logact=loga.protein,state=state,name=name))
      myspecies <- rbind(myspecies,protspecies)
      rownames(myspecies) <- 1:nrow(myspecies)
      ## update the affinity values
      names(protein.affinity) <- ispecies
      a <- c(a,protein.affinity)
      a <- a[-ires]
    }

  }

  # put together return values
  # constant T and P, in Kelvin and bar
  T <- args$T
  P <- args$P
  # the variable names and values
  vars <- args$vars
  vals <- args$vals
  # if T or P is variable it's not constant;
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
  # get Eh
  iEh <- match("Eh", names(args.orig))
  if(!is.na(iEh)) {
    vars[[iEh]] <- "Eh"
    # we have to reconstruct the values used because
    # they got converted to log_a(e-) at an unknown temperature
    Eharg <- args.orig[[iEh]]
    if(length(Eharg) > 3) Ehvals <- Eharg
    else if(length(Eharg) == 3) Ehvals <- seq(Eharg[1], Eharg[2], length.out=Eharg[3])
    else if(length(Eharg) == 2) Ehvals <- seq(Eharg[1], Eharg[2], length.out=256)
    vals[[iEh]] <- Ehvals
  }
  # get pe and pH
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
  # use the variable names for the vals list 20190415
  names(vals) <- vars

  # content of return value depends on buffer request
  if(return.buffer) return(c(tb, list(vars=vars, vals=vals)))
  # for argument recall, include all arguments (except sout) in output 20190117
  allargs <- c(args.orig, list(property=property, exceed.Ttr=exceed.Ttr, exceed.rhomin=exceed.rhomin,
    return.buffer=return.buffer, balance=balance, iprotein=iprotein, loga.protein=loga.protein))
  # add IS value only if it given as an argument 20171101
  # (even if its value is 0, the presence of IS will trigger diagram() to use "m" instead of "a" in axis labels)
  iIS <- match("IS", names(args.orig))
  if(!is.na(iIS)) a <- list(fun="affinity", args=allargs, sout=sout, property=property,
                            basis=mybasis, species=myspecies, T=T, P=P, IS=args$IS, vars=vars, vals=vals, values=a)
  else a <- list(fun="affinity", args=allargs, sout=sout, property=property,
                 basis=mybasis, species=myspecies, T=T, P=P, vars=vars, vals=vals, values=a)
  if(buffer) a <- c(a, list(buffer=tb))
  return(a)
}

