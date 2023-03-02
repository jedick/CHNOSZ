# CHNOSZ/util-affinity.R
# Helper functions for affinity()

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.units.R")

slice.affinity <- function(affinity,d=1,i=1) {
  # Take a slice of affinity along one dimension
  a <- affinity
  for(j in 1:length(a$values)) {
    # Preserve the dimensions (especially: names(mydim))
    # - fix for change in behavior of aperm in R-devel 2015-11-17
    mydim <- dim(a$values[[j]])
    a$values[[j]] <- as.array(slice(a$values[[j]],d=d,i=i))
    # The dimension from which we take the slice vanishes
    dim(a$values[[j]]) <- mydim[-d]
  }
  return(a)
}

### Unexported functions ###

energy <- function(what,vars,vals,lims,T=298.15,P="Psat",IS=0,sout=NULL,exceed.Ttr=FALSE,exceed.rhomin=FALSE,transect=FALSE) {
  # 20090329 Dxtracted from affinity() and made to
  # deal with >2 dimensions (variables)

  # Calculate "what" property
  # logK - logK of formation reactions
  # logact.basis - logarithms of activities of basis species in reactions
  # logact.species - logarithms of activities of species in reactions
  # logQ - logQ of formation reactions
  # A - A/2.303RT of formation reactions
  # "vars": vector of variables (T,P,basisnames,IS)
  # "vals": list of values for each variable
  # "lims": used to get the dimensions (resolution for each variable)

  ### Some argument checking
  if(length(unique(vars)) != length(vars)) stop("please supply unique variable names")
  ## Basis definition / number of basis species 
  mybasis <- basis()
  nbasis <- nrow(mybasis)
  ## Species definition / number of species
  myspecies <- get("thermo", CHNOSZ)$species
  if(is.character(what)) {
    if(is.null(myspecies)) stop('species properties requested, but species have not been defined')
    nspecies <- nrow(myspecies)
    if(!identical(what,"logact.basis")) ispecies <- 1:nspecies
  }
  ## The dimensions of our arrays
  resfun <- function(lim) lim[3]
  mydim <- sapply(lims,resfun)
  if(transect) {
    if(!isTRUE(all.equal(min(mydim),max(mydim)))) stop("variables define a transect but their lengths are not all equal")
    mydim <- mydim[[1]]
  }
  ## The number of dimensions we have
  if(transect) nd <- 1 else nd <- length(vars) # == length(mydim)
  ## Basis names / which vars denote basis species
  basisnames <- rownames(mybasis)
  ibasisvar <- match(vars,basisnames)
  varisbasis <- !is.na(ibasisvar)
  ibasisvar <- ibasisvar[!is.na(ibasisvar)]
  ## Which vars are in P,T,IS
  varissubcrt <- vars %in% c("P","T","IS")
  if(sum(varissubcrt) > 2) stop("only up to 2 of P,T,IS are supported")
  ## Categorize the basis species:
  # 0 - not in the vars; 1 - one of the vars
  ibasis <- 1:nbasis
  ibasis0 <- ibasis[!ibasis %in% ibasisvar]
  ibasis1 <- ibasis[ibasis %in% ibasisvar]
  if(identical(what,"logact.basis")) ispecies <- ibasis
  ## What subcrt variable is used to make a 2-D grid?
  workaround.IS <- FALSE
  if(sum(varissubcrt) > 1 & !transect) {
    grid <- vars[varissubcrt][1]
    if("IS" %in% vars) {
      if(grid != "IS") workaround.IS <- TRUE
      grid <- "IS"
    }
  } else grid <- NULL
  ### Done argument processing

  ### Function to index the variables in a permuted order
  # by swapping ivar for the first
  # e.g. for nd=5, ivars(4)=c(4,2,3,1,5)
  ivars <- function(ivar,iv=NULL) {
    if(nd==0) return(1)
    if(is.null(iv)) iv <- 1:nd
    iv.1 <- iv[1]
    iv[1] <- ivar
    iv[ivar] <- iv.1
    return(iv)
  }

  ### Functions for logact / logQ
  # A basis species not in var
  logact.basis0.fun <- function(ibasis) {
    logact <- mybasis$logact[ibasis]
    # for the case where this basis species is buffered
    if(!can.be.numeric(logact)) logact <- 0
    else logact <- as.numeric(logact)
    return(array(logact,mydim))
  }
  # A basis species in var
  logact.basis1.fun <- function(ivar) {
    dim.in <- dim(vals[[ivar]])
    if(is.null(dim.in)) dim.in <- 0
    # why doesn't this work?
    # if(identical(dim.in,mydim))
    if(all(dim.in==mydim) | transect) return(vals[[ivar]])
    else return(aperm(array(vals[[ivar]],mydim[ivars(ivar)]),ivars(ivar)))
  }
  # Any basis species
  logact.basis.fun <- function(ibasis) {
    if(ibasis %in% ibasis0) return(logact.basis0.fun(ibasis))
    else return(logact.basis1.fun(match(basisnames[ibasis],vars)))
  }
  # All basis species
  logact.basis <- function() lapply(ibasis,logact.basis.fun)
  # logact of a single species
  logact.species.fun <- function(ispecies) array(myspecies$logact[ispecies],mydim)
  ## Contributions to logQ
  # From a single basis species 
  logQ.basis.fun <- function(ibasis,coeff) - coeff * logact.basis.fun(ibasis)
  # From all basis species in a single formation reaction
  logQ.basis.species <- function(ispecies) 
    Reduce("+", mapply(logQ.basis.fun,ibasis,myspecies[ispecies,1:nbasis],SIMPLIFY=FALSE))
  # All basis species in all reactions
  logQ.basis <- function() mapply(logQ.basis.species,1:nspecies,SIMPLIFY=FALSE)
  # By a single species
  logQ.species.fun <- function(ispecies,coeff) coeff * logact.species.fun(ispecies)
  # By all species
  logQ.species <- function() 
    mapply(logQ.species.fun,1:nspecies,1,SIMPLIFY=FALSE)
  # Total logQ of all reactions
  logQ <- function() lsum(logQ.basis(),logQ.species())

  ### Function for calling subcrt
  sout.fun <- function(property="logK") {
    species <- c(mybasis$ispecies,myspecies$ispecies)
    if(!is.null(sout)) {
      # Extract the needed species from a provided sout 20190131
      isout <- match(species, sout$species$ispecies)
      this.sout <- sout
      this.sout$species <- this.sout$species[isout, ]
      this.sout$out <- this.sout$out[isout]
      return(this.sout) 
    } else {
      ## subcrt arguments
      if("T" %in% vars) T <- vals[[which(vars=="T")]]
      if("P" %in% vars) P <- vals[[which(vars=="P")]]
      if("IS" %in% vars) IS <- vals[[which(vars=="IS")]]
      s.args <- list(species=species,property=property,T=T,P=P,IS=IS,grid=grid,convert=FALSE,exceed.Ttr=exceed.Ttr,exceed.rhomin=exceed.rhomin)
      sout <- do.call("subcrt",s.args)
      # 20191210 workaround a limitation in subcrt(): only IS (if present) can be the 'grid' variable
      if(workaround.IS) {
        lenIS <- length(IS)
        # Only T or P has length > 1
        lenTP <- max(length(T), length(P))
        # Make an indexing matrix "byrow" then vectorize it (not byrow) to get permutation indices
        indmat <- matrix(1:(lenIS * lenTP), nrow = lenIS, byrow = TRUE)
        ind <- as.vector(indmat)
        sout$out <- lapply(sout$out, function(x) x[indmat, ])
      }
      # Species indices are updated by subcrt() for minerals with phase transitions
      # e.g. i <- info("chalcocite"); subcrt(i, T=200)$species$ispecies == i + 1
      # so we should keep the original species index to be able to find the species in a provided 'sout'
      # (noted for Mosaic diagram section of anintro.Rmd 20190203)
      sout$species$ispecies <- species
      return(sout)
    }
  }

  ### Functions for logK/subcrt props
  # The logK contribution by any species or basis species
  X.species <- function(ispecies,coeff,X) coeff * sout$out[[ispecies]][,names(sout$out[[ispecies]])==X]
  # The logK contribution by all basis species in a reaction
  X.basis <- function(ispecies,X) Reduce("+", mapply(X.species,ibasis,-myspecies[ispecies,ibasis],X,SIMPLIFY=FALSE))
  # The logK of any reaction
  X.reaction <- function(ispecies,X) X.species((ispecies+nbasis),1,X) + X.basis(ispecies,X)
  # To get logK or subcrt props or other values into the dimensions we are using
  dim.fun <- function(x,idim=NULL) {
    if(is.null(idim)) {
      if(transect) idim <- 1
      else if(is.null(grid)) {
        # One of T,P,IS
        ivar <- which(vars %in% c("T","P","IS"))
        if(length(ivar)==0) ivar <- 1
        idim <- ivars(ivar)
      } else {
        # Two of T,P,IS
        ivar1 <- which(varissubcrt)[1]
        ivar2 <- which(varissubcrt)[2]
        idim <- ivars(ivar2, ivars(ivar1))
      }
    }
    return(aperm(array(x,mydim[idim]),idim))
  }
  # Properties of all species
  X.fun <- function(X) lapply(lapply(ispecies,X.reaction,X),dim.fun)
  logK <- function() lapply(ispecies,X.reaction,"logK")
  # A/2.303RT
  A <- function() {
    out <- mapply(`-`, X.fun("logK"), logQ(), SIMPLIFY = FALSE)
    # Deal with affinities of protein ionization here 20120527
    if("H+" %in% rownames(mybasis)) {
      # Which species are proteins
      isprotein <- grepl("_", myspecies$name)
      if(any(isprotein)) {
        if(get("thermo", CHNOSZ)$opt$ionize.aa) {
            message("affinity: ionizing proteins ...")
            # The rownumbers in thermo()$protein
            ip <- pinfo(myspecies$name[isprotein])
            # Get the affinity of ionization
            iHplus <- match("H+", rownames(mybasis))
            # as.numeric is needed in case the logact column is character mode
            # due to buffer definition (that means we can't do pH buffers)
            pH <- -as.numeric(mybasis$logact[iHplus])
            A.ionization <- A.ionization(ip, vars, vals, T = T, P = P, pH = pH, transect = transect)
            # Add it to the affinities of formation reactions of the non-ionized proteins
            out[isprotein] <- lsum(out[isprotein], A.ionization)
        } else message("affinity: NOT ionizing proteins because thermo()$opt$ionize.aa is FALSE")
      }
    }
    return(out)
  }

  ### Call the necessary functions
  if(!is.character(what)) {
    # Expand numeric values into our dimensions
    # (used by energy.args() for calculating pe=f(Eh,T) )
    # TODO: document that sout here denotes the dimension
    # we're expanding into
    return(dim.fun(what,ivars(sout$out)))
  } else if(what %in% c('G','H','S','Cp','V','E','kT','logK')) {
    # Get subcrt properties for reactions
    sout <- sout.fun(what)
    a <- X.fun(what)
  } else if(length(agrep("species",what)) > 0) {
    # Get subcrt properties for species
    # e.g. what=G.species, Cp.species etc
    mywhat <- s2c(what,sep=".",keep.sep=FALSE)[1]
    sout <- sout.fun(mywhat)
    a <- lapply(mapply(X.species,ispecies+nbasis,1,mywhat,SIMPLIFY=FALSE),dim.fun)
  } else {
    # Get affinities or logact.basis, logact.species or logQ
    if(what=="A") sout <- sout.fun("logK")
    a <- do.call(what,list())
  }

  ### Use species indices as the names 
  if(what=="logact.basis") names(a) <- mybasis$ispecies
  else names(a) <- myspecies$ispecies
  ### Return the results
  return(list(sout=sout,a=a))
}

energy.args <- function(args, transect = FALSE) {
  ## Extracted from affinity() and modified 20090329 jmd
  # Takes a list of arguments which define the dimensions
  # over which to calculate logQ, logK and affinity
  # The names should be T, P, IS and names of basis species
  # (or pH, pe, Eh)
  thermo <- get("thermo", CHNOSZ)
  ## Inputs are like c(T1,T2,res)
  # and outputs are like seq(T1,T2,length.out=res)
  # unless transect: do the variables specify a transect? 20090627
  transect <- any(transect, any(sapply(args,length) > 3))
  ## Convert T, P args and take care of 
  # Single values for T, P or IS
  T <- 298.15
  P <- "Psat"
  IS <- 0
  T.is.var <- P.is.var <- IS.is.var <- FALSE
  arg.is.T <- names(args)=="T"
  arg.is.P <- names(args)=="P"
  arg.is.IS <- names(args)=="IS"
  if(any(arg.is.T)) {
    T <- args[[which(arg.is.T)]]
    if(length(T) > 1) T.is.var <- TRUE
    if(transect) args[[which(arg.is.T)]] <- T <- envert(T,'K')
    else args[[which(arg.is.T)]][1:2] <- T[1:2] <- envert(T[1:2],'K')
    T <- T[!is.na(T)]
  }
  if(any(arg.is.P)) {
    P <- args[[which(arg.is.P)]]
    if(length(P) > 1) P.is.var <- TRUE
    if(!identical(P,"Psat")) {
      if(transect) args[[which(arg.is.P)]] <- P <- envert(P,'bar')
      else args[[which(arg.is.P)]][1:2] <- P[1:2] <- envert(P[1:2],'bar')
    }
    P <- P[!is.na(P)]
  }
  if(any(arg.is.IS)) {
    IS <- args[[which(arg.is.IS)]]
    if(length(IS) > 1) IS.is.var <- TRUE
  }
  # Report non-variables to user
  if(!T.is.var) {
    Tunits <- T.units()
    if(Tunits=="C") Tunits <- "\u00BAC"
    message('affinity: temperature is ', outvert(T, 'K'), ' ', Tunits)
  }
  if(!P.is.var) {
    if(identical(P,"Psat")) message("affinity: pressure is Psat")
    else message('affinity: pressure is ',outvert(P,'bar'),' ',P.units())
  }
  if(!IS.is.var & !identical(IS,0)) message('affinity: ionic strength is ',IS)
  # Default value for resolution
  res <- 256
  # Where we store the output
  what <- "A"
  vars <- character()
  vals <- list(NA)
  # This needs to have 1 as the third component because
  # energy() uses it to build an array with the given dimension
  lims <- list(c(NA, NA, 1))
  # Clean out non-variables
  if(any(arg.is.T) & !T.is.var) args <- args[names(args)!="T"]
  if(any(arg.is.P) & !P.is.var) args <- args[names(args)!="P"]
  if(any(arg.is.IS) & !IS.is.var) args <- args[names(args)!="IS"]
  # The property we're interested in
  if("what" %in% names(args)) {
    what <- args[[names(args)=="what"]]
    args <- args[names(args)!="what"]
  }
  # Assemble the variables
  if(length(args) > 0) {
    for(i in 1:length(args)) {
      nametxt <- names(args)[i]
      if(transect) lims.orig <- c(min(args[[i]]),max(args[[i]]))
      else lims.orig <- args[[i]][1:2]
      if(names(args)[i]=="pH") {
        names(args)[i] <- "H+"
        if(transect) args[[i]] <- -args[[i]]
        else args[[i]][1:2] <- -args[[i]][1:2]
      } 
      if(names(args)[i]=="pe") {
        names(args)[i] <- "e-"
        if(transect) args[[i]] <- -args[[i]]
        else args[[i]][1:2] <- -args[[i]][1:2]
      }
      if(length(args[[i]]) < 3 & !transect) args[[i]] <- c(args[[i]],res)
      vars[length(vars)+1] <- names(args)[i]
      if(transect) {
        vals[[length(vars)]] <- args[[i]]
        lims[[length(vars)]] <- c(lims.orig,length(vals[[i]]))
      } else {
        vals[[length(vars)]] <- seq(args[[i]][1],args[[i]][2],length.out=args[[i]][3])
        lims[[length(vars)]] <- args[[i]]
      }
      names(lims)[length(vars)] <- names(args)[i]
      # Say something about the identities, ranges, and units of the variables
      unittxt <- ""
      # Number of values
      if(transect) n <- length(args[[i]]) else n <- args[[i]][3]
      # Physical state
      ibasis <- match(nametxt, rownames(thermo$basis))
      if(isTRUE(as.logical(ibasis))) {
        if(thermo$basis$state[ibasis]=="gas") nametxt <- paste("log10(f_", nametxt, ")", sep="") 
        else nametxt <- paste("log10(a_", nametxt, ")", sep="") 
      } else {
        # Stop if the argument doesn't correspond to a basis species, T, P, or IS
        if(!nametxt %in% c("T", "P", "IS")) {
          if(! (nametxt=="pH" & 'H+' %in% rownames(thermo$basis) | nametxt %in% c("pe", "Eh") & 'e-' %in% rownames(thermo$basis))) {
            stop(nametxt, " is not one of T, P, or IS, and does not match any basis species")
          }
        }
      }
      # Temperature and pressure and Eh
      if(nametxt=="T") unittxt <- " K"
      if(nametxt=="P") unittxt <- " bar"
      if(nametxt=="Eh") unittxt <- " V"
      message("affinity: variable ", length(vars), " is ", nametxt, 
        " at ", n, " values from ", lims.orig[1], " to ", lims.orig[2], unittxt)
    }
  }
  args <- list(what=what,vars=vars,vals=vals,lims=lims,T=T,P=P,IS=IS,transect=transect)

  # Convert Eh to pe
  if("Eh" %in% args$vars) {
    # Get Eh into our dimensions
    Eh.args <- args
    # What variable is Eh
    Eh.var <- which(args$vars=="Eh")
    Eh.args$what <- args$vals[[Eh.var]]
    Eh.args$sout$out <- Eh.var
    Eh <- do.call("energy",Eh.args)
    # Get temperature into our dimensions
    T.args <- args  
    if("T" %in% args$vars) {
      T.var <- which(args$vars=="T")
      T.args$what <- args$vals[[T.var]]
    } else {
      T.var <- 1
      T.args$what <- T
    }
    T.args$sout$out <- T.var
    T <- do.call("energy",T.args)
    # Do the conversion on vectors
    mydim <- dim(Eh)
    Eh <- as.vector(Eh)
    T <- as.vector(T)
    pe <- convert(Eh,"pe",T=T)
    dim(pe) <- mydim
    # Update the arguments list
    args$vars[Eh.var] <- "e-"
    args$vals[[Eh.var]] <- -pe
  }
  return(args)
}

A.ionization <- function(iprotein, vars, vals, T=298.15, P="Psat", pH=7, transect=FALSE) {
  # A function to build a list of values of A/2.303RT of protein ionization
  # that can be used by energy(); 20120527 jmd
  # Some of the variables might not affect the values (e.g. logfO2)
  # What are the variables that affect the values
  T <- convert(T, "C")
  if(!is.na(iT <- match("T", vars))) T <- convert(vals[[iT]], "C")
  if(!is.na(iP <- match("P", vars))) P <- vals[[iP]]
  if(!is.na(iHplus <- match("H+", vars))) pH <- -vals[[iHplus]]
  # Is it a transect (single points) or a grid?
  if(transect) TPpH <- list(T=T, P=P, pH=pH)
  else {
    # Make a grid of all combinations
    # Put T, P, pH in same order as they show up vars
    egargs <- list(T=T, P=P, pH=pH)
    # order(c(NA, NA, NA)) might segfault in some versions of R (seen in R 2.15.3 on Linux)
    if(is.na(iT) & is.na(iP) & is.na(iHplus)) TPpHorder <- c(1, 2, 3)
    else TPpHorder <- order(c(iT, iP, iHplus))
    egargs <- c(egargs[TPpHorder], list(stringsAsFactors=FALSE))
    TPpH <- do.call(expand.grid, egargs)
    # Figure out the dimensions of T-P-pH (making sure to drop any that aren't in vars)
    TPpHdim <- numeric(3)
    TPpHdim[iT] <- length(T)
    TPpHdim[iP] <- length(P)
    TPpHdim[iHplus] <- length(pH)
    TPpHdim <- TPpHdim[!TPpHdim==0]
    # Figure out the dimensions of the other vars
    othervars <- vars[!vars %in% c("T", "P", "H+")]
    iother <- match(othervars, vars)
    otherdim <- sapply(vals, length)[iother]
    # If Eh was given to energy.args, values of pe were calculated in all dimensions
    # Figure out the original length of the Eh variable
    ieminus <- match("e-", vars[iother])
    if(!is.na(ieminus)) {
      otherdim[ieminus] <- NA
      edim <- dim(vals[[iother[ieminus]]])
      # Loop TPpHdim and otherdim, taking out each one from edim
      for(dim in c(TPpHdim, otherdim)) {
        id <- match(dim, edim)
        edim[id] <- NA
      }
      otherdim[ieminus] <- edim[!is.na(edim)]
    }
    # The permutation vector
    if(length(iother) > 0) allvars <- c(vars[-iother], vars[iother])
    else allvars <- vars
    perm <- match(vars, allvars)
  }
  # Initialize output list
  out <- vector("list", length(iprotein))
  # Get aa from iprotein
  aa <- pinfo(iprotein)
  # Calculate the values of A/2.303RT as a function of T-P-pH
  A <- ionize.aa(aa=aa, property="A", T=TPpH$T, P=TPpH$P, pH=TPpH$pH)
  if(transect) {
    # Just turn the values into a list
    for(i in 1:length(iprotein)) out[[i]] <- A[, i]
  } else for(i in 1:length(iprotein)) {
    # What are the values of ionization affinity for this protein
    thisA <- A[, i]
    # Apply the dimensions of T-P-pH
    tpphdim <- TPpHdim
    if(length(tpphdim)==0) tpphdim <- 1
    thisA <- array(thisA, tpphdim)
    # Grow into the dimensions of all vars
    alldim <- c(TPpHdim, otherdim)
    if(length(alldim)==0) alldim <- 1
    thisA <- array(thisA, alldim)
    # Now permute the array to put dimensions in same order as the variables
    thisA <- aperm(thisA, perm)
    # Store in output list
    out[[i]] <- thisA
  }
  return(out)
}
