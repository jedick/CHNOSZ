# CHNOSZ/util.misc.R
# some utility functions for the CHNOSZ package
# speciate/thermo.R 20051021 jmd

dPdTtr <- function(ispecies, ispecies2 = NULL) {
  # calculate dP/dT for a phase transition
  # (argument is index of the lower-T phase)
  thermo <- get("thermo")
  if(is.null(ispecies2)) ispecies2 <- ispecies + 1
  pars <- info(c(ispecies, ispecies2), check.it=FALSE)
  # if these aren't the same mineral, we shouldn't be here
  if(as.character(pars$name[1]) != as.character(pars$name[2])) stop("different names for species ", ispecies, " and ", ispecies2)
  # the special handling for quartz and coesite interfere with this function,
  # so we convert to uppercase names to prevent cgl() from calling quartz_coesite()
  pars$name <- toupper(pars$name)
  props <- cgl(c("G", "S", "V"), pars, P=0, T=thermo$obigt$z.T[ispecies])
  # we really hope the G's are the same ...
  #if(abs(props[[2]]$G - props[[1]]$G) > 0.1) warning('dP.dT: inconsistent values of G for different phases of ',ispecies,call.=FALSE)
  dP.dT <- convert( ( props[[2]]$S - props[[1]]$S ) / ( props[[2]]$V - props[[1]]$V ), 'cm3bar' )
  return(dP.dT)
}

Ttr <- function(ispecies, ispecies2 = NULL, P = 1, dPdT = NULL) {
  # calculate a phase transition temperature for given P
  TtrPr <- get("thermo")$obigt$z.T[ispecies]
  # the constant slope, dP/dT
  if(is.null(dPdT)) dPdT <- dPdTtr(ispecies, ispecies2)
  Pr <- 1
  TtrPr + (P - Pr) / dPdT
}

GHS_Tr <- function(ispecies, Htr) {
  # calculate G, H, and S at Tr for cr2, cr3, ... phases 20170301
  # Htr: enthalpy(ies) of transition
  # ispecies: the species index for cr (the lowest-T phase)
  thisinfo <- info(ispecies)
  name <- thisinfo$name
  # start from Tr (T=298.15 K)
  Tprev <- 298.15
  # the GHS at T
  Gf <- thisinfo$G
  Hf <- thisinfo$H
  S <- thisinfo$S
  # where to store the calculated GHS at Tr
  Gf_Tr <- Hf_Tr <- S_Tr <- numeric()
  for(i in 1:(length(Htr)+1)) {
    # check that we have the correct one of cr, cr2, cr3, ...
    if(i==1) thiscr <- "cr" else thiscr <- paste0("cr", i)
    if(thisinfo$state!=thiscr | thisinfo$name!=name) stop(paste("species", thisis, "is not", name, thiscr))
    # if we're above cr (lowest-T), calculate the equivalent GHS at Tr
    if(i > 1) {
      # set the starting GHS to 0 (in case they're NA - we only need the increments over temperature)
      thisinfo$G <- thisinfo$H <- thisinfo$S <- 0
      # the HS increments from 298.15 to Ttr
      HSinc <- cgl(c("H", "S"), parameters=thisinfo, T=c(298.15, Ttr))
      Hf_Tr <- c(Hf_Tr, Hf - diff(HSinc[[1]]$H))
      S_Tr <- c(S_Tr, S - diff(HSinc[[1]]$S))
      # plug in the calculated S_Tr to calculate the G increment correctly
      thisinfo$S <- tail(S_Tr, 1)
      Ginc <- cgl("G", parameters=thisinfo, T=c(298.15, Ttr))
      Gf_Tr <- c(Gf_Tr, Gf - diff(Ginc[[1]]$G))
    }
    # the temperature of the next transition
    Ttr <- thisinfo$T
    # the GHS increments from Tprev to Ttr
    GHCinc <- cgl(c("G", "H", "S"), parameters=thisinfo, T=c(Tprev, Ttr))
    # the GHS + transition at Tr
    Gf <- Gf + diff(GHCinc[[1]]$G)
    Hf <- Hf + diff(GHCinc[[1]]$H) + Htr[i]
    S <- S + diff(GHCinc[[1]]$S) + Htr[i] / Ttr
    # prepare next phase
    thisis <- ispecies + i
    thisinfo <- info(thisis)
    Tprev <- Ttr
  }
  list(Gf_Tr=Gf_Tr, Hf_Tr=Hf_Tr, S_Tr=S_Tr)
}

unitize <- function(logact=NULL,length=NULL,logact.tot=0) {
  # scale the logarithms of activities given in loga
  # so that the logarithm of total activity of residues
  # is zero (i.e. total activity of residues is one),
  # or some other value set in loga.tot.
  # length indicates the number of residues in each species.
  # if loga is NULL, take the logarithms of activities from
  # the current species definition. if any of those species
  # are proteins, get their lengths using protein.length.
  thermo <- get("thermo")
  if(is.null(logact)) {
    if(is.null(thermo$species)) stop("loga is NULL and no species are defined")
    ts <- thermo$species
    logact <- ts$logact
    length <- rep(1,length(logact))
    ip <- grep("_",ts$name)
    if(length(ip) > 0) length[ip] <- protein.length(ts$name[ip])
  }
  # the lengths of the species
  if(is.null(length)) length <- 1
  length <- rep(length,length.out=length(logact)) 
  # remove the logarithms
  act <- 10^logact
  # the total activity
  act.tot <- sum(act*length)
  # the target activity
  act.to.get <- 10^logact.tot
  # the factor to apply
  act.fact <- act.to.get/act.tot
  # apply the factor
  act <- act * act.fact
  # take the logarithms
  log10(act)
  # done!
}

### unexported functions ###

# return, in order, which column(s) of species all have non-zero values.
which.balance <- function(species) {
  # find the first basis species that
  # is present in all species of interest
  # ... it can be used to balance the system
  nbasis <- function(species) return(ncol(species)-4)
  ib <- NULL
  nb <- 1
  nbs <- nbasis(species)
  for(i in 1:nbs) {
    coeff <- species[,i]
    if(length(coeff)==length(coeff[coeff!=0])) {
      ib <- c(ib,i)
      nb <- nb + 1
    } else ib <- c(ib,NA)
  }
  return(ib[!is.na(ib)])
}
