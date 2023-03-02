# CHNOSZ/util.misc.R
# Some utility functions for the CHNOSZ package
# speciate/thermo.R 20051021 jmd

dPdTtr <- function(ispecies, ispecies2 = NULL) {
  # Calculate dP/dT for a phase transition
  # (argument is index of the lower-T phase)
  thermo <- get("thermo", CHNOSZ)
  if(is.null(ispecies2)) ispecies2 <- ispecies + 1
  pars <- info(c(ispecies, ispecies2), check.it=FALSE)
  # If these aren't the same mineral, we shouldn't be here
  if(as.character(pars$name[1]) != as.character(pars$name[2])) stop("different names for species ", ispecies, " and ", ispecies2)
  # The special handling for quartz and coesite interfere with this function,
  # so we convert to uppercase names to prevent cgl() from calling quartz_coesite()
  pars$name <- toupper(pars$name)
  props <- cgl(c("G", "S", "V"), pars, P=0, T=thermo$OBIGT$z.T[ispecies])
  # The G's should be the same ...
  #if(abs(props[[2]]$G - props[[1]]$G) > 0.1) warning('dP.dT: inconsistent values of G for different phases of ',ispecies,call.=FALSE)
  dP.dT <- convert( ( props[[2]]$S - props[[1]]$S ) / ( props[[2]]$V - props[[1]]$V ), 'cm3bar' )
  return(dP.dT)
}

Ttr <- function(ispecies, ispecies2 = NULL, P = 1, dPdT = NULL) {
  # Calculate a phase transition temperature for given P
  TtrPr <- get("thermo", CHNOSZ)$OBIGT$z.T[ispecies]
  # The constant slope, dP/dT
  if(is.null(dPdT)) dPdT <- dPdTtr(ispecies, ispecies2)
  Pr <- 1
  TtrPr + (P - Pr) / dPdT
}

GHS_Tr <- function(ispecies, Htr) {
  # Calculate G, H, and S at Tr for cr2, cr3, ... phases 20170301
  # Htr: enthalpy(ies) of transition
  # ispecies: the species index for cr (the lowest-T phase)
  thisinfo <- info(ispecies)
  name <- thisinfo$name
  # Start from Tr (T=298.15 K)
  Tprev <- 298.15
  # The GHS at T
  Gf <- thisinfo$G
  Hf <- thisinfo$H
  S <- thisinfo$S
  # Where to store the calculated GHS at Tr
  Gf_Tr <- Hf_Tr <- S_Tr <- numeric()
  for(i in 1:(length(Htr)+1)) {
    # Check that we have the correct one of cr, cr2, cr3, ...
    if(i==1) thiscr <- "cr" else thiscr <- paste0("cr", i)
    if(thisinfo$state!=thiscr | thisinfo$name!=name) stop(paste("species", thisis, "is not", name, thiscr))
    # If we're above cr (lowest-T), calculate the equivalent GHS at Tr
    if(i > 1) {
      # Set the starting GHS to 0 (in case they're NA - we only need the increments over temperature)
      thisinfo$G <- thisinfo$H <- thisinfo$S <- 0
      # The HS increments from 298.15 to Ttr
      HSinc <- cgl(c("H", "S"), parameters=thisinfo, T=c(298.15, Ttr))
      Hf_Tr <- c(Hf_Tr, Hf - diff(HSinc[[1]]$H))
      S_Tr <- c(S_Tr, S - diff(HSinc[[1]]$S))
      # Plug in the calculated S_Tr to calculate the G increment correctly
      thisinfo$S <- tail(S_Tr, 1)
      Ginc <- cgl("G", parameters=thisinfo, T=c(298.15, Ttr))
      Gf_Tr <- c(Gf_Tr, Gf - diff(Ginc[[1]]$G))
    }
    # The temperature of the next transition
    Ttr <- thisinfo$T
    # The GHS increments from Tprev to Ttr
    GHCinc <- cgl(c("G", "H", "S"), parameters=thisinfo, T=c(Tprev, Ttr))
    # The GHS + transition at Tr
    Gf <- Gf + diff(GHCinc[[1]]$G)
    Hf <- Hf + diff(GHCinc[[1]]$H) + Htr[i]
    S <- S + diff(GHCinc[[1]]$S) + Htr[i] / Ttr
    # Prepare next phase
    thisis <- ispecies + i
    thisinfo <- info(thisis)
    Tprev <- Ttr
  }
  list(Gf_Tr=Gf_Tr, Hf_Tr=Hf_Tr, S_Tr=S_Tr)
}

unitize <- function(logact=NULL,length=NULL,logact.tot=0) {
  # Scale the logarithms of activities given in loga
  # so that the logarithm of total activity of residues
  # is zero (i.e. total activity of residues is one),
  # or some other value set in loga.tot.
  # length indicates the number of residues in each species.
  # If loga is NULL, take the logarithms of activities from the current species definition
  # If any of those species are proteins, get their lengths using protein.length
  thermo <- get("thermo", CHNOSZ)
  if(is.null(logact)) {
    if(is.null(thermo$species)) stop("loga is NULL and no species are defined")
    ts <- thermo$species
    logact <- ts$logact
    length <- rep(1,length(logact))
    ip <- grep("_",ts$name)
    if(length(ip) > 0) length[ip] <- protein.length(ts$name[ip])
  }
  # The lengths of the species
  if(is.null(length)) length <- 1
  length <- rep(length,length.out=length(logact)) 
  # Remove the logarithms
  act <- 10^logact
  # The total activity
  act.tot <- sum(act*length)
  # The target activity
  act.to.get <- 10^logact.tot
  # The factor to apply
  act.fact <- act.to.get/act.tot
  # Apply the factor
  act <- act * act.fact
  # Take the logarithms
  log10(act)
  # Done!
}

### Unexported functions ###

# Return, in order, which column(s) of species all have non-zero values.
which.balance <- function(species) {
  # Find the first basis species that
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
