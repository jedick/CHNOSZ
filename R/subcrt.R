# CHNOSZ/subcrt.R
# Calculate standard molal thermodynamic propertes
# 20060817 jmd

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.args.R")
#source("util.character.R")
#source("info.R")
#source("util.units.R")
#source("util.data.R")
#source("species.R")
#source("AD.R")
#source("nonideal.R")
#source("hkf.R")
#source("cgl.R")

subcrt <- function(species, coeff = 1, state = NULL, property = c("logK", "G", "H", "S", "V", "Cp"),
  T = seq(273.15, 623.15, 25), P = "Psat", grid = NULL, convert = TRUE, exceed.Ttr = FALSE,
  exceed.rhomin = FALSE, logact = NULL, autobalance = TRUE, use.polymorphs = TRUE, IS = 0) {

  # Revise the call if the states are the second argument 
  if(!is.null(coeff[1])) {
    if(is.character(coeff[1])) {
      newstate <- coeff
      # This is missing coeff and T in order that missing values are correctly detected further below 20230621
      newargs <- list(species = species, state = newstate,
        property = property, P = P, grid = grid, convert = convert,
        exceed.Ttr = exceed.Ttr, exceed.rhomin = exceed.rhomin, logact = logact,
        autobalance = autobalance, use.polymorphs = use.polymorphs, IS = IS)
      if(!missing(state)) {
        if(is.numeric(state[1])) newcoeff <- state else stop("If they are both given, one of arguments 2 and 3 should be numeric reaction coefficients")
        newargs <- c(list(coeff = newcoeff), newargs)
      }
      if(!missing(T)) newargs <- c(list(T = T), newargs)
      return(do.call(subcrt, newargs))
    }
  }

  do.reaction <- FALSE
  if(!missing(coeff)) do.reaction <- TRUE

  # Species and states are made the same length
  if(!is.null(state[1])) {
    if(length(state) > length(species)) species <- rep(species, length.out = length(state))
    if(length(species) > length(state)) state <- rep(state, length.out = length(species))
    state <- state.args(state)
  }

  # Allowed properties
  properties <- c("rho", "logK", "G", "H", "S", "Cp", "V", "kT", "E")
  # Property checking
  calcprop <- property
  notprop <- property[!calcprop %in% properties]
  if(length(notprop) == 1) stop(paste("invalid property name:", paste(notprop, collapse = " ")))
  if(length(notprop) > 1) stop(paste("invalid property names:", paste(notprop, collapse = " ")))
  # Length checking
  if(do.reaction & length(species)!=length(coeff)) 
    stop("the length of 'coeff' must equal the number of species")
  if(!is.null(logact)) {
    if(length(logact) != length(species)) stop("the length of 'logact' must equal the number of species")
  }
  # Normalize temperature and pressure units
  if(!missing(T)) {
    if(convert) T <- envert(T, "K")
    else if(!missing(convert) & convert) T <- envert(T, "K")
  }
  if(is.numeric(P[1])) {
    if(convert) P <- envert(P, "bar")
  }

  # Warn for too high temperatures for Psat 20171110
  warnings <- character()
  if(identical(P, "Psat") & any(T > 647.067)) {
    nover <- sum(T > 647.067)
    if(nover==1) vtext <- "value" else vtext <- "values"
    warnings <- c(warnings, paste0("P = 'Psat' undefined for T > Tcritical (", nover, " T ", vtext, ")"))
  }

  # Are we gridding?
  isPsat <- FALSE
  do.grid <- FALSE
  if(!is.null(grid)) if(!is.logical(grid)) do.grid <- TRUE
  newIS <- IS
  if(do.grid) {
    if(grid == "T") {
      newT <- numeric()
      for(i in 1:length(T)) newT <- c(newT, rep(T[i], length(P)))
      newP <- rep(P, length(T))
      T <- newT; P <- newP
    }
    if(grid == "P") {
      newP <- numeric()
      for(i in 1:length(P)) newP <- c(newP, rep(P[i], length(T)))
      newT <- rep(T, length(P))
      T <- newT; P <- newP
    }
    if(grid == "IS") {
      ll <- length(T)
      if(length(P) > 1) ll <- length(P)
      newIS <- numeric()
      for(i in 1:length(IS)) newIS <- c(newIS, rep(IS[i], ll))
      tpargs <- TP.args(T = T, P = P)
      T <- rep(tpargs$T, length.out = length(newIS))
      P <- rep(tpargs$P, length.out = length(newIS))
    }
  } else {
    # For AD, remember if P = "Psat" 20190219
    if(identical(P, "Psat")) isPsat <- TRUE
    # Expansion of Psat and equivalence of argument lengths
    tpargs <- TP.args(T = T,P = P)
    T <- tpargs$T; P <- tpargs$P
    if(length(newIS) > length(T)) T <- rep(T, length.out = length(newIS))
    if(length(newIS) > length(P)) P <- rep(P, length.out = length(newIS))
  }

  # Get species information
  thermo <- get("thermo", CHNOSZ)
  # Before 20110808, we sent numeric species argument through info() to get species name and state(s)
  # But why slow things down if we already have a species index?
  if(is.numeric(species[1])) {
    ispecies <- species
    species <- as.character(thermo$OBIGT$name[ispecies])
    state <- as.character(thermo$OBIGT$state[ispecies])
    newstate <- as.character(thermo$OBIGT$state[ispecies])
  } else {
    # Species are named ... look up the index
    ispecies <- numeric()
    newstate <- character()
    for(i in 1:length(species)) {
      # Get the species index for a named species
      if(!can.be.numeric(species[i])) sindex <- info.character(species[i], state[i])
      else {
        # Check that a coerced-to-numeric argument is a rownumber of thermo()$OBIGT
        sindex <- as.numeric(species[i])
        if(!sindex %in% 1:nrow(thermo$OBIGT)) stop(paste(species[i], "is not a row number of thermo()$OBIGT"))
      }
      # info.character() has the possible side-effect of adding a protein; re-read thermo to use the (possible) additions
      thermo <- get("thermo", CHNOSZ)
      if(is.na(sindex[1])) stop("no info found for ", species[i], " ",state[i])
      if(!is.null(state[i])) is.cr <- state[i]=="cr" else is.cr <- FALSE
      # If we found multiple species, take the first one (useful for minerals with polymorphs)
      if(thermo$OBIGT$state[sindex[1]] == "cr" & (is.null(state[i]) | is.cr)) {
        newstate <- c(newstate, "cr")
        ispecies <- c(ispecies, sindex[1])
      } else {
        newstate <- c(newstate, as.character(thermo$OBIGT$state[sindex[1]]))
        ispecies <- c(ispecies, sindex[1])
      }
    }
  }

  # Stop if species not found
  noname <- is.na(ispecies)
  if(TRUE %in% noname)
    stop(paste("species", species[noname], "not found.\n"))

  # Take care of mineral phases
  state <- as.character(thermo$OBIGT$state[ispecies])
  name <- as.character(thermo$OBIGT$name[ispecies])
  # A counter of all species considered
  # iphases is longer than ispecies if multiple polymorphs (cr, cr2, ...) are present
  # polymorph.species shows which of ispecies correspond to iphases
  # Before 20091114: the success of this depends on there not being duplicated aqueous or other
  #   non-mineral species (i.e., two entries in OBIGT for Cu+ mess this up when running the skarn example).
  # After 20091114: we can deal with duplicated species (aqueous at least)
  iphases <- polymorph.species <- coeff.new <- numeric()
  for(i in 1:length(ispecies)) {
     # Add check for use.polymorphs argument 20230620
     if(newstate[i] == "cr" & use.polymorphs) {
       # Check for available polymorphs in OBIGT
       polymorph.states <- c("cr", "cr2", "cr3", "cr4", "cr5", "cr6", "cr7", "cr8", "cr9") 
       tghs <- thermo$OBIGT[(thermo$OBIGT$name %in% name[i]) & thermo$OBIGT$state %in% polymorph.states, ]
       # Only take the first one if they are duplicated non-mineral species (i.e., not polymorphs)
       if(all(tghs$state == tghs$state[1])) tghs <- thermo$OBIGT[ispecies[i], ]
     } else tghs <- thermo$OBIGT[ispecies[i], ]
     iphases <- c(iphases, as.numeric(rownames(tghs))) 
     polymorph.species <- c(polymorph.species, rep(ispecies[i], nrow(tghs)))
     coeff.new <- c(coeff.new, rep(coeff[i], nrow(tghs)))
  }

  # Where we keep info about the species involved
  # Add model information 20220919
  model <- thermo$OBIGT$model[iphases]
  # Label specific water model;
  # this is also how we record a "wet" reaction
  isH2O <- model == "H2O"
  isH2O[is.na(isH2O)] <- FALSE
  model[isH2O] <- paste0("water.", thermo$opt$water)
  reaction <- data.frame(coeff = coeff.new, name = thermo$OBIGT$name[iphases],
    formula = thermo$OBIGT$formula[iphases], state = thermo$OBIGT$state[iphases],
    ispecies = iphases, model = model, stringsAsFactors = FALSE)
  # Make the rownames readable ... but they have to be unique
  if(length(unique(iphases))==length(iphases)) rownames(reaction) <- as.character(iphases)
  # Which species use models for aqueous species
  # This breaks things if we have state = "aq" and model = "CGL" 20230220
  #isaq <- reaction$state == "aq"
  isaq <- toupper(reaction$model) %in% c("HKF", "AD", "DEW")

  # Produce message about conditions
  if(length(species)==1 & convert==FALSE) {
    # No message produced here (internal calls from other functions)
  } else {
    # Include units here 20190530
    uT <- outvert(T, "K")
    if(identical(grid, "IS")) uT <- unique(uT)
    Tunits <- T.units()
    if(Tunits == "C") Tunits <- "\u00BAC"
    if(length(uT) == 1) T.text <- paste(uT, Tunits) else {
      T.text <- paste0(length(uT), " values of T (", Tunits, ")")
    }
    uP <- outvert(P, "bar")
    if(length(P) == 1) {
      if(can.be.numeric(P)) P.text <- paste(round(as.numeric(uP),2), P.units())
      else P.text <- paste0("P (", P.units(), ")")
    } else P.text <- paste0("P (", P.units(), ")")
    if(identical(P[[1]], "Psat")) P.text <- P
    if(any(c(isH2O, isaq))) P.text <- paste(P.text, " (wet)", sep = "")
    E.text <- paste0("[energy units: ", E.units(), "]")
    message(paste("subcrt:", length(species), "species at", T.text, "and", P.text, E.text))
  }

  # Inform about unbalanced reaction
  if(do.reaction) {
    # The mass balance; should be zero for a balanced reaction
    mss <- makeup(ispecies, coeff, sum = TRUE)
    # Take out very small numbers
    mss[abs(mss) < 1e-7] <- 0
    # Report and try to fix any non-zero mass balance
    if(any(mss != 0)) {
      # The missing composition: the negative of the mass balance
      miss <- -mss
      # Drop elements that are zero
      miss <- miss[miss != 0]
      message("subcrt: reaction is not balanced; it is missing this composition:")
      # We have to do this awkward dance to send a formatted matrix to message
      message(paste(capture.output(print(miss)), collapse = "\n"))
      # Look for basis species that have our compositoin
      tb <- thermo$basis
      if(!is.null(tb) & autobalance) {
        if(all(names(miss) %in% colnames(tb)[1:nrow(tb)])) {
          # The missing composition in terms of the basis species
          bc <- species.basis(species = NULL, mkp = as.matrix(miss))
          # Drop zeroes
          bc.new <- bc[, (bc[1, ] != 0), drop = FALSE]
          # and get the states
          b.state <- as.character(thermo$basis$state)[bc[1, ] != 0]
          bc <- bc.new
          # Special thing for Psat
          if(identical(P[[1]], "Psat")) P <- "Psat"
          else P <- outvert(P, "bar")
          # Add to logact values if present
          if(!is.null(logact)) {
            ila <- match(colnames(bc), rownames(thermo$basis))
            nla <- !(can.be.numeric(thermo$basis$logact[ila]))
            if(any(nla)) warning("subcrt: logact values of basis species",
              c2s(rownames(thermo$basis)[ila]), "are NA.")
            logact <- c(logact, thermo$basis$logact[ila])
          }
          # Warn user and do it!
          ispecies.new <- tb$ispecies[match(colnames(bc),rownames(tb))]
          b.species <- thermo$OBIGT$formula[ispecies.new]
          if(identical(species,b.species) & identical(state,b.state))
            message("subcrt: balanced reaction, but it is a non-reaction; restarting...")
          else message("subcrt: adding missing composition from basis definition and restarting...")
          newspecies <- c(species, tb$ispecies[match(colnames(bc), rownames(tb))])
          newcoeff <- c(coeff, as.numeric(bc[1, ]))
          newstate <- c(state, b.state)
          return(subcrt(species = newspecies, coeff = newcoeff, state = newstate,
            property = property, T = outvert(T, "K"), P = P, grid = grid, convert = convert, logact = logact, exceed.Ttr = FALSE))
        } else warnings <- c(warnings, paste("reaction among", paste(species, collapse = ","), "was unbalanced, missing", as.chemical.formula(miss)))
      } else warnings <- c(warnings, paste("reaction among", paste(species, collapse = ","), "was unbalanced, missing", as.chemical.formula(miss)))
    }
  }

  # Calculate the properties
  # If we want affinities we must have logK; include it in the ouput
  if(!is.null(logact)) if(!"logK" %in% calcprop) calcprop <- c("logK", calcprop)
  # If logK but not G was requested, we need to calculate G
  eosprop <- calcprop
  if("logK" %in% calcprop & ! "G" %in% calcprop) eosprop <- c(eosprop, "G")
  # Also get G if we are dealing with mineral phases
  if(!"G" %in% eosprop & length(iphases) > length(ispecies)) eosprop <- c(eosprop, "G")
  # Don't request logK or rho from the eos ...
  eosprop <- eosprop[!eosprop %in% c("logK", "rho")]
  # The reaction result will go here
  outprops <- list()
  # Aqueous species and H2O properties
  if(TRUE %in% isaq) {
    # 20110808 get species parameters using OBIGT2eos()
    # (this is faster than using info() and is how we get everything in the same units)
    param <- OBIGT2eos(thermo$OBIGT[iphases[isaq], ], "aq", fixGHS = TRUE, toJoules = TRUE)
    # Aqueous species with model = "AD" use the AD model 20210407
    model <- thermo$OBIGT$model[iphases[isaq]]
    model[is.na(model)] <- ""
    isAD <- model == "AD"
    # Always get density
    H2O.props <- "rho"
    # Calculate A_DH and B_DH if we're using the B-dot (Helgeson) equation
    if(any(IS != 0) & thermo$opt$nonideal %in% c("Bdot", "Bdot0", "bgamma", "bgamma0")) H2O.props <- c(H2O.props, "A_DH", "B_DH")
    # Get other properties for H2O only if it's in the reaction
    if(any(isH2O)) H2O.props <- c(H2O.props, eosprop)
    # In case everything is AD, run hkf (for water properties) but exclude all species
    hkfpar <- param
    if(all(isAD)) hkfpar <- param[0, ]
    hkfstuff <- hkf(eosprop, parameters = hkfpar, T = T, P = P, H2O.props = H2O.props)
    p.aq <- hkfstuff$aq
    H2O.PT <- hkfstuff$H2O
    # Set properties to NA for density below 0.35 g/cm3 (a little above the critical isochore, threshold used in SUPCRT92) 20180922
    if(!exceed.rhomin & !all(isAD)) {
      ilowrho <- H2O.PT$rho < 350
      ilowrho[is.na(ilowrho)] <- FALSE
      if(any(ilowrho)) {
        for(i in 1:length(p.aq)) p.aq[[i]][ilowrho, ] <- NA
        if(sum(ilowrho) == 1) ptext <- "pair" else ptext <- "pairs"
        warnings <- c(warnings, paste0("below minimum density for applicability of revised HKF equations (", sum(ilowrho), " T,P ", ptext, ")"))
      }
    }
    # Calculate properties using Akinfiev-Diamond model 20190219
    if(any(isAD)) {
      # get the parameters with the right names
      param <- OBIGT2eos(param[isAD, ], "aq", toJoules = TRUE)
      p.aq[isAD] <- AD(eosprop, parameters = param, T = T, P = P, isPsat = isPsat)
    }
    # Calculate activity coefficients if ionic strength is not zero
    if(any(IS != 0)) {
      if(thermo$opt$nonideal %in% c("Bdot", "Bdot0", "bgamma", "bgamma0")) p.aq <- nonideal(iphases[isaq], p.aq, newIS, T, P, H2O.PT$A_DH, H2O.PT$B_DH)
      else if(thermo$opt$nonideal=="Alberty") p.aq <- nonideal(iphases[isaq], p.aq, newIS, T)
    }
    outprops <- c(outprops, p.aq)
  } else if(any(isH2O)) {
    # We're not using the HKF, but still want water
    H2O.PT <- water(c("rho", eosprop), T = T, P = P)
  }

  # Crystalline, gas, or liquid (except water) species
  iscgl <- reaction$model %in% c("CGL", "Berman")

  if(TRUE %in% iscgl) {
    param <- OBIGT2eos(thermo$OBIGT[iphases[iscgl],], "cgl", fixGHS = TRUE, toJoules = TRUE)
    p.cgl <- cgl(eosprop, parameters = param, T = T, P = P)
    # Replace Gibbs energies with NA where the
    # phases are beyond their temperature range
    if("G" %in% eosprop) {
      # 20080304 This code is weird and hard to read - needs a lot of cleanup!
      # 20120219 Cleaned up somewhat; using exceed.Ttr and NA instead of do.phases and 999999
      # the numbers of the cgl species (becomes 0 for any that aren't cgl)
      ncgl <- iscgl
      ncgl[iscgl] <- 1:nrow(param)
      for(i in 1:length(iscgl)) {
        # Not if we're not cgl
        if(!iscgl[i]) next
        # Name and state
        myname <- reaction$name[i]
        mystate <- reaction$state[i]
        # If we are considering multiple polymorphs, and if this polymorph is cr2 or higher, check if we're below the transition temperature
        if(length(iphases) > length(ispecies) & i > 1) {
          if(!(reaction$state[i] %in% c("liq", "cr", "gas")) & reaction$name[i-1] == reaction$name[i]) {
            # After add.OBIGT("SUPCRT92"), quartz cr and cr2 are not next to each other in thermo()$OBIGT,
            # so use iphases[i-1] here, not iphases[i]-1  20181107
            Ttr <- Ttr(iphases[i-1], iphases[i], P = P, dPdT = dPdTtr(iphases[i-1], iphases[i]))
            if(all(is.na(Ttr))) next
            if(any(T <= Ttr)) {
              status.Ttr <- "(extrapolating G)"
              if(!exceed.Ttr) {
                # put NA into the value of G
                p.cgl[[ncgl[i]]]$G[T <= Ttr] <- NA
                status.Ttr <- "(using NA for G)"
              } 
              #message(paste("subcrt: some points below transition temperature for", myname, mystate, status.Ttr))
            }
          }
        }
        # Check if we're above the temperature limit or transition temperature
        # T limit (or Ttr) from the database
        warn.above <- TRUE
        Ttr <- thermo$OBIGT$z.T[iphases[i]]
        # Calculate Ttr at higher P if a phase transition is present
        if(i < nrow(reaction)) {
          # If the next one is cr2, cr3, etc we have a transition
          if(reaction$state[i+1] %in% c("cr1", "cr2", "cr3", "cr4", "cr5", "cr6", "cr7", "cr8", "cr9") & reaction$name[i+1] == reaction$name[i]) {
            Ttr <- Ttr(iphases[i], iphases[i+1], P = P, dPdT = dPdTtr(iphases[i], iphases[i+1]))
            # We don't warn here about the transition
            warn.above <- FALSE
          }
        }
        if(any(is.na(Ttr))) next
        if(!all(Ttr == 0) & any(T > Ttr)) {
          status.Ttr <- "(extrapolating G)"
          if(!exceed.Ttr) {
            p.cgl[[ncgl[i]]]$G[T > Ttr] <- NA
            status.Ttr <- "(using NA for G)"
          }
          Tmax <- min(T[T > Ttr])
          if(warn.above) message(paste("subcrt: temperature(s) of", Tmax, "K and above exceed limit for", myname, mystate, status.Ttr))
        }
        # Use variable-pressure standard Gibbs energy for gases if varP is TRUE (not the default)
        if(mystate == "gas" & thermo$opt$varP) p.cgl[[ncgl[i]]]$G <- p.cgl[[ncgl[i]]]$G - convert(log10(P), "G", T = T)
      }
    }
    outprops <- c(outprops,p.cgl)
  }

  # Water
  if(any(isH2O)) {
    p.H2O <- H2O.PT[, match(eosprop, colnames(H2O.PT)), drop = FALSE]
    p.H2O <- list(p.H2O)
    outprops <- c(outprops, rep(p.H2O, sum(isH2O == TRUE)))
  }

  # logK
  if("logK" %in% calcprop) {
    for(i in 1:length(outprops)) {
      outprops[[i]] <- cbind(outprops[[i]],data.frame(logK = convert(outprops[[i]]$G, "logK", T = T)))
      colnames(outprops[[i]][ncol(outprops[[i]])]) <- "logK"
    }
  }

  # Ordering the output
  # The indices of the species in outprops thus far
  ns <- 1:nrow(reaction)
  is <- c(ns[isaq],ns[iscgl],ns[isH2O])
  if(length(ns) != length(is)) stop("subcrt: not all species are accounted for.")
  v <- list()
  for(i in 1:length(is)) v[[i]] <- outprops[[match(ns[i], is)]]
  outprops <- v

  # Deal with polymorphs (cr,cr2) here:
  # We have to eliminate rows from outprops, 
  # reaction and values from isaq, iscgl, isH2O
  out.new <- list()
  reaction.new <- reaction
  isaq.new <- logical()
  iscgl.new <- logical()
  isH2O.new <- logical()
  for(i in 1:length(ispecies)) {
    are.polymorphs <- which(ispecies[i]==polymorph.species)
    # Deal with repeated species here
    if(TRUE %in% duplicated(iphases[are.polymorphs])) {
      # Only take the first, not the duplicates
      ndups <- sum(ispecies==ispecies[i])
      npolymorphs <- length(are.polymorphs) / ndups
      are.polymorphs <- are.polymorphs[1:npolymorphs]
    }
    if(length(are.polymorphs) > 1) {
      message(paste("subcrt:", length(are.polymorphs), "polymorphs for", thermo$OBIGT$name[ispecies[i]], "... "), appendLF = FALSE)
      # Assemble the Gibbs energies for each species
      for(j in 1:length(are.polymorphs)) {
        G.this <- outprops[[are.polymorphs[j]]]$G
        if(sum(is.na(G.this)) > 0 & exceed.Ttr) warning(paste("subcrt: NAs found for G of ",
          reaction$name[are.polymorphs[j]], " ", reaction$state[are.polymorphs[j]], " at T-P point(s) ", 
          c2s(which(is.na(G.this)), sep = " "), sep = ""), call. = FALSE)
        if(j == 1) G <- as.data.frame(G.this)
        else G <- cbind(G, as.data.frame(G.this))
      }
      # Find the minimum-energy polymorph at each T-P point
      stable.polymorph <- numeric()
      out.new.entry <- outprops[[are.polymorphs[1]]]
      for(j in 1:nrow(G)) {
        ps <- which.min(as.numeric(G[j, ]))
        if(length(ps)==0) {
          # minimum not found (we have NAs)
          # - no non-NA value of G to begin with (e.g. aegerine) --> probably should use lowest-T phase
          #ps <- 1
          # - above temperature limit for the highest-T phase (subcrt.Rd skarn example) --> use highest-T phase 20171110
          ps <- ncol(G)
          if(exceed.Ttr) warning("subcrt: stable polymorph for ", reaction$name[are.polymorphs[ps]], " at T-P point ", j, 
          " undetermined (using ", reaction$state[are.polymorphs[ps]], ")", call. = FALSE)
        } 
        stable.polymorph <- c(stable.polymorph, ps)
        out.new.entry[j, ] <- outprops[[ are.polymorphs[ps] ]][j, ]
      }

      # Update our objects
      out.new[[i]] <- cbind(out.new.entry, data.frame(polymorph = stable.polymorph))
      reaction.new[i, ] <- reaction[are.polymorphs[stable.polymorph[1]], ]
      # Mark the minerals with multiple polymorphs
      reaction.new$state[i] <- "cr*"
      isaq.new <- c(isaq.new, isaq[are.polymorphs[stable.polymorph[1]]])
      iscgl.new <- c(iscgl.new, iscgl[are.polymorphs[stable.polymorph[1]]])
      isH2O.new <- c(isH2O.new, isH2O[are.polymorphs[stable.polymorph[1]]])
      # Info for the user
      up <- unique(stable.polymorph)
      if(length(up) > 1) { word <- "are"; p.word <- "polymorphs" }
      else { word <- "is"; p.word <- "polymorph" }
      message(paste(p.word, paste(unique(stable.polymorph), collapse = ","), word, "stable"))
    } else {
      # Multiple polymorphs aren't involved ... things stay the same
      out.new[[i]] <- outprops[[are.polymorphs]]
      reaction.new[i, ] <- reaction[are.polymorphs, ]
      reaction.new$state[i] <- reaction$state[are.polymorphs]
      isaq.new <- c(isaq.new, isaq[are.polymorphs])
      iscgl.new <- c(iscgl.new, iscgl[are.polymorphs])
      isH2O.new <- c(isH2O.new, isH2O[are.polymorphs])
    }
  }

  outprops <- out.new
  # Remove the rows that were added to keep track of phase transitions
  reaction <- reaction.new[1:length(ispecies), ]
  # The manipulations above should get the correct species indices and state labels,
  # but if species are (intentionally) repeated, include only the first
  # (and possibly incorrect) reaction coefficients, so use the originals here 20180822
  reaction$coeff <- coeff
  isaq <- isaq.new
  iscgl <- iscgl.new
  isH2O <- isH2O.new

  # Adjust the output order of the properties
  for(i in 1:length(outprops)) {
    # the calculated properties are first
    ipp <- match(calcprop, colnames(outprops[[i]]))
    # move polymorph/loggam columns to end
    if("polymorph" %in% colnames(outprops[[i]])) ipp <- c(ipp, match("polymorph", colnames(outprops[[i]]))) 
    if("loggam" %in% colnames(outprops[[i]])) ipp <- c(ipp, match("loggam", colnames(outprops[[i]]))) 
    outprops[[i]] <- outprops[[i]][, ipp, drop = FALSE]
  }

  # Add up reaction properties
  if(do.reaction) {
    o <- 0
    morphcols <- NULL
    # do our affinity calculations here
    if(!is.null(logact)) {
      logQ <- logK <- rep(0, length(T))
      for(i in 1:length(coeff)) {
        logK <- logK + outprops[[i]]$logK * coeff[i]
        logQ <- logQ + logact[i] * coeff[i]
      }
      reaction <- cbind(reaction, logact)
      A <- logK - logQ
      # convert A/2.303RT (dimensionless) to J mol-1
      # then outvert to the user's units from J mol-1
      A <- outvert(convert(-A, "G", T = T), "J")
    }
    # Loop over reaction coefficients
    for(i in 1:length(coeff)) {
      # Assemble polymorph columns separately
      if("polymorph" %in% colnames(outprops[[i]])) {
         sc <- as.data.frame(outprops[[i]]$polymorph)
         outprops[[i]] <- outprops[[i]][, -match("polymorph", colnames(outprops[[i]]))]
         colnames(sc) <- as.character(reaction$name[i])
         if(is.null(morphcols)) morphcols <- sc
         else morphcols <- cbind(morphcols, sc)
      }
      # Include a zero loggam column if needed (for those species that are ideal)
      o.i <- outprops[[i]]
      if("loggam" %in% colnames(o.i)) if(!"loggam" %in% colnames(o))
        o <- cbind(o, loggam = 0)
      if("loggam" %in% colnames(o)) if(!"loggam" %in% colnames(o.i))
        o.i <- cbind(o.i, loggam = 0)
      # the real addition of properties
      o <- o + o.i * coeff[i]
    }
    # Output for reaction (stack on polymorph columns if exist)
    if(!is.null(morphcols)) OUT <- list(reaction = reaction,out = o,polymorphs = morphcols)
    else OUT <- list(reaction = reaction,out = o)
  } else {
    # Output for species: strip the coeff column from reaction
    reaction <- reaction[,-match("coeff",colnames(reaction))]
    OUT <- c(list(species = reaction),outprops)
  }
  # Append T, P, rho, A, logQ columns and convert units
  for(i in 2:length(OUT)) {
    # affinity and logQ
    if(i==2) if(do.reaction & !is.null(logact)) {
      OUT[[i]] <- cbind(OUT[[i]], data.frame(logQ = logQ, A = A))
    }
    # 20120114 Only prepend T, P, rho columns if we have more than one T
    # 20171020 or if the "property" argument is missing (it's nice to see everything using e.g. subcrt("H2O", T = 150))
    # 20171021 or if the "property" argument is not missing, but is identical to the default (happens when auto-balancing reactions)
    if(length(T) > 1 | missing(property) | identical(property, c("logK", "G", "H", "S", "V", "Cp"))) {
      # 20090329 Added checks for converting T, P units
      if(convert) T.out <- outvert(T, "K") else T.out <- T
      if(convert) P.out <- outvert(P, "bar") else P.out <- P
      # Try to stuff in a column of rho if we have aqueous species
      # watch out! supcrt-ish densities are in g/cc not kg/m3
      if("rho" %in% calcprop | ( (missing(property) | identical(property, c("logK", "G", "H", "S", "V", "Cp"))) &
                                any(c(isaq, isH2O))) & (names(OUT)[i]) != "polymorph") 
        OUT[[i]] <- cbind(data.frame(T = T.out, P = P.out, rho = H2O.PT$rho/1000), OUT[[i]])
      else
        OUT[[i]] <- cbind(data.frame(T = T.out, P = P.out, OUT[[i]]))
    }
  }
  # Put ionic strength next to any loggam columns
  for(i in 2:length(OUT)) {
    if("loggam" %in% colnames(OUT[[i]])) OUT[[i]] <- cbind(OUT[[i]], IS = newIS)
  }
  # More fanagling for species
  if(!do.reaction) {
    OUT <- list(species = OUT$species, out = OUT[2:length(OUT)])
    # add names to the output
    names(OUT$out) <- as.character(reaction$name)
  }
  # Rewritten code to convert energy units 20220325
  if(convert) {
    if(do.reaction) {
      isenergy <- colnames(OUT$out) %in% c("G", "H", "S", "Cp")
      if(any(isenergy)) OUT$out[, isenergy] <- outvert(OUT$out[, isenergy], "J")
    } else {
      isenergy <- colnames(OUT$out[[1]]) %in% c("G", "H", "S", "Cp")
      if(any(isenergy)) {
        for(i in 1:length(OUT$out)) OUT$out[[i]][, isenergy] <- outvert(OUT$out[[i]][, isenergy], "J")
      }
    }
  }
  # Add warnings to output 20180922
  if(length(warnings) > 0) {
    OUT <- c(OUT, list(warnings = warnings))
    for(warn in warnings) warning(warn)
  }
  return(OUT)
}
