# CHNOSZ/subcrt.R
# calculate standard molal thermodynamic propertes
# 20060817 jmd

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.args.R")
#source("util.character.R")
#source("info.R")
#source("util.units.R")
#source("util.data.R")
#source("species.R")
#source("AkDi.R")

subcrt <- function(species, coeff = 1, state = NULL, property = c("logK", "G", "H", "S", "V", "Cp"),
  T = seq(273.15, 623.15, 25), P = "Psat", grid = NULL, convert = TRUE, exceed.Ttr = FALSE,
  exceed.rhomin = FALSE, logact = NULL, action.unbalanced = "warn", IS = 0) {

  # revise the call if the states have 
  # come as the second argument 
  if(!is.null(coeff[1])) {
    if(is.numeric(state[1])) newcoeff <- state else newcoeff <- 1
    if(is.character(coeff[1])) newstate <- coeff else newstate <- NULL
    if(is.character(coeff[1])) {
      if(missing(T)) {
        if(identical(newcoeff,1) & !(identical(newcoeff,state))) 
          return(subcrt(species,state=coeff,property=property,P=P,grid=grid,
            convert=convert,exceed.Ttr=exceed.Ttr,logact=logact))
          else return(subcrt(species,coeff=newcoeff,state=coeff,property=property,
            P=P,grid=grid,convert=convert,exceed.Ttr=exceed.Ttr,logact=logact))
      } else {
        if(identical(newcoeff,1) & !(identical(newcoeff,state))) 
          return(subcrt(species,state=coeff,property=property,T=T,P=P,grid=grid,
            convert=convert,exceed.Ttr=exceed.Ttr,logact=logact))
          else return(subcrt(species,coeff=newcoeff,state=coeff,property=property,
            T=T,P=P,grid=grid,convert=convert,exceed.Ttr=exceed.Ttr,logact=logact))
      }
    }
  }

  do.reaction <- FALSE
  if(!missing(coeff)) do.reaction <- TRUE

  # species and states are made the same length
  if(!is.null(state[1])) {
    if(length(state) > length(species)) species <- rep(species,length.out=length(state))
    if(length(species) > length(state)) state <- rep(state,length.out=length(species))
    state <- state.args(state)
  }

  # allowed properties
  properties <- c("rho", "logK", "G", "H", "S", "Cp", "V", "kT", "E")
  # property checking
  calcprop <- property
  notprop <- property[!calcprop %in% properties]
  if(length(notprop) == 1) stop(paste("invalid property name:", paste(notprop, collapse=" ")))
  if(length(notprop) > 1) stop(paste("invalid property names:", paste(notprop, collapse=" ")))
  # length checking
  if(do.reaction & length(species)!=length(coeff)) 
    stop("the length of 'coeff' must equal the number of species")
  if(!is.null(logact)) {
    if(length(logact)!=length(species)) stop("the length of 'logact' must equal the number of species")
  }
  # normalize temperature and pressure units
  if(!missing(T)) {
    if(convert) T <- envert(T,'K')
    else if(!missing(convert) & convert) T <- envert(T,'K')
  }
  if(is.numeric(P[1])) {
    if(convert) P <- envert(P,'bar')
  }

  # warn for too high temperatures for Psat 20171110
  warnings <- character()
  if(identical(P, "Psat") & any(T > 647.067)) {
    nover <- sum(T > 647.067)
    if(nover==1) vtext <- "value" else vtext <- "values"
    warnings <- c(warnings, paste0("P = 'Psat' undefined for T > Tcritical (", nover, " T ", vtext, ")"))
  }

  # gridding?
  isPsat <- FALSE
  do.grid <- FALSE
  if(!is.null(grid)) if(!is.logical(grid)) do.grid <- TRUE
  newIS <- IS
  if(do.grid) {
    if(grid=='T') {
      newT <- numeric()
      for(i in 1:length(T)) newT <- c(newT,rep(T[i],length(P)))
      newP <- rep(P,length(T))
      T <- newT; P <- newP
    }
    if(grid=='P') {
      newP <- numeric()
      for(i in 1:length(P)) newP <- c(newP,rep(P[i],length(T)))
      newT <- rep(T,length(P))
      T <- newT; P <- newP
    }
    if(grid=='IS') {
      ll <- length(T)
      if(length(P) > 1) ll <- length(P)
      newIS <- numeric()
      for(i in 1:length(IS)) newIS <- c(newIS,rep(IS[i],ll))
      tpargs <- TP.args(T=T,P=P)
      T <- rep(tpargs$T,length.out=length(newIS))
      P <- rep(tpargs$P,length.out=length(newIS))
    }
  } else {
    # for AkDi, remember if P = "Psat" 20190219
    if(identical(P, "Psat")) isPsat <- TRUE
    # expansion of Psat and equivalence of argument lengths
    tpargs <- TP.args(T=T,P=P)
    T <- tpargs$T; P <- tpargs$P
    if(length(newIS) > length(T)) T <- rep(T, length.out=length(newIS))
    if(length(newIS) > length(P)) P <- rep(P, length.out=length(newIS))
  }

  # get species information
  thermo <- get("thermo", CHNOSZ)
  # pre-20110808, we sent numeric species argument through info() to
  # get species name and state(s)
  # but why slow things down if we already have a species index?
  # so now phase species stuff will only work for character species names
  if(is.numeric(species[1])) {
    ispecies <- species
    species <- as.character(thermo$obigt$name[ispecies])
    state <- as.character(thermo$obigt$state[ispecies])
    newstate <- as.character(thermo$obigt$state[ispecies])
  } else {
    # from names, get species indices and states and possibly
    # keep track of phase species (cr,cr2 ...)
    ispecies <- numeric()
    newstate <- character()
    for(i in 1:length(species)) {
      # get the species index for a named species
      if(!can.be.numeric(species[i])) si <- info.character(species[i], state[i])
      else {
        # check that a numeric argument is a rownumber of thermo$obigt
        si <- as.numeric(species[i])
        if(!si %in% 1:nrow(thermo$obigt)) stop(paste(species[i], "is not a row number of thermo$obigt"))
      }
      # that could have the side-effect of adding a protein; re-read thermo
      thermo <- get("thermo", CHNOSZ)
      if(is.na(si[1])) stop('no info found for ',species[i],' ',state[i])
      if(!is.null(state[i])) is.cr <- state[i]=='cr' else is.cr <- FALSE
      if(thermo$obigt$state[si[1]]=='cr' & (is.null(state[i]) | is.cr)) {
        newstate <- c(newstate,'cr')
        ispecies <- c(ispecies,si[1])
      } else {
        newstate <- c(newstate,as.character(thermo$obigt$state[si[1]]))
        ispecies <- c(ispecies,si[1])
      }
    }
  }

  # stop if species not found
  noname <- is.na(ispecies)
  if(TRUE %in% noname)
    stop(paste('species',species[noname],'not found.\n'))

  # take care of mineral phases
  state <- as.character(thermo$obigt$state[ispecies])
  name <- as.character(thermo$obigt$name[ispecies])
  # a counter of all species considered
  # iphases is longer than ispecies if cr,cr2 ... phases are present
  # phasespecies shows which of ispecies correspond to iphases
  # pre-20091114: the success of this depends on there not being duplicated aqueous or other
  # non-mineral-phase species (i.e., two entries in obigt for Cu+ screw this up
  # when running the skarn example).
  # after 20091114: we can deal with duplicated species (aqueous at least)
  iphases <- phasespecies <- coeff.new <- numeric()
  for(i in 1:length(ispecies)) {
     if(newstate[i]=='cr') {
       searchstates <- c('cr','cr2','cr3','cr4','cr5','cr6','cr7','cr8','cr9') 
       tghs <- thermo$obigt[(thermo$obigt$name %in% name[i]) & thermo$obigt$state %in% searchstates,]
       # we only take one if they are in fact duplicated species and not phase species
       if(all(tghs$state==tghs$state[1])) tghs <- thermo$obigt[ispecies[i],]
     } else tghs <- thermo$obigt[ispecies[i],]
     iphases <- c(iphases,as.numeric(rownames(tghs))) 
     phasespecies <- c(phasespecies,rep(ispecies[i],nrow(tghs)))
     coeff.new <- c(coeff.new,rep(coeff[i],nrow(tghs)))
  }

  # where we keep info about the species involved
  reaction <- data.frame(coeff = coeff.new, name = thermo$obigt$name[iphases],
    formula = thermo$obigt$formula[iphases], state = thermo$obigt$state[iphases],
    ispecies = iphases, stringsAsFactors = FALSE)
  # make the rownames readable ... but they have to be unique
  if(length(unique(iphases))==length(iphases)) rownames(reaction) <- as.character(iphases)

  # wetness etc.
  isH2O <- reaction$name=='water' & reaction$state=='liq'
  isaq <- reaction$state=='aq'

  # produce message about conditions
  if(length(species)==1 & convert==FALSE) {
    # no message produced here (internal calls from other functions)
  } else {
    # include units here 20190530
    uT <- outvert(T, "K")
    if(identical(grid,'IS')) uT <- unique(uT)
    if(length(uT)==1) T.text <- paste(uT, T.units()) else {
      T.text <- paste0(length(uT), " values of T (", T.units(), ")")
    }
    uP <- outvert(P, "bar")
    if(length(P)==1) {
      if(can.be.numeric(P)) P.text <- paste(round(as.numeric(uP),2), P.units())
      else P.text <- paste0("P (", P.units(), ")")
    } else P.text <- paste0("P (", P.units(), ")")
    if(identical(P[[1]],'Psat')) P.text <- P
    if(any(c(isH2O,isaq))) P.text <- paste(P.text,' (wet)',sep='')
    E.text <- paste0("[energy units: ", E.units(), "]")
    message(paste("subcrt:", length(species), "species at", T.text, "and", P.text, E.text))
  }

  # inform about unbalanced reaction
  if(do.reaction) {
    # the mass balance ... is zero for a balanced reaction
    mss <- makeup(ispecies, coeff, sum=TRUE)
    # take out very small numbers
    mss[abs(mss) < 1e-7] <- 0
    # report and try to fix any non-zero mass balance
    if(any(mss!=0) & !is.null(action.unbalanced)) {
      # the missing composition: the negative of the mass balance
      miss <- -mss
      # drop elements that are zero
      miss <- miss[miss!=0]
      message("subcrt: reaction is not balanced; it is missing this composition:")
      # we have to do this awkward dance to send a formatted matrix to message
      message(paste(capture.output(print(miss)), collapse="\n"))
      # look for basis species that have our compositoin
      tb <- thermo$basis
      if(!is.null(tb)) {
        if(all(names(miss) %in% colnames(tb)[1:nrow(tb)])) {
          # the missing composition as formula
          ft <- as.chemical.formula(miss)
          # the basis species needed to supply it
          bc <- species.basis(ft)
          # drop zeroes
          bc.new <- bc[,(bc[1,]!=0),drop=FALSE]
          # and get the states
          b.state <- as.character(thermo$basis$state)[bc[1,]!=0]
          bc <- bc.new
          # special thing for Psat
          if(identical(P[[1]], "Psat")) P <- "Psat"
          else P <- outvert(P,"bar")
          # add to logact values if present
          if(!is.null(logact)) {
            ila <- match(colnames(bc),rownames(thermo$basis))
            nla <- !(can.be.numeric(thermo$basis$logact[ila]))
            if(any(nla)) warning('subcrt: logact values of basis species',
              c2s(rownames(thermo$basis)[ila]),'are NA.')
            logact <- c(logact,thermo$basis$logact[ila])
          }
          # warn user and do it!
          ispecies.new <- tb$ispecies[match(colnames(bc),rownames(tb))]
          b.species <- thermo$obigt$formula[ispecies.new]
          if(identical(species,b.species) & identical(state,b.state))
            message("subcrt: balanced reaction, but it is a non-reaction; restarting...")
          else message('subcrt: adding missing composition from basis definition and restarting...')
          newspecies <- c(species, tb$ispecies[match(colnames(bc), rownames(tb))])
          newcoeff <- c(coeff, as.numeric(bc[1, ]))
          newstate <- c(state, b.state)
          return(subcrt(species=newspecies, coeff=newcoeff, state=newstate,
            property=property, T=outvert(T, "K"), P=P, grid=grid, convert=convert, logact=logact, exceed.Ttr=FALSE))
        } else if(identical(action.unbalanced,'warn')) 
            warnings <- c(warnings, paste('reaction was unbalanced, missing', as.chemical.formula(miss)))
      } else {
        if(identical(action.unbalanced,'warn')) 
          warnings <- c(warnings, paste('reaction was unbalanced, missing', as.chemical.formula(miss)))
      }
    }
  }

  # calculate the properties
  # if we want affinities we must have logK; include it in the ouput
  if(!is.null(logact)) if(!'logK' %in% calcprop) calcprop <- c('logK', calcprop)
  # if logK but not G was requested, we need to calculate G
  eosprop <- calcprop
  if('logK' %in% calcprop & ! 'G' %in% calcprop) eosprop <- c(eosprop, 'G')
  # also get G if we are dealing with mineral phases
  if(!'G' %in% eosprop & length(iphases) > length(ispecies)) eosprop <- c(eosprop, 'G')
  # don't request logK or rho from the eos ...
  eosprop <- eosprop[!eosprop %in% c('logK','rho')]
  # the reaction result will go here
  outprops <- list()
  # aqueous species and H2O properties
  if(TRUE %in% isaq) {
    # 20110808 get species parameters using obigt2eos() (faster than using info())
    param <- obigt2eos(thermo$obigt[iphases[isaq],], "aq", fixGHS = TRUE, tocal = TRUE)
    # aqueous species with NA for Z use the AkDi model
    isAkDi <- is.na(param$Z)
    # always get density
    H2O.props <- "rho"
    # calculate A_DH and B_DH if we're using the B-dot (Helgeson) equation
    if(any(IS != 0) & thermo$opt$nonideal %in% c("Bdot", "Bdot0", "bgamma", "bgamma0")) H2O.props <- c(H2O.props, "A_DH", "B_DH")
    # get other properties for H2O only if it's in the reaction
    if(any(isH2O)) H2O.props <- c(H2O.props, eosprop)
    # in case everything is AkDi, run hkf (for water properties) but exclude all species
    hkfpar <- param
    if(all(isAkDi)) hkfpar <- param[0, ]
    hkfstuff <- hkf(eosprop, parameters = hkfpar, T = T, P = P, H2O.props=H2O.props)
    p.aq <- hkfstuff$aq
    H2O.PT <- hkfstuff$H2O
    # set properties to NA for density below 0.35 g/cm3 (a little above the critical isochore, threshold used in SUPCRT92) 20180922
    if(!exceed.rhomin & !all(isAkDi)) {
      ilowrho <- H2O.PT$rho < 350
      ilowrho[is.na(ilowrho)] <- FALSE
      if(any(ilowrho)) {
        for(i in 1:length(p.aq)) p.aq[[i]][ilowrho, ] <- NA
        if(sum(ilowrho)==1) ptext <- "pair" else ptext <- "pairs"
        warnings <- c(warnings, paste0("below minimum density for applicability of revised HKF equations (", sum(ilowrho), " T,P ", ptext, ")"))
      }
    }
    # calculate properties using Akinfiev-Diamond model 20190219
    if(any(isAkDi)) {
      # get the parameters with the right names
      param <- obigt2eos(param[isAkDi, ], "aq")
      p.aq[isAkDi] <- AkDi(eosprop, parameters = param, T = T, P = P, isPsat = isPsat)
    }
    # calculate activity coefficients if ionic strength is not zero
    if(any(IS != 0)) {
      if(thermo$opt$nonideal %in% c("Bdot", "Bdot0", "bgamma", "bgamma0")) p.aq <- nonideal(iphases[isaq], p.aq, newIS, T, P, H2O.PT$A_DH, H2O.PT$B_DH)
      else if(thermo$opt$nonideal=="Alberty") p.aq <- nonideal(iphases[isaq], p.aq, newIS, T)
    }
    outprops <- c(outprops, p.aq)
  } else if(any(isH2O)) {
    # we're not using the HKF, but still want water
    H2O.PT <- water(c("rho", eosprop), T = T, P = P)
  }

  # crystalline, gas, liquid (except water) species
  cglstates <- c("liq", "cr", "gas", "cr2", "cr3", "cr4", "cr5", "cr6", "cr7", "cr8", "cr9")
  iscgl <- reaction$state %in% cglstates & reaction$name != "water"

  if(TRUE %in% iscgl) {
    param <- obigt2eos(thermo$obigt[iphases[iscgl],], "cgl", fixGHS = TRUE)
    p.cgl <- cgl(eosprop, parameters = param, T = T, P = P)
    # replace Gibbs energies with NA where the
    # phases are beyond their temperature range
    if('G' %in% eosprop) {
      # 20080304 this code is weird and hard to read - needs a lot of cleanup!
      # 20120219 cleaned up somewhat; using exceed.Ttr and NA instead of do.phases and 999999
      # the numbers of the cgl species (becomes 0 for any that aren't cgl)
      ncgl <- iscgl
      ncgl[iscgl] <- 1:nrow(param)
      for(i in 1:length(iscgl)) {
        # not if we're not cgl
        if(!iscgl[i]) next
        # name and state
        myname <- reaction$name[i]
        mystate <- reaction$state[i]
        # if we are considering multiple phases, and if this phase is cr2 or higher, check if we're below the transition temperature
        if(length(iphases) > length(ispecies) & i > 1) {
          if(!(reaction$state[i] %in% c('liq','cr','gas')) & reaction$name[i-1] == reaction$name[i]) {
            # after add.obigt("SUPCRT92"), quartz cr and cr2 are not next to each other in thermo$obigt,
            # so use iphases[i-1] here, not iphases[i]-1  20181107
            Ttr <- Ttr(iphases[i-1], iphases[i], P=P, dPdT = dPdTtr(iphases[i-1], iphases[i]))
            if(all(is.na(Ttr))) next
            if(any(T < Ttr)) {
              status.Ttr <- "(extrapolating G)"
              if(!exceed.Ttr) {
                # put NA into the value of G
                p.cgl[[ncgl[i]]]$G[T<Ttr] <- NA
                status.Ttr <- "(using NA for G)"
              } 
              #message(paste('subcrt: some points below transition temperature for',myname, mystate, status.Ttr))
            }
          }
        }
        # check if we're above the temperature limit or transition temperature
        # T limit (or Ttr) from the database
        warn.above <- TRUE
        Ttr <- thermo$obigt$z.T[iphases[i]]
        # calculate Ttr at higher P if a phase transition is present
        if(i < nrow(reaction)) {
          # if the next one is cr2, cr3, etc we have a transition
          if(reaction$state[i+1] %in% c("cr1", "cr2", "cr3", "cr4", "cr5", "cr6", "cr7", "cr8", "cr9") & reaction$name[i+1] == reaction$name[i]) {
            Ttr <- Ttr(iphases[i], iphases[i+1], P = P, dPdT = dPdTtr(iphases[i], iphases[i+1]))
            # we don't warn here about the transition
            warn.above <- FALSE
          }
        }
        if(any(is.na(Ttr))) next
        if(!all(Ttr == 0) & any(T >= Ttr)) {
          status.Ttr <- "(extrapolating G)"
          if(!exceed.Ttr) {
            p.cgl[[ncgl[i]]]$G[T >= Ttr] <- NA
            status.Ttr <- "(using NA for G)"
          }
          Tmax <- min(T[T >= Ttr])
          if(warn.above) message(paste("subcrt: temperature(s) of", Tmax, "K and above exceed limit for", myname, mystate, status.Ttr))
        }
      }
    }
    outprops <- c(outprops,p.cgl)
  }

  # water
  if(any(isH2O)) {
    p.H2O <- H2O.PT[, match(eosprop, colnames(H2O.PT)), drop=FALSE]
    p.H2O <- list(p.H2O)
    outprops <- c(outprops, rep(p.H2O, sum(isH2O == TRUE)))
  }

  # use variable-pressure standard Gibbs energy for gases
  isgas <- reaction$state %in% "gas" 
  if(any(isgas) & "G" %in% eosprop & thermo$opt$varP) {
    for(i in which(isgas)) outprops[[i]]$G <- outprops[[i]]$G - convert(log10(P), "G", T=T)
  }

  # logK
  if('logK' %in% calcprop) {
    for(i in 1:length(outprops)) {
      outprops[[i]] <- cbind(outprops[[i]],data.frame(logK=convert(outprops[[i]]$G,'logK',T=T)))
      colnames(outprops[[i]][ncol(outprops[[i]])]) <- 'logK'
    }
  }

  # ordering the output
  # the indices of the species in outprops thus far
  ns <- 1:nrow(reaction)
  is <- c(ns[isaq],ns[iscgl],ns[isH2O])
  if(length(ns)!=length(is)) stop('subcrt: not all species are accounted for.')
  v <- list()
  for(i in 1:length(is))  v[[i]] <- outprops[[match(ns[i],is)]]
  outprops <- v

  # deal with phases (cr,cr2) here
  # we have to eliminate rows from outprops, 
  # reaction and values from isaq, iscgl, isH2O
  out.new <- list()
  reaction.new <- reaction
  isaq.new <- logical()
  iscgl.new <- logical()
  isH2O.new <- logical()
  for(i in 1:length(ispecies)) {
    arephases <- which(ispecies[i]==phasespecies)
    # deal with repeated species here
    if(TRUE %in% duplicated(iphases[arephases])) {
      # only take the first, not the duplicates
      ndups <- length(which(ispecies==ispecies[i]))
      nphases <- length(arephases) / ndups
      arephases <- arephases[1:nphases]
    }
    if(length(arephases)>1) {
      message(paste('subcrt:',length(arephases),'phases for',thermo$obigt$name[ispecies[i]],'... '), appendLF=FALSE)
      # assemble the Gibbs energies for each species
      for(j in 1:length(arephases)) {
        G.this <- outprops[[arephases[j]]]$G
        if(length(which(is.na(G.this))) > 0 & exceed.Ttr) warning(paste('subcrt: NAs found for G of ',
          reaction$name[arephases[j]],' ',reaction$state[arephases[j]],' at T-P point(s) ',
          c2s(which(is.na(G.this)),sep=' '),sep=''),call.=FALSE)
        if(j==1) G <- as.data.frame(G.this)
        else G <- cbind(G,as.data.frame(G.this))
      }
      # find the minimum-energy phase at each T-P point
      phasestate <- numeric()
      out.new.entry <- outprops[[arephases[1]]]
      for(j in 1:nrow(G)) {
        ps <- which.min(as.numeric(G[j,]))
        if(length(ps)==0) {
          # minimum not found (we have NAs)
          # - no non-NA value of G to begin with, e.g. aegerine) --> probably should use lowest-T phase
          #ps <- 1
          # - above temperature limit for the highest-T phase (subcrt.Rd skarn example) --> use highest-T phase 20171110
          ps <- ncol(G)
          if(exceed.Ttr) warning('subcrt: stable phase for ',reaction$name[arephases[ps]],' at T-P point ',j,
          ' undetermined (using ',reaction$state[arephases[ps]],')',call.=FALSE)
        } 
        phasestate <- c(phasestate,ps)
        out.new.entry[j,] <- outprops[[ arephases[ps] ]][j,]
      }

      # update our objects
      out.new[[i]] <- cbind(out.new.entry,data.frame(polymorph=phasestate))
      reaction.new[i,] <- reaction[arephases[phasestate[1]],]
      # mark the minerals with multiple phases
      reaction.new$state[i] <- "cr*"
      isaq.new <- c(isaq.new,isaq[arephases[phasestate[1]]])
      iscgl.new <- c(iscgl.new,iscgl[arephases[phasestate[1]]])
      isH2O.new <- c(isH2O.new,isH2O[arephases[phasestate[1]]])
      # info for the user
      up <- unique(phasestate)
      if(length(up)>1) { word <- 'are'; p.word <- 'phases' }
      else { word <- 'is'; p.word <- 'phase' }
      message(paste(p.word,paste(unique(phasestate), collapse=","),word,'stable'))
    } else {
      # multiple phases aren't involved ... things stay the same
      out.new[[i]] <- outprops[[arephases]]
      reaction.new[i, ] <- reaction[arephases, ]
      reaction.new$state[i] <- reaction$state[arephases]
      isaq.new <- c(isaq.new,isaq[arephases])
      iscgl.new <- c(iscgl.new,iscgl[arephases])
      isH2O.new <- c(isH2O.new,isH2O[arephases])
    }
  }

  outprops <- out.new
  # remove the rows that were added to keep track of phase transitions
  reaction <- reaction.new[1:length(ispecies),]
  # the manipulations above should get the correct species indices and state labels,
  # but if species are (intentionally) repeated, include only the first
  # (and possibly incorrect) reaction coefficients, so use the originals here 20180822
  reaction$coeff <- coeff
  isaq <- isaq.new
  iscgl <- iscgl.new
  isH2O <- isH2O.new

  # adjust the output order of the properties
  for(i in 1:length(outprops)) {
    # the calculated properties are first
    ipp <- match(calcprop, colnames(outprops[[i]]))
    # move polymorph/loggam columns to end
    if('polymorph' %in% colnames(outprops[[i]])) ipp <- c(ipp,match('polymorph',colnames(outprops[[i]]))) 
    if('loggam' %in% colnames(outprops[[i]])) ipp <- c(ipp,match('loggam',colnames(outprops[[i]]))) 
    outprops[[i]] <- outprops[[i]][,ipp,drop=FALSE]
  }

  # add up reaction properties
  if(do.reaction) {
    o <- 0
    morphcols <- NULL
    # do our affinity calculations here
    if(!is.null(logact)) {
      logQ <- logK <- rep(0,length(T))
      for(i in 1:length(coeff)) {
        logK <- logK + outprops[[i]]$logK * coeff[i]
        logQ <- logQ + logact[i] * coeff[i]
      }
      reaction <- cbind(reaction,logact)
      A <- logK - logQ
      # convert A/2.303RT (no dims) to cal mol-1
      # then to the user's units (outvert) from cal
      A <- outvert(convert(-A,'G',T=T),'cal')
    }
    # loop over reaction coefficients
    for(i in 1:length(coeff)) {
      # assemble polymorph columns separately
      if('polymorph' %in% colnames(outprops[[i]])) {
         sc <- as.data.frame(outprops[[i]]$polymorph)
         outprops[[i]] <- outprops[[i]][,-match('polymorph',colnames(outprops[[i]]))]
         colnames(sc) <- as.character(reaction$name[i])
         if(is.null(morphcols)) morphcols <- sc
         else morphcols <- cbind(morphcols,sc)
      }
      # include a zero loggam column if needed (for those species that are ideal)
      o.i <- outprops[[i]]
      if('loggam' %in% colnames(o.i)) if(!'loggam' %in% colnames(o))
        o <- cbind(o,loggam=0)
      if('loggam' %in% colnames(o)) if(!'loggam' %in% colnames(o.i))
        o.i <- cbind(o.i,loggam=0)
      # the real addition of properties
      o <- o + o.i * coeff[i]
    }
    # output for reaction (stack on polymorph columns if exist)
    if(!is.null(morphcols)) OUT <- list(reaction=reaction,out=o,polymorphs=morphcols)
    else OUT <- list(reaction=reaction,out=o)
  } else {
    # output for species: strip the coeff column from reaction
    reaction <- reaction[,-match('coeff',colnames(reaction))]
    OUT <- c(list(species=reaction),outprops)
  }
  # append T,P,rho, A, logQ columns and convert units
  for(i in 2:length(OUT)) {
    # affinity and logQ
    if(i==2) if(do.reaction & !is.null(logact)) {
      OUT[[i]] <- cbind(OUT[[i]],data.frame(logQ=logQ,A=A))
    }
    # 20120114 only prepend T, P, rho columns if we have more than one T
    # 20171020 or if the 'property' argument is missing (it's nice to see everything using e.g. subcrt("H2O", T=150))
    # 20171021 or if the 'property' argument is not missing, but is identical to the default (happens when auto-balancing reactions)
    if(length(T) > 1 | missing(property) | identical(property, c("logK", "G", "H", "S", "V", "Cp"))) {
      # 20090329 added checks for converting T, P units
      if(convert) T.out <- outvert(T,"K") else T.out <- T
      if(convert) P.out <- outvert(P,"bar") else P.out <- P
      # try to stuff in a column of rho if we have aqueous species
      # watch out! supcrt-ish densities are in g/cc not kg/m3
      if('rho' %in% calcprop | ( (missing(property) | identical(property, c("logK", "G", "H", "S", "V", "Cp"))) &
                                any(c(isaq,isH2O))) & (names(OUT)[i])!='polymorph') 
        OUT[[i]] <- cbind(data.frame(T=T.out,P=P.out,rho=H2O.PT$rho/1000),OUT[[i]])
      else
        OUT[[i]] <- cbind(data.frame(T=T.out,P=P.out,OUT[[i]]))
    }
    if(convert) {
      for(j in 1:ncol(OUT[[i]])) {
        if(colnames(OUT[[i]])[j] %in% c('G','H','S','Cp')) OUT[[i]][,j] <- outvert(OUT[[i]][,j],'cal')
      }
    }
  }
  # put ionic strength next to any loggam columns
  for(i in 2:length(OUT)) {
    if('loggam' %in% colnames(OUT[[i]])) OUT[[i]] <- cbind(OUT[[i]],IS=newIS)
  }
  # more fanagling for species
  if(!do.reaction) {
    OUT <- list(species=OUT$species, out=OUT[2:length(OUT)])
    # add names to the output
    names(OUT$out) <- as.character(reaction$name)
  }
  # add warnings to output 20180922
  if(length(warnings) > 0) {
    OUT <- c(OUT, list(warnings=warnings))
    for(warn in warnings) warning(warn)
  }
  return(OUT)
}

