# CHNOSZ/cgl.R
# Calculate standard thermodynamic properties of non-aqueous species
# 20060729 jmd

cgl <- function(property = NULL, parameters = NULL, T = 298.15, P = 1) {
  # Calculate properties of crystalline, liquid (except H2O) and gas species
  Tr <- 298.15
  Pr <- 1
  # The number of T, P conditions
  ncond <- max(c(length(T), length(P)))
  # Initialize output list
  out <- list()
  # Loop over each species
  for(k in 1:nrow(parameters)) {
    # The parameters for *this* species
    PAR <- parameters[k, ]
    if(PAR$model == "Berman") {
      # Use Berman equations (parameters not in thermo()$OBIGT)
      properties <- Berman(PAR$name, T = T, P = P)
      iprop <- match(property, colnames(properties))
      values <- properties[, iprop, drop = FALSE]
    } else {
      # In CHNOSZ, we have
      # 1 cm^3 bar --> convert(1, "calories") == 0.02390057 cal
      # but REAC92D.F in SUPCRT92 uses
      cm3bar_to_cal <- 0.023901488 # cal
      cm3bar_to_J <- convert(cm3bar_to_cal, "J")
      # Start with NA values
      values <- data.frame(matrix(NA, ncol = length(property), nrow=ncond))
      colnames(values) <- property
      # A test for availability of heat capacity coefficients (a, b, c, d, e, f)
      # based on the column assignments in thermo()$OBIGT
      if(any(!is.na(PAR[, 15:20]))) {
        # We have at least one of the heat capacity coefficients;
        # zero out any NA's in the rest (leave lambda and T of transition (columns 19-20) alone)
        PAR[, 15:20][, is.na(PAR[, 15:20])] <- 0
        # calculate the heat capacity and its integrals
        Cp <- PAR$a + PAR$b*T + PAR$c*T^-2 + PAR$d*T^-0.5 + PAR$e*T^2 + PAR$f*T^PAR$lambda
        intCpdT <- PAR$a*(T - Tr) + PAR$b*(T^2 - Tr^2)/2 + PAR$c*(1/T - 1/Tr)/-1 + PAR$d*(T^0.5 - Tr^0.5)/0.5 + PAR$e*(T^3-Tr^3)/3
        intCpdlnT <- PAR$a*log(T / Tr) + PAR$b*(T - Tr) + PAR$c*(T^-2 - Tr^-2)/-2 + PAR$d*(T^-0.5 - Tr^-0.5)/-0.5  + PAR$e*(T^2 - Tr^2)/2
        # Do we also have the lambda parameter (Cp term with adjustable exponent on T)?
        if(!is.na(PAR$lambda) & !identical(PAR$lambda, 0)) {
           # Equations for lambda adapted from Helgeson et al., 1998 (doi:10.1016/S0016-7037(97)00219-6)
           if(PAR$lambda == -1) intCpdT <- intCpdT + PAR$f*log(T/Tr) 
           else intCpdT <- intCpdT - PAR$f*( T^(PAR$lambda + 1) - Tr^(PAR$lambda + 1) ) / (PAR$lambda + 1)
           intCpdlnT <- intCpdlnT + PAR$f*(T^PAR$lambda - Tr^PAR$lambda) / PAR$lambda
        }
      } else {
        # Use constant heat capacity if the coefficients are not available
        Cp <- PAR$Cp
        intCpdT <- PAR$Cp*(T - Tr)
        intCpdlnT <- PAR$Cp*log(T / Tr)
        # In case Cp is listed as NA, set the integrals to 0 at Tr
        intCpdT[T == Tr] <- 0
        intCpdlnT[T == Tr] <- 0
      }
      # Volume and its integrals
      if(PAR$name %in% c("quartz", "coesite")) {
        # Volume calculations for quartz and coesite
        qtz <- quartz_coesite(PAR, T, P)
        V <- qtz$V
        intVdP <- qtz$intVdP
        intdVdTdP <- qtz$intdVdTdP
      } else {
        # For other minerals, volume is constant (Helgeson et al., 1978)
        V <- rep(PAR$V, ncond)
        # If the volume is NA, set its integrals to zero
        if(is.na(PAR$V)) intVdP <- intdVdTdP <- numeric(ncond)
        else {
          intVdP <- PAR$V*(P - Pr) * cm3bar_to_J
          intdVdTdP <- 0
        }
      }
      # Get the values of each of the requested thermodynamic properties
      for(i in 1:length(property)) {
        if(property[i] == "Cp") values[, i] <- Cp
        if(property[i] == "V") values[, i] <- V
        if(property[i] == "E") values[, i] <- rep(NA, ncond)
        if(property[i] == "kT") values[, i] <- rep(NA, ncond)
        if(property[i] == "G") {
          # Calculate S * (T - Tr), but set it to 0 at Tr (in case S is NA)
          Sterm <- PAR$S*(T - Tr)
          Sterm[T==Tr] <- 0
          values[, i] <- PAR$G - Sterm + intCpdT - T*intCpdlnT + intVdP
        }
        if(property[i] == "H") values[, i] <- PAR$H + intCpdT + intVdP - T*intdVdTdP
        if(property[i] == "S") values[, i] <- PAR$S + intCpdlnT - intdVdTdP
      }
    } # End calculations using parameters from thermo()$OBIGT
    out[[k]] <- values
  } # End loop over species

  return(out)
}

### Unexported function ###

# Calculate GHS and V corrections for quartz and coesite 20170929
# (these are the only mineral phases for which SUPCRT92 uses an inconstant volume)
quartz_coesite <- function(PAR, T, P) {
  # The corrections are 0 for anything other than quartz and coesite
  if(!PAR$name %in% c("quartz", "coesite")) return(list(G = 0, H = 0, S = 0, V = 0))
  ncond <- max(c(length(T), length(P)))
  # Tr, Pr and TtPr (transition temperature at Pr)
  Pr <- 1      # bar
  Tr <- 298.15 # K
  TtPr <- 848  # K
  # Constants from SUP92D.f
  aa <- 549.824
  ba <- 0.65995
  ca <- -0.4973e-4
  VPtTta <- 23.348
  VPrTtb <- 23.72
  Stran <- 0.342
  # Constants from REAC92D.f
  VPrTra <- 22.688 # VPrTr(a-quartz)
  Vdiff <- 2.047   # VPrTr(a-quartz) - VPrTr(coesite)
  k <- 38.5       # dPdTtr(a/b-quartz)
  #k <- 38.45834    # calculated in CHNOSZ: dPdTtr(info("quartz"))
  # Code adapted from REAC92D.f
  qphase <- gsub("cr", "", PAR$state)
  if(qphase == 2) {
    Pstar <- P
    Sstar <- rep(0, ncond)
    V <- rep(VPrTtb, ncond)
  } else {
    Pstar <- Pr + k * (T - TtPr)
    Sstar <- rep(Stran, ncond)
    V <- VPrTra + ca*(P-Pr) + (VPtTta - VPrTra - ca*(P-Pr))*(T-Tr) / (TtPr + (P-Pr)/k - Tr)
  }
  Pstar[T < TtPr] <- Pr
  Sstar[T < TtPr] <- 0
  if(PAR$name == "coesite") {
    VPrTra <- VPrTra - Vdiff
    VPrTtb <- VPrTtb - Vdiff
    V <- V - Vdiff
  }
  cm3bar_to_cal <- 0.023901488
  GVterm <- cm3bar_to_cal * (VPrTra * (P - Pstar) + VPrTtb * (Pstar - Pr) -
    0.5 * ca * (2 * Pr * (P - Pstar) - (P^2 - Pstar^2)) -
    ca * k * (T - Tr) * (P - Pstar) +
    k * (ba + aa * ca * k) * (T - Tr) * log((aa + P/k) / (aa + Pstar/k)))
  SVterm <- cm3bar_to_cal * (-k * (ba + aa * ca * k) *
    log((aa + P/k) / (aa + Pstar/k)) + ca * k * (P - Pstar)) - Sstar
  # Note the minus sign on "SVterm" in order that intdVdTdP has the correct sign
  list(intVdP = convert(GVterm, "J"), intdVdTdP = convert(-SVterm, "J"), V = V)
}
