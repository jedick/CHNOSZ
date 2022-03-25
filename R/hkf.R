# CHNOSZ/hkf.R
# Calculate thermodynamic properties using equations of state
# 11/17/03 jmd

## If this file is interactively sourced, the following is also needed to provide unexported functions:
#source("util.args.R")

hkf <- function(property = NULL, parameters = NULL, T = 298.15, P = 1,
  contrib = c("n", "s", "o"), H2O.props="rho") {
  # Calculate G, H, S, Cp, V, kT, and/or E using the revised HKF equations of state
  # H2O.props - H2O properties needed for subcrt() output
  # Constants
  Tr <- 298.15 # K
  Pr <- 1      # bar
  Theta <- 228 # K
  Psi <- 2600  # bar
  # Make T and P equal length
  if(!identical(P, "Psat")) {
    if(length(P) < length(T)) P <- rep(P, length.out = length(T))
    if(length(T) < length(P)) T <- rep(T, length.out = length(P))
  }
  # Nonsolvation, solvation, and origination contribution
  notcontrib <- ! contrib %in% c("n", "s", "o")
  if(TRUE %in% notcontrib) stop(paste("contrib must be in c('n', 's', 'o); got", c2s(contrib[notcontrib])))
  # Get water properties
  # rho - for subcrt() output and g function
  # Born functions and epsilon - for HKF calculations
  H2O.props <- c(H2O.props, "QBorn", "XBorn", "YBorn", "epsilon")
  thermo <- get("thermo", CHNOSZ)
  if(grepl("SUPCRT", thermo$opt$water)) {
    # Using H2O92D.f from SUPCRT92: alpha, daldT, beta - for partial derivatives of omega (g function)
    H2O.props <- c(H2O.props, "alpha", "daldT", "beta")
  }
  if(grepl("IAPWS", thermo$opt$water)) {
    # Using IAPWS-95: NBorn, UBorn - for compressibility, expansibility
    H2O.props <- c(H2O.props, "NBorn", "UBorn")
  }
  if(grepl("DEW", thermo$opt$water)) {
    # Using DEW model: get beta to calculate dgdP
    H2O.props <- c(H2O.props, "beta")
  }
  H2O <- water(H2O.props, T = c(Tr, T), P = c(Pr, P))
  H2O.PrTr <- H2O[1, ]
  H2O.PT <- H2O[-1, ]
  ZBorn <- -1 / H2O.PT$epsilon
  ZBorn.PrTr <- -1 / H2O.PrTr$epsilon
  # A list to store the result
  aq.out <- list()
  nspecies <- nrow(parameters)
  for(k in seq_len(nspecies)) {
    # Loop over each species
    PAR <- parameters[k, ]
    # Substitute Cp and V for missing EoS parameters
    # Here we assume that the parameters are in the same position as in thermo()$OBIGT
    # We don't need this if we're just looking at solvation properties (Cp_s_var, V_s_var)
    if("n" %in% contrib) {
      # Put the heat capacity in for c1 if both c1 and c2 are missing
      if(all(is.na(PAR[, 18:19]))) PAR[, 18] <- PAR$Cp
      # Put the volume in for a1 if a1, a2, a3 and a4 are missing
      if(all(is.na(PAR[, 14:17]))) PAR[, 14] <- convert(PAR$V, "joules")
      # Test for availability of the EoS parameters
      hasEOS <- any(!is.na(PAR[, 14:21]))
      # If at least one of the EoS parameters is available, zero out any NA's in the rest
      if(hasEOS) PAR[, 14:21][, is.na(PAR[, 14:21])] <- 0
    }
    # Compute values of omega(P,T) from those of omega(Pr,Tr)
    # Using g function etc. (Shock et al., 1992 and others)
    omega <- PAR$omega  # omega.PrTr
    # The derivatives are zero unless the g function kicks in
    dwdP <- dwdT <- d2wdT2 <- numeric(length(T))
    Z <- PAR$Z
    omega.PT <- rep(PAR$omega, length(T))
    if(!identical(Z, 0) & !is.na(Z) & !identical(PAR$name, "H+")) {
      # Compute derivatives of omega: g and f functions (Shock et al., 1992; Johnson et al., 1992)
      rhohat <- H2O.PT$rho/1000  # just converting kg/m3 to g/cm3
      g <- gfun(rhohat, convert(T, "C"), P, H2O.PT$alpha, H2O.PT$daldT, H2O.PT$beta)
      # After SUPCRT92/reac92.f
      # eta needs to be converted to Joules! 20220325
      eta <- convert(1.66027E5, "J")
      reref <- Z^2 / (omega/eta + Z/(3.082 + 0))
      re <- reref + abs(Z) * g$g
      omega.PT <- eta * (Z^2/re - Z/(3.082 + g$g))
      Z3 <- abs(Z^3)/re^2 - Z/(3.082 + g$g)^2
      Z4 <- abs(Z^4)/re^3 - Z/(3.082 + g$g)^3
      dwdP <- (-eta * Z3 * g$dgdP)
      dwdT <- (-eta * Z3 * g$dgdT)
      d2wdT2 <- (2 * eta * Z4 * g$dgdT^2 - eta * Z3 * g$d2gdT2)
    }
    # Loop over each property
    w <- NULL
    for(i in 1:length(property)) {
      PROP <- property[i]
      # Over nonsolvation, solvation, or origination contributions
      hkf.p <- numeric(length(T))
      for(icontrib in contrib) {
        # Various contributions to the properties
        if(icontrib == "n") {
          # Nonsolvation ghs equations
          if(PROP == "H") {
            p.c <- PAR$c1*(T-Tr) - PAR$c2*(1/(T-Theta)-1/(Tr-Theta))
            p.a <- PAR$a1*(P-Pr) + PAR$a2*log((Psi+P)/(Psi+Pr)) + 
              ((2*T-Theta)/(T-Theta)^2)*(PAR$a3*(P-Pr)+PAR$a4*log((Psi+P)/(Psi+Pr)))
            p <- p.c + p.a
          } else if(PROP == "S") {
            p.c <- PAR$c1*log(T/Tr) - 
              (PAR$c2/Theta)*( 1/(T-Theta)-1/(Tr-Theta) + 
              log( (Tr*(T-Theta))/(T*(Tr-Theta)) )/Theta )
            p.a <- (T-Theta)^(-2)*(PAR$a3*(P-Pr)+PAR$a4*log((Psi+P)/(Psi+Pr)))
            p <- p.c + p.a
          } else if(PROP == "G") {
            p.c <- -PAR$c1*(T*log(T/Tr)-T+Tr) - 
              PAR$c2*( (1/(T-Theta)-1/(Tr-Theta))*((Theta-T)/Theta) - 
              (T/Theta^2)*log((Tr*(T-Theta))/(T*(Tr-Theta))) )
            p.a <- PAR$a1*(P-Pr) + PAR$a2*log((Psi+P)/(Psi+Pr)) + 
              (PAR$a3*(P-Pr) + PAR$a4*log((Psi+P)/(Psi+Pr)))/(T-Theta)
            p <- p.c + p.a
            # If the origination contribution is not NA at Tr,Pr, ensure the solvation contribution is 0, not NA
            if(!is.na(PAR$G)) p[T==Tr & P==Pr] <- 0
          # Nonsolvation cp v kt e equations
          } else if(PROP == "Cp") {
            p <- PAR$c1 + PAR$c2 * ( T - Theta ) ^ (-2)        
          } else if(PROP == "V") {
            p <- convert(PAR$a1, "cm3bar") + 
              convert(PAR$a2, "cm3bar") / (Psi + P) +
              (convert(PAR$a3, "cm3bar") + convert(PAR$a4, "cm3bar") / (Psi + P)) / (T - Theta)
          } else if(PROP == "kT") {
            p <- (convert(PAR$a2, "cm3bar") + 
              convert(PAR$a4, "cm3bar") / (T - Theta)) * (Psi + P) ^ (-2)
          } else if(PROP == "E") {
            p <- convert( - (PAR$a3 + PAR$a4 / convert((Psi + P), "joules")) * 
              (T - Theta) ^ (-2), "cm3bar")
          }
        }
        if(icontrib == "s") {
          # Solvation ghs equations
          if(PROP == "G") {
            p <- -omega.PT*(ZBorn+1) + omega*(ZBorn.PrTr+1) + omega*H2O.PrTr$YBorn*(T-Tr)
            # If the origination contribution is not NA at Tr,Pr, ensure the solvation contribution is 0, not NA
            if(!is.na(PAR$G)) p[T==Tr & P==Pr] <- 0
          }
          if(PROP == "H") 
            p <- -omega.PT*(ZBorn+1) + omega.PT*T*H2O.PT$YBorn + T*(ZBorn+1)*dwdT +
                   omega*(ZBorn.PrTr+1) - omega*Tr*H2O.PrTr$YBorn
          if(PROP == "S") 
            p <- omega.PT*H2O.PT$YBorn + (ZBorn+1)*dwdT - omega*H2O.PrTr$YBorn
          # Solvation cp v kt e equations
          if(PROP == "Cp") p <- omega.PT*T*H2O.PT$XBorn + 2*T*H2O.PT$YBorn*dwdT + 
            T*(ZBorn+1)*d2wdT2
          if(PROP == "V") p <- -convert(omega.PT, "cm3bar") * 
            H2O.PT$QBorn + convert(dwdP, "cm3bar") * (-ZBorn - 1)
          # TODO: the partial derivatives of omega are not included here here for kt and e
          # (to do it, see p. 820 of SOJ+92 ... but kt requires d2wdP2 which we don"t have yet)
          if(PROP == "kT") p <- convert(omega, "cm3bar") * H2O.PT$NBorn
          if(PROP == "E") p <- -convert(omega, "cm3bar") * H2O.PT$UBorn
        }
        if(icontrib == "o") {
          # Origination ghs equations
          if(PROP == "G") {
            p <- PAR$G - PAR$S * (T-Tr)
            # Don't inherit NA from PAR$S at Tr
            p[T==Tr] <- PAR$G
          }
          else if(PROP == "H") p <- PAR$H
          else if(PROP == "S") p <- PAR$S
          # Origination eos equations (Cp, V, kT, E): senseless
          else p <- 0 * T
        }
        # Accumulate the contribution
        hkf.p <- hkf.p + p
      }
      wnew <- data.frame(hkf.p)
      if(i > 1) w <- cbind(w, wnew) else w <- wnew
    }
    colnames(w) <- property
    aq.out[[k]] <- w
  }

  return(list(aq = aq.out, H2O = H2O.PT))
}

### Unexported functions ###

gfun <- function(rhohat, Tc, P, alpha, daldT, beta) {
  ## g and f functions for describing effective electrostatic radii of ions
  ## split from hkf() 20120123 jmd      
  ## based on equations in
  ## Shock EL, Oelkers EH, Johnson JW, Sverjensky DA, Helgeson HC, 1992
  ## Calculation of the Thermodynamic Properties of Aqueous Species at High Pressures 
  ## and Temperatures: Effective Electrostatic Radii, Dissociation Constants and 
  ## Standard Partial Molal Properties to 1000 degrees C and 5 kbar
  ## J. Chem. Soc. Faraday Trans., 88(6), 803-826  doi:10.1039/FT9928800803
  # rhohat - density of water in g/cm3
  # Tc - temperature in degrees Celsius
  # P - pressure in bars
  # start with an output list of zeros
  out0 <- numeric(length(rhohat))
  out <- list(g=out0, dgdT=out0, d2gdT2=out0, dgdP=out0)
  # only rhohat less than 1 will give results other than zero
  idoit <- rhohat < 1 & !is.na(rhohat)
  rhohat <- rhohat[idoit]
  Tc <- Tc[idoit]
  P <- P[idoit]
  alpha <- alpha[idoit]
  daldT <- daldT[idoit]
  beta <- beta[idoit]
  # eta in Eq. 1
  eta <- 1.66027E5
  # Table 3
  ag1 <- -2.037662
  ag2 <- 5.747000E-3
  ag3 <- -6.557892E-6
  bg1 <- 6.107361
  bg2 <- -1.074377E-2
  bg3 <- 1.268348E-5
  # Eq. 25
  ag <- ag1 + ag2 * Tc + ag3 * Tc ^ 2
  # Eq. 26
  bg <- bg1 + bg2 * Tc + bg3 * Tc ^ 2
  # Eq. 24
  g <- ag * (1 - rhohat) ^ bg
  # Table 4
  af1 <- 0.3666666E2
  af2 <- -0.1504956E-9
  af3 <- 0.5017997E-13
  # Eq. 33
  f <- 
    ( ((Tc - 155) / 300) ^ 4.8 + af1 * ((Tc - 155) / 300) ^ 16 ) *
    ( af2 * (1000 - P) ^ 3 + af3 * (1000 - P) ^ 4 ) 
  # limits of the f function (region II of Fig. 6)
  ifg <- Tc > 155 & P < 1000 & Tc < 355
  # in case any T or P are NA
  ifg <- ifg & !is.na(ifg)
  # Eq. 32
  g[ifg] <- g[ifg] - f[ifg]
  # at P > 6000 bar (in DEW calculations), g is zero 20170926
  g[P > 6000] <- 0
  ## now we have g at P, T
  # put the results in their right place (where rhohat < 1)
  out$g[idoit] <- g
  ## the rest is to get its partial derivatives with pressure and temperature
  ## after Johnson et al., 1992
  # alpha - coefficient of isobaric expansivity (K^-1)
  # daldT - temperature derivative of coefficient of isobaric expansivity (K^-2)
  # beta - coefficient of isothermal compressibility (bar^-1)
  # if these are NULL or NA (for IAPWS-95 and DEW), we skip the calculation
  if(is.null(alpha)) alpha <- NA
  if(is.null(daldT)) daldT <- NA
  if(is.null(beta)) beta <- NA
  # Eqn. 76
  d2fdT2 <- (0.0608/300*((Tc-155)/300)^2.8 + af1/375*((Tc-155)/300)^14) * (af2*(1000-P)^3 + af3*(1000-P)^4)
  # Eqn. 75
  dfdT <- (0.016*((Tc-155)/300)^3.8 + 16*af1/300*((Tc-155)/300)^15) * (af2*(1000-P)^3 + af3*(1000-P)^4)
  # Eqn. 74
  dfdP <- -(((Tc-155)/300)^4.8 + af1*((Tc-155)/300)^16) * (3*af2*(1000-P)^2 + 4*af3*(1000-P)^3)
  d2bdT2 <- 2 * bg3  # Eqn. 73
  d2adT2 <- 2 * ag3  # Eqn. 72
  dbdT <- bg2 + 2*bg3*Tc  # Eqn. 71
  dadT <- ag2 + 2*ag3*Tc  # Eqn. 70
  if(!all(is.na(alpha)) & !all(is.na(daldT))) {
    # Eqn. 69
    dgadT <- bg*rhohat*alpha*(1-rhohat)^(bg-1) + log(1-rhohat)*g/ag*dbdT  
    D <- rhohat
    # transcribed from SUPCRT92/reac92.f
    dDdT <- -D * alpha
    #dDdP <- D * beta
    dDdTT <- -D * (daldT - alpha^2)
    Db <- (1-D)^bg
    dDbdT <- -bg*(1-D)^(bg-1)*dDdT + log(1-D)*Db*dbdT
    dDbdTT <- -(bg*(1-D)^(bg-1)*dDdTT + (1-D)^(bg-1)*dDdT*dbdT + 
      bg*dDdT*(-(bg-1)*(1-D)^(bg-2)*dDdT + log(1-D)*(1-D)^(bg-1)*dbdT)) +
      log(1-D)*(1-D)^bg*d2bdT2 - (1-D)^bg*dbdT*dDdT/(1-D) + log(1-D)*dbdT*dDbdT
    d2gdT2 <- ag*dDbdTT + 2*dDbdT*dadT + Db*d2adT2
    d2gdT2[ifg] <- d2gdT2[ifg] - d2fdT2[ifg]
    dgdT <- g/ag*dadT + ag*dgadT  # Eqn. 67
    dgdT[ifg] <- dgdT[ifg] - dfdT[ifg]
    # phew! done with those derivatives
    out$dgdT[idoit] <- dgdT
    out$d2gdT2[idoit] <- d2gdT2
  }
  if(!all(is.na(beta))) {
    dgdP <- -bg*rhohat*beta*g*(1-rhohat)^-1  # Eqn. 66
    dgdP[ifg] <- dgdP[ifg] - dfdP[ifg]
    out$dgdP[idoit] <- dgdP
  }
  return(out)
}

