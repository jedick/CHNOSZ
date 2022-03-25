# CHNOSZ/Berman.R 20170930
# Calculate thermodynamic properties of minerals using equations from:
#   Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals
#      in the system Na2O-K2O-CaO-MgO-FeO-Fe2O3-Al2O3-SiO2-TiO2-H2O-CO2.
#      J. Petrol. 29, 445-522. https://doi.org/10.1093/petrology/29.2.445

Berman <- function(name, T = 298.15, P = 1, check.G=FALSE, calc.transition=TRUE, calc.disorder=TRUE) {
  # Reference temperature and pressure
  Pr <- 1
  Tr <- 298.15
  # Make T and P the same length
  ncond <- max(length(T), length(P))
  T <- rep(T, length.out=ncond)
  P <- rep(P, length.out=ncond)

  # Get parameters in the Berman equations
  # Start with thermodynamic parameters provided with CHNOSZ
  dat <- thermo()$Berman
  # Is there a user-supplied data file?
  userfile <- get("thermo", CHNOSZ)$opt$Berman
  userfileexists <- FALSE
  if(!is.na(userfile)) {
    if(userfile!="") {
      if(file.exists(userfile)) {
        userfileexists <- TRUE
        BDat_user <- read.csv(userfile, as.is=TRUE)
        dat <- rbind(BDat_user, dat)
      } else stop("the file named in thermo()$opt$Berman (", userfile, ") does not exist")
    } 
  }
  # Remove duplicates (only the first, i.e. most recent entry is kept)
  dat <- dat[!duplicated(dat$name), ]
  # Remove the multipliers on volume parameters
  vcols <- 13:16 # columns with v1, v2, v3, v4
  multexp <- c(5, 5, 5, 8)
  dat[, vcols] <- t(t(dat[, vcols]) / 10^multexp)
  # If name is missing, return the entire data frame (used in test-Berman.R)
  if(missing(name)) return(dat) else {
    # Which row has data for this mineral?
    irow <- which(dat$name == name)
    if(length(irow)==0) {
      if(userfileexists) stop("Data for ", name, " not available. Please add it to ", userfile)
      if(!userfileexists) stop("Data for ", name, " not available. Please add it to your_data_file.csv and run thermo('opt$Berman' = 'path/to/your_data_file.csv')")
    }
    dat <- dat[irow, ]
  }

  # The function works fine with just the following assign() call,
  # but an explicit dummy assignment here is used to avoid "Undefined global functions or variables" in R CMD check
  GfPrTr <- HfPrTr <- SPrTr <- Tlambda <- Tmax <- Tmin <- Tref <- VPrTr <-
    d0 <- d1 <- d2 <- d3 <- d4 <- Vad <- dTdP <- k0 <- k1 <- k2 <- k3 <-
    k4 <- k5 <- k6 <- l1 <- l2 <- v1 <- v2 <- v3 <- v4 <- NA
  # Assign values to the variables used below
  for(i in 1:ncol(dat)) assign(colnames(dat)[i], dat[, i])
  # Get the entropy of the elements using the chemical formula in thermo()$OBIGT
  OBIGT <- thermo()$OBIGT
  formula <- OBIGT$formula[match(name, OBIGT$name)]
  SPrTr_elements <- entropy(formula)
  # Check that G in data file is the G of formation from the elements --> Benson-Helgeson convention (DG = DH - T*DS)
  if(check.G) {
    GfPrTr_calc <- HfPrTr - Tr * (SPrTr - SPrTr_elements)
    Gdiff <- GfPrTr_calc - GfPrTr
    #if(is.na(GfPrTr)) warning(paste0(name, ": GfPrTr(table) is NA"), call.=FALSE)
    if(!is.na(GfPrTr)) if(abs(Gdiff) >= 1000) warning(paste0(name, ": GfPrTr(calc) - GfPrTr(table) is too big! == ",
                                          round(GfPrTr_calc - GfPrTr), " J/mol"), call.=FALSE)
    # (the tabulated GfPrTr is unused below)
  }

  ### Thermodynamic properties ###
  # Calculate Cp and V (Berman, 1988 Eqs. 4 and 5)
  # k4, k5, k6 terms from winTWQ documentation (doi:10.4095/223425)
  Cp <- k0 + k1 * T^-0.5 + k2 * T^-2 + k3 * T^-3 + k4 * T^-1 + k5 * T + k6 * T^2
  P_Pr <- P - Pr
  T_Tr <- T - Tr
  V <- VPrTr * (1 + v1 * T_Tr + v2 * T_Tr^2 + v3 * P_Pr + v4 * P_Pr^2)
  ## Calculate Ga (Ber88 Eq. 6) (superseded 20180328 as it does not include k4, k5, k6)
  #Ga <- HfPrTr - T * SPrTr + k0 * ( (T - Tr) - T * (log(T) - log(Tr)) ) +
  #  2 * k1 * ( (T^0.5 - Tr^0.5) + T*(T^-0.5 - Tr^-0.5) ) -
  #  k2 * ( (T^-1 - Tr^-1) - T / 2 * (T^-2 - Tr^-2) ) -
  #  k3 * ( (T^-2 - Tr^-2) / 2 - T / 3 * (T^-3 - Tr^-3) ) +
  #  VPrTr * ( (v3 / 2 - v4) * (P^2 - Pr^2) + v4 / 3 * (P^3 - Pr^3) +
  #    (1 - v3 + v4 + v1 * (T - Tr) + v2 * (T - Tr)^2) * (P - Pr) )
  # Calculate Ha (symbolically integrated using sympy - expressions not simplified)
  intCp <- T*k0 - Tr*k0 + k2/Tr - k2/T + k3/(2*Tr^2) - k3/(2*T^2) + 2.0*k1*T^0.5 - 2.0*k1*Tr^0.5 + 
    k4*log(T) - k4*log(Tr) + k5*T^2/2 - k5*Tr^2/2 - k6*Tr^3/3 + k6*T^3/3
  intVminusTdVdT <- -VPrTr + P*(VPrTr + VPrTr*v4 - VPrTr*v3 - Tr*VPrTr*v1 + VPrTr*v2*Tr^2 - VPrTr*v2*T^2) +
    P^2*(VPrTr*v3/2 - VPrTr*v4) + VPrTr*v3/2 - VPrTr*v4/3 + Tr*VPrTr*v1 + VPrTr*v2*T^2 - VPrTr*v2*Tr^2 + VPrTr*v4*P^3/3
  Ha <- HfPrTr + intCp + intVminusTdVdT
  # Calculate S (also symbolically integrated)
  intCpoverT <- k0*log(T) - k0*log(Tr) - k3/(3*T^3) + k3/(3*Tr^3) + k2/(2*Tr^2) - k2/(2*T^2) + 2.0*k1*Tr^-0.5 - 2.0*k1*T^-0.5 +
    k4/Tr - k4/T + T*k5 - Tr*k5 + k6*T**2/2 - k6*Tr**2/2
  intdVdT <- -VPrTr*(v1 + v2*(-2*Tr + 2*T)) + P*VPrTr*(v1 + v2*(-2*Tr + 2*T))
  S <- SPrTr + intCpoverT - intdVdT
  # Calculate Ga --> Berman-Brown convention (DG = DH - T*S, no S(element))
  Ga <- Ha - T * S

  ### Polymorphic transition properties ***
  if(!is.na(Tlambda) & !is.na(Tref) & any(T > Tref) & calc.transition) {
    # Starting transition contributions are 0
    Cptr <- Htr <- Str <- Gtr <- numeric(ncond)
    ## Ber88 Eq. 8: Cp at 1 bar
    #Cplambda_1bar <- T * (l1 + l2 * T)^2
    # Eq. 9: Tlambda at P
    Tlambda_P <- Tlambda + dTdP * (P - 1)
    # Eq. 8a: Cp at P
    Td <- Tlambda - Tlambda_P
    Tprime <- T + Td
    # With the condition that Tref < Tprime < Tlambda(1bar)
    iTprime <- Tref < Tprime & Tprime < Tlambda
    # Handle NA values (arising from NA in input P values e.g. Psat above Tcritical) 20180925
    iTprime[is.na(iTprime)] <- FALSE
    Tprime <- Tprime[iTprime]
    Cptr[iTprime] <- Tprime * (l1 + l2 * Tprime)^2
    # We got Cp, now calculate the integrations for H and S
    # The lower integration limit is Tref
    iTtr <- T > Tref
    Ttr <- T[iTtr]
    Tlambda_P <- Tlambda_P[iTtr]
    Td <- Td[iTtr]
    # Handle NA values 20180925
    Tlambda_P[is.na(Tlambda_P)] <- Inf
    # The upper integration limit is Tlambda_P
    Ttr[Ttr >= Tlambda_P] <- Tlambda_P[Ttr >= Tlambda_P]
    # Derived variables
    tref <- Tref - Td
    x1 <- l1^2 * Td + 2 * l1 * l2 * Td^2 + l2^2 * Td^3
    x2 <- l1^2 + 4 * l1 * l2 * Td + 3 * l2^2 * Td^2
    x3 <- 2 * l1 * l2 + 3 * l2^2 * Td
    x4 <- l2 ^ 2
    # Eqs. 10, 11, 12
    Htr[iTtr] <- x1 * (Ttr - tref) + x2 / 2 * (Ttr^2 - tref^2) + x3 / 3 * (Ttr^3 - tref^3) + x4 / 4 * (Ttr^4 - tref^4)
    Str[iTtr] <- x1 * (log(Ttr) - log(tref)) + x2 * (Ttr - tref) + x3 / 2 * (Ttr^2 - tref^2) + x4 / 3 * (Ttr^3 - tref^3)
    Gtr <- Htr - T * Str
    # Apply the transition contributions
    Ga <- Ga + Gtr
    Ha <- Ha + Htr
    S <- S + Str
    Cp <- Cp + Cptr
  }

  ### Disorder thermodynamic properties ###
  if(!is.na(Tmin) & !is.na(Tmax) & any(T > Tmin) & calc.disorder) {
    # Starting disorder contributions are 0
    Cpds <- Hds <- Sds <- Vds <- Gds <- numeric(ncond)
    # The lower integration limit is Tmin
    iTds <- T > Tmin
    Tds <- T[iTds]
    # The upper integration limit is Tmax
    Tds[Tds > Tmax] <- Tmax
    # Ber88 Eqs. 15, 16, 17
    Cpds[iTds] <- d0 + d1*Tds^-0.5 + d2*Tds^-2 + d3*Tds + d4*Tds^2
    Hds[iTds] <- d0*(Tds - Tmin) + d1*(Tds^0.5 - Tmin^0.5)/0.5 +
      d2*(Tds^-1 - Tmin^-1)/-1 + d3*(Tds^2 - Tmin^2)/2 + d4*(Tds^3 - Tmin^3)/3
    Sds[iTds] <- d0*(log(Tds) - log(Tmin)) + d1*(Tds^-0.5 - Tmin^-0.5)/-0.5 +
      d2*(Tds^-2 - Tmin^-2)/-2 + d3*(Tds - Tmin) + d4*(Tds^2 - Tmin^2)/2
    # "d5 is a constant computed in such as way as to scale the disordring enthalpy to the volume of disorder" (Berman, 1988)
    # 20180331: however, having a "d5" that isn't a coefficient in the same equation as d0, d1, d2, d3, d4 is confusing notation
    # therefore, CHNOSZ now uses "Vad" for this variable, following the notation in the Theriak-Domino manual
    # Eq. 18; we can't do this if Vad == 0 (dolomite and gehlenite)
    if(Vad != 0) Vds <- Hds / Vad
    # Berman puts the Vds term directly into Eq. 19 (commented below), but that necessarily makes Gds != Hds - T * Sds
    #Gds <- Hds - T * Sds + Vds * (P - Pr)
    # Instead, we include the Vds term with Hds
    Hds <- Hds + Vds * (P - Pr)
    # Disordering properties above Tmax (Eq. 20)
    ihigh <- T > Tmax
    # Again, Berman put the Sds term (for T > Tmax) into Eq. 20 for Gds (commented below), which would also make Gds != Hds - T * Sds
    #Gds[ihigh] <- Gds[ihigh] - (T[ihigh] - Tmax) * Sds[ihigh]
    # Instead, we add the Sds[ihigh] term to Hds
    Hds[ihigh] <- Hds[ihigh] - (T[ihigh] - Tmax) * Sds[ihigh]
    # By writing Gds = Hds - T * Sds, the above two changes w.r.t. Berman's
    # equations affect the computed values only for Hds, not Gds
    Gds <- Hds - T * Sds
    # Apply the disorder contributions
    Ga <- Ga + Gds
    Ha <- Ha + Hds
    S <- S + Sds
    V <- V + Vds
    Cp <- Cp + Cpds
  }

  ### (for testing) Use G = H - TS to check that integrals for H and S are written correctly
  Ga_fromHminusTS <- Ha - T * S
  if(!isTRUE(all.equal(Ga_fromHminusTS, Ga))) stop(paste0(name, ": incorrect integrals detected using DG = DH - T*S"))

  ### Thermodynamic and unit conventions used in SUPCRT ###
  # Use entropy of the elements in calculation of G --> Benson-Helgeson convention (DG = DH - T*DS)
  Gf <- Ga + Tr * SPrTr_elements
  # The output will just have "G" and "H"
  G <- Gf
  H <- Ha
  # Convert J/bar to cm^3/mol
  V <- V * 10

  data.frame(T = T, P = P, G = G, H = H, S = S, Cp = Cp, V = V)
}
