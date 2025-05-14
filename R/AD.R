# CHNOSZ/AD.R
# Akinfiev-Diamond model for aqueous species
# 20190219 First version
# 20220206 Calculate S, Cp, and V

AD <- function(property = NULL, parameters = NULL, T = 298.15, P = 1, isPsat = TRUE) {

  # Some constants (from Akinfiev and Diamond, 2004 doi:10.1016/j.fluid.2004.06.010)
  MW <- 18.0153 # g mol-1
  NW <- 1000 / MW # mol kg-1
  #R <- 8.31441 # J K-1 mol-1  20190219
  R <- 8.314463  # https://physics.nist.gov/cgi-bin/cuu/Value?r 20230630
  # R expressed in volume units
  RV <- 10 * R # cm3 bar K-1 mol-1

  # Calculate H2O fugacity and derivatives of density
  # These calculations are done through (unexported) functions
  # to be able to test their output in test-AD.R 20220206
  f1 <- .f1(T, P, isPsat)
  drho1_dT <- mapply(.drho1_dT, T, P, MoreArgs = list(isPsat = isPsat))
  drho1_dP <- mapply(.drho1_dP, T, P, MoreArgs = list(isPsat = isPsat))
  d2rho1_dT2 <- mapply(.d2rho1_dT2, T, P, MoreArgs = list(isPsat = isPsat))

  # Calculate other properties of H2O solvent
  waterTP <- water(c("rho", "S", "Cp", "V"), T = T, P = P)
  # Density (g cm-3)
  rho1 <- waterTP$rho / 1000
  # Entropy (dimensionless)
  S1 <- waterTP$S / R
  # Heat capacity (dimensionless)
  Cp1 <- waterTP$Cp / R
  # Volume (cm3 mol-1)
  V1 <- waterTP$V
  # Calculate properties of ideal H2O gas
  S1_g <- sapply(T, .S1_g)
  Cp1_g <- sapply(T, .Cp1_g)

  # Initialize a list for the output
  out <- list()
  # Loop over species
  nspecies <- nrow(parameters)
  for(i in seq_len(nspecies)) {

    # Get thermodynamic parameters for the gas and calculate properties at T, P
    PAR <- parameters[i, ]
    gasprops <- subcrt(PAR$name, "gas", T = T, P = P, convert = FALSE)$out[[1]]
    # Send a message
    message("AD: Akinfiev-Diamond model for ", PAR$name, " gas to aq")
    # Start with an NA-filled data frame
    myprops <- as.data.frame(matrix(NA, ncol = length(property), nrow = length(T)))
    colnames(myprops) <- property

    # Loop over properties
    for(j in seq_along(property)) {

      if(property[[j]] == "G") {
        # Get gas properties (J mol-1)
        G_gas <- gasprops$G
        # Calculate G_hyd (J mol-1)
        G_hyd <- R*T * ( -log(NW) + (1 - PAR$xi) * log(f1) + PAR$xi * log(RV * T * rho1 / MW) + rho1 * (PAR$a + PAR$b * (1000/T)^0.5) )
        # Calculate the chemical potential (J mol-1)
        G <- G_gas + G_hyd
        # Insert into data frame of properties
        myprops$G <- G
      }

      if(property[[j]] == "S") {
        # Get S_gas
        S_gas <- gasprops$S
        # Calculate S_hyd
        S_hyd <- R * (
            (1 - PAR$xi) * (S1 - S1_g)
          + log(NW)
          - (PAR$xi + PAR$xi * log(RV * T / MW) + PAR$xi * log(rho1) + PAR$xi * T / rho1 * drho1_dT)
          - (
               PAR$a * (rho1 + T * drho1_dT)
             + PAR$b * (0.5 * 10^1.5 * T^-0.5 * rho1 + 10^1.5 * T^0.5 * drho1_dT)
          )
        )
        S <- S_gas + S_hyd
        myprops$S <- S
      }

      if(property[[j]] == "Cp") {
        # Get Cp_gas
        Cp_gas <- gasprops$Cp
        # Calculate Cp_hyd
        Cp_hyd <- R * (
            (1 - PAR$xi) * (Cp1 - Cp1_g)
          - (PAR$xi + 2 * PAR$xi * T / rho1 * drho1_dT - PAR$xi * T^2 / rho1^2 * drho1_dT^2 + PAR$xi * T^2 / rho1 * d2rho1_dT2)
        ) - R*T * (
            PAR$a * (2 * drho1_dT + T * d2rho1_dT2)
          + PAR$b * (-0.25 * 10^1.5 * T^-1.5 * rho1 + 10^1.5 * T^-0.5 * drho1_dT + 10^1.5 * T^0.5 * d2rho1_dT2)
        )
        Cp <- Cp_gas + Cp_hyd
        myprops$Cp <- Cp
      }

      if(property[[j]] == "V") {
        # Get V_gas
        V_gas <- 0
        # Calculate V_hyd
        V_hyd <- V1 * (1 - PAR$xi) + PAR$xi * RV * T / rho1 * drho1_dP + RV * T * drho1_dP * (PAR$a + PAR$b * (1000/T)^0.5)
        V <- V_gas + V_hyd
        myprops$V <- V
      }

    }
    # Calculate enthalpy 20220206
    myprops$H <- myprops$G - 298.15 * entropy(PAR$formula) + T * myprops$S
    out[[i]] <- myprops
  }
  out
}

### UNEXPORTED FUNCTIONS ###

.f1 <- function(T, P, isPsat) {
  # Get H2O fugacity (bar)
  GH2O_P <- water("G", T = T, P = P)$G
  GH2O_1 <- water("G", T = T, P = 1)$G
  f1 <- exp ( (GH2O_P - GH2O_1) / (8.31441 * T) )
  # For Psat, calculate the real liquid-vapor curve (not floored at 1 bar)
  if(isPsat) {
    P <- water("Psat", T = T, P = "Psat", Psat_floor = NULL)$Psat
    f1[P < 1] <- P[P < 1]
  }
  f1
}

.rho1 <- function(T, P) {
  # Density of H2O (g cm-3)
  water("rho", T = T, P = P)$rho / 1000
}

.drho1_dT <- function(T, P, isPsat) {
  # Partial derivative of density with respect to temperature at constant pressure 20220206
  dT <- 0.1
  T1 <- T - dT
  T2 <- T + dT
  # Add 1 bar to P so the derivative doesn't blow up when P = Psat at T > 100 degC
  # TODO: Is there a better way?
  if(isPsat) P <- P +1
  rho1 <- .rho1(c(T1, T2), P)
  diff(rho1) / (T2 - T1)
}

.drho1_dP <- function(T, P, isPsat) {
  # Partial derivative of density with respect to pressure at constant temperature 20220206
  dP <- 0.1
  P1 <- P - dP
  P2 <- P + dP
  # Subtract 1 degC from T so the derivative doesn't blow up when P = Psat at T > 100 degC
  # TODO: Is there a better way?
  if(isPsat) T <- T - 1
  rho1 <- .rho1(T, c(P1, P2))
  diff(rho1) / (P2 - P1)
}

.d2rho1_dT2 <- function(T, P, isPsat) {
  # Second partial derivative of density with respect to temperature at constant pressure 20220206
  # NOTE: dT, Tval, and P <- P + 1 are chosen to produce demo/AD.R;
  # these may not be the best settings for other T-P ranges  20220207
  dT <- 0.2
  # Calculate density at seven temperature values
  Tval <- seq(T - 3 * dT, T + 3 * dT, dT)
  # TODO: Is there a better way to calculate the partial derivative for P = Psat?
  if(isPsat) P <- P + 2 else P <- P + 1
  rho1 <- .rho1(Tval, P)
  # At P = 281 bar there are identical values of rho1 between T = 683.15 and 683.35 K
  # Don't allow duplicates because they produce artifacts in the second derivative 20220207
  if(any(duplicated(rho1))) {
    message(paste("AD: detected identical values of rho1 in second derivative calculation; returning NA at", T, "K and", P, "bar"))
    return(NA)
  }
  # https://stackoverflow.com/questions/11081069/calculate-the-derivative-of-a-data-function-in-r
  spl <- smooth.spline(Tval, rho1)
  # The second derivative of the fitted spline function at the fourth point (i.e., T)
  predict(spl, deriv = 2)$y[4]
}

.S1_g <- function(T) {
  # Entropy of the ideal gas (dimensionless)
  .Fortran(C_ideal2, T, 0, 0, 0, 0, 0, 0, 0)[[4]]
}

.Cp1_g <- function(T) {
  # Heat capacity of the ideal gas (dimensionless)
  .Fortran(C_ideal2, T, 0, 0, 0, 0, 0, 0, 0)[[8]]
}
