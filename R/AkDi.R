# CHNOSZ/AkDi.R
# Akinfiev-Diamond model for aqueous species
# 20190219 First version
# 20220206 Calculate S, Cp, and V

AkDi <- function(property = NULL, parameters = NULL, T = 298.15, P = 1, isPsat = TRUE) {

  # Some constants (from Akinfiev and Diamond, 2004 doi:10.1016/j.fluid.2004.06.010)
  MW <- 18.0153 # g mol-1
  NW <- 1000 / MW # mol kg-1
  R <- 8.31441 # J K-1 mol-1
  # R expressed in volume units
  RV <- 10 * R # cm3 bar K-1 mol-1

  # Calculate properties of H2O
  # Moved these calculations to unexported functions
  # to be able to test their output in test-AkDi.R 20220206
  f1 <- .f1(T, P, isPsat)
  rho1 <- .rho1(T, P)
  drho1_dT <- mapply(.drho1_dT, T, P)
  drho1_dP <- mapply(.drho1_dP, T, P)
  d2rho1_dT2 <- mapply(.d2rho1_dT2, T, P)
  S1_g <- sapply(T, .S1_g)
  S1 <- mapply(.S1, T, P, MoreArgs = list(isPsat = isPsat))

#ADstuff <<- list(f1 = f1, drho1_dT = drho1_dT, drho1_dP = drho1_dP, d2rho1_dT2 = d2rho1_dT2, S1_g = S1_g, S1 = S1)

  # Initialize a list for the output
  out <- list()
  # Loop over species
  nspecies <- nrow(parameters)
  for(i in seq_len(nspecies)) {
    PAR <- parameters[i, ]
    # Send a message
    message("AkDi: Akinfiev-Diamond model for ", PAR$name, " gas to aq")
    # Start with an NA-filled data frame
    myprops <- as.data.frame(matrix(NA, ncol = length(property), nrow = length(T)))
    colnames(myprops) <- property

    # Loop over properties
    for(j in seq_along(property)) {

      if(property[[j]] == "G") {
        # Get gas properties (J mol-1)
        G_gas <- subcrt(PAR$name, "gas", T = T, P = P, convert = FALSE)$out[[1]]$G
        # TODO: Does this work if E.units is cal or J?
        G_gas <- convert(G_gas, "J", T = T)
        # Calculate G_hyd (J mol-1)
        G_hyd <- R*T * ( -log(NW) + (1 - PAR$xi) * log(f1) + PAR$xi * log(RV * T * rho1 / MW) + rho1 * (PAR$a + PAR$b * (1000/T)^0.5) )
        # Calculate the chemical potential (J mol-1)
        G <- G_gas + G_hyd
        # Convert J to cal
        G <- convert(G, "cal", T = T)
        # Insert into data frame of properties
        myprops$G <- G
      }

      if(property[[j]] == "S") {
        # Get S_gas
        S_gas <- subcrt(PAR$name, "gas", T = T, P = P, convert = FALSE)$out[[1]]$S
        S_gas <- convert(S_gas, "J", T = T)
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
        S <- convert(S, "cal", T = T)
        myprops$S <- S
      }

    }
    out[[i]] <- myprops
  }
  out
}

### UNEXPORTED FUNCTIONS ###

.f1 <- function(T, P, isPsat) {
  # Get H2O fugacity (bar)
  GH2O_P <- water("G", T = T, P = P)$G
  GH2O_1 <- water("G", T = T, P = 1)$G
  f1 <- exp ( (GH2O_P - GH2O_1) / (1.9872 * T) )
  # For Psat, calculate the real liquid-vapor curve (not 1 bar below 100 degC)
  if(isPsat) {
    P <- water("Psat", T = T, P = "Psat", P1 = FALSE)$Psat
    f1[P < 1] <- P[P < 1]
  }
  f1
}

.rho1 <- function(T, P) {
  # Density of H2O (g cm-3)
  water("rho", T = T, P = P)$rho / 1000
}

.drho1_dT <- function(T, P) {
  # Partial derivative of density with respect to temperature at constant pressure 20220206
  dT <- 0.1
  T1 <- T - dT
  T2 <- T + dT
  # Add 1 bar to P so the derivative doesn't blow up when P = Psat at T > 100 degC
  # TODO: Is there a better way?
  rho1 <- .rho1(c(T1, T2), P + 1)
  diff(rho1) / (T2 - T1)
}

.drho1_dP <- function(T, P) {
  # Partial derivative of density with respect to pressure at constant temperature 20220206
  dP <- 0.1
  P1 <- P - dP
  P2 <- P + dP
  rho1 <- .rho1(T, c(P1, P2))
  diff(rho1) / (P2 - P1)
}

.d2rho1_dT2 <- function(T, P) {
  # Second partial derivative of density with respect to temperature at constant pressure 20220206
  dT <- 0.1
  # Calculate density at five temperature values
  Tval <- seq(T - 2 * dT, T + 2 * dT, dT)
  rho1 <- .rho1(Tval, P)
  # https://stackoverflow.com/questions/11081069/calculate-the-derivative-of-a-data-function-in-r
  spl <- smooth.spline(Tval, rho1)
  # The second derivative of the fitted spline function at the third point (i.e., T)
  predict(spl, deriv = 2)$y[3]
}

.S1_g <- function(T) {
  # Entropy of the ideal gas (dimensionless)
  .Fortran(C_ideal2, T, 0, 0, 0, 0, 0, 0, 0)[[4]]
}

.S1 <- function(T, P, isPsat) {
  # Entropy of solvent H2O (dimensionless)
  if(isPsat) S1 <- water("S", T = T, P = "Psat")$S / 1.9872
  else S1 <- water("S", T = T, P = P)$S / 1.9872
  S1
}
