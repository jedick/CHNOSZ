# CHNOSZ/logB.to.OBIGT.R
# Get thermodynamic parameters from formation constants (logB) as a function of temperature
# 20220324 v1 Regress three parameters (G, S, and Cp)
# 20221202 v2 Regress HKF parameters (assume constant pressure and optimize omega parameter for charged species)

logB.to.OBIGT <- function(logB, species, coeffs, T, P, npar = 3, optimize.omega = FALSE, tolerance = 0.05) {

  # We need at least five temperature values to optimize omega
  if(optimize.omega & length(unique(T)) < 5) {
    message("logB.to.OBIGT: forcing optimize.omega = FALSE because there are < 5 T values")
    optimize.omega <- FALSE
  }
  # We need five parameters to optimize omega 20221209
  if(optimize.omega & npar != 5) {
    message("logB.to.OBIGT: forcing optimize.omega = FALSE because npar != 5")
    optimize.omega <- FALSE
  }

  ## Get the formula of the formed species (used as the species name in OBIGT)
  ### NOTE: the formed species has to be *last* in this list
  name <- tail(species, 1)
  # Add the formed species to the thermodynamic database with ΔG°f = 0 at all temperatures
  # Be sure to zap properties of a formed species that is already in the database
  suppressMessages(mod.OBIGT(name, formula = name, E_units = E.units(), G = 0, S = 0, Cp = 0, zap = TRUE))
  # Calculate logK of the formation reaction with ΔG°f = 0 for the formed species
  logK0 <- suppressMessages(subcrt(species, coeffs, T = T, P = P)$out$logK)

  ## Get Gibbs energy of species from logB of reaction
  # Calculate T in Kelvin
  TK <- convert(T, "K")
  # logB gives values for ΔG°r of the reaction
  Gr <- convert(logB, "G", TK)
  # logK0 gives values for ΔG°r of the reaction with ΔG°f = 0 for the formed species
  Gr0 <- convert(logK0, "G", TK)
  # Calculate ΔG°f of the formed species
  Gf <- Gr - Gr0

  ## Fit to HKF 20221202

  # Get blank set of HKF parameters
  PAR <- info(info("H+"))
  # DON'T keep the name H+, becuase it tells hkf() to not calculate g function
  PAR$name <- name
  # Insert charge of formed species (for calculating g function and variable omega)
  Z <- makeup(c(name, "Z0"), sum = TRUE)["Z"]
  # Force charge to be 0 when optimize.omega is FALSE 20221208
  if(!optimize.omega) Z <- 0
  PAR$Z <- Z
  # Get values of T and P (so "Psat" becomes numeric) 20221208
  tpargs <- TP.args(T = TK, P = P)
  # Calculate variable terms for S, c1, c2, and omega
  # NOTE: Units of Kelvin for HKF
  S_var <- hkf("G", T = tpargs$T, P = tpargs$P, parameters = transform(PAR, S = 1))$aq[[1]]$G
  c1_var <- hkf("G", T = tpargs$T, P = tpargs$P, parameters = transform(PAR, c1 = 1))$aq[[1]]$G
  c2_var <- hkf("G", T = tpargs$T, P = tpargs$P, parameters = transform(PAR, c2 = 1))$aq[[1]]$G
  omega_var <- hkf("G", T = tpargs$T, P = tpargs$P, parameters = transform(PAR, omega = 1))$aq[[1]]$G

  if(!optimize.omega) {

    # Default values
    omega <- NA
    c2 <- NA
    c1 <- NA
    S <- NA
    if(npar == 5) {
      # Perform linear regression with constant omega
      Glm <- lm(Gf ~ S_var + c1_var + c2_var + omega_var)
      omega <- Glm$coefficients["omega_var"]
      c2 <- Glm$coefficients["c2_var"]
      c1 <- Glm$coefficients["c1_var"]
      S <- Glm$coefficients["S_var"]
      G <- Glm$coefficients["(Intercept)"]
    } else if(npar == 4) {
      # Now with fewer parameters
      Glm <- lm(Gf ~ S_var + c1_var + c2_var)
      c2 <- Glm$coefficients["c2_var"]
      c1 <- Glm$coefficients["c1_var"]
      S <- Glm$coefficients["S_var"]
      G <- Glm$coefficients["(Intercept)"]
    } else if(npar == 3) {
      Glm <- lm(Gf ~ S_var + c1_var)
      c1 <- Glm$coefficients["c1_var"]
      S <- Glm$coefficients["S_var"]
      G <- Glm$coefficients["(Intercept)"]
    } else if(npar == 2) {
      Glm <- lm(Gf ~ S_var)
      S <- Glm$coefficients["S_var"]
      G <- Glm$coefficients["(Intercept)"]
    } else if(npar == 1) {
      G <- mean(Gf)
    } else stop("invalid value for npar (should be from 1 to 5)")

  } else {

    # Function to perform linear regression for a given omega(Pr,Tr)
    Gfun <- function(omega) {
      PAR$omega <- omega
      omega_var <- hkf("G", T = tpargs$T, P = tpargs$P, parameters = PAR)$aq[[1]]$G
      # Subtract omega term from Gf
      Gf_omega <- Gf - omega_var
      # Perform linear regression with all terms except omega
      lm(Gf_omega ~ S_var + c1_var + c2_var)
    }
    # Calculate the sum of squares of residuals for given omega(Pr,Tr)
    Sqfun <- function(omega) sum(Gfun(omega)$residuals ^ 2)
    # Find the omega(Pr,Tr) that minimizes the sum of squares of residuals
    omega <- optimize(Sqfun, c(-1e10, 1e10))$minimum
    # Construct the linear model with this omega
    Glm <- Gfun(omega)
    # Get coefficients other than omega
    G <- Glm$coefficients["(Intercept)"]
    S <- Glm$coefficients["S_var"]
    c1 <- Glm$coefficients["c1_var"]
    c2 <- Glm$coefficients["c2_var"]

  }

  # NAs propagate as zero in the HKF equations
  if(is.na(S)) S <- 0
  if(is.na(c1)) c1 <- 0
  if(is.na(c2)) c2 <- 0
  if(is.na(omega)) omega <- 0
  # Calculate Cp at 25 °C (not used in HKF - just for info() and check.EOS())
  PAR$c1 <- c1
  PAR$c2 <- c2
  PAR$omega <- omega
  Cp <- hkf("Cp", parameters = PAR)$aq[[1]]$Cp

  ## Update the thermodynamic parameters of the formed species
  # NOTE: we have to apply OOM scaling
  ispecies <- do.call(mod.OBIGT, list(name = name, G = G, S = S, Cp = Cp, c1 = c1, c2 = c2 * 1e-4, omega = omega * 1e-5, z = Z))
  # Calculate logK of the formation reaction with "real" ΔG°f for the formed species
  logK <- suppressMessages(subcrt(species, coeffs, T = T, P = P)$out$logK)
  # Calculate the mean absolute difference
  mad <- mean(abs(logK - logB))
  message(paste("logB.to.OBIGT: mean absolute difference between logB (experimental) and logK (calculated) is", round(mad, 4)))
  # Check that calculated values are close to input values
  stopifnot(all.equal(logK, logB, tolerance = tolerance, scale = 1))
  # Return the species index in OBIGT
  ispecies

}
