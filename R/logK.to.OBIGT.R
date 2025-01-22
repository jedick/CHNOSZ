# CHNOSZ/logK.to.OBIGT.R
# Get thermodynamic parameters from formation constants (logK) as a function of temperature
# 20220324 v1 Regress three parameters (G, S, and Cp)
# 20221202 v2 Regress HKF parameters (assume constant pressure and optimize omega parameter for charged species)

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.args.R")
#source("hkf.R")
#source("cgl.R")

logK.to.OBIGT <- function(logK, species, coeffs, T, P, name = NULL, state = "aq", V = NA, npar = 3, optimize.omega = FALSE, tolerance = 0.05) {

  # We need at least five temperature values to optimize omega
  if(optimize.omega & length(unique(T)) < 5) {
    message("logK.to.OBIGT: forcing optimize.omega = FALSE because there are < 5 T values")
    optimize.omega <- FALSE
  }
  # We need five parameters to optimize omega 20221209
  if(optimize.omega & npar != 5) {
    message("logK.to.OBIGT: forcing optimize.omega = FALSE because npar != 5")
    optimize.omega <- FALSE
  }

  ## Get the formula of the formed species (used as the species name in OBIGT)
  ### NOTE: the formed species has to be *last* in this list
  formula <- tail(species, 1)
  if(is.null(name)) name <- formula
  # Add the formed species to the thermodynamic database with ΔG°f = 0 at all temperatures
  # - but with a non-zero G if V != 0 (for minerals)
  # - zap properties of a formed species that is already in the database
  suppressMessages(inew <- mod.OBIGT(name, formula = formula, state = state, E_units = E.units(), G = 0, S = 0, Cp = 0, V = V, zap = TRUE))
  # Use species indices in case formed species has same formula and state (but different name) as a species in the database
  ispecies <- info(species)
  ispecies[length(ispecies)] <- inew
  # Calculate logK of the formation reaction with ΔG°f = 0 for the formed species
  logK0 <- suppressMessages(subcrt(ispecies, coeffs, T = T, P = P)$out$logK)

  ## Get Gibbs energy of species from logK of reaction
  # Calculate T in Kelvin
  TK <- convert(T, "K")
  # logK gives values for ΔG°r of the reaction
  Gr <- convert(logK, "G", TK)
  # logK0 gives values for ΔG°r of the reaction with ΔG°f = 0 for the formed species
  Gr0 <- convert(logK0, "G", TK)
  # Calculate ΔG°f of the formed species
  Gf <- Gr - Gr0

  if(state == "aq") {

    ## Fit to HKF 20221202

    # Get blank set of HKF parameters
    PAR <- info(info("H+"))
    # DON'T keep the name H+, becuase it tells hkf() to not calculate g function
    PAR$name <- name
    # Insert charge of formed species (for calculating g function and variable omega)
    Z <- makeup(c(formula, "Z0"), sum = TRUE)["Z"]
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
    # NOTE: we have to apply OOM scaling for HKF
    do.call(mod.OBIGT, list(name = name, formula = formula, state = state, G = G, S = S, Cp = Cp, c1 = c1, c2 = c2 * 1e-4, omega = omega * 1e-5, z = Z))

  } else if(state == "cr") {

    ## Fit to CGL (Maier-Kelley Cp equation) 20250122

    # Get blank set of CGL parameters
    PAR <- data.frame(name = NA, abbrv = NA, formula = NA, state = "cr", ref1 = NA, ref2 = NA, date = NA, model = "CGL", E_units = "J",
      G = 0, H = 0, S = 0, Cp = 0, V = 0, a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, lambda = 0, T = 0)
    # Get values of T and P (so "Psat" becomes numeric) 20221208
    tpargs <- TP.args(T = TK, P = P)
    # Calculate variable terms for S, a, b, and c
    S_var <- cgl("G", T = tpargs$T, P = tpargs$P, parameters = transform(PAR, S = 1))[[1]]$G
    a_var <- cgl("G", T = tpargs$T, P = tpargs$P, parameters = transform(PAR, a = 1))[[1]]$G
    b_var <- cgl("G", T = tpargs$T, P = tpargs$P, parameters = transform(PAR, b = 1))[[1]]$G
    c_var <- cgl("G", T = tpargs$T, P = tpargs$P, parameters = transform(PAR, c = 1))[[1]]$G

    # Default values
    c <- 0
    b <- 0
    a <- 0
    S <- 0
    if(npar == 5) {
      # Perform linear regression with all parameters
      Glm <- lm(Gf ~ S_var + a_var + b_var + c_var)
      c <- Glm$coefficients["c_var"]
      b <- Glm$coefficients["b_var"]
      a <- Glm$coefficients["a_var"]
      S <- Glm$coefficients["S_var"]
      G <- Glm$coefficients["(Intercept)"]
    } else if(npar == 4) {
      # Now with fewer parameters
      Glm <- lm(Gf ~ S_var + a_var + b_var)
      b <- Glm$coefficients["b_var"]
      a <- Glm$coefficients["a_var"]
      S <- Glm$coefficients["S_var"]
      G <- Glm$coefficients["(Intercept)"]
    } else if(npar == 3) {
      Glm <- lm(Gf ~ S_var + a_var)
      a <- Glm$coefficients["a_var"]
      S <- Glm$coefficients["S_var"]
      G <- Glm$coefficients["(Intercept)"]
    } else if(npar == 2) {
      Glm <- lm(Gf ~ S_var)
      S <- Glm$coefficients["S_var"]
      G <- Glm$coefficients["(Intercept)"]
    } else if(npar == 1) {
      G <- mean(Gf)
    } else stop("invalid value for npar (should be from 1 to 5)")

    # Calculate Cp at 25 °C (for info() and check.EOS())
    PAR$a <- a
    PAR$b <- b
    PAR$c <- c
    print(PAR)
    Cp <- cgl("Cp", parameters = PAR)[[1]]$Cp

    ## Update the thermodynamic parameters of the formed species
    # NOTE: do NOT apply OOM scaling for CGL
    # Use maximum provided T as T limit for Cp equation
    do.call(mod.OBIGT, list(name = name, formula = formula, state = state, G = G, S = S, Cp = Cp, V = V,
                            a = a, b = b, c = c, d = 0, e = 0, f = 0, lambda = 0, T = max(TK)))

  } else {
    stop(paste0("Unrecognized state (", state, "); should be aq or cr"))
  }

  # Calculate logK of the formation reaction with "real" ΔG°f for the formed species
  logK_calc <- suppressMessages(subcrt(ispecies, coeffs, T = T, P = P)$out$logK)
  print(inew)
  # Calculate the mean absolute difference
  mad <- mean(abs(logK_calc - logK))
  message(paste("logK.to.OBIGT: mean absolute difference between experimental and calculated logK is", round(mad, 4)))
  # Check that calculated values are close to input values
  stopifnot(all.equal(logK_calc, logK, tolerance = tolerance, scale = 1))
  # Return the index of the new species in OBIGT
  inew

}
