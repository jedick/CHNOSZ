# CHNOSZ/logB_to_OBIGT.R
# Fit G, S, and Cp to formation constants (logB) as a function of temperature
# 20220324 jmd

logB_to_OBIGT <- function(logB, species, coeff, T, P, tolerance = 0.05) {

  # We need at least three temperature values
  if(length(unique(T)) < 3) stop("at least 3 unique temperature values are needed")

  ## Get the formula of the formed species (used as the species name in OBIGT)
  ### NOTE: the formed species has to be *last* in this list
  name <- tail(species, 1)
  # Add the formed species to the thermodynamic database with ΔG°f = 0 at all temperatures
  # Be sure to zap properties of a formed species that is already in the database
  suppressMessages(mod.OBIGT(name, formula = name, E_units = E.units(), G = 0, S = 0, Cp = 0, zap = TRUE))
  # Calculate logK of the formation reaction with ΔG°f = 0 for the formed species
  logK0 <- suppressMessages(subcrt(species, coeff, T = T, P = P)$out$logK)

  ## Get Gibbs energy of species from logB of reaction
  # Calculate T in Kelvin
  TK <- convert(T, "K")
  # logB gives values for ΔG°r of the reaction
  Gr <- convert(logB, "G", TK)
  # logK0 gives values for ΔG°r of the reaction with ΔG°f = 0 for the formed species
  Gr0 <- convert(logK0, "G", TK)
  # Calculate ΔG°f of the formed species
  Gf <- Gr - Gr0

  ## Solve for G, S, and Cp
  # Make an 'lm' model object for given Cp
  Gfun <- function(Cp = 0) {
    Tr <- 298.15
    TTr <- TK - Tr
    # Subtract Cp term from Gf
    GCp <- Cp * (TK - Tr - TK * log(TK / Tr))
    GCp[is.na(GCp)] <- 0
    GfCp <- Gf - GCp
    # Write linear model in Ttr -- slope is -S
    lm(GfCp ~ TTr)
  }
  # Calculate the sum of squares of residuals for given Cp
  Sqfun <- function(Cp) sum(Gfun(Cp)$residuals ^ 2)
  # Find the Cp with minimum sum of squares of residuals
  Cp <- optimize(Sqfun, c(-100, 100))$minimum
  # Calculate the fitted G and S for this Cp
  G <- Gfun(Cp)$coefficients[[1]]
  S <- - Gfun(Cp)$coefficients[[2]]

  ## Update the thermodynamic parameters of the formed species
  ispecies <- do.call(mod.OBIGT, list(name = name, G = G, S = S, Cp = Cp))
  # Calculate logK of the formation reaction with "real" ΔG°f for the formed species
  logK <- suppressMessages(subcrt(species, coeff, T = T, P = P)$out$logK)
  # Calculate the mean absolute difference
  mad <- mean(abs(logK - logB))
  message(paste("logB_to_OBIGT: mean absolute difference between logB (experimental) and logK (calculated) is", round(mad, 4)))
  # Check that calculated values are close to input values
  stopifnot(all.equal(logK, logB, tolerance = tolerance, scale = 1))
  # Return the species index in OBIGT
  ispecies

}
