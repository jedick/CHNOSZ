# CHNOSZ/NaCl.R
# Calculate ionic strength and molalities of species given a total molality of NaCl
# 20181102 First version (jmd)
# 20181106 Use activity coefficients of Na+ and Nacl
# 20221122 Make it work with m_NaCl = 0
# 20221210 Rewritten to include pH dependence;
#   uses affinity() and equilibrate() instead of algebraic equations

NaCl <- function(m_NaCl = 1, T = 25, P = "Psat", pH = NA, attenuate = FALSE) {

  # Store existing thermo data frame
  thermo <- get("thermo", CHNOSZ)
  # Get length of longest variable
  nTP <- max(length(T), length(P), length(pH))
  pH.arg <- pH
  pH <- rep(pH, length.out = nTP)

  # Start with complete dissociation into Na+ and Cl-,
  # so ionic strength and molality of Na+ are equal to m_NaCl
  m_Naplus <- IS <- m_NaCl
  # Make them same length as T and P
  IS <- rep(IS, length.out = nTP)
  m_Naplus <- rep(m_Naplus, length.out = nTP)
  # Set tolerance for convergence to 1/100th of m_NaCl
  tolerance <- m_NaCl / 100

  # If m_NaCl is 0, return 0 for all variables 20221122
  zeros <- rep(0, nTP)
  if(m_NaCl == 0) return(list(IS = zeros, m_Naplus = zeros, m_Clminus = zeros, m_NaCl0 = zeros, m_HCl0 = zeros))

  maxiter <- 100
  for(i in 1:maxiter) {
    # Setup chemical system and calculate affinities
    if(identical(pH.arg, NA)) {
      basis(c("Na+", "Cl-", "e-"))
      species(c("Cl-", "NaCl"))
      a <- suppressMessages(affinity(T = T, P = P, "Na+" = log10(m_Naplus), IS = IS, transect = TRUE))
    } else {
      basis(c("Na+", "Cl-", "H+", "e-"))
      species(c("Cl-", "NaCl", "HCl"))
      a <- suppressMessages(affinity(T = T, P = P, pH = pH, "Na+" = log10(m_Naplus), IS = IS, transect = TRUE))
    }
    # Speciate Cl-
    e <- suppressMessages(equilibrate(a, loga.balance = log10(m_NaCl)))
    # Get molality of each Cl-bearing species
    m_Clminus <- 10^e$loga.equil[[1]]
    m_NaCl0 <- 10^e$loga.equil[[2]]
    if(identical(pH.arg, NA)) m_HCl0 <- NA else m_HCl0 <- 10^e$loga.equil[[3]]
    # Store previous ionic strength and molality of Na+
    IS_prev <- IS
    m_Naplus_prev <- m_Naplus
    # Calculate new molality of Na+ and deviation
    m_Naplus <- m_NaCl - m_NaCl0
    # Only go halfway to avoid overshoot
    if(attenuate) m_Naplus <- (m_Naplus_prev + m_Naplus) / 2
    dm_Naplus <- m_Naplus - m_Naplus_prev
    # Calculate ionic strength and deviation
    IS <- (m_Naplus + m_Clminus) / 2
    dIS <- IS - IS_prev
    # Keep going until the deviations in ionic strength and molality of Na+ at all temperatures are less than tolerance
    converged <- abs(dIS) < tolerance & abs(dm_Naplus) < tolerance
    if(all(converged)) {
      # Add one step without attenuating the deviation of m_Naplus
      if(attenuate) attenuate <- FALSE else break
    }
  }
  if(i == maxiter) {
    if(attenuate) stop(paste("reached", maxiter, "iterations without converging"))
    else stop(paste("reached", maxiter, "iterations without converging (try setting attenuate = TRUE)"))
  }

  # Restore thermo data frame
  assign("thermo", thermo, CHNOSZ)
  # Return the calculated values
  list(IS = IS, m_Naplus = m_Naplus, m_Clminus = m_Clminus, m_NaCl0 = m_NaCl0, m_HCl0 = m_HCl0)

}
