# NaCl.R
# Calculate ion activities and ionic strength
# for a total molality of NaCl,
# taking account of ion association: Na+ + Cl- = NaCl(aq)
# 20181102 jmd first version
# 20181105 use activity coefficient of Na+
# 20181106 use activity coefficient of NaCl
# 20221122 make it work with m_tot = 0

NaCl <- function(T = seq(100, 500, 100), P = 1000, m_tot = 2, ...) {
  # Define a function for the reaction quotient
  logQ <- function(m_Cl, gam_NaCl, gam_Na, gam_Cl) {
    # Starting with Q = a_NaCl / (a_Na+ * a_Cl-),
    # substitute m_tot = m_NaCl + m_Cl and m_Cl = m_Na
    # to write:
    log10( (m_tot - m_Cl) * gam_NaCl / (m_Cl ^ 2 * gam_Na * gam_Cl) )
  }
  # Define a function for affinity = log(K / Q)
  A <- function(m_Cl, gam_NaCl, gam_Na, gam_Cl, logK) logK - logQ(m_Cl, gam_NaCl, gam_Na, gam_Cl)
  # Calculate equilibrium constant at all temperatures (standard conditions: IS = 0)
  logK <- subcrt(c("Na+", "Cl-", "NaCl"), c(-1, -1, 1), T = T, P = P, ...)$out$logK
  # Calculate Debye-Huckel parameters at all temperatures
  wout <- water(c("A_DH", "B_DH"), T = convert(T, "K"), P = P)
  # Initialize output variables
  N <- length(T)
  ISout <- a_Cl <- numeric(N)
  # Initial guess for m_Cl and ionic strength assuming complete dissociation of NaCl
  IS <- m_Cl <- rep(m_tot, N)
  # The corresponding total molality of dissolved species (NaCl + Cl- + Na+)
  m_star <- (m_tot - m_Cl) + 2*m_Cl
  # The species indices for Na+, Cl-, and NaCl(aq)
  ispecies <- info(c("Na+", "Cl-", "NaCl"))
  # Data frame needed for nonideal()
  speciesprops <- rep(list(data.frame(G = numeric(N))), length(ispecies))
  # Only try to find root for non-zero NaCl concentration 20221122
  if(m_tot > 0) {
    # Start by doing calculations for all temperatures
    doit <- !logical(N)
    while(any(doit)) {
      # Calculate activity coefficient at given ionic strength
      gammas <- suppressMessages(nonideal(ispecies, speciesprops, IS = IS, T = convert(T, "K"), P = P, A_DH = wout$A_DH, B_DH = wout$B_DH, m_star = m_star))
      gam_Na <- 10^gammas[[1]]$loggam
      gam_Cl <- 10^gammas[[2]]$loggam
      gam_NaCl <- 10^gammas[[3]]$loggam
      # In case Setchenow calculations are turned off 20191209
      if(length(gam_NaCl) == 0) gam_NaCl <- rep(1, length(T))
      # Solve for m_Cl
      for(i in which(doit)) m_Cl[i] <- uniroot(A, c(0, m_tot), gam_NaCl = gam_NaCl[i], gam_Na = gam_Na[i], gam_Cl = gam_Cl[i], logK = logK[i])$root
      # Calculate new total molality
      m_star <- (m_tot - m_Cl) + 2*m_Cl
      # Calculate new ionic strength and deviation
      ISnew <- m_Cl
      dIS <- ISnew - IS
      # Set new ionic strength
      IS <- ISnew
      # Keep going until the deviation in ionic strength at any temperature is less than 0.01
      doit <- abs(dIS) > 0.01
    }
  }
  # Assemble final molality of Cl- and gammas
  gammas <- suppressMessages(nonideal(ispecies, speciesprops, IS = IS, T = convert(T, "K"), P = P, A_DH = wout$A_DH, B_DH = wout$B_DH))
  gam_Na <- 10^gammas[[1]]$loggam
  gam_Cl <- 10^gammas[[2]]$loggam
  gam_NaCl <- 10^gammas[[3]]$loggam
  # Return the calculated values
  list(IS = IS, m_Cl = m_Cl, gam_Na = gam_Na, gam_Cl = gam_Cl, gam_NaCl = gam_NaCl)
}
