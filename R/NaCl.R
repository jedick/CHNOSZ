# NaCl.R
# calculate ion activities and ionic strength
# given a total molality of NaCl
# taking account of ion association: Na+ + Cl- = NaCl(aq)
# 20181102 jmd first version
# 20181105 use activity coefficient of Na+
# 20181106 use activity coefficient of NaCl

NaCl <- function(T=seq(100, 500, 100), P=1000, m_tot=2, ...) {
  # define a function for the reaction quotient
  logQ <- function(m_Cl, gam_NaCl, gam_Na, gam_Cl) {
    # starting with Q = a_NaCl / (a_Na+ * a_Cl-),
    # substitute m_tot = m_NaCl + m_Cl and m_Cl = m_Na
    # to write:
    log10( (m_tot - m_Cl) * gam_NaCl / (m_Cl ^ 2 * gam_Na * gam_Cl) )
  }
  # define a function for affinity = log(K / Q)
  A <- function(m_Cl, gam_NaCl, gam_Na, gam_Cl, logK) logK - logQ(m_Cl, gam_NaCl, gam_Na, gam_Cl)
  # calculate equilibrium constant at all temperatures (standard conditions: IS = 0)
  logK <- subcrt(c("Na+", "Cl-", "NaCl"), c(-1, -1, 1), T = T, P = P, ...)$out$logK
  # calculate Debye-Huckel parameters at all temperatures
  wout <- water(c("A_DH", "B_DH"), T = convert(T, "K"), P = P)
  # initialize output variables
  N <- length(T)
  ISout <- a_Cl <- numeric(N)
  # initial guess for m_Cl and ionic strength assuming complete dissociation of NaCl
  IS <- m_Cl <- rep(m_tot, N)
  # the corresponding total molality of dissolved species (NaCl + Cl- + Na+)
  m_star <- (m_tot - m_Cl) + 2*m_Cl
  # the species indices for Na+, Cl-, and NaCl(aq)
  ispecies <- info(c("Na+", "Cl-", "NaCl"))
  # we start by doing calculations for all temperatures
  doit <- !logical(N)
  while(any(doit)) {
    # calculate activity coefficient at given ionic strength
    speciesprops <- rep(list(data.frame(G=numeric(N))), length(ispecies))
    gammas <- suppressMessages(nonideal(ispecies, speciesprops, IS=IS, T=convert(T, "K"), P=P, A_DH=wout$A_DH, B_DH=wout$B_DH, m_star=m_star))
    gam_Na <- 10^gammas[[1]]$loggam
    gam_Cl <- 10^gammas[[2]]$loggam
    gam_NaCl <- 10^gammas[[3]]$loggam
    # solve for m_Cl
    for(i in which(doit)) m_Cl[i] <- uniroot(A, c(0, m_tot), gam_NaCl=gam_NaCl[i], gam_Na=gam_Na[i], gam_Cl=gam_Cl[i], logK=logK[i])$root
    # calculate new total molality
    m_star <- (m_tot - m_Cl) + 2*m_Cl
    # calculate new ionic strength and deviation
    ISnew <- m_Cl
    dIS <- ISnew - IS
    # set new ionic strength
    IS <- ISnew
    # keep going until the deviation in ionic strength at any temperature is less than 0.01
    doit <- abs(dIS) > 0.01
  }
  # assemble final molality of Cl- and gammas
  gammas <- suppressMessages(nonideal(ispecies, speciesprops, IS=IS, T=convert(T, "K"), P=P, A_DH=wout$A_DH, B_DH=wout$B_DH))
  gam_Na <- 10^gammas[[1]]$loggam
  gam_Cl <- 10^gammas[[2]]$loggam
  # return the calculated values
  list(IS=IS, m_Cl=m_Cl, gam_Na=gam_Na, gam_Cl=gam_Cl, gam_NaCl=gam_NaCl)
}
