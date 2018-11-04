# NaCl.R
# calculate ion activities and ionic strength
# given a total molality of NaCl
# taking account of ion association: Na+ + Cl- = NaCl(aq)
# 20181102 jmd first version

NaCl <- function(T=seq(100, 500, 100), P=1000, m_tot=2) {
  # define a function for the reaction quotient
  logQ <- function(m_Cl, gamma) {
    # starting with Q = a_NaCl / (a_Na+ * a_Cl-),
    # substitute gam_NaCl = 0, m_NaCl + m_Cl = m_tot, m_Cl = m_Na, gam_Cl = gam_Na = gamma
    # to write:
    log10( (m_tot - m_Cl) / (m_Cl * gamma) ^ 2 )
  }
  # define a function for affinity = log(K / Q)
  A <- function(m_Cl, gamma, logK) logK - logQ(m_Cl, gamma)
  # calculate equilibrium constant at all temperatures (standard conditions: IS = 0)
  logK <- subcrt(c("Na+", "Cl-", "NaCl"), c(-1, -1, 1), T = T, P = P)$out$logK
  # calculate Debye-Huckel parameters at all temperatures
  wout <- water(c("A_DH", "B_DH"), T = convert(T, "K"), P = P)
  # initialize output variables
  N <- length(T)
  ISout <- a_Cl <- numeric(N)
  # initial guess for m_Cl and ionic strength assuming complete dissociation of NaCl
  IS <- m_Cl <- rep(m_tot, N)
  # the species index for Cl-
  iCl <- info("Cl-")
  # we start by doing calculations for all temperatures
  doit <- !logical(N)
  while(any(doit)) {
    # calculate activity coefficient at given ionic strength
    gamma <- suppressMessages(10^nonideal(iCl, list(data.frame(G=numeric(N))), IS=IS, T=convert(T, "K"), P=P, A_DH=wout$A_DH, B_DH=wout$B_DH)[[1]]$loggam)
    # solve for m_Cl
    for(i in which(doit)) m_Cl[i] <- uniroot(A, c(0, m_tot), gamma=gamma[i], logK=logK[i])$root
    # calculate new ionic strength and deviation
    ISnew <- m_Cl
    dIS <- ISnew - IS
    # set net ionic strength
    IS <- ISnew
    # keep going until the deviation in ionic strength at any temperature is less than 0.01
    doit <- abs(dIS) > 0.01
  }
  # assemble final gamma and activity of Cl-
  gamma <- suppressMessages(10^nonideal(iCl, list(data.frame(G=numeric(N))), IS=IS, T=convert(T, "K"), P=P, A_DH=wout$A_DH, B_DH=wout$B_DH)[[1]]$loggam)
  a_Cl <- m_Cl * gamma
  # return the calculated values
  list(IS=IS, gamma=gamma, a_Cl=a_Cl)
}
