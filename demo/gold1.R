# CHNOSZ/demo/gold.R: Au solubility calculations
# 20181101 jmd first version
# 20181109 add calculation of K+ molality
# 20190127 update Au species in OBIGT, not here

library(CHNOSZ)

# Set up system
# Use H2S here: it's the predominant species at the pH of the KMQ buffer -- see sulfur()
basis(c("Au", "Al2O3", "quartz", "Fe", "K+", "Cl-", "H2S", "H2O", "oxygen", "H+"))
# set molality of K+ in completely dissociated 0.5 molal KCl
# NOTE: This value is used only for making the legend;
# activities corrected for ionic strength are computed below
basis("K+", log10(0.5))

# Create a pH buffer
mod.buffer("KMQ", c("quartz", "muscovite", "K-feldspar"), "cr", 0)

# Define colors for Au(HS)2-, AuHS, AuOH, AuCl2-
# after Williams-Jones et al., 2009
# (doi:10.2113/gselements.5.5.281)
col <- c("#ED4037", "#F58645", "#0F9DE2", "#22CC88")

# Sulfur logfO2-pH diagrams showing redox and pH buffers at four temperatures 20181031
sulfur <- function() {
  species(c("H2S", "HS-", "S3-", "SO2", "HSO4-", "SO4-2"))
  T <- c(200, 300, 400, 500)
  P <- 1000
  O2min <- c(-50, -40, -30, -25)
  O2max <- c(-30, -20, -20, -10)
  par(mfrow=c(2, 2))
  for(i in 1:4) {
    a <- affinity(pH = c(0, 14), O2 = c(O2min[i], O2max[i]), T = T[i], P = 1000)
    diagram(a)
    basis("H+", "KMQ")
    pH_KMQ <- -affinity(T = T[i], P = P, return.buffer = TRUE)$`H+`
    abline(v = pH_KMQ, lty = 2)
    basis("O2", "HM")
    O2_HM <- affinity(T = T[i], P = P, return.buffer = TRUE)$O2
    abline(h = O2_HM, lty = 2, col = "blue")
    text(12, O2_HM, "HM", adj = c(0, -0.5), col = "blue")
    basis("O2", "PPM")
    O2_PPM <- affinity(T = T[i], P = P, return.buffer = TRUE)$O2
    abline(h = O2_PPM, lty = 2, col = "blue")
    text(12, O2_PPM, "PPM", adj = c(0, -0.5), col = "blue")
    # remove the buffers for next plot
    basis("O2", 0)
    basis("pH", 0)
  }
}

# log(m_Au)-pH diagram like Fig. 7 of Akinfiev and Zotov, 2001
# (http://pleiades.online/cgi-perl/search.pl/?type=abstract&name=geochem&number=10&year=1&page=990)
Au_pH1 <- function() {
  # Apply PPM buffer for fO2 and aH2S
  basis("O2", "PPM")
  basis("H2S", "PPM")
  # Calculate solubility of gold
  species("Au")
  iaq <- info(c("Au(HS)2-", "AuHS", "AuOH"))
  # (set IS = 0 for diagram to show "log m" instead of "log a")
  s <- solubility(iaq, pH = c(3, 8), T = 300, P = 1000, IS = 0)
  # Make diagram and show total log molality
  diagram(s, type = "loga.equil", ylim = c(-10, -5), col = col, lwd = 2, lty = 1)
  diagram(s, add = TRUE, lwd = 3, lty = 2)
  # Add neutral pH line
  pH <- -subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = 300, P = 1000)$out$logK/2
  abline(v = pH, lty = 3)
  # Make legend and title
  dprop <- describe.property(c("T", "P", "IS"), c(300, 1000, 0))
  legend("topleft", dprop, bty = "n")
  dbasis <- describe.basis(c(9, 7))
  legend("bottomright", dbasis, bty = "n")
  title(main=("After Akinfiev and Zotov, 2001, Fig. 7"), font.main = 1)
}

# log(m_Au)-pH diagram similar to Fig. 12b of Stefansson and Seward, 2004
# (doi:10.1016/j.gca.2004.04.006)
Au_pH2 <- function() {
  # Apply PPM buffer for fO2 and aH2S
  basis("O2", "PPM")
  basis("H2S", "PPM")
  # Calculate solubility of gold
  # (set IS = 0 for diagram to show "log m" instead of "log a")
  species("Au")
  iaq <- info(c("Au(HS)2-", "AuHS", "AuOH", "AuCl2-"))
  s <- solubility(iaq, pH = c(3, 8), T = 450, P = 1000, IS = 0)
  # Make diagram and show total log molality
  diagram(s, type = "loga.equil", ylim = c(-8, -3), col = col, lwd = 2, lty = 1)
  diagram(s, add = TRUE, lwd = 3, lty = 2)
  # Add neutral pH line
  pH <- -subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = 450, P = 1000)$out$logK/2
  abline(v = pH, lty = 3)
  # Make legend and title
  dprop <- describe.property(c("T", "P", "IS"), c(450, 1000, 0))
  legend("topleft", dprop, bty = "n")
  dbasis <- describe.basis(c(6, 9, 7))
  legend("topright", dbasis, bty = "n")
  title(main=("After Stef\u00e1nsson and Seward, 2004, Fig. 12b"), font.main = 1, cex.main = 1.1)
}

# Estimate the Cl- molality and ionic strength for a hypothetical 
# NaCl solution with total chloride equal to specified NaCl + KCl solution,
# then estimate the molality of K+ in that solution 20181109
chloride <- function(T, P, m_NaCl, m_KCl) {
  NaCl <- NaCl(m_NaCl = m_NaCl + m_KCl, T = T, P = P)
  # Calculate logK of K+ + Cl- = KCl, adjusted for ionic strength
  logKadj <- subcrt(c("K+", "Cl-", "KCl"), c(-1, -1, 1), T = T, P = P, IS = NaCl$IS)$out$logK
  # What is the molality of K+ from 0.5 mol KCl in solution with 2 mol total Cl
  m_Kplus <- m_KCl / (10^logKadj * NaCl$m_Clminus + 1)
  list(IS = NaCl$IS, m_Clminus = NaCl$m_Clminus, m_Kplus = m_Kplus)
}

# log(m_Au)-T diagram like Fig. 2B of Williams-Jones et al., 2009
# (doi:10.2113/gselements.5.5.281)
Au_T1 <- function() {
  # Apply PPM buffer for fO2 and aH2S
  basis("O2", "PPM")
  basis("H2S", "PPM")
  # Apply KMQ buffer for pH
  basis("H+", "KMQ")
  # Estimate solution composition for 1.5 m NaCl and 0.5 m KCl
  chl <- chloride(T = seq(150, 550, 10), P = 1000, m_NaCl = 1.5, m_KCl = 0.5)
  # Calculate solubility of gold
  species("Au")
  iaq <- info(c("Au(HS)2-", "AuHS", "AuOH", "AuCl2-"))
  s <- solubility(iaq, T = seq(150, 550, 10), `Cl-` = log10(chl$m_Clminus), `K+` = log10(chl$m_Kplus), P = 1000, IS = chl$IS)
  # Make diagram and show total log molality
  diagram(s, type = "loga.equil", ylim = c(-10, -3), col = col, lwd = 2, lty = 1)
  diagram(s, add = TRUE, lwd = 3, lty = 2)
  # Make legend and title
  dP <- describe.property("P", 1000)
  dNaCl <- expression(italic(m)[NaCl] == 1.5)
  dKCl <- expression(italic(m)[KCl] == 0.5)
  legend("topleft", c(dP, dNaCl, dKCl), bty = "n")
  dH2S <- describe.basis(7, molality=TRUE)
  dO2 <- describe.basis(9)
  dpH <- describe.basis(10)
  legend(300, -3, c(dH2S, dO2, dpH), bty = "n")
  title(main=("After Williams-Jones et al., 2009, Fig. 2B"), font.main = 1)
}

# log(m_Au)-T diagram like Fig. 2A of Williams-Jones et al., 2009 and Fig. 8a of Pokrovski et al., 2014
# (doi:10.2113/gselements.5.5.281)
# (doi:10.1144/SP402.4)
Au_T2 <- function() {
  # Total S = 0.01 m
  basis("H2S", -2)
  # Apply HM buffer for fO2
  basis("O2", "HM")
  # Apply KMQ buffer for pH
  basis("H+", "KMQ")
  # Estimate solution composition for 1.5 m NaCl and 0.5 m KCl
  chl <- chloride(T = seq(150, 550, 10), P = 1000, m_NaCl = 1.5, m_KCl = 0.5)
  # Calculate solubility of gold
  species("Au")
  iaq <- info(c("Au(HS)2-", "AuHS", "AuOH", "AuCl2-"))
  s <- solubility(iaq, T = seq(150, 550, 10), `Cl-` = log10(chl$m_Clminus), `K+` = log10(chl$m_Kplus), P = 1000, IS = chl$IS)
#  # Uncomment to calculate solubility considering speciation of sulfur
#  bases <- c("H2S", "HS-", "S3-", "SO2", "HSO4-", "SO4-2")
#  s <- solubility(iaq, bases = bases, T = seq(150, 550, 10), `Cl-` = log10(chl$m_Clminus), `K+` = log10(chl$m_Kplus), P = 1000, IS = chl$IS)
  # Make diagram and show total log molality
  diagram(s, type = "loga.equil", ylim = c(-10, -3), col = col, lwd = 2, lty = 1)
  diagram(s, add = TRUE, lwd = 3, lty = 2)
  # Make legend and title
  dP <- describe.property("P", 1000)
  dNaCl <- expression(italic(m)[NaCl] == 1.5)
  dKCl <- expression(italic(m)[KCl] == 0.5)
  legend("topleft", c(dP, dNaCl, dKCl), bty = "n")
  dH2S <- expr.species("H2S", value = 0.01, molality = TRUE)
  dO2 <- describe.basis(9)
  dpH <- describe.basis(10)
  legend(300, -3, c(dH2S, dO2, dpH), bty = "n")
  title(main=("After Williams-Jones et al., 2009, Fig. 2A"), font.main = 1)
}

# Make plots
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
Au_pH1()
Au_pH2()
Au_T1()
Au_T2()
par(opar)
