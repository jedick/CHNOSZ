# CHNOSZ/demo/gold.R: Au solubility calculations
# 20181101 jmd first version

## additions to OBIGT:
# Au(HS) from Akinfiev and Zotov, 2010
# (doi:10.1134/S0016702910070074)
# corrected H taken from Pokrovski et al., 2014
# (doi:10.1144/SP402.4)
mod.obigt("Au(HS)", formula = "Au(HS)", state = "aq", ref1 = "AZ10", date = today(),
  G = 8344, H = 13193, S = 50.86, Cp = 1.8, V = 56.5,
  a1 = 9.4965, a2 = 15.4057, a3 = -0.3052, a4 = -3.1459,
  c1 = -38.1356, c2 = 19.6524, omega = 0, z = 0)
# AuOH from Pokrovski et al., 2014
mod.obigt("AuOH", formula = "AuOH", state = "aq", ref1 = "PAB+14", date = today(),
  G = -32716, H = -41533, S = 21.89, Cp = -11.1, V = 32.4,
  a1 = 6.1937, a2 = 7.3415, a3 = 2.8644, a4 = -3.0825,
  c1 = -3.0370, c2 = -3.9635, omega = 0, z = 0)

## modifications to OBIGT:
# AuCl2- from Akinfiev and Zotov, 2001 (reported in AZ10)
# (http://pleiades.online/cgi-perl/search.pl/?type=abstract&name=geochem&number=10&year=1&page=990)
mod.obigt("AuCl2-", formula = "AuCl2-", state = "aq", ref1 = "AZ01", date = today(),
  G = -36795, H = -46664, S = 47.16, Cp = -26.4, V = 68.6,
  a1 = 11.4774, a2 = 20.2425, a3 = -2.2063, a4 = -3.6158,
  c1 = 27.0677, c2 = -22.240, omega = 0.8623, z = -1)
# Au(HS)2- from Pokrovski et al., 2014
mod.obigt("Au(HS)2-", G = 3487, H = 4703, S = 77.46, Cp = 3.3, V = 75.1,
  a1 = 12.3373, a2 = 22.3421, a3 = 3.0317, a4 = -3.7026,
  c1 = -53.6010, c2 = 31.4030, omega = 0.7673, z = -1)

# set up system
# use H2S here: it's the predominant species at the pH of the QMK buffer -- see sulfur()
basis(c("Al2O3", "quartz", "Fe", "Au", "K+", "Cl-", "H2S", "H2O", "oxygen", "H+"))
# set molality of K+ in completely dissociated 0.5 molal KCl
# NOTE: This value is used only for making the legend;
# activities corrected for ionic strength are computed below
basis("K+", log10(0.5))

# create a pH buffer
mod.buffer("QMK", c("quartz", "muscovite", "K-feldspar"), "cr", 0)

# define colors for Au(HS)2-, Au(HS), AuOH, AuCl2-
# after Williams-Jones et al., 2009
# (doi:10.2113/gselements.5.5.281)
col <- c("#ED4037", "#F58645", "#0F9DE2", "#22CC88")

# sulfur logfO2-pH diagrams showing redox and pH buffers at four temperatures 20181031
sulfur <- function() {
  species(delete = TRUE)
  species(c("H2S", "HS-", "HSO4-", "SO4-2"))
  T <- c(200, 300, 400, 500)
  P <- 1000
  O2min <- c(-50, -40, -30, -25)
  O2max <- c(-30, -20, -20, -10)
  par(mfrow=c(2, 2))
  for(i in 1:4) {
    a <- affinity(pH = c(0, 14), O2 = c(O2min[i], O2max[i]), T = T[i], P = 1000)
    diagram(a)
    basis("H+", "QMK")
    pH_QMK <- -affinity(T = T[i], P = P, return.buffer = TRUE)$`H+`
    abline(v = pH_QMK, lty = 2)
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
  species(c("Au(HS)2-", "Au(HS)", "AuOH"))
  # apply PPM buffer for fO2 and aH2S
  basis("O2", "PPM")
  basis("H2S", "PPM")
  # calculate affinity and solubility
  # (set IS = 0 for diagram to show "log m" instead of "log a")
  a <- affinity(pH = c(3, 8), T = 300, P = 1000, IS = 0)
  s <- solubility(a)
  # make diagram and show total log molality
  diagram(s, ylim = c(-10, -5), col = col, lwd = 2, lty = 1)
  diagram(s, add = TRUE, type = "loga.balance", lwd = 3, lty = 2)
  # add neutral pH line
  pH <- -subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = 300, P = 1000)$out$logK/2
  abline(v = pH, lty = 3)
  # make legend and title
  dprop <- describe.property(c("T", "P", "IS"), c(300, 1000, 0))
  legend("topleft", dprop, bty = "n")
  dbasis <- describe.basis(ibasis = c(9, 7))
  legend("bottomright", dbasis, bty = "n")
  title(main=("After Akinfiev and Zotov, 2001, Fig. 7"), font.main = 1)
}

# log(m_Au)-pH diagram similar to Fig. 12b of Stefansson and Seward, 2004
# (doi:10.1016/j.gca.2004.04.006)
Au_pH2 <- function() {
  species(c("Au(HS)2-", "Au(HS)", "AuOH", "AuCl2-"))
  # apply PPM buffer for fO2 and aH2S
  basis("O2", "PPM")
  basis("H2S", "PPM")
  # calculate affinity and solubility
  # (set IS = 0 for diagram to show "log m" instead of "log a")
  a <- affinity(pH = c(3, 8), T = 450, P = 1000, IS = 0)
  s <- solubility(a)
  # make diagram and show total log molality
  diagram(s, ylim = c(-8, -3), col = col, lwd = 2, lty = 1)
  diagram(s, add = TRUE, type = "loga.balance", lwd = 3, lty = 2)
  # add neutral pH line
  pH <- -subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = 450, P = 1000)$out$logK/2
  abline(v = pH, lty = 3)
  # make legend and title
  dprop <- describe.property(c("T", "P", "IS"), c(450, 1000, 0))
  legend("topleft", dprop, bty = "n")
  dbasis <- describe.basis(ibasis = c(6, 9, 7))
  legend("topright", dbasis, bty = "n")
  title(main=("After Stef\u00e1nsson and Seward, 2004, Fig. 12b"), font.main = 1, cex.main = 1.1)
}

# log(m_Au)-T diagram like Fig. 2B of Williams-Jones et al., 2009
# (doi:10.2113/gselements.5.5.281)
Au_T1 <- function() {
  species(c("Au(HS)2-", "Au(HS)", "AuOH", "AuCl2-"))
  # apply PPM buffer for fO2 and aH2S
  basis("O2", "PPM")
  basis("H2S", "PPM")
  # apply QMK buffer for pH
  basis("H+", "QMK")
  # calculate solution composition for 2 mol/kg NaCl
  NaCl <- NaCl(T = seq(150, 550, 10), P = 1000, m_tot=2)
  a_Cl <- NaCl$m_Cl * NaCl$gam_Cl
  # using this ionic strength, calculate the activity of K+
  # assuming complete dissociation of 0.5 mol/kg KCl
  gam_K <- 10^subcrt("K+", T = seq(150, 550, 10), P = 1000, IS=NaCl$IS)$out$`K+`$loggam
  a_K <- 0.5 * gam_K
  # calculate affinity and solubility
  a <- affinity(T = seq(150, 550, 10), `Cl-` = log10(a_Cl), `K+` = log10(a_K), P = 1000, IS = NaCl$IS)
  s <- solubility(a)
  # make diagram and show total log molality
  diagram(s, ylim = c(-10, -4), col = col, lwd = 2, lty = 1)
  diagram(s, add = TRUE, type = "loga.balance", lwd = 3, lty = 2)
  # make legend and title
  dP <- describe.property("P", 1000)
  dNaCl <- expression(NaCl == 2~mol~kg^-1)
  dK <- describe.basis(ibasis=5, use.molality=TRUE)
  legend("topleft", c(dP, dNaCl, dK), bty = "n")
  dbasis <- describe.basis(ibasis = c(9, 7, 10))
  legend("topright", dbasis, bty = "n")
  title(main=("After Williams-Jones et al., 2009, Fig. 2B"), font.main = 1)
}

# log(m_Au)-T diagram like Fig. 2A of Williams-Jones et al., 2009 and Fig. 8a of Pokrovski et al., 2014
# (doi:10.2113/gselements.5.5.281)
# (doi:10.1144/SP402.4)
Au_T2 <- function() {
  species(c("Au(HS)2-", "Au(HS)", "AuOH", "AuCl2-"))
  # approximate activity of H2S for total S = 0.01 m
  basis("H2S", -2)
  # apply HM buffer for fO2
  basis("O2", "HM")
  # apply QMK buffer for pH
  basis("H+", "QMK")
  # calculate solution composition for 2 mol/kg NaCl
  NaCl <- NaCl(T = seq(150, 550, 10), P = 1000, m_tot=2)
  a_Cl <- NaCl$m_Cl * NaCl$gam_Cl
  # using this ionic strength, calculate the activity of K+
  # assuming complete dissociation of 0.5 mol/kg KCl
  gam_K <- 10^subcrt("K+", T = seq(150, 550, 10), P = 1000, IS=NaCl$IS)$out$`K+`$loggam
  a_K <- 0.5 * gam_K
  # calculate affinity and solubility
  a <- affinity(T = seq(150, 550, 10), `Cl-` = log10(a_Cl), `K+` = log10(a_K), P = 1000, IS = NaCl$IS)
  s <- solubility(a)
  # make diagram and show total log molality
  diagram(s, ylim = c(-10, -2), col = col, lwd = 2, lty = 1)
  diagram(s, add = TRUE, type = "loga.balance", lwd = 3, lty = 2)
  # make legend and title
  dP <- describe.property("P", 1000)
  dNaCl <- expression(NaCl == 2~mol~kg^-1)
  dK <- describe.basis(ibasis=5, use.molality=TRUE)
  legend("topleft", c(dP, dNaCl, dK), bty = "n")
  dbasis <- describe.basis(ibasis = c(9, 7, 10))
  legend("topright", dbasis, bty = "n")
  title(main=("After Williams-Jones et al., 2009, Fig. 2A"), font.main = 1)
}

# make plots
par(mfrow=c(2, 2))
Au_pH1()
Au_pH2()
Au_T1()
Au_T2()