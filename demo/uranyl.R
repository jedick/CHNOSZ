# CHNOSZ/demo/uranyl.R
# Total (carbonate|sulfate)-pH diagrams for uranyl species, after Migdisov et al., 2024
# 20241116 jmd
library(CHNOSZ)

# Conditions
logm_U <- log10(3.16e-5)
m_tot <- 1  # mol NaCl / kg H2O
T <- 200
P <- "Psat"
pH_lim <- c(2, 10)
CS_lim <- c(-4, 1)
res <- 500

# Calculations for NaCl
NaCl <- NaCl(m_tot = m_tot, T = T, P = P)
IS <- NaCl$IS
logm_Na <- log10(NaCl$m_Na)
logm_Cl <- log10(NaCl$m_Cl)

# Total carbonate-pH
iaq <- retrieve("U", ligands = c("C", "O", "H", "Cl", "Na"), state = "aq")
icr <- retrieve("U", ligands = c("C", "O", "H", "Cl", "Na"), state = "cr")
basis(c("UO2+2", "CO3-2", "Na+", "Cl-", "H+", "H2O", "O2"))
basis(c("Na+", "Cl-"), c(logm_Na, logm_Cl))
species(iaq, logm_U)
species(icr, add = TRUE)
bases <- c("CO3-2", "HCO3-", "CO2")
m <- mosaic(bases, pH = c(pH_lim, res), "CO3-2" = c(CS_lim, res), T = T, P = P, IS = IS)
diagram(m$A.species)
diagram(m$A.bases, add = TRUE, col = 2, lty = 2, col.names = 2)
title("Uranyl-carbonate complexation at 200 \u00b0C, after Migdisov et al., 2024", font.main = 1)

# Total sulfate-pH
iaq <- retrieve("U", ligands = c("S", "O", "H", "Cl", "Na"), state = "aq")
icr <- retrieve("U", ligands = c("S", "O", "H", "Cl", "Na"), state = "cr")
basis(c("UO2+2", "SO4-2", "Na+", "Cl-", "H+", "H2O", "O2"))
basis(c("Na+", "Cl-"), c(logm_Na, logm_Cl))
species(iaq, logm_U)
species(icr, add = TRUE)
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
m <- mosaic(bases, pH = c(pH_lim, res), "SO4-2" = c(CS_lim, res), T = T, P = P, IS = IS)
diagram(m$A.species)
diagram(m$A.bases, add = TRUE, col = 2, lty = 2, col.names = 2)
title("Uranyl-sulfate complexation at 200 \u00b0C, after Migdisov et al., 2024", font.main = 1)
