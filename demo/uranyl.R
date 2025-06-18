# CHNOSZ/demo/uranyl.R
# Total (carbonate|sulfate)-pH diagrams for uranyl species, after Migdisov et al., 2024
# 20241116 first version

library(CHNOSZ)

# Conditions
logm_U <- log10(3.16e-5)
m_NaCl <- 1  # mol NaCl / kg H2O
T <- 200
P <- "Psat"
pH_lim <- c(2, 10)
CS_lim <- c(-4, 1)
res <- 500

# Calculations for NaCl
NaCl <- NaCl(m_NaCl = m_NaCl, T = T, P = P)
IS <- NaCl$IS
logm_Naplus <- log10(NaCl$m_Naplus)
logm_Clminus <- log10(NaCl$m_Clminus)

# Total carbonate-pH
iaq <- retrieve("U", ligands = c("C", "O", "H", "Cl", "Na"), state = "aq")
# Only get minerals that have non-NA G at the indicated temperature
icr <- retrieve("U", ligands = c("C", "O", "H", "Cl", "Na"), state = "cr", T = T)
basis(c("UO2+2", "CO3-2", "Na+", "Cl-", "H+", "H2O", "O2"))
basis(c("Na+", "Cl-"), c(logm_Naplus, logm_Clminus))
species(iaq, logm_U)
species(icr, add = TRUE)
bases <- c("CO3-2", "HCO3-", "CO2")
# Suppress warnings about being above T limit of Cp equation for beta-UO2(OH)2, UO2.25, and beta-UO2.3333
m <- suppressWarnings(mosaic(bases, pH = c(pH_lim, res), "CO3-2" = c(CS_lim, res), T = T, P = P, IS = IS))
diagram(m$A.species)
diagram(m$A.bases, add = TRUE, col = 8, lty = 2, col.names = 8, italic = TRUE)
title("Uranyl-carbonate complexation at 200 \u00b0C, after Migdisov et al., 2024", font.main = 1)

# Total sulfate-pH
iaq <- retrieve("U", ligands = c("S", "O", "H", "Cl", "Na"), state = "aq")
icr <- retrieve("U", ligands = c("S", "O", "H", "Cl", "Na"), state = "cr", T = T)
basis(c("UO2+2", "SO4-2", "Na+", "Cl-", "H+", "H2O", "O2"))
basis(c("Na+", "Cl-"), c(logm_Naplus, logm_Clminus))
species(iaq, logm_U)
species(icr, add = TRUE)
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
m <- suppressWarnings(mosaic(bases, pH = c(pH_lim, res), "SO4-2" = c(CS_lim, res), T = T, P = P, IS = IS))
diagram(m$A.species)
diagram(m$A.bases, add = TRUE, col = 8, lty = 2, col.names = 8, italic = TRUE)
title("Uranyl-sulfate complexation at 200 \u00b0C, after Migdisov et al., 2024", font.main = 1)
