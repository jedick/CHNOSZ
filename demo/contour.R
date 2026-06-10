# CHNOSZ/demo/contour.R
# Gold solubility contours on logfO2-pH diagram
# 20181107 initial version
# 20190415 cleanup for demo
# - Diagram based on Williams-Jones et al. (2009) doi:10.2113/gselements.5.5.281
# 20250215 adjust conditions to reproduce plot from Ding et al.
# - After Ding et al. (2023) doi:10.1016/j.oregeorev.2022.105231
# 20260610 use IGEM data
# - Check diagram against Rubtsova et al. (2026) doi:10.1134/S1075701525600690
#   (especially narrow Au2S2-2 field at high pH)

library(CHNOSZ)

# Define plot resolution
res <- 300
# Define temperature (degrees C), pressure (bar)
T <- 250
P <- 300
# Define molality of NaCl and total S
m_NaCl <- 1
sum_S <- 0.01
# Define ranges of pH and logfO2
pH <- c(1, 10)
O2 <- c(-40, -30)
# Make smooth (TRUE) or sharp (FALSE) transitions between basis species 
blend <- TRUE

# Set up system
basis(c("Au", "Cl-", "H2S", "H2O", "oxygen", "H+"))
iaq <- retrieve("Au", c("H", "O", "S", "Cl"), state = "aq")
species(iaq)
basis("H2S", log10(sum_S))
# Calculate solution composition for given NaCl molality
NaCl <- NaCl(m_NaCl = m_NaCl, T = T, P = P)
basis("Cl-", log10(NaCl$m_Clminus))
# Calculate affinity with changing basis species across diagram
bases <- c("H2S", "HS-", "S3-", "SO2", "HSO4-", "SO4-2")
m <- mosaic(bases, pH = c(pH, res), O2 = c(O2, res), T = T, P = P, IS = NaCl$IS, blend = blend)
# Show predominance fields for S-species
diagram(m$A.bases, col = 8, col.names = 8, lty = 3, italic = TRUE)
# Show predominance fields for Au-species
# NOTE: balance = 1 is used to correctly handle *aqueous* multinuclear species ()
d <- diagram(m$A.species, col = 4, col.names = 4, lty = 2, bold = TRUE, add = TRUE, balance = 1)
# Add dot-dash line for water stability limit
water.lines(d, lty = 4)

# Calculate and plot solubility of Au (use named 'bases' argument to trigger mosaic calculation)
species("Au")
s <- solubility(iaq, bases = bases, pH = c(pH, res), O2 = c(O2, res), T = T, P = P, IS = NaCl$IS, blend = blend)
# Convert to ppb
s <- convert(s, "ppb")
diagram(s, levels = c(1, 10, 100, 1000), col = 1, add = TRUE)
# Add legend and title
dP <- describe.property(c("T", "P"), c(T, P))
lexpr <- c(dP, lNaCl(m_NaCl), lS(sum_S))
legend("topright", legend = lexpr, bty = "n")
title(main = ("Au solubility (ppb), after Ding et al. (2023) and Rubtsova et al. (2026)"), font.main = 1)
