# Solubility and speciation of Au with polynuclear species (Au2S2-2)
# After Tagirov et al. (2025) doi:10.1016/j.gca.2024.08.022
# 20260611 initial version

library(CHNOSZ)

par(mfrow = c(2, 1))

# Define temperature (degrees C), pressure (bar)
T <- 300
P <- 500
# Define molality of NaCl and total S
# NOTE: for simplicity m_NaCl is used for Cl- and ionic strength;
# see other demos for using NaCl() to calculate these more accurately
m_NaCl <- 1
S_tot <- 0.05
# Define pH range
pH <- c(1, 12)

# Define O2 buffer
mod.buffer("NNO", c("nickel", "bunsenite"), state = "cr", logact = 0)
O2_buffer <- "NNO"
dlogfO2 <- 1.5

# Set up system
# NOTE: Au must be first otherwise we get this cryptic error message for solubility():
# Error in basis(ispecies, loga.numeric) : species names are not unique
basis(c("Au", "Ni", "Cl-", "H2S", "H2O", "oxygen", "H+"))
basis("H2S", log10(S_tot))
basis("Cl-", log10(m_NaCl))

# Define mineral and solutes
species("Au")
iaq <- retrieve("Au", c("H", "O", "S", "Cl"), state = "aq")

# Calculate logfO2 of NNO buffer at T
basis("O2", O2_buffer)
logfO2 <- affinity(return.buffer = TRUE, T = T, P = P)$O2
basis("O2", logfO2 + dlogfO2)

### Solubility calculation

# The candidate S-bearing basis species
S_species <- c("H2S", "HS-", "S3-", "SO2", "HSO4-", "SO4-2")
s <- solubility(iaq, bases = S_species, pH = pH, T = T, P = P, IS = m_NaCl)
# Plot partial solubility then total solubility
d_Au <- diagram(s, type="loga.equil", ylim = c(-10, -3), dy = 0.4)
diagram(s, type="loga.balance", col = 4, lwd = 2, add = TRUE)
title("Use solubility() to reproduce/extend Fig. 7 of Tagirov et al. (2025)", font.main = 1)
# Add T,P legend
legend("topleft", legend = lTP(T, P), bty = "n")
# Add composition legend
l_NaCl <- bquote(italic(m)*NaCl == .(m_NaCl))
l_S <- bquote(italic(m)*S[tot] == .(S_tot))
l_O2 <- bquote(log~italic(f)*O[2] == .(O2_buffer) + .(dlogfO2))
legend("topright", legend = c(l_NaCl, l_S, l_O2), bty = "n")

# Extract the pH range where Au2S2-2 predominates
species_index <- which(names(iaq) == "Au2S2-2")
species_predominant <- d_Au$predominant == species_index
species_pH <- d_Au$vals$pH[species_predominant]
info <- "CORRECT: solubility() has a small pH range for Au2S2-2 predominance"
# Becuase the demo isn't run by a test runner, wrap the test in stopifnot(isTRUE(...)))
stopifnot(isTRUE(expect_equal(round(range(species_pH), 1), c(9.2, 9.9), info = info)))

### Speciation calculation

# Define the total activity of Au
Au_tot <- -4
# Load the Au-bearing species (no Au(cr))
species(iaq)
# Calculate affinities over mosaic of S-bearing basis species
m <- mosaic(S_species, pH = pH, T = T, P = P, IS = m_NaCl)

# Distribute Au species with fixed total activity
e1 <- equilibrate(m$A.species, balance = 1, loga.balance = Au_tot)
d1 <- diagram(e1, ylim = c(-8, -3))
title("Use equilibrate() with balance=1 for speciation at fixed total Au", line=1.5, font.main = 1)
subtext <- bquote("NOTE: balance=NULL (the default) overestimates"~.(expr.species("Au2S2-2"))~"predominance")
mtext(subtext, line = 0, cex = par("cex"))
species_predominant <- d1$predominant == species_index
species_pH <- d1$vals$pH[species_predominant]
info <- "CORRECT: equilibrate() with balance=1 has a small pH range for Au2S2-2 predominance"
stopifnot(isTRUE(expect_equal(round(range(species_pH), 1), c(9.2, 9.9), info = info)))

### Incorrect speciation calculation (not plotted)

e0 <- equilibrate(m$A.species, balance = NULL, loga.balance = Au_tot)
d0 <- diagram(e0, ylim = c(-8, -3), plot.it = FALSE)
species_predominant <- d0$predominant == species_index
species_pH <- d0$vals$pH[species_predominant]
info <- "INCORRECT: equilibrate() with balance=NULL overestimates Au2S2-2 predominance"
stopifnot(isTRUE(expect_equal(round(range(species_pH), 1), c(8.2, 11.8), info = info)))

