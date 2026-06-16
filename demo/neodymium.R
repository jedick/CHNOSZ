# Solubility of multiple Nd minerals
# After Migdisov et al. (2016) doi:10.1016/j.chemgeo.2016.06.005
# 20260614 first version

library(CHNOSZ)

# Setup database
reset()
# Use Haas et al. (1995) data for REE complexes
add.OBIGT("SLOP98b")

# Fit logK for NdSO4+ and Nd(SO4)2-
# Data from Migdisov and Williams-Jones (2008)
T_logK <- c(25, 100, 150, 200, 250)
logK_NdSO4_1 <- c(3.65, 4.32, 4.86, 5.41, 6.06)
logK_NdSO4_2 <- c(5.15, 6.75, 8.17, 10.17, 12.40)
# Fit logK values and add to OBIGT
logK.to.OBIGT(logK_NdSO4_1, c("Nd+3", "SO4-2", "NdSO4+"), c(-1, -1, 1), T = T_logK, P = "Psat")
logK.to.OBIGT(logK_NdSO4_2, c("Nd+3", "SO4-2", "Nd(SO4)2-"), c(-1, -2, 1), T = T_logK, P = "Psat")

# Define substrates and solutes for solubility calculation
icr <- info(c("NdF3", "Nd(OH)3"))
iaq <- retrieve("Nd", c("H", "O", "S", "F", "Cl"), state = "aq")
# Remove Nd(OH)4- which is not shown by Migdisov et al.
iaq <- iaq[names(iaq) != "NdO2-"]

# Define temperature (degrees C), pressure (bar)
T <- 300
P <- 500
# Compositions from the figure caption
wt_pct_NaCl <- 10
wt_pct_Na2SO4 <- 2
# Fig. 21 of Migdisov et al. has 500 ppm HF,
# but that produces a larger stability field for NdF3
ppm_HF <- 100
ppm_Nd <- 200
# Convert to moles of NaCl, total S and F, and maximum Nd
# nb. wt_pct * 10 and ppm * 1e-3 both give parts per thousand, which is ca.
# g solute / kg H2O; divide this by molecular weight to get molality
m_NaCl <- wt_pct_NaCl * 10 / mass("NaCl")
S_tot <- wt_pct_Na2SO4 * 10 / mass("Na2SO4")
F_tot <- ppm_HF * 1e-3 / mass("HF")
Nd_max <- ppm_Nd * 1e-3 / mass("Nd")
# Define ranges of pH and logfO2
pH <- c(1, 10)
O2 <- c(-40, -20)

# Set up system
# NOTE: Nd must be first otherwise we get this cryptic error message for solubility():
# Error in basis(ispecies, loga.numeric) : species names are not unique
basis(c("Nd+3", "F-", "Cl-", "H2S", "H2O", "oxygen", "H+"))
basis("Cl-", log10(m_NaCl))
basis("H2S", log10(S_tot))
basis("F-", log10(F_tot))
basis("Nd+3", log10(Nd_max))

# Define candidate basis species to swap through for mosaic
# (the first species in each set must be present in the basis() definition above)
bases <- list(
  c("H2S", "HS-", "SO2", "HSO4-", "SO4-2"),
  c("F-", "HF")
)

### logfO2-pH diagram

if(FALSE) {
  species(icr)
  species(iaq, log10(Nd_max), add = TRUE)
  m <- mosaic(bases = bases, pH = pH, O2 = O2, T = T, P = P, IS = m_NaCl)
  diagram(m$A.species)
  diagram(m$A.bases[[1]], col = 8, col.names = 8, lty = 2, italic = TRUE, add = TRUE)
  diagram(m$A.bases[[2]], col = 2, col.names = 2, lty = 2, italic = TRUE, add = TRUE)
}

# Based on this we use:
# - logfO2 = -30
# - NdCl+2 as the predominant species at low pH

### Solubility calculation

basis("O2", -30)
# Dissolve three substrates (two minerals and the predominant aqueous species at low pH)
species(icr)
species("NdCl+2", log10(Nd_max), add = TRUE)

# Calculate solubility 
s <- solubility(iaq, bases = bases, pH = pH, T = T, P = P, IS = m_NaCl)
# Plot total solubility then aqueous speciation
diagram(s$aqueous, col = 4, lwd = 2, ylim = c(-9, -1))
diagram(s$aqueous, type="loga.equil", dy = 0.1, add = TRUE)

# Get stability ranges of substrates
d <- diagram(s$substrate, type = "loga.equil", plot.it = FALSE)
for(i in 1:nrow(d$species)) {
  # Get the mean pH for the stability range of this substrate
  mean_pH <- mean(d$vals$pH[d$predominant == i])
  text(mean_pH, -2, expr.species(d$species$name[i]), cex = 1.5)
}
# Add vertical lines at stability boundaries
abline(v = d$vals$pH[which(diff(d$predominant) != 0)], lty = 2, col = 8)

# Add legend and title
TP_text <- lTP(T, P)
composition_text <- bquote(.(wt_pct_NaCl)~"wt% NaCl,"~.(wt_pct_Na2SO4)~"wt%"~
                           .(expr.species("Na2SO4"))~.(ppm_HF)~"ppm HF, "~.(ppm_Nd)~" ppm Nd")
legend_text <- bquote(list(.(TP_text), .(composition_text)))
legend("top", legend = legend_text, bty = "n")
title("Solubility of Nd minerals, after Migdisov et al. (2016)", font.main = 1)

# Add note about lower solubility at high pH
arrows(9, -5.1, 9, -6.2)
text(9, -6.5, "New data\nmake this lower")
