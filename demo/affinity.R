# CHNOSZ/demo/affinity.R
## Affinities of metabolic reactions
## After Amend and Shock, 2001, Fig. 7
##  Amend, J. P. and Shock, E. L. (2001) Energetics of overall metabolic reactions of thermophilic and hyperthermophilic Archaea and Bacteria.
##  FEMS Microbiol. Rev. 25, 175--243. https://doi.org/10.1016/S0168-6445(00)00062-0
library(CHNOSZ)

# Use aq state for all basis species (including O2)
basis(c("CO2", "H2", "NH3", "O2", "H2S", "H+"), "aq")
# We're going to make H2O
species("H2O")
# A function to create the plots
doplot <- function(T) {
  res <- 20
  # calculate affinity/2.303RT as a function of loga(H2) and loga(O2)
  a <- affinity(H2 = c(-10, 0, res), O2 = c(-10, 0, res), T = T)
  # Temperature in Kelvin
  T.K <- convert(T, "K")
  # Convert dimensionless affinity (A/2.303RT) to Gibbs energy (J/mol)
  GJ <- convert(a$values[[1]], "G", T.K)
  # Convert J/mol to kJ/mol
  GkJ <- GJ / 1000
  # Now contour the values
  xyvals <- seq(-10, 0, length.out = res)
  contour(x = xyvals, y = xyvals, z = t(GkJ), levels = seq(-150, -250, -20),
    labcex = 1, xlab = axis.label("H2"), ylab = axis.label("O2"))
  # Show the temperature
  legend("topleft", bg = "white", cex = 1,
    legend = describe.property("T", T, digits = 0, ret.val = TRUE) )
}
# Plot layout with space for title at top
opar <- par(no.readonly = TRUE)
layout(matrix(c(1, 1, 2, 3, 4, 5), ncol=2, byrow = TRUE), heights = c(1, 4, 4))

par(mar = c(0, 0, 0, 0))
plot.new()
# We use subcrt() to generate a reaction for titling the plot
rxnexpr <- describe.reaction(subcrt("H2O", 1)$reaction, states = "all")
# Also in the title is the property with its units
Gexpr <- axis.label("DGr", prefix="k")[[2]]
text(0.5, 0.6, substitute(paste(G~"(kJ/mol) for"~r), list(G = Gexpr, r = rxnexpr)), cex = 2)
text(0.5, 0.2, "after Amend and Shock, 2001 Figure 7", cex = 2)
# Now make the plots
par(mar = c(3, 3, 0.5, 0.5), cex = 1.3, mgp = c(2, 1, 0))
sapply(c(25, 55, 100, 150), doplot)
# affinity() can handle the three dimensions simultaneously
print(affinity(H2 = c(-10, 0, 3), O2 = c(-10, 0, 3), T = c(25, 150, 4))$values)

# Reset plot settings
layout(matrix(1))
par(opar)

## Amino acid synthesis at low and high temperatures, based on:
##  Amend, J. P. and Shock, E. L. (1998) Energetics of amino acid synthesis in hydrothermal ecosystems.
##  Science 281, 1659--1662. https://doi.org/10.1126/science.281.5383.1659
# Select the basis species and species of interest
# and set their activities, first for the 18 degree C case
basis(c("H2O", "CO2", "NH4+", "H2", "H+", "H2S"),
  log10(c(1, 1e-4, 5e-8, 2e-9, 5e-9, 1e-15)))
species(sort(aminoacids("Z")),
  log10(c(3.9, 0.7, 1.1, 3.3, 0.5, 3.8, 1.0, 5.8, 1.2, 0.7,
  0.8, 1.0, 2.8, 0.5, 0.5, 4.6, 5.8, 0.6, 0.9, 2.8)/1e9))
T <- 18
TK <- convert(T, "K")
# Calculate A/2.303RT (dimensionless), convert to G of reaction (J/mol)
a <- affinity(T = T)
G.18.J <- convert(unlist(a$values), "G", T = TK)
# Convert to kJ/mol
G.18.kJ <- G.18.J / 1000
# The 100 degree C case
basis(c("H2O", "CO2", "NH4+", "H2", "H+", "H2S"),
  log10(c(1, 2.2e-3, 2.9e-6, 3.4e-4, 1.9e-6, 1.6e-3)))
species(1:20, log10(c(2.8e-9, 5.0e-10, 7.9e-10, 2.4e-9, 3.6e-10,
  2.7e-9, 7.2e-10, 4.2e-9, 8.6e-10, 5.0e-10, 5.7e-10, 7.2e-10, 2.0e-9,
  3.6e-10,3.6e-10, 3.3e-9, 4.2e-9, 4.3e-10, 6.5e-10, 2.0e-9)))
T <- 100
TK <- convert(T, "K")
a <- affinity(T = T)
G.100.J <- convert(unlist(a$values), "G", T = TK)
G.100.kJ <- G.100.J / 1000
# Rhe average oxidation states of carbon
Z.C <- ZC(thermo()$OBIGT$formula[thermo()$species$ispecies])
# Put everything together like Table 3 in the paper
print(out <- data.frame(G.18 = G.18.kJ, G.100 = G.100.kJ, Z.C = Z.C))
# Make a plot; set units to get correct label
plot(out$Z.C, out$G.18, pch = 20, xlim = c(-1.1, 1.1), ylim = c(-200, 500), 
  xlab = axis.label("ZC"), ylab = axis.label("DGr", prefix = "k"))
points(out$Z.C, out$G.100, col = "red", pch = 20)
legend("topleft", pch = c(20, 20), col = c("black", "red"),
  legend = describe.property(c("T", "T"), c(18, 100)))
title(main = "Amino acid synthesis, after Amend and Shock, 1998")
# 9 amino acids have negative delta Gr under hydrothermal conditions
# (cf. AS98 with 11; we are using more recent thermodynamic data)
stopifnot(sum(out$G.100 < 0) == 9)
