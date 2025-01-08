# CHNOSZ/demo/minsol.R
# Make solubility diagram with multiple minerals
# 20190530 jmd first version (plot_Zn.R)
# 20201008 combine solubility contours for different minerals
# 20201014 added to CHNOSZ
# 20220715 Allow to change metals (renamed from zinc.R to minsol.R)

library(CHNOSZ)
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))

# System variables
metal <- "Zn"
res <- 300
T <- 100
P <- "Psat"
Stot <- 1e-3
pH <- c(0, 14, res)
O2 <- c(-62, -40, res)
# Mass fraction NaCl in saturated solution at 100 degC, from CRC handbook
wNaCl <- 0.2805  

# Set up basis species
basis(c(metal, "H2S", "Cl-", "oxygen", "H2O", "H+"))
basis("H2S", log10(Stot))
# Molality of NaCl
mNaCl <- 1000 * wNaCl / (mass("NaCl") * (1 - wNaCl))
# Estimate ionic strength and molality of Cl-
NaCl <- NaCl(m_NaCl = mNaCl, T = T)
basis("Cl-", log10(NaCl$m_Clminus))

# Add minerals and aqueous species
icr <- retrieve(metal, c("Cl", "S", "O"), state = "cr")
iaq <- retrieve(metal, c("Cl", "S", "O"), state = "aq")
logm_metal <- -3
species(icr)
species(iaq, logm_metal, add = TRUE)

# Calculate affinities and make diagram
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
m <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS)
d <- diagram(m$A.species, bold = TRUE)
diagram(m$A.bases, add = TRUE, col = 8, col.names = 8, lty = 3, italic = TRUE)
title(bquote(log * italic(m)[.(metal)*"(aq) species"] == .(logm_metal)))
label.figure("A")

# Add legend
plot.new()
l <- c(
  lTP(T, P),
  lNaCl(mNaCl),
  lS(Stot)
)
legend("topleft", legend = lex(l), bty = "n", cex = 1.5)
# Describe steps
par(xpd = NA)
legend("bottomleft", c("Predominance diagram: molality of aqueous", "species defines one solubility contour.",
  "Take away aqueous species to see", "all possible minerals.",
  "Calculate solubility for each mineral separately", "then find the minimum to plot solubilities", "of stable minerals across the diagram."),
       pch = c("A", "", "B", "", "C", "", ""), inset = c(-0.1, 0), cex = 0.95)
par(xpd = FALSE)

# Make diagram for minerals only 20201007
species(icr)
mcr <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS)
diagram(mcr$A.species, col = 2)
label.figure("B")

# Calculate *minimum* solubility among all the minerals 20201008
# (i.e. saturation condition for the solution)
# Use solubility() 20210303
s <- solubility(iaq, bases = bases, pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS, in.terms.of = metal)
# Specify contour levels
levels <- seq(-12, 9, 3)
diagram(s, levels = levels, contour.method = "flattest")

# Show the mineral stability boundaries
diagram(mcr$A.species, names = NA, add = TRUE, lty = 2, col = 2)
title(paste("Solubilities of", length(icr), "minerals"), font.main = 1, line = 1.5)
title(bquote(log[10]~"moles of"~.(metal)~"in solution"), line = 0.7)
label.figure("C")

par(opar)
