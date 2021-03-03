# CHNOSZ/demo/zinc.R
# Make Zn solubility diagram with multiple minerals
# 20190530 jmd first version (plot_Zn.R)
# 20201008 combine solubility contours for different minerals
# 20201014 added to CHNOSZ

library(CHNOSZ)
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))

# System variables
res <- 300
T <- 100
P <- "Psat"
Stot <- 1e-3
pH <- c(0, 14, res)
O2 <- c(-62, -40, res)
# Mass percent NaCl in saturated solution at 100 degC, from CRC handbook
w2 <- 0.2805  

# Set up basis species
basis(c("ZnO", "H2S", "Cl-", "oxygen", "H2O", "H+"))
basis("H2S", log10(Stot))
# Molality of NaCl in saturated solution
m2 <- 1000 * w2 / (mass("NaCl") * (1 - w2))
# Estimate ionic strength and molality of Cl-
sat <- NaCl(T = 100, m_tot = m2)
basis("Cl-", log10(sat$m_Cl))

# Add minerals and aqueous species
icr <- (c("sphalerite", "zincite"))
iaq <- (c("ZnCl4-2", "ZnO2-2"))
logm_Zn <- -3
species(c(icr, iaq), c(0, 0, logm_Zn, logm_Zn))

# Calculate affinities and make diagram
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
m <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = sat$IS)
fill <- c("lightgoldenrod1", "goldenrod1", "skyblue", "skyblue")
col.names <- c("red", "red", "blue", "blue")
d <- diagram(m$A.species, fill = fill, col.names = col.names, bold = TRUE)
diagram(m$A.bases, add = TRUE, col = "slategray", lwd = 2, lty = 3, names = NA)
title(bquote(log * italic(m)["Zn(aq) species"] == .(logm_Zn)))
label.figure("A")

# Add legend
plot.new()
l <- c(
  lTP(T, P),
  lNaCl(m2),
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

# Make diagram for Zn minerals only 20201007
if(packageVersion("CHNOSZ") <= "1.3.6") species(delete = TRUE)
species(icr)
mcr <- mosaic(bases, pH = pH, O2 = O2, T = T, P = P, IS = sat$IS)
diagram(mcr$A.species, col = 2)
label.figure("B")

# Calculate *minimum* solubility among all the minerals 20201008
# (i.e. saturation condition for the solution)
# Use solubilities() 20210303
# FIXME: why do we need to set dissociation = FALSE here?
s <- solubilities(mcr, iaq, in.terms.of = "Zn", dissociation = FALSE)
# Specify contour levels
levels <- seq(-12, 9, 3)
diagram(s, type = "loga.balance", levels = levels, contour.method = "flattest")

# Show the mineral stability boundaries
diagram(mcr$A.species, names = NA, add = TRUE, lty = 2, col = 2)
title("Solubilities of 2 minerals", font.main = 1, line = 1.5)
title(bquote(log[10]~"moles of Zn in solution"), line = 0.7)
label.figure("C")

par(opar)
