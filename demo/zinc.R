# CHNOSZ/demo/zinc.R
# make Zn solubility diagram with multiple minerals
# 20190530 jmd first version (plot_Zn.R)
# 20201008 combine solubility contours for different minerals
# 20201014 added to CHNOSZ

library(CHNOSZ)
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))

# system variables
res <- 300
T <- 100
P <- "Psat"
Stot <- 1e-3
pH <- c(0, 14, res)
O2 <- c(-62, -40, res)
# mass percent NaCl in saturated solution at 100 degC, from CRC handbook
w2 <- 0.2805  

# set up basis species
basis(c("ZnO", "H2S", "Cl-", "oxygen", "H2O", "H+"))
basis("H2S", log10(Stot))
# molality of NaCl in saturated solution
m2 <- 1000 * w2 / (mass("NaCl") * (1 - w2))
# estimate ionic strength and molality of Cl-
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

## Solubility plot 20201008

# Make a list to store the calculated solubilities for each mineral
slist <- list()
# Loop over minerals
minerals <- c("zincite", "sphalerite")
for(i in seq_along(minerals)) {
  # Define basis species with mineral to dissolve
  basis(c(minerals[i], "H2S", "Cl-", "oxygen", "H2O", "H+"))
  basis("H2S", log10(Stot))
  basis("Cl-", log10(sat$m_Cl))
  # Add aqueous species (no need to define activities here - they will be calculated)
  species(iaq)
  # Calculate affinities of formation reactions, using mosaic() to speciate the S-bearing basis species
  m <- mosaic(bases, pH = pH, O2 = O2, T = T, IS = sat$IS)
  # Calculate solubility of this mineral
  s <- solubility(m$A.species, in.terms.of = "Zn", dissociation = FALSE)
  # Store the solubilities in the list
  slist[[i]] <- s$loga.balance
}

# The overall solubility is the *minimum* among all the minerals
smin <- do.call(pmin, slist)
# Put this into the last-computed 'solubility' object
s$loga.balance <- smin
# Specify contour levels
levels <- seq(-12, 9, 3)
diagram(s, type = "loga.balance", levels = levels, contour.method = "flattest")
# Show the mineral stability boundaries
diagram(mcr$A.species, names = NA, add = TRUE, lty = 2, col = 2)
title("Solubilities of 2 minerals", font.main = 1, line = 1.5)
title(bquote(log[10]~"moles of Zn in solution"), line = 0.7)
label.figure("C")

par(opar)
