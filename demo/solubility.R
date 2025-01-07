# CHNOSZ/demo/solubility.R: solubility of CO2 and calcite
# 20150306 jmd first version; used uniroot() to find zero affinity
# 20181031 use new vectorized, non-uniroot solubility(); add T-pH plots
library(CHNOSZ)

# For comparison with published CO2 solubility plot, see Fig. 4.5 in
#   Stumm and Morgan, 1996, Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters
#   (New York: John Wiley & Sons), 3rd edition

# For comparison with published calcite solubility plot, see Fig. 4A in
#   Manning et al., 2013, Reviews in Mineralogy & Geochemistry, v. 75, pp. 109-148
#   (doi: 10.2138/rmg.2013.75.5)

opar <- par(no.readonly = TRUE)
layout(matrix(1:4, nrow = 2))

# Set pH and T range and resolution, constant temperature and ionic strength
pH <- c(0, 14)
T <- c(0, 300)
res <- 100
T1 <- 25
IS <- 0

# Start with CO2
basis(c("CO2", "H2O", "O2", "H+"))
# This is ca. atmospheric PCO2
species("carbon dioxide", -3.5)
iaq <- info(c("CO2", "HCO3-", "CO3-2"))
s <- solubility(iaq, pH = c(pH, res), T = T1, IS = IS)
# First plot total activity line
diagram(s, ylim = c(-10, 4), lwd = 4, col = "green2")
# Plot activities of species
diagram(s, type = "loga.equil", ylim=c(-10, 4), add = TRUE, dy = 1)
# Add legend
lexpr <- as.expression(c("total", expr.species("CO2", state = "aq"),
  expr.species("HCO3-"), expr.species("CO3-2")))
legend("topleft", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
  col = c("green2", rep("black", 3)), legend = lexpr)
title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
  list(what = expr.species("CO2"), T = T1)), line = 1.6)
mtext("cf. Fig. 4.5 of Stumm and Morgan, 1996")
# Check the endpoints
stopifnot(round(s$loga.balance[c(1, res)])==c(-5, 6))

# CO2 T-pH plot
s <- solubility(iaq, pH = c(pH, res), T = c(T, res), IS = IS)
diagram(s)
title(main = substitute("Solubility of"~what, list(what = expr.species("CO2"))))

# Now do calcite
basis(c("CO2", "Ca+2", "H2O", "O2", "H+"))
species("calcite")
iaq <- info(c("CO2", "HCO3-", "CO3-2"))
# Optional: use dissociate = 2 to get straight lines like Fig. 4A of Manning et al., 2013
s <- solubility(iaq, pH = c(pH, res), T = T1, IS = IS, dissociate = TRUE)
diagram(s, ylim = c(-10, 4), lwd = 4, col = "green2")
diagram(s, type = "loga.equil", add = TRUE, dy = 1)
legend("topright", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
  col = c("green2", rep("black", 3)), legend = lexpr)
title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
  list(what = "calcite", T = T1)), line = 1.6)
mtext("cf. Fig. 4A of Manning et al., 2013")
# Check the endpoints
stopifnot(round(s$loga.balance[c(1, res)])==c(4, -4))

# Calcite T-pH plot
s <- solubility(iaq, pH = c(pH, res), T = c(T, res), IS = IS, dissociate = TRUE)
diagram(s)
title(main = "Solubility of calcite", font.main = 1)

layout(matrix(1))
par(opar)
