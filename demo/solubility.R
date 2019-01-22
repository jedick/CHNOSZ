# CHNOSZ/demo/solubility.R: solubility of CO2 and calcite
# 20150306 jmd first version; used uniroot() to find zero affinity
# 20181031 use new vectorized, non-uniroot solubility(); add T-pH plots

# for comparison with published CO2 solubility plot, see Fig. 4.5 in
# Stumm and Morgan, 1996, Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters
# (New York: John Wiley & Sons), 3rd edition

# for comparison with published calcite solubility plot, see Fig. 4A in
# Manning et al., 2013, Reviews in Mineralogy & Geochemistry, v. 75, pp. 109-148
# (doi: 10.2138/rmg.2013.75.5)

opar <- par(no.readonly = TRUE)
layout(matrix(1:4, nrow = 2))

# set pH and T range and resolution, constant temperature and ionic strength
pH <- c(0, 14)
T <- c(0, 300)
res <- 100
T1 <- 25
IS <- 0

# start with CO2
basis(c("carbon dioxide", "H2O", "O2", "H+"))
# ca. atmospheric PCO2
basis("CO2", -3.5)
species(c("CO2", "HCO3-", "CO3-2"))
a <- affinity(pH = c(pH, res), T = T1, IS = IS)
s <- solubility(a)
# first plot total activity line
diagram(s, ylim = c(-10, 4), type = "loga.balance", lwd = 4, col = "green2")
# add activities of species
diagram(s, ylim=c(-10, 4), add = TRUE, dy = 1)
# add legend
lexpr <- as.expression(c("total", expr.species("CO2", state = "aq"),
  expr.species("HCO3-"), expr.species("CO3-2")))
legend("topleft", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
  col = c("green2", rep("black", 3)), legend = lexpr)
title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
  list(what = expr.species("CO2"), T = T1)), line = 1.6)
mtext("cf. Fig. 4.5 of Stumm and Morgan, 1996")
# check the endpoints
stopifnot(round(s$loga.balance[c(1, res)])==c(-5, 6))

# CO2 T-pH plot
a <- affinity(pH = c(pH, res), T = c(T, res), IS = IS)
s <- solubility(a)
diagram(s, type = "loga.balance")
title(main = substitute("Solubility of"~what, list(what = expr.species("CO2"))))

# now do calcite
basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
species(c("CO2", "HCO3-", "CO3-2"))
a <- affinity(pH = c(pH, res), T = T1, IS = IS)
s <- solubility(a)
diagram(s, ylim = c(-10, 4), type = "loga.balance", lwd = 4, col = "green2")
diagram(s, add = TRUE, dy = 1)
legend("topright", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
  col = c("green2", rep("black", 3)), legend = lexpr)
title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
  list(what = "calcite", T = T1)), line = 1.6)
mtext("cf. Fig. 4A of Manning et al., 2013")
# check the endpoints
stopifnot(round(s$loga.balance[c(1, res)])==c(4, -4))

# calcite T-pH plot
a <- affinity(pH = c(pH, res), T = c(T, res), IS = IS)
s <- solubility(a)
diagram(s, type = "loga.balance")
title(main = "Solubility of calcite", font.main = 1)

par(opar)
