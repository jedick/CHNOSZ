# CHNOSZ/demo/mosaic.R
# 20141221 first version

# Fe-minerals and aqueous species in Fe-S-O-H-C system
# after Garrels and Christ, 1965 Figure 7.21
# to reproduce their diagram as closely as posssible, use their thermodynamic data (from Appendix 2)
mod.OBIGT(c("Fe+2", "Fe+3"), G = c(-20300, -2520))
mod.OBIGT(c("hematite", "magnetite", "pyrrhotite", "pyrite", "siderite"), G = c(-177100, -242400, -23320, -36000, -161060))
mod.OBIGT(c("SO4-2", "HS-", "H2S", "HSO4-"), G = c(-177340, 3010, -6540, -179940))
mod.OBIGT(c("CO2", "HCO3-", "CO3-2"), G = c(-92310, -140310, -126220))
# conditions and system definition
pH <- c(0, 14, 500)
Eh <- c(-1, 1, 500)
T <- 25
basis(c("FeO", "SO4-2", "H2O", "H+", "e-", "CO3-2"))
basis("SO4-2", -6)
basis("CO3-2", 0)
species(c("Fe+2", "Fe+3"), -6)
species(c("pyrrhotite", "pyrite", "hematite", "magnetite", "siderite"), add = TRUE)
# two sets of changing basis species:
# speciate SO4-2, HSO4-, HS-, H2S as a function of Eh and pH
# speciate CO3-2, HCO3-, CO2 as a function of pH
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
bases2 <- c("CO3-2", "HCO3-", "CO2")
# calculate affinities using the relative abundances of different basis species
# (using default blend = TRUE)
# note curved lines, particularly at the boundaries with siderite
m1 <- mosaic(bases, bases2, pH = pH, Eh = Eh, T = T)
# make a diagram and add water stability lines
diagram(m1$A.species, lwd = 2)
water.lines(m1$A.species, col = "seagreen", lwd = 1.5)
# show lines for Fe(aq) = 10^-4 M
species(c("Fe+2", "Fe+3"), -4)
m2 <- mosaic(bases, bases2, pH = pH, Eh = Eh, T = T)
diagram(m2$A.species, add = TRUE, names = FALSE)
# overlay the sulfur and carbonate basis species predominance fields
d <- diagram(m1$A.bases, add = TRUE, col = "red", col.names = "red", lty = 3, limit.water = FALSE)
d <- diagram(m1$A.bases2, add = TRUE, col = "blue", names = FALSE, lty = 3, limit.water = FALSE)
text(d$namesx, -0.8, as.expression(sapply(m1$A.bases2$species$name, expr.species)), col = "blue")
# add legend and title
dP <- describe.property(c("T", "P"), c(25, 1))
legend("top", dP, bty = "n")
dS <- expression(sum(S)*"(aq)" == 10^-6~italic(m))
dC <- expression(sum(C)*"(aq)" == 1~italic(m))
legend("topright", c(dS, dC), bty = "n")
title(main=paste("Iron minerals, sulfur, and carbonate in water,",
  "after Garrels and Christ, 1965, Figure 7.21", sep = "\n"), font.main = 1)
# reset the database, as it was changed in this example
reset()
