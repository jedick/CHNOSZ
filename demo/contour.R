# CHNOSZ/demo/contour.R
# gold solubility contours on logfO2-pH diagram
# 20181107 initial version
# 20190415 cleanup for demo

# After Williams-Jones et al., 2009, Fig. 3
# doi:10.2113/gselements.5.5.281

# define temperature (degrees C), pressure (bar), grid resolution
T <- 250
P <- 500
res <- 600
# make smooth (TRUE) or sharp (FALSE) transitions between basis species 
blend <- TRUE

# set up system
basis(c("Au", "Cl-", "H2S", "H2O", "oxygen", "H+"))
species(c("Au(HS)2-", "AuHS", "AuOH", "AuCl2-"))
# this get us close to total S = 0.01 m
basis("H2S", -2)
# calculate solution composition for 1 mol/kg NaCl
NaCl <- NaCl(T = T, P = P, m_tot=1)
basis("Cl-", log10(NaCl$m_Cl))
# calculate affinity with changing basis species
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
m <- mosaic(bases, pH = c(2, 10, res), O2 = c(-41, -29, res), T = T, P = P, IS = NaCl$IS, blend = blend)
# calculate and plot solubility
s <- solubility(m$A.species)
# convert to ppb
s <- convert(s, "ppb")
diagram(s, type="loga.balance", levels = c(1, 10, 100, 1000))
# show predominance fields
diagram(m$A.bases, add=TRUE, col = "red", col.names = "red", limit.water = FALSE, lty = 2, italic = TRUE)
diagram(m$A.species, add=TRUE, col = "blue", col.names = "blue", limit.water = FALSE, lwd = 2, bold = TRUE)
# add legend and title
dP <- describe.property(c("T", "P"), c(250, 500))
legend("top", dP, bty = "n", inset = c(0, 0.06))
dNaCl <- expression(NaCl == 1~mol~kg^-1)
dS <- expression(sum(S) ~"in basis" == 0.01~mol~kg^-1)
legend("topright", c(dNaCl, dS), bty = "n", inset = c(0, 0.05))
title(main = ("Solubility of gold (ppb), after Williams-Jones et al., 2009, Fig. 3"), font.main = 1)
