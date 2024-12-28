# CHNOSZ/demo/aluminum.R 
# 20171018 comparisons with calculations of SUPCRTBL
# 20190415 add diagrams from Tutolo et al., 2014
library(CHNOSZ)

## Set up plotting area
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))

###########
### Plot 1: boehmite - kaolinite equilibrium
###########
# After Zhu and Lu, 2009 (doi:10.1016/j.gca.2009.03.015)
# Experimental data from Table 1 of Hemley et al., 1980 (doi:10.2113/gsecongeo.75.2.210)
xT <- c(200, 200, 200, 200, 250, 250, 265, 300, 300, 300, 300)
xlogaSiO2 <- -c(2.54, 2.59, 2.65, 2.77, 2.21, 2.32, 2.12, 1.90, 1.95, 1.94, 1.90)
## Set up basis species so that axis.label shows activity of SiO2
basis(c("Al2O3","SiO2", "H2O", "O2"))
T <- 125:350
thermo.plot.new(xlim = range(T), ylim = c(-3.5, -1.5), xlab = axis.label("T"), ylab = axis.label("SiO2"))
points(xT, xlogaSiO2)
basis(delete = TRUE)
## First calculation: CHNOSZ default
# kaolinite from Berman, 1988
# boehmite from Hemingway et al., 1991
r1 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T = T, P = 1000, exceed.Ttr = TRUE) 
## Second calculation: get SiO2(aq) from Apps and Spycher, 2004
add.OBIGT("SiO2")
r2 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T = T, P = 1000, exceed.Ttr = TRUE) 
reset()
## Third calculation: get Si(OH)4 from Akinfiev and Plyasunov, 2014
add.OBIGT("AD")
r3 <- subcrt(c("boehmite", "Si(OH)4", "H2O", "kaolinite"), c(-1, -1, 1.5, 0.5), T = T, P = 1000, exceed.Ttr = TRUE) 
reset()
## Fourth calculation: minerals as in SUPCRT92
add.OBIGT("SUPCRT92") # gets kaolinite and boehmite from HDNB78
# We need exceed.Ttr = TRUE because the T limit for boehmite is 500 K (Helgeson et al., 1978)
r4 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T = T, P = 1000, exceed.Ttr = TRUE) 
reset()
## log activity of SiO2 is -ve logK
lines(T, -r1$out$logK, lwd = 1.5)
lines(T, -r2$out$logK, col = "red", lty = 2)
lines(T, -r3$out$logK, col = "purple", lty = 3)
lines(T, -r4$out$logK, col = "blue1", lty = 4)
## Add points calculated using the SUPCRTBL package
points(seq(125, 350, 25), -c(3.489, 3.217, 2.967, 2.734, 2.517, 2.314, 2.124, 1.946, 1.781, 1.628), pch = 4, col = "red")
## Add legend and title
title(main = describe.reaction(r4$reaction), cex.main = 1.1)
legend("bottomright", lty = c(0, 0, 1, 2, 3, 0), pch = c(1, 4, NA, NA, NA, NA), lwd = c(1, 1, 1.5, 1, 1, 0),
       col = c("black", "red", "black", "red", "purple", NA), bty = "n", cex = 0.9,
       legend = c("Hemley et al., 1980", "SUPCRTBL", "CHNOSZ", 'add.OBIGT("SiO2")', 'add.OBIGT("AD")', ""))
legend("bottomright", lty = 4, pch = NA, lwd = 1, col = "blue", legend = 'add.OBIGT("SUPCRT92")', bty = "n", cex = 0.9)
legend("topleft", c("Boehmite - Kaolinite", "After Zhu and Lu, 2009 Fig. A1"), bty = "n")
# Helgeson et al., 1978 (HDNB78): http://www.worldcat.org/oclc/13594862
# Shock et al., 1989 (SHS89): doi:10.1016/0016-7037(89)90341-4
# Berman, 1988 (Ber88): doi:10.1093/petrology/29.2.445
# Holland and Powell, 2011 (HP11): 10.1111/j.1525-1314.2010.00923.x
# Hemingway et al., 1991 (HRA91): https://pubs.usgs.gov/publication/70016664
# Apps and Spycher, 2004 (AS04): Bechtel SAIC Company, LLC ANL-NBS-HS-000043 REV 00 (DOC.20041118.0004)

###########
### Plot 2: dawsonite solubility
###########
# After Zimmer et al., 2016 (doi:10.1016/j.cageo.2016.02.013)
# Dxperimental data from Benezeth et al., 2007 Table 5 (doi:10.1016/j.gca.2007.07.003)
# (averages for each temperature in a single run)
T <- c(100.1, 100.1, 150.1, 100.1, 150.1, 99.8, 99.8, 200.7, 99.8, 50.1, 75.1, 100.3, 150.1)
logK <- -c(14.825, 14.735, 13.625, 14.79, 13.665, 14.725, 14.1775, 12.74, 14.4925, 16.8625, 15.61, 14.51, 13.455)
thermo.plot.new(c(25, 250), c(-18, -10), axis.label("T"), axis.label("logK"))
points(T, logK)
# Calculation 1: CHNOSZ default
T <- 0:250
species <- c("dawsonite", "H2O", "Al(OH)4-", "HCO3-", "Na+", "H+")
coeffs <- c(-1, -2, 1, 1, 1, 1)
Daw1 <- subcrt(species, coeffs, T = T)
lines(T, Daw1$out$logK, lwd = 1.5)
# Calculation 2: dawsonite with Cp = 0
mod.OBIGT("dawsonite", Cp = 0, a = 0, b = 0, c = 0)
Daw2 <- subcrt(species, coeffs, T = T)
lines(T, Daw2$out$logK, col = "red", lty = 2)
## Add points calculated using the SUPCRTBL package
#points(seq(25, 250, 25), c(-17.829, -16.523, -15.402, -14.425, -13.568, -12.815, -12.154, -11.581, -11.094, -10.699), pch=4, col="red")
## 20190417: recalculated using the SUPCRTBL package (timestamp: 20190309)
##   with a locally updated data file that includes heat capacity coefficients of dawsonite
##   from Robie and Hemingway, 1995, with typos corrected in Tutolo et al., 2014
points(seq(25, 250, 25), c(-17.829, -16.546, -15.485, -14.599, -13.856, -13.236, -12.724, -12.312, -11.997, -11.782), pch=4, col="red")
## Add legend and title
title(main = describe.reaction(Daw1$reaction), cex.main = 0.95)
legend("bottomright", lty = c(0, 0, 0, 1, 2), pch = c(1, 4, NA, NA, NA), col = c("black", "red", NA, "black", "red"), lwd = c(1, 1, 0, 1.5, 1),
       bty = "n", cex = 0.9, legend = c("Ben\u00e9z\u00e9th et al., 2007", "SUPCRTBL with Cp", "  coefficients for dawsonite", "CHNOSZ", "Cp(dawsonite) = 0"))
legend("topleft", c("Dawsonite solubility", "After Zimmer et al., 2016 Fig. 2"), bty = "n")
reset()

###########
### Plot 3: kaolinite solubility
###########
# After Tutolo et al., 2014, Fig. 2 (doi:10.1016/j.gca.2014.02.036)
dat <- read.csv(system.file("extdata/misc/TKSS14_Fig2.csv", package = "CHNOSZ"))
thermo.plot.new(c(3.5, 1.5), c(-2, 14), quote(1000 / italic(T)*"(K)"), quote(p*italic(K)))
points(dat)
# Plot line: default database
invTK <- seq(3.5, 1.6, -0.02)
T <- 1000/invTK - 273.15
sres <- subcrt(c("kaolinite", "OH-", "H2O", "Al(OH)4-", "SiO2"), c(-1, -2, -1, 2, 2), T = T)
pK <- -sres$out$logK
lines(invTK, pK, lwd = 1.5)
# Plot line: SiO2 from Apps and Spycher, 2004
add.OBIGT("SiO2")
sres <- subcrt(c("kaolinite", "OH-", "H2O", "Al(OH)4-", "SiO2"), c(-1, -2, -1, 2, 2), T = T)
pK <- -sres$out$logK
lines(invTK, pK, col = "red", lty = 2)
reset()
# Plot line: Si(OH)4 from Akinfiev and Plyasunov, 2014
add.OBIGT("AD")
sres <- subcrt(c("kaolinite", "OH-", "H2O", "Al(OH)4-", "Si(OH)4"), c(-1, -2, -5, 2, 2), T = T)
pK <- -sres$out$logK
lines(invTK, pK, col = "purple", lty = 3)
reset()
# Plot line: SUPCRT92
add.OBIGT("SUPCRT92")
sres <- subcrt(c("kaolinite", "OH-", "H2O", "Al(OH)4-", "SiO2"), c(-1, -2, -1, 2, 2), T = T)
pK <- -sres$out$logK
lines(invTK, pK, col = "blue", lty = 4)
reset()
# Add points calculated using the SUPCRTBL package
T <- seq(25, 300, 25)
invTK <- 1000/(T + 273.15)
points(invTK, c(12.621, 11.441, 10.383, 9.402, 8.477, 7.597, 6.756, 5.948, 5.171, 4.422, 3.703, 3.023), pch = 4, col = "red")
# Add title and legend
par(xpd = NA)
title(main = describe.reaction(sres$reaction), cex.main = 1.1)
par(xpd = FALSE)
legend("topright", c("Kaolinite solubility", "After Tutolo et al., 2014 Fig. 2"), bty = "n")
legend("bottomleft", lty = c(0, 0, 0, 1, 2, 3, 4), pch = c(1, NA, 4, NA, NA, NA, NA),
       lwd = c(1, 1, 1, 1.5, 1, 1, 1), col = c("black", "black", "red", "black", "red", "purple", "blue"),
       legend = c("Various sources \u2013", "  see Tutolo et al., 2014", "SUPCRTBL", "CHNOSZ", 'add.OBIGT("SiO2")', 'add.OBIGT("AD")', 'add.OBIGT("SUPCRT92")'),
       bty = "n", cex = 0.9)

###########
### Plot 4: albite - K-feldspar exchange
###########
# After Tutolo et al., 2014, Fig. 5 (doi:10.1016/j.gca.2014.02.036)
# Experimental data from Merino, 1975, Table 4 (doi:10.1016/0016-7037(75)90085-X)
# Plot line calculated using default database
basis(c("Al2O3", "SiO2", "K+", "Na+", "O2", "H2O", "H+"))
species(c("albite", "K-feldspar"))
T <- 100
P <- 150
a <- affinity("K+" = c(4, 7), "Na+" = c(6, 9), T = T, P = P)
diagram(a, lwd = 1.5, xlab = ratlab("K+"), ylab = ratlab("Na+"), names = FALSE)
# Plot experimental data
dat <- read.csv(system.file("extdata/misc/Mer75_Table4.csv", package = "CHNOSZ"))
points(dat$log.aK..aH.., dat$log.aNa..aH..)
# Plot line calculated using SUPCRT92 data
add.OBIGT("SUPCRT92")
a <- affinity("K+" = c(4, 7), "Na+" = c(6, 9), T = 100, P = 150)
diagram(a, col = "blue", lty = 4, add = TRUE, names = FALSE)
# Add SUPCRTBL calculation
logK_BL <- 2.092
logaK <- seq(4, 7, 0.5)
logaNa <- logaK + logK_BL
points(logaK, logaNa, pch = 4, col = "red")
# Add title and legend
sres <- subcrt(c("albite", "K+", "K-feldspar", "Na+"), c(-1, -1, 1, 1))
title(main = describe.reaction(sres$reaction), cex.main = 1.1)
legend("topleft", c("Albite - K-feldspar", "After Tutolo et al., 2014 Fig. 5"), bty = "n", cex = 0.9)
legend("bottomright", lty = c(0, 0, 1, 0), pch = c(1, 4, NA, NA), lwd = c(1, 1, 1.5, 0), col = c("black", "red", "black", NA),
       legend = c("Merino, 1975", "SUPCRTBL", "CHNOSZ", ""), bty = "n", cex = 0.9)
legend("bottomright", lty = 4, pch = NA, lwd = 1, col = "blue", legend = 'add.OBIGT("SUPCRT92")', bty = "n", cex = 0.9)
legend("left", describe.property(c("T", "P"), c(T, P)), bty = "n")
reset()

par(opar)
