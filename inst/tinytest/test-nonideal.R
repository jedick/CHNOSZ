# Load default settings for CHNOSZ
reset()

# 20161219
info <- "loggam and logK values are consistent"
oldnon <- nonideal("Alberty")
rxn1 <- subcrt(c("anhydrite", "Ca+2", "SO4-2"), c(-1, 1, 1), P = 1, T = 25, I = 0)
rxn2 <- subcrt(c("anhydrite", "Ca+2", "SO4-2"), c(-1, 1, 1), P = 1, T = 25, I = 0.7)
expect_equal(rxn1$out$logK, rxn2$out$loggam + rxn2$out$logK, info = info)
nonideal(oldnon)

# The maximum absolute pairwise difference between x and y
maxdiff <- function(x, y) max(abs(y - x))

# 20171011
info <- "A and B parameters are calculated correctly"
## Psat
# values from Helgeson, 1967 "Solution chemistry and metamorphism"
# (chapter in http://www.worldcat.org/oclc/152725534)
T <- c(25, 50, 100, 200, 300, 350)
A <- c(0.5095, 0.5354, 0.6019, 0.8127, 1.2979, 1.9936)
B <- c(0.3284, 0.3329, 0.3425, 0.3659, 0.4010, 0.4300)
SUP <- water.SUPCRT92(c("A_DH", "B_DH"), T = convert(T, "K"), P = "Psat")
IAP <- water.IAPWS95(c("A_DH", "B_DH"), T = convert(T, "K"), P = "Psat")
expect_true(maxdiff(SUP$A_DH, A) < 0.18, info = info)
expect_true(maxdiff(IAP$A_DH, A) < 0.11, info = info)
expect_true(maxdiff(SUP$B_DH / 1e8, B) < 0.012, info = info)
expect_true(maxdiff(IAP$B_DH / 1e8, B) < 0.008, info = info)

# values digitized from Fig. 10 of Manning et al., 2013
# doi: 10.2138/rmg.2013.75.5
## 5 kbar
T5 <- c(25, seq(100, 1000, 100))
A5 <- c(0.441, 0.49, 0.583, 0.685, 0.817, 0.983, 1.164, 1.409, 1.673, 1.938, 2.187)
B5 <- c(0.328, 0.336, 0.350, 0.363, 0.377, 0.391, 0.405, 0.421, 0.434, 0.445, 0.453)
SUP5 <- water.SUPCRT92(c("A_DH", "B_DH"), T = convert(T5, "K"), P = rep(5000, 11))
IAP5 <- water.IAPWS95(c("A_DH", "B_DH"), T = convert(T5, "K"), P = rep(5000, 11))
DEW5 <- water.DEW(c("A_DH", "B_DH"), T = convert(T5, "K"), P = rep(5000, 11))
# DEW is the winner here
expect_true(maxdiff(SUP5$A_DH, A5) < 0.47, info = info)
expect_true(maxdiff(IAP5$A_DH, A5) < 0.26, info = info)
expect_true(maxdiff(DEW5$A_DH, A5) < 0.14, info = info)
expect_true(maxdiff(SUP5$B_DH / 1e8, B5) < 0.036, info = info)
expect_true(maxdiff(IAP5$B_DH / 1e8, B5) < 0.021, info = info)
expect_true(maxdiff(DEW5$B_DH / 1e8, B5) < 0.013, info = info)

## 30 kbar
T30 <- seq(700, 1000, 100)
A30 <- c(0.625, 0.703, 0.766, 0.815)
B30 <- c(0.386, 0.400, 0.409, 0.416)
DEW30 <- water.DEW(c("A_DH", "B_DH"), T = convert(T30, "K"), P = rep(30000, 4))
expect_true(maxdiff(DEW30$A_DH, A30) < 0.06, info = info)
expect_true(maxdiff(DEW30$B_DH / 1e8, B30) < 0.024, info = info)

info <- "Affinity transect incorporates IS correctly"
basis("CHNOS+")
species("acetate")
# calculations at single combinations of logfO2 and IS
basis("O2", -80); a80_0 <- affinity()
basis("O2", -60); a60_1 <- affinity(IS = 1)
# calculations on a transect with those endpoints
a <- affinity(O2 = seq(-80, -60, length.out = 4), IS = seq(0, 1, length.out = 4))
expect_equal(a$values[[1]][1], a80_0$values[[1]][1], info = info)
expect_equal(a$values[[1]][4], a60_1$values[[1]][1], info = info)
# 20171013: that was working fine, but how about a more complicated case involving T?
a25_0 <- affinity()
a50_1 <- affinity(T = 50, IS = 1)
a <- affinity(T = seq(25, 50, length.out = 4), IS = seq(0, 1, length.out = 4))
expect_equal(a$values[[1]][1], a25_0$values[[1]][1], info = info)
expect_equal(a$values[[1]][4], a50_1$values[[1]][1], info = info)

# 20171221
info <- "Nonideality calculations work for Zn"
# nonideal() had a bug where both Z and Zn were identified as the charge
# in the species formula, producing an error in this calculation
expect_true(is.list(subcrt(c("Zn+2", "Cl-", "ZnCl+"), c(-1, -1, 1), T = 200, P = 16, IS = 0.05)), info = info)

# 20181105
info <- "Activity coefficients are similar to those from HCh"
# ionic strength of solution and activity coefficients of Na+ and Cl-
# calculated with HCh version 3.7 (Shvarov and Bastrakov, 1999) at 1000 bar,
# 100, 200, and 300 degress C, and 1 to 6 molal NaCl,
# using the default "B-dot" activity coefficient model (Helgeson, 1969)
# and the default setting for the Setchenow equation,
# for which the only non-zero term is the mole fraction to molality conversion factor
IS.HCh <- list(`100` = c(0.992, 1.969, 2.926, 3.858, 4.758, 5.619),
               `300` = c(0.807, 1.499, 2.136, 2.739, 3.317, 3.875),
               `500` = c(0.311, 0.590, 0.861, 1.125, 1.385, 1.642))
gamCl.HCh <- list(`100` = c(0.565, 0.545, 0.551, 0.567, 0.589, 0.615),
                  `300` = c(0.366, 0.307, 0.275, 0.254, 0.238, 0.224),
                  `500` = c(0.19, 0.137, 0.111, 0.096, 0.085, 0.077))
gamNa.HCh <- list(`100` = c(0.620, 0.616, 0.635, 0.662, 0.695, 0.730),
                  `300` = c(0.421, 0.368, 0.339, 0.318, 0.302, 0.288),
                  `500` = c(0.233, 0.180, 0.155, 0.138, 0.126, 0.117))
gamNaCl.HCh <- list(`100` = c(0.965, 0.933, 0.904, 0.876, 0.850, 0.827),
                    `300` = c(0.968, 0.941, 0.915, 0.892, 0.870, 0.849),
                    `500` = c(0.977, 0.955, 0.935, 0.915, 0.897, 0.879))
# calculate activity coefficent of Cl- at each temperature
gamCl.100 <- 10^subcrt("Cl-", T = 100, P = 1000, IS = IS.HCh$`100`)$out$`Cl-`$loggam
gamCl.300 <- 10^subcrt("Cl-", T = 300, P = 1000, IS = IS.HCh$`300`)$out$`Cl-`$loggam
gamCl.500 <- 10^subcrt("Cl-", T = 500, P = 1000, IS = IS.HCh$`500`)$out$`Cl-`$loggam
expect_true(maxdiff(gamCl.100, gamCl.HCh$`100`) < 0.07, info = info)
expect_true(maxdiff(gamCl.300, gamCl.HCh$`300`) < 0.03, info = info)
expect_true(maxdiff(gamCl.500, gamCl.HCh$`500`) < 0.009, info = info)
# calculate activity coefficent of Na+ at each temperature
gamNa.100 <- 10^subcrt("Na+", T = 100, P = 1000, IS = IS.HCh$`100`)$out$`Na+`$loggam
gamNa.300 <- 10^subcrt("Na+", T = 300, P = 1000, IS = IS.HCh$`300`)$out$`Na+`$loggam
gamNa.500 <- 10^subcrt("Na+", T = 500, P = 1000, IS = IS.HCh$`500`)$out$`Na+`$loggam
expect_true(maxdiff(gamNa.100, gamNa.HCh$`100`) < 0.08, info = info)
expect_true(maxdiff(gamNa.300, gamNa.HCh$`300`) < 0.03, info = info)
expect_true(maxdiff(gamNa.500, gamNa.HCh$`500`) < 0.013, info = info)
# calculate activity coefficent of NaCl at each temperature
gamNaCl.100 <- 10^subcrt("NaCl", T = 100, P = 1000, IS = IS.HCh$`100`)$out$`NaCl`$loggam
gamNaCl.300 <- 10^subcrt("NaCl", T = 300, P = 1000, IS = IS.HCh$`300`)$out$`NaCl`$loggam
gamNaCl.500 <- 10^subcrt("NaCl", T = 500, P = 1000, IS = IS.HCh$`500`)$out$`NaCl`$loggam
expect_true(maxdiff(gamNaCl.100, gamNaCl.HCh$`100`) < 0.09, info = info)
expect_true(maxdiff(gamNaCl.300, gamNaCl.HCh$`300`) < 0.09, info = info)
expect_true(maxdiff(gamNaCl.500, gamNaCl.HCh$`500`) < 0.10, info = info)

# 20190603
info <- "G, H, S, and Cp corrections are calculated consistently"
nonideal("Alberty")

# first test: DG = DH - T * DS
T <- seq(25, 100, 25)
GHSstd <- subcrt("Na+", T = T)$out[[1]][, c("G", "H", "S")]
GHSadj <- subcrt("Na+", T = T, IS = 0.25)$out[[1]][, c("G", "H", "S")]
D <- GHSadj - GHSstd
TK <- convert(T, "K")
expect_equal(D$G, D$H - TK * D$S, info = info)

# second test: Cp = dH/dT
# calculate the derivative numerically at 50 degC
T <- c(49.9, 50.1)
# first do it for standard properties (at I = 0) to make sure the test is set up correctly
Cpstd <- subcrt("Na+", T = 50)$out[[1]][, "Cp"]
Hstd <- subcrt("Na+", T = T)$out[[1]][, "H"]
expect_equal(Cpstd, diff(Hstd) / diff(T), tol = 1e-6, info = info)
# now at I = 0.25 mol kg-1
Cpadj <- subcrt("Na+", T = 50, IS = 0.25)$out[[1]][, "Cp"]
Hadj <- subcrt("Na+", T = T, IS = 0.25)$out[[1]][, "H"]
expect_equal(Cpadj, diff(Hadj) / diff(T), tol = 1e-6, info = info)

# another way to test things: compare our calculations with values of the Debye-Huckel limiting slopes (Alberty, 2001, Table 1)
Glim <- c(2.56494, 2.70073, 2.84196, 2.91482, 2.98934, 3.14349)
Hlim <- c(1.075, 1.213, 1.3845, 1.4775, 1.5775, 1.800)
Cplim <- c(13.255, 15.41, 17.90, 19.27, 20.725, 23.885)
T <- c(0, 10, 20, 25, 30, 40)
E.units("J")
IS <- 0.25
GHSCpstd <- subcrt("Na+", T = T)$out[[1]][, c("G", "H", "S", "Cp")]
GHSCpadj <- subcrt("Na+", T = T, IS = IS)$out[[1]][, c("G", "H", "S", "Cp")]
D <- GHSCpadj - GHSCpstd
# now make the checks
IBterm <- IS ^ 0.5 / (1 + 1.6 * IS ^ 0.5)
expect_equal(Glim * 1000, - D$G / IBterm, tol = 1e-4, info = info)
expect_equal(Hlim * 1000, D$H / IBterm, tol = 1e-3, info = info)
expect_equal(Cplim, D$Cp / IBterm, tol = 1e-4, info = info)

info <- "Affinity has IS-T and T-IS calculations"
basis("CHNOS+")
species(c("H2S", "HS-", "HSO4-", "SO4-2", "SO2"))
# reference values: affinity as a function of IS at T = 250 degC
a250 <- affinity(IS = seq(0, 1, 0.2), T = 250)
# increasing ionic strength should raise the affinity of HS-
expect_true(all(diff(a250$values[[2]]) > 0), info = info)
# T-IS grid; use different resolutions to help identify indexing bugs
a_T.IS <- affinity(IS = c(0, 1, 6), T = c(190, 310, 5))
expect_equivalent(a_T.IS$values[[2]][, 3], as.numeric(a250$values[[2]]), info = info)
# IS-T grid
a_IS.T <- affinity(T = c(190, 310, 5), IS = c(0, 1, 6))
expect_equivalent(a_IS.T$values[[2]][3, ], as.numeric(a250$values[[2]]), info = info)
