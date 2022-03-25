# Load default settings for CHNOSZ
reset()

# Get properties of water from DEW implementation in CHNOSZ
water("DEW")
# Use DEW species parameters in OBIGT database
add.OBIGT("DEW")

info <- "Density of water is calculated correctly"
pressure <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
temperature <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
# Density from R translation of DEW macro functions
RDensity <- calculateDensity(pressure, temperature)
# Density from DEW spreadsheet
DEWDensity <- c(1.108200, 0.597623, 1.196591, 0.798331, 1.321050, 1.000735, 1.578116, 1.287663)
expect_equal(RDensity, DEWDensity, tolerance = 1e-6, info = info)

info <- "Gibbs energy of water is calculated correctly"
pressure <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
temperature <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
# Gibbs energies from R translation of DEW macro functions
RGibbs <- calculateGibbsOfWater(pressure, temperature)
# Gibbs energies from DEW spreadsheet
DEWGibbs <- c(-56019.85419280258, -84262.028821198, -54155.004480575895, -81210.38766217149,
              -50735.122222685815, -76433.07602205424, -41823.26077175943, -65187.48113532527)
expect_equal(RGibbs, DEWGibbs, info = info)

info <- "Dielectric constant of water is calculated correctly"
pressure <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
temperature <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
# epsilon from R translation of DEW macro functions
Repsilon <- calculateEpsilon(calculateDensity(pressure, temperature), temperature)
# epsilon from DEW spreadsheet
DEWepsilon <- c(65.63571, 6.10465, 72.40050, 8.97800, 82.16244, 12.13131, 103.12897, 16.97266)
expect_equal(Repsilon, DEWepsilon, tolerance = 1e-7, info = info)

info <- "Born coefficient Q is calculated correctly"
pressure <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
temperature <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
# Q from R translation of DEW macro functions
RQ <- calculateQ(calculateDensity(pressure, temperature), temperature)
# Q from DEW spreadsheet
DEWQ <- c(0.32319817, 14.50286092, 0.19453478, 3.12650897,
          0.10918151, 0.87729257,  0.05068788, 0.20640645) / 1e6
expect_equal(RQ, DEWQ, info = info)

info <- "g function is calculated correctly"
pressure <- c(1000, 1000, 5000, 5000, 10000)
temperature <- c(100, 1000, 100, 1000, 100)
# note that values returned for alpha, daldT, beta are NA
w <- water(c("rho", "alpha", "daldT", "beta", "Psat"), T=convert(temperature, "K"), P=pressure)
# g from CHNOSZ functions
Rg <- CHNOSZ:::gfun(w$rho/1000, temperature, pressure, w$alpha, w$daldT, w$beta)$g
# g from R translation of DEW macro functions (not used in CHNOSZ)
DEWg <- CHNOSZ:::calculateG(pressure, temperature, w$rho/1000)
expect_equal(Rg, DEWg, info = info)

## The following tests use reference values in calories
E.units("cal")

info <- "Gibbs energies of species are calculated correctly"
P <- c(5000, 5000, 10000, 10000, 20000, 20000, 50000, 50000)
T <- c(100, 1000, 100, 1000, 100, 1000, 100, 1000)
RG_HCl <- subcrt("HCl", P=P, T=T)$out[[1]]$G
DEWG_HCl <- c(-28784.99, -58496.85, -26520.94, -55276.92, -21928.89, -50337.19, -8014.34, -36746.87)
expect_equal(RG_HCl, DEWG_HCl, tolerance = 1e-5, info = info)
RG_Cl <- subcrt("Cl-", P=P, T=T)$out[[1]]$G
DEWG_Cl <- c(-30054.59, -22839.35, -27910.68, -28094.07, -23568.45, -27959.67, -10443.07, -18744.93)
expect_equal(RG_Cl, DEWG_Cl, tolerance = 1e-7, info = info)

info <- "Delta G, logK, and Delta V of reactions are calculated correctly"
# These are reactions corresponding to Fig. 1b of Sverjensky et al., 2014 (Nat. Geosci.).
# The properties are calculated using parameters from the DEW spreadsheet,
# which are not necessarily identical those that were used for the paper.
T <- 600
P <- c(5000, 50000)
R1 <- subcrt(c("H2O", "CO2", "HCO3-", "H+"), c(-1, -1, 1, 1), T=T, P=P)$out
R2 <- subcrt(c("HCO3-", "CO3-2", "H+"), c(-1, 1, 1), T=T, P=P)$out
R3 <- subcrt(c("acetic acid", "acetate", "H+"), c(-1, 1, 1), T=T, P=P)$out
R4 <- subcrt(c("H2O", "CO2", "acetic acid", "oxygen"), c(-2, -2, 1, 2), T=T, P=P)$out
R5 <- subcrt(c("oxygen", "CH4", "acetic acid", "H2O"), c(-2, -2, 1, 2), T=T, P=P)$out

# Delta G values calculated using the DEW spreadsheet (May 2017 version)
DEW_DG <- c(38267.404507442814, 14893.170655988564,   # R1
            41407.05995898576,  21347.599525026497,   # R2
            28109.598104640143, 16112.928184675075,   # R3
            186960.22705581048, 133640.9631638353,    # R4
            -141552.60059404257, -134279.54670605875) # R5
# the aqueous-only reactions
expect_equal(c(R1$G, R2$G, R3$G), DEW_DG[1:6], tolerance = 1e-7, info = info)
# note that there is a small error for oxygen in the DEW spreadsheet (wrong c parameter),
# but even so, these tests have a lower tolerance because of the larger magnitude of the difference
expect_equal(c(R4$G, R5$G), DEW_DG[7:10], tolerance = 1e-9, info = info)

# logK values calculated using the DEW spreadsheet
DEW_logK <- c(-9.58455651110442, -3.7301833667366027,
              -10.370923015131565, -5.346776889042665,
              -7.040405143911882, -4.035687100632826, 
              -46.826558649850625, -33.47207316851283,
              35.45364304557209, 33.632014510923746) 
expect_equal(c(R1$logK, R2$logK, R3$logK, R4$logK, R5$logK), DEW_logK, tolerance = 1e-3, info = info)

# Delta V values calculated using the DEW spreadsheet
DEW_DV <- c(-45.26925983499276, -14.640599169742725,
            -47.95180846799733, -9.469432706749927, 
            -27.042087540311922, -6.836267057722694,
            -18.1550937649195, 5.513800665649967,
            -37.37077435045512, -45.08570662275889)
# for the aqueous species we're getting very close results
# (at P=5000 bar this depends on calculating drhodP -> beta -> dgdP -> dwdP -> V correctly, which is not tested above)
expect_equal(c(R1$V, R2$V, R3$V), DEW_DV[1:6], tolerance = 1e-15, info = info)
# TODO: why does DEW spreadsheet use V (O2,g) == 24.465?
#expect_equal(c(R4$V, R5$V), DEW_DV[7:10])

info <- "Calculated logK values are consistent with Extended Deep Earth Water paper"
# Reference logK values are from Appendix D of Huang and Sverjensky, 2019 (doi:10.1016/j.gca.2019.03.027)
# Select T and P for comparisons
T <- c(300, 1100)
P <- c(10000, 60000)
# Adjust calculated logKs for different values of the
# gas constant used in the DEW spreadsheet and CHNOSZ
RoverR <- 1.9858775 / 1.9872
# Calculate logK for each reaction
logK1 <- subcrt(c("H2O", "CO2", "H2CO3"), c(-1, -1, 1), T = T, P = P)$out$logK / RoverR
logK2 <- subcrt(c("AlO2(SiO2)-", "H+", "Al+3", "H2O", "SiO2"), c(-1, -4, 1, 2, 1), T = T, P = P)$out$logK / RoverR
logK3 <- subcrt(c("Ca(HCO3)+", "Ca+2", "HCO3-"), c(-1, 1, 1), T = T, P = P)$out$logK / RoverR
logK4 <- subcrt(c("Ca(HCOO)+", "Ca+2", "HCOO-"), c(-1, 1, 1), T = T, P = P)$out$logK / RoverR
logK5 <- subcrt(c("Ca(HSiO3)+", "H+", "Ca+2", "H2O", "SiO2"), c(-1, -1, 1, 1, 1), T = T, P = P)$out$logK / RoverR
logK6 <- subcrt(c("Fe(HCOO)+", "Fe+2", "HCOO-"), c(-1, 1, 1), T = T, P = P)$out$logK / RoverR
logK7 <- subcrt(c("Fe(HSiO3)+", "H+", "Fe+2", "H2O", "SiO2"), c(-1, -1, 1, 1, 1), T = T, P = P)$out$logK / RoverR
logK8 <- subcrt(c("H2O", "SiO2", "HSiO3-", "H+"), c(-1, -1, 1, 1), T = T, P = P)$out$logK / RoverR
logK9 <- subcrt(c("MgO", "H+", "Mg+2", "H2O"), c(-1, -2, 1, 1), T = T, P = P)$out$logK / RoverR
logK10 <- subcrt(c("Mg(SiO2)(HCO3)+", "Mg+2", "SiO2", "HCO3-"), c(-1, 1, 1, 1), T = T, P = P)$out$logK / RoverR
logK11 <- subcrt(c("NaHCO3", "Na+", "HCO3-"), c(-1, 1, 1), T = T, P = P)$out$logK / RoverR
logK12 <- subcrt(c("Si3O6", "SiO2"), c(-1, 3), T = T, P = P)$out$logK / RoverR
# Make the comparisons:
# - Reference logK values are from Appendix D of Huang and Sverjensky, 2019 (doi:10.1016/j.gca.2019.03.027)
# - Tolerance (tol) is added if there are differences between the reference and calculated values after rounding
#   (scale = 1 is used to compare absolute differences)
expect_equal(round(logK1,  2), c(-0.76, 1.34), info = info)
expect_equal(round(logK2,  4), c( 6.7755,  0.9413), tol = 0.003,  scale = 1, info = info)
expect_equal(round(logK3,  4), c(-2.3759, -5.6840), info = info)
expect_equal(round(logK4,  4), c(-2.2837, -5.5822), tol = 0.0001, scale = 1, info = info)
expect_equal(round(logK5,  4), c( 4.0349,  0.1797), tol = 0.002,  scale = 1, info = info)
expect_equal(round(logK6,  4), c(-7.5354, -8.0238), info = info)
expect_equal(round(logK7,  4), c(-0.6883, -1.6363), tol = 0.002,  scale = 1, info = info)
expect_equal(round(logK8,  4), c(-7.0651, -5.5067), tol = 0.002,  scale = 1, info = info)
expect_equal(round(logK9,  4), c( 8.2759,  4.3493), tol = 0.002,  scale = 1, info = info)
expect_equal(round(logK10, 4), c(-6.8106, -8.5888), tol = 0.0001, scale = 1, info = info)
expect_equal(round(logK11, 4), c(-0.2447, -2.9235), info = info)
expect_equal(round(logK12, 4), c( 3.3283,  0.4527), info = info)
