# Load default settings for CHNOSZ
reset()
# Set subcrt() to output values in Joules
E.units("J")

# 20220206
info <- "Database is set up correctly"
# Activate the AD model for aqueous nonelectrolytes
iAD <- add.OBIGT("AD")
# Get the indices of the modified aqueous species
iaq <- iAD[info(iAD)$state == "aq"]
# Each aqueous species should be associated with the AD model
expect_equal(unique(info(iaq)$abbrv), "AD", info = info)
# Each aqueous species should have a gaseous counterpart
igas <- info(info(iAD)$name, "gas")
expect_true(!any(is.na(igas)), info = info)

# First version of tests used a circular reference (values previously calculated in CHNOSZ) 20190220
# Tests now use values calculated with the AD_full program provided by N. Akinfiev and E. Bastrakov 20210206
info <- "AD produces correct values for CO2 along saturation curve"
P <- "Psat"
T <- c(50, 150, 250, 350)
# J mol-1
G_ref <- c(-389128.8, -405167.1, -425423.6, -450572.9)
# J K-1 mol-1
S_ref <- c(135.082, 183.295, 226.32, 366.626)
Cp_ref <- c(189.791, 178.452, 273.662, 7231.723)
V_ref <- c(32.57, 38.125, 58.269, 298.659)
# Calculate values using AD model in CHNOSZ
sout <- subcrt("CO2", T = T, P = P)$out[[1]]
expect_equal(sout$G, G_ref, tolerance = 13, scale = 1, info = info)
expect_equal(sout$S, S_ref, tolerance = 0.8, scale = 1, info = info)
expect_equal(sout$Cp[1:3], Cp_ref[1:3], tolerance = 3, scale = 1, info = info)
expect_equal(sout$V[1:3], V_ref[1:3], tolerance = 0.11, scale = 1, info = info)
# Cp and V get much larger, and so do the differences, near the critical point
expect_equal(sout$Cp[4], Cp_ref[4], tolerance = 800, scale = 1, info = info)
expect_equal(sout$V[4], V_ref[4], tolerance = 24, scale = 1, info = info)

info <- "AD produces correct values for CO2 at 1000 bar"
P <- 1000
T <- c(300, 400, 500, 600)
G_ref <- c(-431672.1, -455065.1, -480692, -507005.3)
S_ref <- c(221.255, 246.139, 263.759, 259.812)
Cp_ref <- c(150.613, 157.014, 54.489, -76.384)
V_ref <- c(45.831, 67.408, 107.656, 129.264)
# Calculate values using AD model in CHNOSZ
sout <- subcrt("CO2", T = T, P = P)$out[[1]]
expect_equal(sout$G, G_ref, tolerance = 11, scale = 1, info = info)
expect_equal(sout$S, S_ref, tolerance = 0.1, scale = 1, info = info)
expect_equal(sout$Cp, Cp_ref, tolerance = 14, scale = 1, info = info)
expect_equal(sout$V, V_ref, tolerance = 0.4, scale = 1, info = info)

info <- "AD gives consistent values of G, H, and S"
T_in_Kelvin <- convert(T, "K")
S_of_elements_in_Joules <- convert(entropy("CO2"), "J")
expect_equal(sout$H - T_in_Kelvin * sout$S + 298.15 * S_of_elements_in_Joules, sout$G)

# 20220206
info <- "Fugacity, density, and density derivatives of H2O are close to values in Akinfiev and Diamond (2003)"
# This is the value of fugacity given in the paper
expect_equal(log(CHNOSZ:::.f1(298.15, 1, TRUE)), -3.4773, tolerance = 0.026, scale = 1, info = info)
# This is the value of fugacity output by AD_full
expect_equal(log(CHNOSZ:::.f1(298.15, 1, TRUE)), -3.4524, tolerance = 0.0007, scale = 1, info = info)
# The following values are from the paper
# g / cm3
expect_equal(CHNOSZ:::.rho1(298.15, 1),        0.9970,      tolerance = 0.0001, scale = 1, info = info)
# g / cm3 / K
expect_equal(CHNOSZ:::.drho1_dT(298.15, 1, FALSE),   -0.0002571,   tolerance = 0.000002, scale = 1, info = info)
# g / cm3 / bar
expect_equal(CHNOSZ:::.drho1_dP(298.15, 1, FALSE),    0.00004511,  tolerance = 0.0000001, scale = 1, info = info)
# g / cm3 / K^2
expect_equal(CHNOSZ:::.d2rho1_dT2(298.15, 1, FALSE), -0.000009503, tolerance = 0.000007, scale = 1, info = info)

# 20190220 Compare Gibbs energies at 25 degrees calculated with AD model to default OBIGT database
info <- "Gibbs energies at 25 degree C are comparable between AD and default OBIGT database"
# Remove some hydroxides that aren't in the default database
iaq <- iaq[!info(iaq)$name %in% c("Si(OH)4", "As(OH)3")]
# This would produce an error if any calculations failed
# (e.g. because gases corresponding to any aqueous species were unavailable)
sAD <- subcrt(iaq, T = 25)
GAD <- do.call(rbind, sAD$out)$G
# Now calculate the parameters using default OBIGT database
OBIGT()
sOBIGT <- subcrt(iaq, T = 25)
GOBIGT <- do.call(rbind, sOBIGT$out)$G
# Mean absolute differences are less than 280 J / mol
expect_equal(GAD, GOBIGT, tolerance = 280, scale = 1, info = info)
# The largest differences are for HCl, ethane, and B(OH)3
expect_equal(sort(info(iaq[abs(GAD - GOBIGT) > 900])$name), sort(c("HCl", "ethane", "B(OH)3")))
