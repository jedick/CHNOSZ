# Load default settings for CHNOSZ
reset()
# Set subcrt() to output values in Joules
E.units("J")

# 20220206
info <- "Database is set up correctly"
# Activate the AkDi model for aqueous nonelectrolytes
iAkDi <- add.OBIGT("AkDi")
# Get the indices of the modified aqueous species
iaq <- iAkDi[info(iAkDi)$state == "aq"]
# Each aqueous species should be associated with the AkDi model
expect_equal(unique(info(iaq)$abbrv), "AkDi", info = info)
# Each aqueous species should have a gaseous counterpart
igas <- info(info(iAkDi)$name, "gas")
expect_true(!any(is.na(igas)), info = info)

# First version of tests used a circular reference (values previously calculated in CHNOSZ) 20190220
# Tests now use values calculated with the AD_full program provided by N. Akinfiev and E. Bastrakov 20210206
info <- "AkDi produces correct values for CO2 along saturation curve"
P <- "Psat"
T <- seq(50, 350, 100)
# J mol-1
G_ref <- c(-389128.8, -405167.1, -425423.6, -450572.9)
# J K-1 mol-1
S_ref <- c(135.082, 183.295, 226.32, 366.626)
# Calculate values using AkDi model in CHNOSZ
sout <- subcrt("CO2", T = T, P = P)$out[[1]]
expect_equal(sout$G, G_ref, tolerance = 13, scale = 1, info = info)
expect_equal(sout$S, S_ref, tolerance = 0.8, scale = 1, info = info)

info <- "AkDi produces correct values for CO2 at 1000 bar"
P <- 1000
T <- c(300, 400, 500, 600)
G_ref <- c(-431672.1, -455065.1, -480692, -507005.3)
S_ref <- c(221.255, 246.139, 263.759, 259.812)
# Calculate values using AkDi model in CHNOSZ
sout <- subcrt("CO2", T = T, P = P)$out[[1]]
expect_equal(sout$G, G_ref, tolerance = 11, scale = 1, info = info)
expect_equal(sout$S, S_ref, tolerance = 0.09, scale = 1, info = info)

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
expect_equal(CHNOSZ:::.drho1_dT(298.15, 1),   -0.0002571,   tolerance = 0.000002, scale = 1, info = info)
# g / cm3 / bar
expect_equal(CHNOSZ:::.drho1_dP(298.15, 1),    0.00004511,  tolerance = 0.00000003, scale = 1, info = info)
# g / cm3 / K^2
expect_equal(CHNOSZ:::.d2rho1_dT2(298.15, 1), -0.000009503, tolerance = 0.0000014, scale = 1, info = info)

# 20190220 Compare Gibbs energies at 25 degrees calculated with AkDi model to default OBIGT database
info <- "Gibbs energies at 25 degree C are comparable between AkDi and default OBIGT database"
# Remove some hydroxides that aren't in the default database
iaq <- iaq[!info(iaq)$name %in% c("Si(OH)4", "As(OH)3")]
# This would produce an error if any calculations failed
# (e.g. because gases corresponding to any aqueous species were unavailable)
sAkDi <- subcrt(iaq, T = 25)
GAkDi <- do.call(rbind, sAkDi$out)$G
# Now calculate the parameters using default OBIGT database
OBIGT()
sOBIGT <- subcrt(iaq, T = 25)
GOBIGT <- do.call(rbind, sOBIGT$out)$G
# Mean absolute differences are less than 280 J / mol
expect_equal(GAkDi, GOBIGT, tolerance = 280, scale = 1, info = info)
# The largest differences are for HCl, ethane, and B(OH)3
expect_equal(sort(info(iaq[abs(GAkDi - GOBIGT) > 900])$name), sort(c("HCl", "ethane", "B(OH)3")))
