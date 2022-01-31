# Load default settings for CHNOSZ
reset()

info <- "water.SUPCRT92() gives expected erros and warnings"
expect_error(water.SUPCRT92("X"), "not available: X", info = info)
expect_error(water.SUPCRT92("Psat"), "please set P='Psat'", info = info)

info <- "water.SUPCRT92() gives expected values for E and kT"
# E = V * alpha, kT = V * beta
expect_equal(water.SUPCRT92("E"), water.SUPCRT92("V") * water.SUPCRT92("alpha"), check.attributes = FALSE, info = info)
expect_equal(water.SUPCRT92("kT"), water.SUPCRT92("V") * water.SUPCRT92("beta"), check.attributes = FALSE, info = info)
# now compare with some real data from Tables V and VI of Fine and Millero, 1973
T <- convert(c(25, 100), "K")
P <- c(100, 1000)
expect_equal(water.SUPCRT92("beta", T, P)[, 1] * 1e6, c(44.100, 37.002), tolerance = 1e-2, info = info)
expect_equal(water.SUPCRT92("alpha", T, P)[, 1] * 1e6, c(268.06, 625.55), tolerance = 1e-2, info = info)

info <- "water.DEW() gives expected erros and messages"
expect_error(water.DEW("G", P = "Psat"), "Psat isn\'t supported", info = info)
expect_message(water.DEW("G", T = c(298.15, 373.15), P = c(1, 1000)), "using SUPCRT calculations", info = info)

# reference

# Fine, R. A. and Millero, F. J. (1973)
#   Compressibility of water as a function of temperature and pressure
#   J. Chem. Phys. 59, 5529--5536. https://doi.org/10.1063/1.1679903
