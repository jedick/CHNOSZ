# Load default settings for CHNOSZ
reset()

info <- "[P|T|E].units() do not accept invalid units"
expect_error(P.units("X"), "units of pressure must be either bar or MPa", info = info)
expect_error(T.units("X"), "units of temperature must be either C or K", info = info)
expect_error(E.units("X"), "units of energy must be either cal or J", info = info)

# Tests added on 20260111

info <- "Error with non-character units specification"
expect_error(convert(25, 1), "please specify a character argument", info = info)

info <- "Expected results with solubility() output"

# Solubility of gaseous SO2
basis(c("sulfur", "oxygen", "H2O", "H+"))
basis("O2", -56)
basis("pH", 6)
species("sulfur dioxide", -20)
iaq <- retrieve("S", c("O", "H"), "aq")
s <- solubility(iaq, T = 125, in.terms.of = "S")

expect_equal(round(convert(s, "ppt")$loga.balance[[1]]), 3, info = info)
expect_equal(round(convert(s, "ppm")$loga.balance[[1]]), 3054, info = info)
expect_equal(round(convert(s, "ppb")$loga.balance[[1]]), 3054344, info = info)
expect_equal(round(convert(s, "logppb")$loga.balance[[1]]), round(log10(3054344)), info = info)
expect_error(convert(s, 1), "please specify a character argument", info = info)
expect_error(convert(s, "ppx"), "units ppx not available", info = info)
expect_error(convert(s[-1], "ppb"), "is not the output from solubility()", info = info)
