# Test added 20250122

## Mineral species
info <- "Mineral species: parameters are fit correctly"
species <- c("iron", "sulfur", "FeS2")
coeffs <- c(-1, -2, 1)
P <- 1000
T <- seq(120, 400, 10)
# Use calculated logK for formation of pyrite
logK <- subcrt(c("iron", "sulfur", "pyrite"), c(-1, -2, 1), T = T, P = P)$out$logK
# Use V for pyrite to reproduce value of G at 25 degC
inew <- logK.to.OBIGT(logK, species, coeffs, T, P, name = "newpyrite", state = "cr", V = info(info("pyrite"))$V, npar = 5)
# Get parameters of pyrite (from OBIGT database) and newpyrite (created by logK.to.OBIGT)
oldpar <- info(info("pyrite"))
newpar <- info(info("newpyrite"))
expect_equal(convert(oldpar$G, "J"), newpar$G, info = info)
expect_equal(convert(oldpar$S, "J"), newpar$S, info = info)
expect_equal(convert(oldpar$a, "J"), newpar$a, info = info)
expect_equal(convert(oldpar$b, "J"), newpar$b, info = info)
expect_equal(convert(oldpar$c, "J"), newpar$c, info = info)

# Tests added on 20260108

info <- "Mineral species: no error for npar = 2 to 4"
expect_silent(logK.to.OBIGT(logK, species, coeffs, T, P, name = "newpyrite", state = "cr", V = info(info("pyrite"))$V, npar = 4), info = info)
expect_silent(logK.to.OBIGT(logK, species, coeffs, T, P, name = "newpyrite", state = "cr", V = info(info("pyrite"))$V, npar = 3), info = info)
expect_silent(logK.to.OBIGT(logK, species, coeffs, T, P, name = "newpyrite", state = "cr", V = info(info("pyrite"))$V, npar = 2), info = info)

info <- "Mineral species: error for npar = 1"
expect_error(logK.to.OBIGT(logK, species, coeffs, T, P, name = "newpyrite", state = "cr", V = info(info("pyrite"))$V, npar = 1), info = info)

## Aqueous species
# CoHS+ from Migdisov et al. (2011)
logK <- c(6.24, 6.02, 5.84, 5.97, 6.52)
T <- c(120, 150, 200, 250, 300)
P <- "Psat"
species <- c("Co+2", "HS-", "CoHS+")
coeffs <- c(-1, -1, 1)

info <- "Aqueous species: optimize.omega is forced to FALSE for insufficient T or npar values"
expect_message(logK.to.OBIGT(logK[1:4], species, coeffs, T[1:4], P, npar = 5, optimize.omega = TRUE), "forcing optimize.omega = FALSE", info = info)
expect_message(logK.to.OBIGT(logK, species, coeffs, T, P, npar = 4, optimize.omega = TRUE), "forcing optimize.omega = FALSE", info = info)

info <- "No error for npar = 4 and optimize.omega = FALSE"
expect_silent(logK.to.OBIGT(logK, species, coeffs, T, P, npar = 3, optimize.omega = FALSE), info = info)

info <- "Error for exceeding default tolerance"
expect_error(logK.to.OBIGT(logK, species, coeffs, T, P, npar = 2, optimize.omega = FALSE), info = info)
expect_error(logK.to.OBIGT(logK, species, coeffs, T, P, npar = 1, optimize.omega = FALSE), info = info)

info <- "No error for optimize.omega = TRUE"
expect_silent(inew <- logK.to.OBIGT(logK, species, coeffs, T, P, npar = 5, optimize.omega = TRUE), info = info)

# NOTE: optimize.omega technically works for this small dataste, but is massively overfitting:
## Plot experimental logK
#plot(T, logK, "n", c(100, 320), c(5.8, 6.8), xlab = axis.label("T"), ylab = quote(log~beta))
#points(T, logK, pch = 19, cex = 2)
## Plot calculated values
#Tfit <- seq(100, 320, 10)
#sres <- subcrt(species, coeffs, T = Tfit)
#lines(sres$out$T, sres$out$logK, col = 4)

info <- "Errors for other conditions"
expect_error(logK.to.OBIGT(logK, species, coeffs, T, P, npar = 6, state = "aq"), "invalid value for npar", info = info)
expect_error(logK.to.OBIGT(logK, species, coeffs, T, P, npar = 6, state = "cr"), "invalid value for npar", info = info)
expect_error(logK.to.OBIGT(logK, species, coeffs, T, P, npar = 6, state = "xx"), "Unrecognized state", info = info)
