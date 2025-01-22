# Test added 20250122
info <- "Parameters are fit correctly for a mineral"
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

