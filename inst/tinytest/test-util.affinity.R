# Load default settings for CHNOSZ
reset()

info <- "A.ionization() builds arrays with correct dimensions"
iprotein <- pinfo(c("LYSC_CHICK", "RNAS1_BOVIN"))
vals <- list(T = convert(c(25, 50, 100), "K"), O2 = c(-80, -70), pH = c(5, 6, 7, 8))
vars <- names(vals)
A <- CHNOSZ:::A.ionization(iprotein, vars, vals)
# For two proteins, a list of length 2
expect_equal(length(A), 2, info = info)
# For 3 variables with the values given to energy.args()
expect_equal(dim(A[[1]]), c(3, 2, 4), check.names = FALSE, info = info)

info <- "A.ionization() handles arguments generated by energy.args()"
# Specifically, when the user asks for Eh, energy.args() produces
# values of pe that are specified in all dimensions
# (because the conversion from Eh is a function of temperature)
# Here, the original length of the Eh variable is 3, so that should be retained by A.ionization
# We need a basis() for the function to poke around in 
basis("CHNOSe")
ea <- CHNOSZ:::energy.args(list(CO2 = c(-5, 0, 6), pH = c(0, 14, 2), Eh = c(-1, 0, 3), T = c(25, 28, 4)))
A <- CHNOSZ:::A.ionization(1, ea$vars, ea$vals)
expect_equal(dim(A[[1]]), c(6, 2, 3, 4), info = info)
# A simpler set of arguments, for testing dimensions with identical extents
# (the reason for the "for(dim in c(TPpHdim, otherdim))" in A.ionization())
ea <- CHNOSZ:::energy.args(list(pH = c(0, 14, 2), Eh = c(-1, 0, 2)))
A <- CHNOSZ:::A.ionization(1, ea$vars, ea$vals)
expect_equal(dim(A[[1]]), c(2, 2), info = info)
# An even simpler set of arguments: be able to handle a single point
A <- CHNOSZ:::A.ionization(1, vars = character(), vals = list())
expect_equal(dim(A[[1]]), 1, info = info)

info <- "energy() includes ionization affinities of proteins when needed"
# Get some results for a nonionized protein
basis("CHNOS")
species(c("alanine","LYSC_CHICK"))
# A 2x3 grid
ea <- CHNOSZ:::energy.args(list(T = c(25, 50, 2), O2 = c(-80, -60, 3)))
Anonion <- do.call(CHNOSZ:::energy, ea)
# Get some results for an ionized protein at pH 11
basis("CHNOS+")
basis("pH", 11)
species(c("alanine", "LYSC_CHICK"))
ea <- CHNOSZ:::energy.args(list(T = c(25, 50, 2), O2 = c(-80, -60, 3)))
Aion.pH11 <- do.call(CHNOSZ:::energy, ea)
# Alanine should not change
expect_equal(Anonion$a[[1]], Aion.pH11$a[[1]], info = info)
# What was the difference for LYSC_CHICK
A.ionization.LYSC_CHICK <- Aion.pH11$a[[2]] - Anonion$a[[2]]
# There should be 2 unique values (logfO2 doesn't matter)
# and they should be equal to the "manual" calculation using ionize.aa
A.ionization.2550 <- ionize.aa(pinfo(pinfo("LYSC_CHICK")), property = "A", T = c(25, 50), pH = 11)
expect_equal(A.ionization.LYSC_CHICK[, 1], A.ionization.2550[, 1], check.attributes = FALSE, info = info)
expect_equal(A.ionization.LYSC_CHICK[, 3], A.ionization.2550[, 1], check.attributes = FALSE, info = info)
