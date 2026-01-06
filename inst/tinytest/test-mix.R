# Load default settings for CHNOSZ
reset()

info <- "mix() and mosaic() yield same affinities for a bimetallic mineral"

# Test that we get the same values for affinity of a bimetallic mineral formed from different single-metal basis species using:
#   1. affinity() of the bimetallic mineral minus mix() for two single-metal systems in correct proportions
#   2. mosaic() with two sets of changing basis species corresponding to the single-metal systems
# 20200724

plot.it <- FALSE
res <- 50
## Uncomment to show plots
#plot.it <- TRUE
#par(mfrow = c(2, 2))

# Set up system
pH <- c(0, 7, res)
Eh <- c(-1, 1, res)
basis(c("copper", "iron", "H2S", "H2O", "e-", "H+"))
iFe.cr <- info(c("iron", "ferrous-oxide", "magnetite", "hematite"))
iFe.aq <- info(c("Fe+2", "FeOH+"))
iCu.cr <- info(c("copper", "cuprite", "tenorite"))
iCu.aq <- info(c("Cu+", "CuOH"))

# METHOD 1 (mix)

# Fe-O-H diagram
# Use logact = 0 for everything so we can compare results with mosaic()
species(c(iFe.aq, iFe.cr), 0)
# Increase the temperature so we get a ferrous oxide field
aFe <- affinity(pH = pH, Eh = Eh, T = 300)
dFe <- diagram(aFe, plot.it = plot.it)

# Cu-O-H diagram
species(c(iCu.aq, iCu.cr), 0)
aCu <- affinity(aFe)  # argument recall
dCu <- diagram(aCu, plot.it = plot.it)

# Combine the diagrams for a 1:5 mixture of Fe:Cu
aFeCu15 <- mix(dFe, dCu, parts = c(1, 5))
dFeCu15 <- diagram(aFeCu15, min.area = 0.01, plot.it = FALSE)
# Calculate affinity of bornite (Cu5FeS4)
species("bornite")
abornite <- affinity(aFe)  # argument recall
# Because of the 1) equal stoichiometry and 2) same basis species
# used for bornite and the mix() diagram, subtracting these makes sense
# 20201221 But now we have to multiply by 5 because predominant.values is divided by the balancing coefficients
abornite_vs_predominant.values <- abornite$values[[1]] - 5 * dFeCu15$predominant.values
if(plot.it) image(abornite_vs_predominant.values)

# METHOD 2 (mosaic)

# Make a mosaic diagram to calculate affinity of bornite with changing basis species
# (which correspond to the Fe-O-H and Cu-O-H diagrams)
iFe <- c(iFe.cr, iFe.aq)
iCu <- c(iCu.cr, iCu.aq)
# TODO: allow numeric values for bases in mosaic()
#mbornite <- mosaic(list(iFe, iCu), pH = pH, Eh = Eh, T = 300, predominant = list(dFe$predominant, dCu$predominant))
Fe <- info(iFe, check.it = FALSE)$name
Cu <- info(iCu, check.it = FALSE)$name
mbornite <- mosaic(list(Fe, Cu), pH = pH, Eh = Eh, T = 300, blend = FALSE)
if(plot.it) image(mbornite$A.species$values[[1]])

# The test: values calculated both ways should be equal
expect_equal(abornite_vs_predominant.values, mbornite$A.species$values[[1]], tol = 1e-5, scale = 1, info = info)

# Test added on 20250106
info <- "Using d3 argument to add bimetallic species produces no errors"
species("bornite")
a3 <- affinity(aFe)
d3 <- diagram(a3, plot.it = FALSE)
# 1:1 mixture (Fe:Cu)
expect_silent(a11 <- mix(dFe, dCu, d3, c(1, 1)), info = info)
# We could visualize it like this:
#diagram(a11)

# Test added on 20250106
info <- "maxh() and rebalance() produce no error"
basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
species(c("pyrite", "pyrrhotite", "magnetite", "hematite"))
aFe <- affinity("Fe+2" = c(0, 12, 20), O2 = c(-40, -16, 20))
dFe <- diagram(aFe)
species(c("covellite", "chalcocite", "tenorite", "cuprite"))
aCu <- affinity(aFe)  # argument recall
dCu <- diagram(aCu)
expect_silent(ac <- mash(dFe, dCu), info = info)
expect_silent(diagram(ac), info = info)
expect_silent(ad <- rebalance(dFe, dCu), info = info)
expect_silent(diagram(ad, balance = 1), info = info)

