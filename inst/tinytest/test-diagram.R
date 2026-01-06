# Load default settings for CHNOSZ
reset()

# Save all plots to a temporary PNG file
pngfile <- tempfile(fileext = ".png")
png(pngfile)

info <- "Expected errors are produced for inconsistent arguments"
expect_error(diagram(list()), "'eout' is not the output from", info = info)
basis("CHNOS")
species(c("glycine", "alanine"))
a <- affinity()
expect_message(diagram(a), "balance: on moles of CO2 in formation reactions", info = info)
e <- equilibrate(a)
expect_error(diagram(e, "Z"), "Z is not a valid diagram type", info = info)

info <- "Expected messages, errors and results arise using output from affinity()"
basis("CHNOS+")
# Fugacity of O2 is buffered here
basis("O2", "CO2-AC")
species(c("formic acid", "formate", "acetic acid", "acetate"))
# 0-D
a <- affinity()
# Equilibrium activities are not possible here
expect_error(diagram(a, "loga.equil"), "'eout' is not the output from equil\\(\\)", info = info)
# We can't calculate the equilibrium activity of a basis species if it's externally buffered
expect_error(diagram(a, "O2"), "is not numeric - was a buffer selected\\?", info = info)
# This one works - a barplot of A/2.303RT
expect_message(diagram(a), "balance: on moles of CO2 in formation reactions", info = info)
# If we're plotting A/2.303RT the values can be divided by balancing coefficient or not
expect_silent(d.1 <- diagram(a, balance = 1), info = info)
expect_silent(d.CO2 <- diagram(a), info = info)
expect_equal(as.numeric(d.CO2$plotvals), as.numeric(d.1$plotvals)/c(1, 1, 2, 2), info = info)
# Now run the calculation over a range of O2 
basis("O2", -90)
# 1-D
a <- affinity(O2 = c(-80, -70))
# Ask for the equilibrium activity of CO2
expect_error(diagram(a, "CO2", groups = list(1:2, 3:4)), "can't plot equilibrium activities of basis species for grouped species", info = info)
expect_error(diagram(a, "CO2", alpha = TRUE), "equilibrium activities of basis species not available with alpha = TRUE", info = info)
expect_silent(d <- diagram(a, "CO2"), info = info)
# Test that the result does in fact correspond to zero affinity of formation, how about for acetate?
a <- affinity(O2 = d$vals[[1]], CO2 = d$plotvals[[4]])
expect_equal(a$values[[4]], array(numeric(256)), info = info)

info <- "'groups' and 'alpha' work as expected"
basis("CHNOS+")
species(c("formic acid", "formate", "acetic acid", "acetate"))
# 1-D
a <- affinity(O2 = c(-80, -60))
e <- equilibrate(a)
# Group the species together
d <- diagram(e, groups = list(1:2, 3:4))
# We should find that their activities have been multiplied by the balance coefficients and summed
n.balance <- CHNOSZ:::balance(a)$n.balance
expect_equal(d$plotvals[[1]], log10(n.balance[1]*10^e$loga.equil[[1]] + n.balance[2]*10^e$loga.equil[[2]]), info = info)
expect_equal(d$plotvals[[2]], log10(n.balance[3]*10^e$loga.equil[[3]] + n.balance[4]*10^e$loga.equil[[4]]), info = info)
# Ask for degrees of formation instead of logarithms of activities
d <- diagram(e, alpha = TRUE)
# We should find that the sum of alphas is one
expect_equal(Reduce("+", d$plotvals), array(rep(1, 256)), check.attributes = FALSE, info = info)

info <- "'normalize' and 'as.residue' work as expected"
basis("CHNOS")
species(c("LYSC_CHICK", "MYG_PHYCA", "RNAS1_BOVIN", "CYC_BOVIN"))
# 1-D
a <- affinity(O2 = c(-80, -70))
expect_error(diagram(a, normalize = TRUE), "can be TRUE only for a 2-D \\(predominance\\) diagram", info = info)
# 2-D
a <- affinity(H2O = c(-10, 0), O2 = c(-80, -70))
d1 <- diagram(a, normalize = TRUE)
e <- equilibrate(a, normalize = TRUE)
expect_silent(d2 <- diagram(e), info = info)
expect_equal(d1$predominant, d2$predominant, info = info)
expect_error(diagram(e, normalize = TRUE), "can be TRUE only if 'eout' is the output from affinity\\(\\)", info = info)
expect_silent(d3 <- diagram(a, as.residue = TRUE), info = info)
expect_silent(e <- equilibrate(a, as.residue = TRUE), info = info)
expect_silent(d4 <- diagram(e), info = info)
expect_equal(d3$predominant, d4$predominant, info = info)

info <- "NaN values from equilibrate() are preserved (as NA in predominance calculation)"
# Example provided by Grayson Boyer 20170411
basis(c("H2", "O2", "CO2"), c(-7.19, -60, -2.65))
species(c("n-hexadecanol", "n-hexadecanoic acid", "n-octadecanol", "n-octadecanoic acid"), c("liq", "liq", "liq", "liq"))
a <- affinity("H2" = c(-12, 0), "O2" = c(-90, -50), T = 30)
e <- equilibrate(a, balance = 1)
d <- diagram(e)
# equilibrate() here with default "boltzmann" method produces
# NaN at very high O2 + low H2 or very low O2 + high H2 
expect_equal(d$predominant[1, 256], as.numeric(NA), info = info)
expect_equal(d$predominant[256, 1], as.numeric(NA), info = info)

## TODO: Exclude this test for now because plot.it = FALSE doesn't produce values for namesx 20190223
#info <- "labels are dropped outside of xlim and ylim ranges"
#basis(c("Fe", "O2", "S2"))
#species(c("iron", "ferrous-oxide", "magnetite",
#  "hematite", "pyrite", "pyrrhotite"))
#a <- affinity(S2 = c(-50, 0), O2 = c(-90, -10), T = 200)
## total range: all species are present
#d <- diagram(a, fill = "heat", xlim = NULL, ylim = NULL)
#expect_equal(sum(is.na(d$namesx)), 0, info = info)
## reduce y-range to exclude hematite
#d <- diagram(a, fill = "heat", xlim = NULL, ylim = c(-90, -50))
#expect_equal(sum(is.na(d$namesx)), 1, info = info)
## reduce x-range to exclude pyrite
#d <- diagram(a, fill = "heat", xlim = c(-50, -20), ylim = c(-90, -50))
#expect_equal(sum(is.na(d$namesx)), 2, info = info)

info <- "P-T diagram has expected geometry"
# Modified from kayanite-sillimanite-andalusite example in ?diagram 20200811
basis(c("corundum", "quartz", "oxygen"))
species(c("kyanite", "sillimanite", "andalusite"))
a <- affinity(T = c(200, 900, 50), P = c(0, 9000, 51), exceed.Ttr = TRUE)
expect_silent(d <- diagram(a), info = info)
expect_equal(species()$name[d$predominant[1, 1]], "andalusite", info = info)
expect_equal(species()$name[d$predominant[1, 51]], "kyanite", info = info)
expect_equal(species()$name[d$predominant[50, 51]], "sillimanite", info = info)
# The location of the triple point - the algorithm gives an ambiguous location within a certain range
tp <- find.tp(d$predominant)
expect_equal(range(tp[, 1]), c(30, 32))
expect_equal(range(tp[, 2]), c(22, 23))

info <- "diagram(type = .) and affinity(return.buffer = TRUE) give the same results"
# Extracted from ?buffer 20200811
O2 <- c(-85, -70, 4)
T <- c(25, 100, 4)
basis("CHNOS")
basis("CO2", 999)
species("acetic acid", -3)
species(1, -10)
a <- affinity(O2 = O2, T = T)
d <- diagram(a, type = "CO2")
and <- as.numeric(d$plotvals[[1]])
# Now do the calculation with affinity(return.buffer = TRUE)
basis("CO2", "AC")
mod.buffer("AC", logact = -10)
a.buffer <- affinity(O2 = O2, T = T, return.buffer = TRUE)
ana <- as.numeric(unlist(a.buffer[[1]]))
expect_equal(ana, and, info = info)

# Tests added on 20260106
info <- "Function works with 'dotted', 'lty.aq', and 'lty.cr' arguments"
basis(c("Fe", "H2O", "H+", "e-"))
species(c("Fe+2", "Fe+3", "magnetite", "hematite"))
a <- affinity(pH = c(0, 12, 50), Eh = c(-1, 1, 50))
expect_silent(diagram(a, dotted = 3), info = info)
expect_silent(diagram(a, lty.cr = 2, lty.aq = 3), info = info)

# Close the graphics device and remove the temporary PNG file
dev.off()
file.remove(pngfile)
