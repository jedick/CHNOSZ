# Load default settings for CHNOSZ
reset()

# initial setup
suppressMessages({
  # for numerical reproducibility, use the now-superseded properties of glycine 20190208
  mod.OBIGT("glycine", G = -90950, H = -124780, S = 39.29)
  basis("CHNOS")
  basis("O2", -65)
  species(c("leucine", "glycine", "glutamic acid"))
  # a 0-D case
  a <- affinity()
  e0 <- equilibrate(a)
  # a 1-D case
  a <- affinity(O2 = c(-75, -65))
  e1 <- equilibrate(a)
  # a 2-D case
  a <- affinity(O2 = c(-75, -65), NH3 = c(-4, 100, 3))
  e2 <- equilibrate(a)
  # a 3-D case
  a <- affinity(O2 = c(-75, -65, 4), H2O = c(-8, 0, 3), NH3 = c(-6, -4, 2))
  e3 <- equilibrate(a)
})

info <- "Inconsistent arguments produce an error"
expect_error(CHNOSZ:::get.objfun("affinity"), "affinity is not a function with an attribute named 'optimum'", info = info)
expect_error(revisit(list(1, 1), plot.it = TRUE), "can't make a plot if 'eout' is not the output from equilibrate\\(\\)", info = info)
expect_error(revisit(list(1, c(1, 2))), "the list provided in 'eout' is not the output from equilibrate\\(\\)", info = info)
expect_error(revisit(e1, "RMSD"), "loga2 must be supplied for RMSD", info = info)
expect_error(revisit(e1, "RMSD", list(1, 1)), "loga2 has different length \\(2\\) than list in eout \\(3\\)", info = info)
# commented because a plot is still initialized ...
#expect_error(revisit(e2, "CV", style.2D = "xxx"), "2D plot style xxx not one of 'contour' or 'image'", info = info)

info <- "0-D, 1-D, 2-D and 3-D calculations give identical results at the same conditions"
r0.qqr <- revisit(e0, "qqr", plot.it = FALSE)
r1.qqr <- revisit(e1, "qqr", plot.it = FALSE)
r2.qqr <- revisit(e2, "qqr", plot.it = FALSE)
r3.qqr <- revisit(e3, "qqr", plot.it = FALSE)
# check that we get the same values
expect_equal(c(r0.qqr$H), c(tail(r1.qqr$H, 1)), info = info)
expect_equal(c(r0.qqr$H), r3.qqr$H[4, 3, 2], info = info)
# check that we get the same index and same optimum
expect_equal(r1.qqr$ixopt, r2.qqr$ixopt, check.attributes = FALSE, info = info)
expect_equal(r1.qqr$optimum, r2.qqr$optimum, info = info)

info <- "Non-referenced objectives give expected results"
# the non-referenced objectives use only the logarithms of activities in eout (loga1)
r1.cv <- revisit(e1, "CV", plot.it = FALSE)
r1.sd <- revisit(e1, "SD", plot.it = FALSE)
r1.shannon <- revisit(e1, "shannon", plot.it = FALSE)
r1.qqr <- revisit(e1, "qqr", plot.it = FALSE)
# the tests will alert us to significant numerical changes
# but so far haven't been independently verified
expect_equal(r1.cv$optimum, 0.29927, tolerance = 1e-5, info = info) 
expect_equal(r1.sd$optimum, 0.00027671, tolerance = 1e-5, info = info) 
expect_equal(r1.shannon$optimum, 1.067228, tolerance = 1e-5, info = info)
expect_equal(r1.qqr$optimum, 0.9999584, tolerance = 1e-5, info = info)

info <- "Referenced objectives give expected results"
# the referenced objectives compare the logarithms of activities (loga1) to reference values (loga2)
# the spearman correlation coefficient
r1.spearman <- revisit(e1, "spearman", c(1, 2, 3), plot.it = FALSE)
expect_equal(c(head(r1.spearman$H, 1)), -1, info = info) # perfect anti-rank correlation
expect_equal(max(r1.spearman$H), 1, info = info)  # perfect rank correlation
# where logarithm of activity of the 3rd species (glutamic acid) maximizes
r1.logact <- revisit(e1, "logact", 3, plot.it = FALSE)
expect_equal(r1.logact$ixopt, 141, info = info)

info <- "DGtr objective gives zero at equilibrium and >0 not at equilibrium"
# let's use n-alkanes
basis(c("CH4", "H2"), c("gas", "gas"))
species(c("methane", "ethane", "propane", "butane"), "liq")
# calculate equilibrium distribution over a range of logaH2
a1 <- affinity(H2 = c(-10, -5, 101), exceed.Ttr = TRUE)
e1 <- equilibrate(a1)
# take the equilibrium distribution at logfH2  =  -7.5 as the reference distribution
loga2 <- list2array(e1$loga.equil)[51, ]
# calculate the DGtr/RT relative to the reference distribution
r1 <- revisit(e1, "DGtr", loga2 = loga2, plot.it = FALSE)
# we should find a minimum of zero at logfH2 = -7.5
expect_equal(min(r1$H), 0, info = info)
expect_equal(r1$xopt, -7.5, info = info)

# we can even go into 2 dimensions
# (it's a slightly longer test, so don't run it on CRAN)
if(at_home()) {
  a2 <- affinity(H2 = c(-10, -5, 101), T = c(0, 100, 101), exceed.Ttr = TRUE)
  e2 <- equilibrate(a2)
  r2 <- revisit(e2, "DGtr", loga2 = loga2, plot.it = FALSE)
  # we should DGtr = 0 at the temperature of the reference distribution (25 degC)
  expect_equal(min(r2$H), 0, info = info)
  expect_equal(r2$yopt, 25, info = info)
}
