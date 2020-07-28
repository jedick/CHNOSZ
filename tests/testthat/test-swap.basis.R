context("swap.basis")

# clear out any previous basis definition or database alterations
suppressMessages(reset())

test_that("swap.basis raises errors when needed", {
  expect_error(swap.basis(c("CO2", "H2O")), "requires an existing basis definition")
  basis("CHNOS+")
  expect_error(swap.basis(c("CO2", "H2O")), "two species must be identified")
  expect_error(swap.basis(c("CO2", "H2O"), c("HCO3-", "H2O")), "can only swap one species for one species")
  expect_error(swap.basis("CH4", "C2H5OH"), "basis species .* is not defined")
  expect_error(swap.basis("CO2", "C60"), "is not available")
  expect_error(swap.basis("CO2", "H2"), "the number of basis species is greater than the number of elements and charge")
})

test_that("basis.logact only accepts defined elements", {
  # setup basis species with two elements: C and H
  basis(c("graphite", "H2"), c("cr", "gas"))
  # we can't get basis activities with one element
  expect_error(basis.logact(c(C=1)), "number of elements in 'emu' is less than those in basis")
})

# 20181111
test_that("swapping works with a buffer (no recalculation of activities)", {
  basis("FeCHNOS+")
  oldb <- basis("O2", "PPM")
  # before version 1.1.3-57, this gave Error in value * (-log(10) * R * T) : non-numeric argument to binary operator
  newb <- swap.basis("O2", "hydrogen")
  # note: logact includes "PPM" for O2 (old) and H2 (new)
  expect_identical(oldb$logact, newb$logact)
})

# 20200728 moved from swap-basis.Rd
test_that("swapping doesn't affect affinities of formation reactions of species", {
  ## swapping basis species while species are defined
  ## and using numeric species indices
  basis("MgCHNOPS+") 
  # load some Mg-ATP species
  species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"))
  # swap in CO2(g) for CO2(aq)
  swap.basis("CO2", "carbon dioxide")
  a1 <- affinity()
  # swap in CH4(g) for CO2(g)
  swap.basis("carbon dioxide", "methane")
  a2 <- affinity()
  # the equilibrium fugacity of CH4 is *very* low
  # swap in CO2(aq) for CH4(g)
  swap.basis("methane", "CO2")
  a3 <- affinity()
  # swapping the basis species didn't affect the affinities
  # of the formation reactions of the species, since
  # the chemical potentials of the elements were unchanged
  expect_equal(a1$values, a2$values)
  expect_equal(a1$values, a3$values)
})
