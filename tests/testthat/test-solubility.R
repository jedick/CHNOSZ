context("solubility")

test_that("solubility() produces stable conditions (affinity = 0)", {
  # set pH range and resolution, constant temperature and ionic strength
  pH <- c(0, 14)
  res <- 100
  T <- 25
  IS <- 0

  ## start with CO2
  basis(c("carbon dioxide", "H2O", "O2", "H+"))
  # ca. atmospheric PCO2
  basis("CO2", -3.5)
  species(c("CO2", "HCO3-", "CO3-2"))
  a <- affinity(pH = c(pH, res), T = T, IS = IS)
  s <- solubility(a)
  # a function to check for stable conditions (affinity = 0)
  # do this by setting activities in species() then calculating the affintiy
  checkfun <- function(i) {
    logact <- sapply(s$loga.equil, "[", i)
    species(1:3, logact)
    basis("pH", s$vals[[1]][i])
    affinity(T = T, IS = IS)
  }
  # check any 'i' here - let's just take two
  expect_equal(max(abs(unlist(checkfun(33)$values))), 0)
  expect_equal(max(abs(unlist(checkfun(99)$values))), 0)

  # now do calcite
  basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
  species(c("CO2", "HCO3-", "CO3-2"))
  a <- affinity(pH = c(pH, res), T = T, IS = IS)
  s <- solubility(a, dissociate = TRUE)
  # here we need to also set the activity of Ca+2
  checkfun <- function(i) {
    logact <- sapply(s$loga.equil, "[", i)
    species(1:3, logact)
    basis("pH", s$vals[[1]][i])
    basis("Ca+2", s$loga.balance[i])
    affinity(T = T, IS = IS)
  }
  expect_equal(max(abs(unlist(checkfun(33)$values))), 0)
  expect_equal(max(abs(unlist(checkfun(99)$values))), 0)
})

test_that("backward compatible and new calling styles produce identical results", {
  # Test added 20210319
  # Calculate solubility of a single substance:
  # Gaseous S2 with a given fugacity

  # Define basis species (any S-bearing basis species is allowed)
  basis(c("HS-", "oxygen", "H2O", "H+"))
  basis("pH", 6)
  # Load the substances (minerals or gases) to be dissolved
  species("S2", -20)
  # List the formed aqueous species
  i_aq <- info(c("SO4-2", "HS-"))
  # Place arguments for affinity() or mosaic() after the first argument of solubility()
  s1 <- solubility(i_aq, O2 = c(-55, -40), T = 125, in.terms.of = "SO4-2")

  # Backward-compatible method limited to dissolving one species:
  # Include S2(g) in the basis species
  basis(c("S2", "oxygen", "H2O", "H+"))
  basis("pH", 6)
  basis("S2", -20)
  # Calculate affinities for formation of aqueous species
  species(c("SO4-2", "HS-"))
  a <- affinity(O2 = c(-55, -40), T = 125)
  s_old <- solubility(a, in.terms.of = "SO4-2")

  expect_identical(s1, s_old)
})
