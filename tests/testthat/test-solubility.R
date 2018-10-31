context("solubility")

test_that("solubility() produces stable conditions (affinity = 0)", {
  # set pH range and resolution, constant temperature and ionic strength
  pH <- c(0, 14)
  res <- 100
  T <- 25
  IS <- 0

  # start with CO2
  basis(c("carbon dioxide", "H2O", "O2", "H+"))
  # ca. atmospheric PCO2
  basis("CO2", -3.5)
  species(c("CO2", "HCO3-", "CO3-2"))
  a <- affinity(pH = c(pH, res), T = T, IS = IS)
  e <- equilibrate(a)
  s <- solubility(e)

  # check for stable conditions (affinity = 0)
  species(1:3, 0)
  atest <- affinity(pH = s$vals[[1]], T = T, IS = IS)
  expect_true(all(sapply(unlist(atest$values) - unlist(s$loga.equil), all.equal, 0)))

  # now do calcite
  basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
  species(c("CO2", "HCO3-", "CO3-2"))
  a <- affinity(pH = c(pH, res), T = T, IS = IS)
  e <- equilibrate(a)
  s <- solubility(e, exp = 2)

  # check for stable conditions (affinity = 0)
  species(1:3, 0)
  atest <- affinity(pH = s$vals[[1]], `Ca+2` = s$loga.balance, T = T, IS = IS)
  expect_true(all(sapply(unlist(atest$values) - unlist(s$loga.equil), all.equal, 0)))
})
