context("mosaic")

# this is a long test ... skip it if we're on R CMD check --as-cran
if(!any(grepl("R_CHECK_TIMINGS", names(Sys.getenv())))) {

test_that("results are consistent with affinity()", {
  basis(c("CO2", "H2O", "NH3", "O2"), c(0, 0, 0, 0))
  species(c("alanine", "glycine"))
  a25 <- affinity()
  # this is a degenerate case because we only allow NH3 to swap for NH3, and CO2 for CO2;
  # however it still exercises the affinity scaling and summing code
  m1_25 <- mosaic("NH3", "CO2")
  # this failed before we divided by loga.tot to get _relative_ abundances of basis species in mosaic.R
  expect_equal(a25$values, m1_25$A.species$values)
  # the next call failed when which.pmax(), called by diagram(), choked on a list of length one
  m2_25 <- mosaic("NH3", "CO2", blend = FALSE)
  expect_equal(a25$values, m2_25$A.species$values)
  # make sure the function works when all affinities are NA
  a500 <- suppressWarnings(affinity(T=500))
  # using blend=TRUE was failing prior to version 1.1.3-37
  m1_500 <- suppressWarnings(mosaic("NH3", "CO2", T=500))
  expect_equal(a500$values, m1_500$A.species$values)
  m2_500 <- suppressWarnings(mosaic("NH3", "CO2", blend = FALSE, T=500))
  expect_equal(a500$values, m2_500$A.species$values)
})

test_that("blend=TRUE produces reasonable values", {
  # a more rigorous test than above. this was failing because loga.tot (actually, a.tot)
  # was computed incorrectly, by sum()ing an unlist()ed list (the affinities of basis species)
  # to produce a single value; corrected by using Reduce for addition of vectors/arrays in the list.
  # example adapted from ?mosaic
  basis(c("FeO", "SO4-2", "H2O", "H+", "e-"))
  basis("SO4-2", -6)
  basis("Eh", -0.15)
  species(c("hematite", "magnetite"))
  # the basis species we'll swap through
  bases <- c("SO4-2", "HSO4-", "HS-", "H2S")         
  # calculate affinities using the predominant basis species
  pH <- c(0, 14, 29)
  m1 <- mosaic(bases, pH = pH, blend = FALSE)
  # calculate affinities with smooth transitions between basis species
  m2 <- mosaic(bases, pH = pH)
  # these species have no S so the results should be similar,
  expect_equivalent(m2$A.species$values[[1]], m1$A.species$values[[1]])
  # now with S-bearing species ...
  species(c("pyrrhotite", "pyrite"))
  m3 <- mosaic(bases, pH = pH, blend = FALSE)
  m4 <- mosaic(bases, pH = pH)
  # the results are different ...
  expect_equal(sapply(m3$A.species$values, "[", 13), sapply(m4$A.species$values, "[", 13), tol=1e-1)
  # but more similar at extreme pH values
  expect_equal(sapply(m3$A.species$values, "[", 1), sapply(m4$A.species$values, "[", 1), tol=1e-7)
  expect_equal(sapply(m3$A.species$values, "[", 29), sapply(m4$A.species$values, "[", 29), tol=1e-13)
})

test_that("mosaic() - equilibrate() produces equilibrium activities", {
  # test added 20190505, based on a calculation sent by Kirt Robinson
  basis(c("CO2", "NH3", "O2", "H2O", "H+"))
  species(c("acetamide", "acetic acid", "acetate"))
  m <- mosaic(c("NH3", "NH4+"), pH = c(0, 14))
  e <- equilibrate(m$A.species)
  # calculate logK for form acetamide from predominant species at low pH
  s1 <- subcrt(c("acetic acid", "NH4+", "acetamide", "water", "H+"), c(-1, -1, 1, 1, 1), T = 25)
  logK1 <- s1$out$logK
  # values of activities
  loga_acetic <- e$loga.equil[[2]]
  loga_NH4 <- m$E.bases[[1]]$loga.equil[[2]]
  loga_acetamide <- e$loga.equil[[1]]
  loga_H2O <- m$E.bases[[1]]$basis$logact[[4]]
  loga_Hplus <- - m$E.bases[[1]]$vals$pH
  logQ1 <- - loga_acetic - loga_NH4 + loga_acetamide + loga_H2O + loga_Hplus
  A1 <- logQ1 - logK1
  ## in CHNOSZ versions before 1.3.2-5 (20190505), the affinity was zero at the pH extremes,
  ## but peaked with a value of 0.3 (log10(2)) at pH 9.2 (equal activities of NH3 and NH4+)
  #plot(m$E.bases[[1]]$vals$pH, A1, type = "l")
  #title(main = describe.reaction(s1$reaction))
  expect_equivalent(A1, rep(0, length(A1)))
})

test_that("mosaic() - solubility() produces equilibrium activities", {
  # test added 20190505, simplified from demo/contour.R with varying pH at constant logfO2
  # define temperature and pressure
  T <- 250
  P <- 500
  # set up system
  basis(c("Au", "Cl-", "H2S", "H2O", "oxygen", "H+"))
  species(c("Au(HS)2-", "AuHS", "AuOH", "AuCl2-"))
  # this get us close to total S = 0.01 m
  basis("H2S", -2)
  # set a low logfO2 to get into H2S - HS- fields
  basis("O2", -40)
  # calculate solution composition for 1 mol/kg NaCl
  NaCl <- NaCl(T = T, P = P, m_tot = 1)
  basis("Cl-", log10(NaCl$m_Cl))
  # calculate affinity with changing basis species
  bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
  m <- mosaic(bases, pH = c(2, 10), T = 250, P = 500, IS = NaCl$IS)
  # calculate solubility
  s <- solubility(m$A.species)
  # calculate logK to form Au(HS)2- in H2S stability region
  # include IS here to compute adjusted logK
  s1 <- subcrt(c("Au", "H2S", "oxygen", "Au(HS)2-", "H2O", "H+"), c(-1, -2, -0.25, 1, 0.5, 1), T = T, P = P, IS = NaCl$IS)
  logK1 <- s1$out$logK
  # calculate logQ with the given or computed activities
  loga_Au <- m$A.bases$basis$logact[[1]]
  loga_H2S <- m$E.bases[[1]]$loga.equil[[1]]
  logf_O2 <- m$A.bases$basis$logact[[5]]
  loga_AuHS2minus <- s$loga.equil[[1]]
  loga_H2O <- m$A.bases$basis$logact[[4]]
  loga_Hplus <- - m$A.bases$vals$pH
  logQ1 <- - 1 * loga_Au - 2 * loga_H2S - 0.25 * logf_O2 + 1 * loga_AuHS2minus + 0.5 * loga_H2O + 1 * loga_Hplus
  # calculate affinity - should be zero!
  A1 <- logQ1 - logK1
  #plot(m$A.bases$vals$pH, A1, type = "l")
  #title(main = describe.reaction(s1$reaction))
  expect_equivalent(A1, rep(0, length(A1)))
})

test_that("mosaic() - equilibrate() produces equilibrium activities that are consistent with
          stability differences of minerals and multiple groups of changing basis species", {
  # test added 20190505, adapted from demo/mosaic.R:
  #   select a constant pH close to equal activities of CO2 - HCO3-,
  #   and a range of Eh that crosses the upper and lower boundaries
  #   of pyrite with siderite (including the H2S - SO4-2 transition)
  basis(c("FeO", "SO4-2", "H2O", "H+", "e-", "CO3-2"))
  basis("SO4-2", -6)
  basis("CO3-2", 0)
  basis("pH", 6.3)
  species(c("pyrrhotite", "pyrite", "hematite", "magnetite", "siderite"))
  # two sets of changing basis species:
  # speciate SO4-2, HSO4-, HS-, H2S as a function of Eh and pH
  # speciate CO3-2, HCO3-, CO2 as a function of pH
  bases <- list(c("SO4-2", "HSO4-", "HS-", "H2S"),
                c("CO3-2", "HCO3-", "CO2"))
  # calculate affinities using the relative abundances of different basis species
  m <- mosaic(bases, Eh = c(-0.5, 0))
  # calculate logK for pyrite-siderite reaction using arbitrarily chosen basis species
  s1 <- subcrt(c("pyrite", "CO2", "H2O", "H+", "e-", "siderite", "H2S"), c(-1, -1, -1, -2, -2, 1, 2), T = 25)
  logK <- s1$out$logK
  # get activities of minerals, water, and H+
  loga_pyrite <- loga_siderite <- loga_H2O <- 0
  loga_Hplus <- m$A.bases[[1]]$basis$logact[[4]]
  # get activities of basis species
  loga_eminus <- - convert(m$A.bases[[1]]$vals$Eh, "pe")
  loga_H2S <- m$E.bases[[1]]$loga.equil[[4]]
  loga_CO2 <- m$E.bases[[2]]$loga.equil[[3]]
  # calculate affinity
  logQ <- -1 * loga_pyrite - 1 * loga_CO2 - 1 * loga_H2O - 2 * loga_Hplus - 2 * loga_eminus + 1 * loga_siderite + 2 * loga_H2S
  A <- logQ - logK
  # the "hand-calculated" value and the affinity calculated by the function should be equal
  Adiff <- A - (m$A.species$values[[2]] - m$A.species$values[[5]])
  #par(mfrow = c(2, 1))
  #diagram(m$A.species)
  #plot(m$A.bases[[1]]$vals$Eh, Adiff, type = "l")
  #title(main = "A(single basis species) - A(all basis species)")
  #legend("topleft", legend = describe.reaction(s1$reaction))
  expect_equivalent(Adiff, rep(0, length(Adiff)))
})

}
