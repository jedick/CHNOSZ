context("nonideal")

# 20161219
test_that("loggam and logK values are consistent", {
  oldnon <- nonideal("Alberty")
  rxn1 <- subcrt(c("anhydrite", "Ca+2", "SO4-2"), c(-1, 1, 1), P=1, T=25, I=0)
  rxn2 <- subcrt(c("anhydrite", "Ca+2", "SO4-2"), c(-1, 1, 1), P=1, T=25, I=0.7)
  expect_equal(rxn1$out$logK, rxn2$out$loggam + rxn2$out$logK)
  nonideal(oldnon)
})

# 20171011
test_that("A and B parameters are calculated correctly", {
  ## Psat
  # values from Helgeson, 1967 "Solution chemistry and metamorphism"
  # (chapter in http://www.worldcat.org/oclc/152725534)
  T <- c(25, 50, 100, 200, 300, 350)
  A <- c(0.5095, 0.5354, 0.6019, 0.8127, 1.2979, 1.9936)
  B <- c(0.3284, 0.3329, 0.3425, 0.3659, 0.4010, 0.4300)
  SUP <- water.SUPCRT92(c("A_DH", "B_DH"), T=convert(T, "K"), P="Psat")
  IAP <- water.IAPWS95(c("A_DH", "B_DH"), T=convert(T, "K"), P="Psat")
  expect_maxdiff(SUP$A_DH, A, 0.18)
  expect_maxdiff(IAP$A_DH, A, 0.11)
  expect_maxdiff(SUP$B_DH / 1e8, B, 0.012)
  expect_maxdiff(IAP$B_DH / 1e8, B, 0.008)

  # values digitized from Fig. 10 of Manning et al., 2013
  # doi: 10.2138/rmg.2013.75.5
  ## 5 kbar
  T5 <- c(25, seq(100, 1000, 100))
  A5 <- c(0.441, 0.49, 0.583, 0.685, 0.817, 0.983, 1.164, 1.409, 1.673, 1.938, 2.187)
  B5 <- c(0.328, 0.336, 0.350, 0.363, 0.377, 0.391, 0.405, 0.421, 0.434, 0.445, 0.453)
  SUP5 <- water.SUPCRT92(c("A_DH", "B_DH"), T=convert(T5, "K"), P=rep(5000, 11))
  IAP5 <- water.IAPWS95(c("A_DH", "B_DH"), T=convert(T5, "K"), P=rep(5000, 11))
  DEW5 <- water.DEW(c("A_DH", "B_DH"), T=convert(T5, "K"), P=rep(5000, 11))
  # DEW is the winner here
  expect_maxdiff(SUP5$A_DH, A5, 0.47)
  expect_maxdiff(IAP5$A_DH, A5, 0.26)
  expect_maxdiff(DEW5$A_DH, A5, 0.14)
  expect_maxdiff(SUP5$B_DH / 1e8, B5, 0.036)
  expect_maxdiff(IAP5$B_DH / 1e8, B5, 0.021)
  expect_maxdiff(DEW5$B_DH / 1e8, B5, 0.013)

  ## 30 kbar
  T30 <- seq(700, 1000, 100)
  A30 <- c(0.625, 0.703, 0.766, 0.815)
  B30 <- c(0.386, 0.400, 0.409, 0.416)
  DEW30 <- water.DEW(c("A_DH", "B_DH"), T=convert(T30, "K"), P=rep(30000, 4))
  expect_maxdiff(DEW30$A_DH, A30, 0.06)
  expect_maxdiff(DEW30$B_DH / 1e8, B30, 0.024)
})

test_that("affinity transect incorporates IS correctly", {
  basis("CHNOS+")
  species("acetate")
  # calculations at single combinations of logfO2 and IS
  basis("O2", -80); a80_0 <- affinity()
  basis("O2", -60); a60_1 <- affinity(IS=1)
  # calculations on a transect with those endpoints
  a <- affinity(O2=seq(-80, -60, length.out=4), IS=seq(0, 1, length.out=4))
  expect_equal(a$values[[1]][1], a80_0$values[[1]][1])
  expect_equal(a$values[[1]][4], a60_1$values[[1]][1])
  # 20171013: that was working fine, but how about a more complicated case involving T?
  a25_0 <- affinity()
  a50_1 <- affinity(T=50, IS=1)
  a <- affinity(T=seq(25, 50, length.out=4), IS=seq(0, 1, length.out=4))
  expect_equal(a$values[[1]][1], a25_0$values[[1]][1])
  expect_equal(a$values[[1]][4], a50_1$values[[1]][1])
})

# 20171221
test_that("nonideality calculations work for Zn", {
  # nonideal() had a bug where both Z and Zn were identified as the charge
  # in the species formula, producing an error in this calculation
  expect_type(subcrt(c("Zn+2", "Cl-", "ZnCl+"), c(-1, -1, 1), T=200, P=16, IS=0.05), "list")   
})

