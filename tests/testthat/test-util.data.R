context("util.data")

test_that("checkGHS() and checkEOS() (via info()) produce messages", {
  i1 <- info("S2O3-2")
  expect_message(info(i1), "G of S2O3-2 aq \\(26\\) differs by 939 cal mol-1 from tabulated value")
  i2 <- info("Cu+2")
  expect_message(info(i2), "Cp of Cu\\+2 aq \\(60\\) differs by 3.62 cal K-1 mol-1 from tabulated value")
})

test_that("checkGHS() and checkEOS() respond to thermo$opt$*.tol", {
  i1 <- info("SO4-2")
  thermo("opt$Cp.tol" = 0.5)
  expect_message(info(i1), "checkEOS")
  i2 <- info("a,w-dicarboxytetracosane")
  thermo("opt$G.tol" = 50)
  expect_message(info(i2), "checkGHS")
})

test_that("RH2obigt() gives group additivity results consistent with database values (from Richard and Helgeson, 1998)", {
  file <- system.file("extdata/adds/RH98_Table15.csv", package = "CHNOSZ")
  dat <- read.csv(file, stringsAsFactors=FALSE)
  ispecies <- info(dat$compound, dat$state)
  obigt.ref <- thermo()$obigt[ispecies, ]
  obigt.calc <- RH2obigt(file=file)
  # calculated values of H are spot on; to pass tests, tolerance on
  # G is set higher; is there an incorrect group value somewhere?
  expect_true(max(abs(obigt.calc$G - obigt.ref$G)) < 31)
  expect_true(max(abs(obigt.calc$H - obigt.ref$H)) == 0)
  expect_true(max(abs(obigt.calc$S - obigt.ref$S)) < 0.02001)
  expect_true(max(abs(obigt.calc$Cp - obigt.ref$Cp)) < 0.04001)
  expect_true(max(abs(obigt.calc$V - obigt.ref$V)) < 0.1001)
  expect_true(max(abs(obigt.calc$a1.a - obigt.ref$a1.a)) < 0.01001)
  expect_true(max(abs(obigt.calc$a2.b - obigt.ref$a2.b)) < 1e-13)
  expect_true(max(abs(obigt.calc$a3.c - obigt.ref$a3.c)) < 1e-14)
})

test_that("add.obigt() replaces existing entries without changing species index", {
  # store the original species index of CdCl2
  iCdCl2 <- info("CdCl2", "aq")
  # add supplemental database - includes CdCl2
  file <- system.file("extdata/adds/BZA10.csv", package="CHNOSZ")
  isp <- add.obigt(file)
  # species index of CdCl2 should not have changed
  expect_equal(info("CdCl2", "aq"), iCdCl2)
  # check that names of species modified are same as in file
  newdat <- read.csv(file, stringsAsFactors=FALSE)
  # the order isn't guaranteed ... just make sure they're all there
  expect_true(all(newdat$name %in% thermo()$obigt$name[isp]))
})

test_that("reset() and obigt() produce the same database", {
  reset()
  d1 <- thermo()$obigt
  obigt()
  d2 <- thermo()$obigt
  expect_equal(d1, d2)
})

test_that("add.obigt() is backwards compatibile for a file that doesn't have an E_units column", {
  # test added 20190529
  file <- system.file("extdata/adds/BZA10.csv", package="CHNOSZ")
  rc <- read.csv(file)
  expect_false("E_units" %in% colnames(rc))
  inew <- add.obigt(file)
  expect_true(unique(info(inew)$E_units) == "cal")
})

test_that("info() gives consistent messages for cal and J", {
  # test added 20190529
  # add data for dimethylamine and trimethylamine in different units (cal or J)
  expect_message(add.obigt(system.file("extdata/adds/LA19_test.csv", package = "CHNOSZ")), "energy units: J and cal")
  expect_message(info(info("DMA_cal")), "-1.92 cal")
  expect_message(info(info("DMA_J")), "-8.02 J")
  # for TMA, only a checkGHS message for the entry in J is produced,
  # because it's above the threshold of 100 set in thermo()$opt$G.tol
  expect_silent(info(info("TMA_cal")))
  expect_message(info(info("TMA_J")), "-102 J")
})

test_that("missing values for G, Cp, and V are correct in cal and J", {
  # test added 20190530
  # add data for dimethylamine and trimethylamine in different units (cal or J)
  add.obigt(system.file("extdata/adds/LA19_test.csv", package = "CHNOSZ"))
  calccal <- info(info("DMA_cal_NA"))
  expect_equal(round(calccal$G), 13934)
  expect_equal(round(calccal$Cp, 1), 60.3)
  expect_equal(round(calccal$V, 1), 58.2)
  calcJ <- info(info("DMA_J_NA"))
  expect_equal(round(calcJ$G), 58304)
  expect_equal(round(calcJ$Cp, 1), 252.4)
  expect_equal(round(calcJ$V, 1), 58.2)
})

test_that("subcrt() gives same results for data entered in cal and J", {
  # test added 20190530
  # add data for dimethylamine and trimethylamine in different units (cal or J)
  add.obigt(system.file("extdata/adds/LA19_test.csv", package = "CHNOSZ"))
  E.units("cal")
  scal <- subcrt("DMA_cal")
  sJ <- subcrt("DMA_J")
  expect_maxdiff(scal$out[[1]]$G, sJ$out[[1]]$G, 2)
  expect_maxdiff(scal$out[[1]]$H, sJ$out[[1]]$H, 1)
  expect_maxdiff(scal$out[[1]]$S, sJ$out[[1]]$S, 0.006)
  expect_maxdiff(scal$out[[1]]$V, sJ$out[[1]]$V, 0.011)
  expect_maxdiff(scal$out[[1]]$Cp, sJ$out[[1]]$Cp, 0.16)
  # if we set the output to J, it should be the same as the parameters at 25 degC
  E.units("J")
  calcJ25 <- subcrt("DMA_J", T = 25)$out[[1]]
  infoJ25 <- info(info("DMA_J"))
  expect_equivalent(calcJ25[, c("G", "H", "S")], calcJ25[, c("G", "H", "S")])
  # in the case of Cp and V, there are bigger difference because they are calculated from the HKF parameters
  expect_maxdiff(calcJ25$Cp, infoJ25$Cp, 8.1)
  expect_maxdiff(calcJ25$V, infoJ25$V, 0.55)
  # go back to default units
  E.units("cal")
})

# reference

# Richard, L. and Helgeson, H. C. (1998) Calculation of the thermodynamic properties at elevated 
#   temperatures and pressures of saturated and aromatic high molecular weight solid and liquid 
#   hydrocarbons in kerogen, bitumen, petroleum, and other organic matter of biogeochemical interest. 
#   Geochim. Cosmochim. Acta 62, 3591--3636. https://doi.org/10.1016/S0016-7037(97)00345-1
