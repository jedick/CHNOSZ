context("recalculate")

# SK95: Shock and Korestky, 1995
test_that("recalculated values from SK95 are correctly enetered in OBIGT", {
  # test added 20190206
  # thermodynamic data entries for amino acid glycinate and alanate complexes
  # (no longer in default database)
  data(thermo)
  add.obigt("OldAA")
  iaa <- info(c("alanate", "glycinate"))
  iGly <- grep("\\(Gly\\)\\+$", thermo$obigt$name)
  iGly2 <- grep("\\(Gly\\)2$", thermo$obigt$name)
  iAlan <- grep("\\(Alan\\)\\+$", thermo$obigt$name)
  iAlan2 <- grep("\\(Alan\\)2$", thermo$obigt$name)
  # get values used in OBIGT
  aa_GHS_OBIGT <- thermo$obigt[iaa, c("G", "H", "S")]
  Gly_GHS_OBIGT <- thermo$obigt[iGly, c("G", "H", "S")]
  Gly2_GHS_OBIGT <- thermo$obigt[iGly2, c("G", "H", "S")]
  Alan_GHS_OBIGT <- thermo$obigt[iAlan, c("G", "H", "S")]
  Alan2_GHS_OBIGT <- thermo$obigt[iAlan2, c("G", "H", "S")]
  # get values from SK95
  SK95 <- system.file("extdata/thermo/SK95.csv", package="CHNOSZ")
  add.obigt(SK95)
  aa_GHS_SK95 <- thermo$obigt[iaa, c("G", "H", "S")]
  Gly_GHS_SK95 <- thermo$obigt[iGly, c("G", "H", "S")]
  Gly2_GHS_SK95 <- thermo$obigt[iGly2, c("G", "H", "S")]
  Alan_GHS_SK95 <- thermo$obigt[iAlan, c("G", "H", "S")]
  Alan2_GHS_SK95 <- thermo$obigt[iAlan2, c("G", "H", "S")]
  # calculate differences for alanate and glycinate
  # nb. round values to avoid floating-point difficulties with unique() etc.
  Alan_GHS_Delta <- round(aa_GHS_SK95[1, ] - aa_GHS_OBIGT[1, ], 2)
  Gly_GHS_Delta <- round(aa_GHS_SK95[2, ] - aa_GHS_OBIGT[2, ], 2)
  # test that the differences are the same in the corresponding complexes
  expect_equivalent(Gly_GHS_Delta, unique(round(Gly_GHS_SK95 - Gly_GHS_OBIGT, 2)))
  expect_equivalent(Gly_GHS_Delta, unique(round((Gly2_GHS_SK95 - Gly2_GHS_OBIGT)/2, 2)))
  expect_equivalent(Alan_GHS_Delta, unique(round(Alan_GHS_SK95 - Alan_GHS_OBIGT, 2)))
  expect_equivalent(Alan_GHS_Delta, unique(round((Alan2_GHS_SK95 - Alan2_GHS_OBIGT)/2, 2)))
})
