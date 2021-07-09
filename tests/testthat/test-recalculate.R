context("recalculate")

# SSH97: Sverjensky et al., 1997
# doi:10.1016/S0016-7037(97)00009-4
test_that("recalculated values for HSiO3- match those in original reference", {
  iSiO2 <- info("SiO2")
  iHSiO3 <- info("HSiO3-")
  # get values used in OBIGT
  reset()
  SiO2_GHS_OBIGT <- thermo()$OBIGT[iSiO2, c("G", "H", "S")]
  HSiO3_GHS_OBIGT <- thermo()$OBIGT[iHSiO3, c("G", "H", "S")]
  # get values from SSH97
  add.OBIGT("SLOP98")
  SiO2_GHS_SSH97 <- thermo()$OBIGT[iSiO2, c("G", "H", "S")]
  HSiO3_GHS_SSH97 <- thermo()$OBIGT[iHSiO3, c("G", "H", "S")]
  # calculate differences
  OBIGT_Delta <- HSiO3_GHS_OBIGT - SiO2_GHS_OBIGT
  SSH97_Delta <- HSiO3_GHS_SSH97 - SiO2_GHS_SSH97
  # perform test
  expect_equal(OBIGT_Delta, SSH97_Delta)
  # clean up
  reset()
})
