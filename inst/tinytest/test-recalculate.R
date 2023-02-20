# Load default settings for CHNOSZ
reset()

# SSH97: Sverjensky et al., 1997
# doi:10.1016/S0016-7037(97)00009-4
info <- "Recalculated values for HSiO3- match those in original reference"
iSiO2 <- info("SiO2")
iHSiO3 <- info("HSiO3-")
# Get values used in OBIGT
SiO2_GHS_OBIGT <- thermo()$OBIGT[iSiO2, c("G", "H", "S")]
HSiO3_GHS_OBIGT <- thermo()$OBIGT[iHSiO3, c("G", "H", "S")]
# Get values from SSH97
add.OBIGT("SLOP98")
SiO2_GHS_SSH97 <- thermo()$OBIGT[iSiO2, c("G", "H", "S")]
HSiO3_GHS_SSH97 <- thermo()$OBIGT[iHSiO3, c("G", "H", "S")]
# Calculate differences
OBIGT_Delta <- HSiO3_GHS_OBIGT - SiO2_GHS_OBIGT
SSH97_Delta <- HSiO3_GHS_SSH97 - SiO2_GHS_SSH97
# Perform test
expect_equal(OBIGT_Delta, SSH97_Delta, info = info)
