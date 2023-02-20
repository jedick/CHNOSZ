# Load default settings for CHNOSZ
reset()

info <- "NAs in thermo()$OBIGT propagate to subcrt()"
# First of all, water is in thermo()$OBIGT but its properties
# are actually calculated using water() so it has NAs for some parameters
expect_equal(info(1)$a, as.numeric(NA), info = info)
# Get the existing value of c for [Ala](cr) (it's 0)
expect_equal(c.Ala <- info(info("[Ala]", "cr"))$c, 0, info = info)
# When we make a protein, its G depends on temperature
expect_true(all(diff(subcrt("LYSC_CHICK", "cr")$out[[1]]$G) < 0), info = info)
# Turn the values of G and S for [Ala](cr) into NA
mod.OBIGT(name = "[Ala]", state = "cr", G = NA, S = NA)
# Now when we make a protein(cr), its G is NA
expect_true(all(is.na(subcrt("RNAS1_BOVIN", "cr")$out[[1]]$G)), info = info)
# Also check propagation of NA for aqueous species
mod.OBIGT(name = "[Ala]", state = "aq", G = NA, S = NA)
expect_true(all(is.na(subcrt("[Ala]", "aq")$out[[1]]$G)), info = info)
