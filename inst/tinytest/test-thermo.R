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


### Tests added 20230310 for changes to thermo() argument handling more like par()

info <- "Alternative indexing styles give the same result"
E1 <- thermo()$opt$E.units  # The "old" way
E2 <- thermo("opt$E.units") # The "new" way in CHNOSZ 2.0.0
expect_equal(E2, E1, info = info)

info <- "Value retrieved for opt$E.units is J"
expect_equal(E1, "J", info = info)

## Assign a strange value to opt$E.units
#oldthermo <- thermo("opt$E.units" = 1234)
#info <- "Restoring old parameter values is possible"
#expect_silent(E3 <- thermo(oldthermo), info = info)
#info <- "Old values are restored correctly"
#E4 <- thermo("opt$E.units")
#expect_equal(E4, E1, info = info)

info <- "Parameters can be selected using c() or argument list"
BS1 <- thermo("basis", "species")
BS2 <- thermo(c("basis", "species"))
expect_equal(BS1, BS2, info = info)

info <- "Single parameter gives atomic vector"
expect_null(names(E1), info = info)
info <- "Two more parameters give named list"
expect_equal(names(BS1), c("basis", "species"), info = info)
