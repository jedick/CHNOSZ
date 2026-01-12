# Load default settings for CHNOSZ
reset()

info <- "pinfo() deals with underscore in first argument"
expect_equal(pinfo("LYSC_CHICK"), 6, info = info)

info <- "pinfo() returns NA for unmatched protein_organism pairs"
expect_equal(pinfo(c("LYSC", "zzzz"), "CHICK"), c(6, NA), info = info)

info <- "pinfo() returns NA if no organism is given"
expect_equal(pinfo(c("LYSC_CHICK", "MYGPHYCA")), c(6, NA), info = info)

# Tests added on 20260112

info <- "protein.basis() works as expected"
basis("CHNOS+")
# Formation reaction of unionized protein using species()
s <- species("LYSC_CHICK")
# Basis species to form ionized protein at low and high pH
basis("pH", 0)
pb0 <- protein.basis("LYSC_CHICK")
basis("pH", 14)
pb14 <- protein.basis("LYSC_CHICK")
# Numbers of CO2, H2O, NH3, and H2S should be equal
expect_equivalent(as.numeric(s[, 1:5]), pb0[, 1:5], info = info)
expect_equivalent(as.numeric(s[, 1:5]), pb14[, 1:5], info = info)
# Protein is more positively charged at low pH
expect_true(pb0[, "H+"] > pb14[, "H+"], info = info)
# Divide by protein length to get normalized basis composition
pb_norm <- protein.basis("LYSC_CHICK", normalize = TRUE)
expect_equal(pb_norm, pb14/129, info = info)
