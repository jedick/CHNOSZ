# Load default settings for CHNOSZ
reset()

info <- "pinfo() deals with underscore in first argument"
expect_equal(pinfo("LYSC_CHICK"), 6, info = info)

info <- "pinfo() returns NA for unmatched protein_organism pairs"
expect_equal(pinfo(c("LYSC", "zzzz"), "CHICK"), c(6, NA), info = info)

info <- "pinfo() returns NA if no organism is given"
expect_equal(pinfo(c("LYSC_CHICK", "MYGPHYCA")), c(6, NA), info = info)

