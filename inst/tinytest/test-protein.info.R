# Load default settings for CHNOSZ
reset()

info <- "pinfo() deals with underscore in first argument"
expect_equal(pinfo("LYSC_CHICK"), 6, info = info)

info <- "pinfo() returns NA for unmatched protein_organism pairs"
expect_equal(pinfo(c("LYSC", "zzzz"), "CHICK"), c(6, NA), info = info)

info <- "pinfo() returns NA if no organism is given"
expect_equal(pinfo(c("LYSC_CHICK", "MYGPHYCA")), c(6, NA), info = info)

info <- "protein.equil() reports values consistent with Dick and Shock (2011)"
protein <- pinfo(c("CSG_METVO", "CSG_METJA"))
# To reproduce the calculations in the paper, use superseded properties of [Met], [Gly], and [UPBB]
mod.OBIGT("[Met]", G = -35245, H = -59310, S = 40.38, E_units = "cal")
mod.OBIGT("[Gly]", G = -6075, H = -5570, S = 17.31, E_units = "cal")
mod.OBIGT("[UPBB]", G = -21436, H = -45220, S = 1.62, E_units = "cal")
basis("CHNOS+")
suppressMessages(swap.basis("O2", "H2"))
pequil <- protein.equil(protein, loga.protein=-3)

# Before 20230630:
# With R = 8.314445 (= 1.9872 * 4.184) in util.units(), ionize.aa(pinfo(protein)) gives:
#            20        18
#[1,] -56.06509 -55.87025
# Astar/RT in the paragraph following Eq. 23, p. 6 of DS11
# (truncated because of rounding)
# expect_true(any(grepl(c("0\\.435.*1\\.36"), pequil)), info = info)

# After 20230630:
# With R = 8.314463 in util.units(), ionize.aa(pinfo(protein)) gives:
#            20        18
#[1,] -56.06511 -55.87027
# The differences propagate up, so the test was changed on 20230630:
expect_true(any(grepl(c("0\\.437.*1\\.36"), pequil)), info = info)

# log10 activities of the proteins in the left-hand column of the same page
expect_true(any(grepl(c("-3\\.256.*-2\\.834"), pequil)), info = info)

# Reference

# Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins 
#   as a function of temperature and redox chemistry in a hot spring. 
#   PLoS ONE 6, e22782. https://doi.org/10.1371/journal.pone.0022782
