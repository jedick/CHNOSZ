# Load default settings for CHNOSZ
reset()

# Test added on 20260112

info <- "MP90.cp() gives expected results"
# Reference values from Fig. 11 of Dick et al. (2006)
expect_equal(round(convert(MP90.cp("AMYA_PYRFU", 100), "cal") / 1e3), 44, info = info)
expect_equal(round(convert(MP90.cp("RNH_ECOLI", 100), "cal") / 1e3, 1), 9.8, info = info)
expect_equal(round(convert(MP90.cp("RNH_THET8", 100), "cal") / 1e3, 1), 10.3, info = info)
