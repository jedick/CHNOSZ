# Load default settings for CHNOSZ
reset()

# Test added on 20240211
# Incorrect result was causing quartz tests in test-subcrt.R to fail
info <- "dPdTtr gives expected result"
add.OBIGT("SUPCRT92")
dPdTtr.calc <- round(dPdTtr(info("quartz", "cr"), info("quartz", "cr2")), 5)
# The reference value was calculated with CHNOSZ_1.4.3
# (prior to bug in dPdTtr introduced by switch to Joules)
expect_equal(dPdTtr.calc, 38.45834, info = info)

# Tests added on 20260109

info <- "GHS_Tr() runs without error"
Htr <- c(326.0, 215.0, 165.0)
iiron <- info("iron")
expect_silent(result <- GHS_Tr(iiron, Htr), info = info)
info <- "GHS_Tr() gives expected output"
expect_equal(round(result$Gf_Tr, 2), c(-618.24, 1293.53, 979.29), info = info)
expect_equal(round(result$Hf_Tr, 2), c(-1768.42, 1865.34, 1252.89), info = info)
expect_equal(round(result$S_Tr, 2), c(2.66, 8.44, 7.44), info = info)

info <- "unitize() gives expected output"
# The proteins have unequal activity, but the activities of residues still add up to one
logact.tot <- 0
logact <- c(-3, -2)
length <- c(100, 200)
loga <- unitize(logact, length, logact.tot)
expect_equal(loga[2] - loga[1], 1)
expect_equal(sum(10^loga * length), 1)
