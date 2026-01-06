# Load default settings for CHNOSZ
reset()

# Tests added on 2026-01-06

info <- "Functions for expressions used in legends give expected output"
expect_equal(lNaCl(0.1), quote(italic(m) * NaCl == 0.1 ~ mol ~ kg^-1), info = info)
expect_equal(lS(0.1), quote(sum(S) == 0.1 ~ mol ~ kg^-1), info = info)
expect_equal(lT(25), quote(25 ~ degree * C), info = info)
expect_equal(lP(1000), quote(1000 ~ bar), info = info)
expect_equal(lTP(100, 2000), quote(list(100 ~ degree * C, 2000 ~ bar)), info = info)
