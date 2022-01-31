# Load default settings for CHNOSZ
reset()

# 20200728 moved from objective.Rd
info <- "Calculations give expected results"
loga1 <- t(-4:-1)
loga2 <- loga1 + 1
expect_true(qqr(loga1) < 1, info = info)
expect_equal(RMSD(loga1, loga1), 0, info = info)
expect_equal(RMSD(loga1, loga2), 1, info = info)
expect_equal(CVRMSD(loga1, loga2), -0.4, info = info) # 1 / mean(-4:-1)
expect_equal(spearman(loga1, loga2), 1, info = info)
expect_equal(spearman(loga1, rev(loga2)), -1, info = info)
# less statistical, more thermodynamical...
expect_equal(DGmix(loga1), -0.1234, info = info) # as expected for decimal logarithms
expect_equal(DDGmix(loga1, loga2), 0.0004, info = info)
