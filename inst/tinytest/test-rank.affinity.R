# Load default settings for CHNOSZ
reset()

# Test added 20250522
basis("CHNOSe")
species(c("SO4-2", "H2S", "HS-"))
aout <- affinity(pH = c(0, 14, 10), Eh = c(-1, 1, 10))
groups <- list(oxidized = 1, reduced = 2:3)
# The tests only pass with rescale = TRUE (the default)
arank <- rank.affinity(aout, groups)
info <- "Different-sized groups have same range of average ranks"
expect_equal(range(arank$values[[1]]), c(1, 3), info = info)
expect_equal(range(arank$values[[2]]), c(1, 3), info = info)
info <- "The sum of average ranks from both groups is 1 + the number of species"
sum_average_ranks <- unique(as.integer(arank$values[[1]] + arank$values[[2]]))
expect_equal(sum_average_ranks, 1 + nrow(species()), info = info)

# Test added 20250527
info <- "Empty groups get removed from output"
groups$oxidized <- numeric()
arank <- rank.affinity(aout, groups)
expect_message(arank <- rank.affinity(aout, groups), "removing empty groups: oxidized", info = info)
expect_equal(arank$species$name, "reduced", info = info)
expect_length(arank$values, 1, info = info)
