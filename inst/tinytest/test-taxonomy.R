# Load default settings for CHNOSZ
reset()

# Tests added on 2026-01-06

info <- "Can climb levels and extract taxonomic names"

# Get information about Homo sapiens from the packaged taxonomy files
taxdir <- system.file("extdata/taxonomy", package = "CHNOSZ")

# Start with a species (taxonomic id for H. sapiens)
id1 <- 9606
expect_equal(getrank(id1, taxdir), "species", info = info)

# Go up one level
id2 <- parent(id1, taxdir)
expect_equal(getrank(id2, taxdir), "genus", info = info)
expect_equal(sciname(id2, taxdir), "Homo", info = info)

# Go to an arbitrary level
id3 <- parent(id1, taxdir, "phylum")
expect_equal(sciname(id3, taxdir), "Chordata", info = info)

# H. sapiens' complete taxonomy
id4 <- allparents(id1, taxdir)
expect_length(sciname(id4, taxdir), 31, info = info)

info <- "Retrieves both 'species' and 'no rank'"
id5 <- c(83333, 4932, 9606, 186497, 243232)
expect_equal(unique(getrank(id5, taxdir)[2:3]), "species", info = info)
expect_equal(unique(getrank(id5, taxdir)[c(1, 4, 5)]), "no rank", info = info)
# Some of these have strain names, e.g. K-12 or DSM 2661
expect_equal(range(sapply(strsplit(sciname(id5, taxdir), " "), length)), c(2, 4), info = info)

info <- "Can go from 'no rank' to 'species'"
id6 <- sapply(id5, function(x) parent(x, taxdir = taxdir, rank = "species"))
expect_equal(unique(getrank(id6, taxdir)), "species", info = info)
# All the speices names are binomials
expect_equal(unique(sapply(strsplit(sciname(id6, taxdir), " "), length)), 2, info = info)

info <- "Use getnodes to keep nodes in memory"
taxdir <- system.file("extdata/taxonomy", package = "CHNOSZ")
nodes <- getnodes(taxdir = taxdir)
expect_equal(dim(nodes), c(63, 3), info = info)
