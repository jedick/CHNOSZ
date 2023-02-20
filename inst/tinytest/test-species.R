# Load default settings for CHNOSZ
reset()

info <- "Species not contained by basis cause errors"
expect_error(species("H2O"), "basis species are not defined", info = info)
expect_error(CHNOSZ:::species.basis("H2O"), "basis species are not defined", info = info)
basis("CHNOS")
expect_error(CHNOSZ:::species.basis("U"), "element\\(s\\) not in the basis\\: U", info = info)
expect_error(species("fayalite"), "element\\(s\\) not in the basis\\: Fe Si", info = info)

info <- "For one or more species, species.basis() keeps track of zeroes and puts elements in order of thermo()$basis"
basis("CHNOS")
test0 <- count.elements("OHN0")
test1 <- count.elements("HN0O")
expect_equal(CHNOSZ:::species.basis(test0), CHNOSZ:::species.basis(test1), info = info)
# We can send multiple species to species.basis() but the argument has to be constructed correctly
expect_equal(unique(as.numeric(CHNOSZ:::species.basis(makeup(c("C", "CCN"))))), 0, info = info)
expect_equal(CHNOSZ:::species.basis(makeup(c("C", "CCN"), count.zero = TRUE))[2, , drop = FALSE], CHNOSZ:::species.basis(makeup("CCN")), info = info)

info <- "Deleting nonexistent species causes error or warning"
expect_error(species("CO2", delete = TRUE), "nonexistent species definition", info = info)
species("H2O")
expect_warning(species("CO2", delete = TRUE), "not present, so can not be deleted", info = info)
expect_null(species("water", delete = TRUE), info = info)
# We should also get NULL if *all* species are deleted
species("H2O")
expect_null(species(delete = TRUE), info = info)

info <- "Non-available species cause error, and species can be added or modified"
basis("CHNOS")
expect_error(species("wate"), "species not available", info = info)
# Add CO2, aq
sdef <- species("CO2")
# We can't add the same species twice
expect_equal(nrow(species("CO2")), 1, info = info)
# Change it to gas
expect_equal(species(1, "gas")$state, "gas", info = info)
# Change its log fugacity to -5
expect_equal(species(1, -5)$logact, -5, info = info)
# Add CO2, aq
expect_equal(nrow(species("CO2", add = TRUE)), 2, info = info)
# Add alanine by index in thermo()$OBIGT
expect_equal(nrow(species(info("alanine"), add = TRUE)), 3, info = info)
# If we just use an index, get only that species
expect_equal(species(3)$name, "alanine", info = info)
# We can add a species with the same name but different state
expect_equal(nrow(species("alanine", "cr", add = TRUE)), 4, info = info)
# We can modify the logact of a named species (only first match)
expect_equal(species("alanine", -7)$logact[3], -7, info = info)

info <- "index_return provides indices for touched species"
basis("CHNOS")
expect_equal(species("CO2", index.return = TRUE), 1, info = info)
# Here it's "touched" (but not added or modified)
expect_equal(species("CO2", index.return = TRUE), 1, info = info)
expect_equal(species(c("H2O", "NH3"), index.return = TRUE, add = TRUE), c(2, 3), info = info)
expect_equal(species(1, "gas", index.return = TRUE), 1, info = info)
