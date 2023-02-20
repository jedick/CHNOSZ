# Load default settings for CHNOSZ
reset()

info <- "Invalid basis definitions cause an error"
expect_error(basis(character()), "argument is empty", info = info)
expect_error(basis(c("CO2", "CO2")), "names are not unique", info = info)
expect_error(basis(c("CO2", "H2O")), "the number of basis species is less than the number of elements", info = info)
expect_error(basis(c("H2O", "O2", "H2")), "the number of basis species is greater than the number of elements", info = info)
expect_error(basis(c("HCN", "H2O", "O2", "H2")), "singular", info = info)
expect_error(basis(c("CN", "H2O", "O2", "H2")), "species not available", info = info)
expect_error(basis(c("CN")), "species not available", info = info)
ina <- nrow(thermo()$OBIGT) + 1
expect_error(basis(ina), "species not available", info = info)
expect_error(CHNOSZ:::preset.basis(c("CN")), "is not a keyword", info = info)
# After all that, the basis should still be undefined
expect_null(basis(), info = info)

info <- "Invalid basis modification requests cause an error"
basis(delete = TRUE)
expect_error(CHNOSZ:::mod.basis("CH4", "gas"), "basis is not defined", info = info)
b <- basis("CHNOS+")
expect_error(CHNOSZ:::mod.basis("CH4", "gas"), "is not a formula of one of the basis species", info = info)
iCH4 <- info("CH4")
expect_error(CHNOSZ:::mod.basis(iCH4, "gas"), "is not a species index of one of the basis species", info = info)
expect_error(CHNOSZ:::mod.basis("CO2", "PPM"), "the elements .* in buffer .* are not in the basis", info = info)
expect_error(CHNOSZ:::mod.basis("CO2", "liq"), "state .* not found", info = info)
# After all that, the basis should be unchanged
expect_equal(basis(), b, info = info)

info <- "Modifying states of basis species is possible"
b1 <- basis(c("copper", "chalcocite"))
b2 <- basis("Cu2S", "cr2")
# We went from chalcocite cr to cr2, which is the next row in the database
expect_equal(sum(b2$ispecies - b1$ispecies), 1, info = info)
expect_error(basis("Cu2S", "cr4"), "state or buffer 'cr4' not found for chalcocite", info = info)
# Can we go from CO2(aq) to CO2(gas) back to CO2(aq)?
basis("CHNOS+")  # first basis species is CO2(aq)
expect_equal(basis("CO2", "gas")$state[1], "gas", info = info)
expect_equal(basis("CO2", "aq")$state[1], "aq", info = info)

# 20220208
info <- "Adding basis species to an existing definition works"
basis(c("iron", "oxygen", "H2O", "H+"), c(0, -80, 1, -7))
species(c("Fe", "Fe+2", "Fe+2"), c(0, -6, -6))
a1 <- affinity(pH = c(0, 14))
expect_error(basis(c("iron"), add = TRUE), "this species is already in the basis definition", info = info)
expect_error(basis(c("iron", "oxygen", "H2O", "H+"), add = TRUE), "these species are already in the basis definition", info = info)
expect_error(basis("Fe+2", add = TRUE), "the number of basis species is greater than the number of elements and charge", info = info)
expect_silent(newbasis <- basis(c("H2S", "Cl-"), c(-5, -3), add = TRUE), info = info)
expect_equal(newbasis$logact[5:6], c(-5, -3), info = info)
newspecies <- species()
expect_equal(colnames(newspecies)[5:6], c("H2S", "Cl-"), info = info)
expect_equal(newspecies$H2S, numeric(3), info = info)
a2 <- affinity(pH = c(0, 14))
expect_identical(a1$values, a2$values, info = info)
