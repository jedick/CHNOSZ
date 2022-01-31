# Load default settings for CHNOSZ
reset()

info <- "info.character() produces expected results and messages"
expect_equal(CHNOSZ:::info.character("acetate", "cr"), NA, info = info)
expect_message(CHNOSZ:::info.character("acetate", "cr"), "only 'aq' is available", info = info)
expect_message(CHNOSZ:::info.character("methane", "cr"), "only 'gas' 'liq' are available", info = info)
expect_message(CHNOSZ:::info.character("methane"), "also available in liq", info = info)
expect_message(CHNOSZ:::info.character("SiO2", "cr"), "also available in.*quartz", info = info)
expect_message(CHNOSZ:::info.character("chalcocite"), "found chalcocite\\(cr\\) with 2 phase transitions", info = info)
# H2O is a special case
expect_equal(CHNOSZ:::info.character("H2O", "aq"), CHNOSZ:::info.character("H2O", "liq"), info = info)

info <- "info.numeric() produces expected errors and messages"
expect_error(CHNOSZ:::info.numeric(9999), "species index 9999 not found in thermo\\(\\)\\$OBIGT", info = info)
iargon <- info("argon", "gas")
expect_message(CHNOSZ:::info.numeric(iargon), "Cp of argon\\(gas\\) is NA; set by EOS parameters to 4.97", info = info)
iMgSO4 <- info("MgSO4")
expect_message(CHNOSZ:::info.numeric(iMgSO4), "V of MgSO4\\(aq\\) is NA; set by EOS parameters to 1.34", info = info)

info <- "info.approx() produces expected messages"
expect_message(CHNOSZ:::info.approx("lactic"), "is similar to lactic acid", info = info)
expect_message(CHNOSZ:::info.approx("lactic acid"), "is ambiguous", info = info)
# note though that info("lactic acid") finds a match because info.character is used first...
expect_equal(info("lactic acid"), grep("lactic acid", thermo()$OBIGT$name), info = info)
# looking in optional databases 20190127
expect_message(info("dickite"), "is in an optional database", info = info)

info <- "info() can be used for cr and aq descriptions of the same species and proteins"
i2 <- info("LYSC_CHICK", c("cr", "aq")) 
expect_equal(thermo()$OBIGT$state[i2], c("cr", "aq"), info = info)
expect_equal(info(i2)[1, ], info(i2[1]), check.attributes = FALSE, info = info)

info <- "info() gives correct column names for species using the AkDi model"
# add an aqueous species conforming to the AkDi model
iCO2 <- mod.OBIGT("CO2", abbrv = "AkDi", a = -8.8321, b = 11.2684, c = -0.0850)
params <- info(iCO2)
expect_equal(params$a, -8.8321, info = info)
expect_equal(params$b, 11.2684, info = info)
expect_equal(params$xi, -0.0850, info = info)
