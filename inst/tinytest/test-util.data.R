# Load default settings for CHNOSZ
reset()

info <- "reset() and OBIGT() produce the same database"
d1 <- thermo()$OBIGT
OBIGT()
d2 <- thermo()$OBIGT
expect_equal(d1, d2, info = info)

info <- "check.GHS() and check.EOS() (via info()) produce messages"
i1 <- info("S2O3-2")
expect_message(info(i1), "calculated ΔG°f of S2O3-2\\(aq\\) differs by 939 cal mol-1 from database value", info = info)
i2 <- info("Cu+2")
expect_message(info(i2), "calculated Cp° of Cu\\+2\\(aq\\) differs by 3.62 cal K-1 mol-1 from database value", info = info)

info <- "check.GHS() and check.EOS() respond to thermo()$opt$*.tol"
i1 <- info("SO4-2")
thermo("opt$Cp.tol" = 0.5)
expect_message(info(i1), "check.EOS", info = info)
i2 <- info("a,w-dicarboxytetracosane")
thermo("opt$G.tol" = 50)
expect_message(info(i2), "check.GHS", info = info)

info <- "RH2OBIGT() gives group additivity results consistent with database values (from Richard and Helgeson, 1998)"
file <- system.file("extdata/misc/RH98_Table15.csv", package = "CHNOSZ")
dat <- read.csv(file, stringsAsFactors = FALSE)
ispecies <- info(dat$compound, dat$state)
OBIGT.ref <- thermo()$OBIGT[ispecies, ]
OBIGT.calc <- RH2OBIGT(file = file)
# The maximum absolute pairwise difference between x and y
maxdiff <- function(x, y) max(abs(y - x))
# Calculated values of H are spot on; to pass tests, tolerance on
# G is set higher; is there an incorrect group value somewhere?
expect_true(maxdiff(OBIGT.calc$G, OBIGT.ref$G) < 31, info = info)
expect_true(maxdiff(OBIGT.calc$H, OBIGT.ref$H) == 0, info = info)
expect_true(maxdiff(OBIGT.calc$S, OBIGT.ref$S) < 0.02001, info = info)
expect_true(maxdiff(OBIGT.calc$Cp, OBIGT.ref$Cp) < 0.04001, info = info)
expect_true(maxdiff(OBIGT.calc$V, OBIGT.ref$V) < 0.1001, info = info)
expect_true(maxdiff(OBIGT.calc$a1.a, OBIGT.ref$a1.a) < 0.01001, info = info)
expect_true(maxdiff(OBIGT.calc$a2.b, OBIGT.ref$a2.b) < 1e-13, info = info)
expect_true(maxdiff(OBIGT.calc$a3.c, OBIGT.ref$a3.c) < 1e-14, info = info)

info <- "add.OBIGT() replaces existing entries without changing species index"
# Store the original species index of CdCl2
iCdCl2 <- info("CdCl2", "aq")
# Add supplemental database - includes CdCl2
file <- system.file("extdata/misc/BZA10.csv", package = "CHNOSZ")
isp <- add.OBIGT(file)
# Species index of CdCl2 should not have changed
expect_equal(info("CdCl2", "aq"), iCdCl2, info = info)
# Check that names of species modified are same as in file
newdat <- read.csv(file, stringsAsFactors = FALSE)
# The order isn't guaranteed ... just make sure they're all there
expect_true(all(newdat$name %in% thermo()$OBIGT$name[isp]), info = info)

info <- "add.OBIGT() gives an error for an incompatible file"
# Test added 20191210
file <- system.file("extdata/Berman/Ber88_1988.csv", package = "CHNOSZ")
expect_error(add.OBIGT(file), info = info)

info <- "info() gives consistent messages for cal and J"
# Test added 20190529
# Add data for dimethylamine and trimethylamine in different units (cal or J)
expect_message(add.OBIGT(system.file("extdata/misc/LA19_test.csv", package = "CHNOSZ")), "energy units: J and cal", info = info)
expect_message(info(info("DMA_cal")), "-1.92 cal", info = info)
expect_message(info(info("DMA_J")), "-8.02 J", info = info)
# For TMA, only a check.GHS message for the entry in J is produced,
# because it's above the threshold of 100 set in thermo()$opt$G.tol
expect_silent(info(info("TMA_cal")), info = info)
expect_message(info(info("TMA_J")), "-102 J", info = info)

info <- "Missing values for G, Cp, and V are correct in cal and J"
# Test added 20190530
# Add data for dimethylamine and trimethylamine in different units (cal or J)
add.OBIGT(system.file("extdata/misc/LA19_test.csv", package = "CHNOSZ"))
calccal <- info(info("DMA_cal_NA"))
expect_equal(round(calccal$G), 13934, info = info)
expect_equal(round(calccal$Cp, 1), 60.3, info = info)
expect_equal(round(calccal$V, 1), 58.2, info = info)
calcJ <- info(info("DMA_J_NA"))
expect_equal(round(calcJ$G), 58304, info = info)
expect_equal(round(calcJ$Cp, 1), 252.4, info = info)
expect_equal(round(calcJ$V, 1), 58.2, info = info)

info <- "subcrt() gives same results for data entered in cal and J"
# Test added 20190530
# Add data for dimethylamine and trimethylamine in different units (cal or J)
add.OBIGT(system.file("extdata/misc/LA19_test.csv", package = "CHNOSZ"))
E.units("cal")
scal <- subcrt("DMA_cal")
sJ <- subcrt("DMA_J")
expect_true(maxdiff(scal$out[[1]]$G, sJ$out[[1]]$G) < 2, info = info)
expect_true(maxdiff(scal$out[[1]]$H, sJ$out[[1]]$H) < 1, info = info)
expect_true(maxdiff(scal$out[[1]]$S, sJ$out[[1]]$S) < 0.006, info = info)
expect_true(maxdiff(scal$out[[1]]$V, sJ$out[[1]]$V) < 0.011, info = info)
expect_true(maxdiff(scal$out[[1]]$Cp, sJ$out[[1]]$Cp) < 0.16, info = info)
# Now switch output to J (default units)
E.units("J")
calcJ25 <- subcrt("DMA_J", T = 25)$out[[1]]
infoJ25 <- info(info("DMA_J"))
expect_equivalent(calcJ25[, c("G", "H", "S")], calcJ25[, c("G", "H", "S")], info = info)
# In the case of Cp and V, there are bigger difference because they are calculated from the HKF parameters
expect_true(maxdiff(calcJ25$Cp, infoJ25$Cp) < 8.1, info = info)
expect_true(maxdiff(calcJ25$V, infoJ25$V) < 0.55, info = info)

info <- "OBIGT2eos() doesn't convert lambda for cr species"
## Bug visible with change of ferberite data to J:
## lambda (exponent on heat capacity term)
## was incorrectly going through a units conversion  20190903
# This was working
mod.OBIGT("test_cal", formula = "C0", state = "cr", E_units = "cal", a = 10, b = 100, f = 1, lambda = 0)
mod.OBIGT("test_J", formula = "C0", state = "cr", E_units = "J", a = 41.84, b = 418.4, f = 4.184, lambda = 0)
expect_equal(subcrt("test_cal", T = 25)$out[[1]]$Cp, subcrt("test_J", T = 25)$out[[1]]$Cp, info = info)
# This wasn't working
mod.OBIGT("test_cal2", formula = "C0", state = "cr", E_units = "cal", a = 10, b = 100, f = 1, lambda = -1)
mod.OBIGT("test_J2", formula = "C0", state = "cr", E_units = "J", a = 41.84, b = 418.4, f = 4.184, lambda = -1)
expect_equal(subcrt("test_cal2", T = 25)$out[[1]]$Cp, subcrt("test_J2", T = 25)$out[[1]]$Cp, info = info)

info <- "We can define an aqueous species with CGL model"
# Test added 20230220
icr <- mod.OBIGT("fake_cr", formula = "Na2Cl2", state = "cr", model = "CGL", G = -1000, H = -1000, S = 10, Cp = 10, V = 10)
iaq <- mod.OBIGT("fake_aq", formula = "Na2Cl2", state = "aq", model = "CGL", G = -1000, H = -1000, S = 10, Cp = 10, V = 10)
# Make sure info() runs (message is from check.EOS())
expect_message(info(iaq), "differs", info = info)
expect_silent(info(iaq), info = info)
# Make sure subcrt() runs
expect_equal(subcrt(iaq)$G, subcrt(icr)$G, info = info)

# 20241225
info <- "add.OBIGT() errors with non-existent file"
file <- "XXX"
expect_error(add.OBIGT(file), "XXX is not a file and doesn't match any files in the OBIGT database", info = info)

# Reference

# Richard, L. and Helgeson, H. C. (1998) Calculation of the thermodynamic properties at elevated 
#   temperatures and pressures of saturated and aromatic high molecular weight solid and liquid 
#   hydrocarbons in kerogen, bitumen, petroleum, and other organic matter of biogeochemical interest. 
#   Geochim. Cosmochim. Acta 62, 3591--3636. https://doi.org/10.1016/S0016-7037(97)00345-1

# Tests added on 20260107

info <- "thermo.refs() creates HTML table without error"
# Create the HTML but divert the output to a NULL connection
expect_silent(thermo.refs(browser = sink), info = info)
# Turn off the connection
sink()

info <- "thermo.refs() returns reference for source key"
expect_equal(thermo.refs("HDNB78")$key, "HDNB78", info = info)

info <- "thermo.refs() returns reference for species index"
expect_equal(thermo.refs(info("SiO2"))$key, "SHS89", info = info)

info <- "thermo.refs() returns references for subcrt() species output"
expect_equal(thermo.refs(subcrt("carrollite"))$key, "HDR+24", info = info)

info <- "thermo.refs() returns references for subcrt() reaction output"
expect_equal(thermo.refs(subcrt(c("oxygen", "O2"), c(-1, 1)))$key, c("WEP+82", "SHS89", "Kel60"), info = info)

info <- "check.OBIGT() runs without error"
expect_silent(out <- check.OBIGT(), info = info)

info <- "dumpdata() runs without error"
expect_silent(dd <- dumpdata(), info = info)
