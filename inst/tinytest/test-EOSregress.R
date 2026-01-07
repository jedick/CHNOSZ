# Load default settings for CHNOSZ
reset()

info <- "EOSvar stops with unknown variables"
if(exists("TX")) rm("TX")
expect_error(EOSvar("TX", T = 25, P = 1), "can't find a variable named TX", info = info)
# TOOD: Why aren't these variables visible to the test runner?
#TX <- 2
#expect_error(EOSvar("TX", T = 25, P = 1), "an object named TX is not a function", info = info)
#TX <- function(T) 2
#expect_error(EOSvar("TX", T = 25, P = 1), "the arguments of TX\\(\\) do not contain T and P", info = info)

info <- "Regressions return known HKF parameters (neutral species)"
# Regress computed values of heat capacity and volume of CH4(aq)
# calculated from HKF parameters on a T-P grid
T <- convert(seq(0, 350, 25), "K")
P <- seq(200, 1000, 100)
# We use calories to compare with OBIGT data for CH4(aq) 20220325
T.units("K")
E.units("cal")
CH4.prop <- subcrt("CH4", T = T, P = P, grid = "T")$out[[1]]
# Terms in the HKF equations for Cp
Cp.var <- c("invTTheta2", "TXBorn")
# Get coefficients in Cp regression
Cp.lm <- EOSregress(CH4.prop[, c("T", "P", "Cp")], Cp.var)
Cp.coeff <- Cp.lm$coefficients
# Terms in the HKF equations for V
V.var <- c("invPPsi", "invTTheta", "invPPsiTTheta", "QBorn")
# Get coefficients in V regression
V.lm <- EOSregress(CH4.prop[, c("T", "P", "V")], V.var)
# Use same units as HKF: convert from cm3.bar to joules (divide by 10) then to calories
V.coeff <- convert(convert(V.lm$coefficients, "joules"), "cal")
## The tests: did we get the HKF parameters that are in the database?
CH4.par <- info(info("CH4"))
# c1 and c2
expect_equivalent(Cp.coeff[1], CH4.par$c1, info = info)
expect_equivalent(Cp.coeff[2], CH4.par$c2, info = info)
# omega (from Cp)
expect_equivalent(Cp.coeff[3], CH4.par$omega, info = info)
# a1, a2, a3 and a4
expect_equivalent(V.coeff[1], CH4.par$a1, info = info)
expect_equivalent(V.coeff[2], CH4.par$a2, info = info)
expect_equivalent(V.coeff[3], CH4.par$a3, info = info)
expect_equivalent(V.coeff[4], CH4.par$a4, info = info)
# omega (from V)
expect_equivalent(V.coeff[5], CH4.par$omega, info = info)

# Tests added on 20260107

info <- "Cp_s_var and V_s_var give expected values"
# This should work for any omega.PrTr != 0
expect_equal(Cp_s_var(omega.PrTr = 1), water("XBorn")[[1]] * 298.15, tolerance = 1e-7, scale = 1)
expect_equal(V_s_var(omega.PrTr = 1), -water("QBorn")[[1]], tolerance = 1e-9, scale = 1)

info <- "EOScalc() and EOScoeffs() run without error"

# Calculate the Cp at some temperature and pressure
expect_silent(EOScalc(Cp.lm$coefficients, 298.15, 1), info = info)
# Get the database values of c1, c2 and omega for CH4(aq)
expect_silent(Cp_coeffs <- EOScoeffs("CH4", "Cp"), info = info)
# Volume parameters
expect_silent(V_coeffs <- EOScoeffs("CH4", "V"), info = info)

info <- "EOSplot() runs without error"

## Make plots comparing the regressions
## with the accepted EOS parameters of CH4
# Read the data from Hnedkovsky and Wood, 1997
f <- system.file("extdata/misc/HW97_Cp.csv", package = "CHNOSZ")
d <- read.csv(f)
# Use data for CH4
d <- d[d$species == "CH4", ]
d <- d[, -1]
# Convert J to cal and MPa to bar
d$Cp <- convert(d$Cp, "cal")
d$P <- convert(d$P, "bar")

# Save plot to a temporary PNG file
pngfile <- tempfile(fileext = ".png")
png(pngfile)
expect_silent(EOSplot(d, T.max = 600), info = info)
# Close the graphics device and remove the temporary PNG file
dev.off()
file.remove(pngfile)
