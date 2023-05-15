# Load default settings for CHNOSZ
reset()

info <- "EOSvar stops with unknown variables"
expect_error(EOSvar("TX", T = 25, P = 1), "can't find a variable named TX", info = info)
# Why can't the test find these?
#TX <- 2
#expect_error(EOSvar("TX", T = 25, P = 1), "an object named TX is not a function")
#TX <- function(T) 2
#expect_error(EOSvar("TX", T = 25, P = 1), "the arguments of TX\\(\\) are not T, P")

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
