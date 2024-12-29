# Load default settings for CHNOSZ
reset()

# The maximum absolute pairwise difference between x and y
maxdiff <- function(x, y) max(abs(y - x))

info <- "Unbalanced reactions give a warning or are balanced given sufficient basis species"
expect_warning(subcrt(c("glucose", "ethanol"), c(-1, 3)), "reaction among glucose,ethanol was unbalanced, missing H-6O3", info = info)
basis("CHNOS")
s <- subcrt(c("malic acid", "citric acid"), c(-1, 1))
expect_equal(s$reaction$coeff, c(-1, 1, -2, -1, 1.5), info = info)
expect_equal(s$reaction$name, c("malic acid", "citric acid", "CO2", "water", "oxygen"), info = info)

info <- "Standard Gibbs energies (kJ/mol) of reactions involving aqueous species are consistent with the literature"
# From Amend and Shock, 2001 [AS01] Table 3
T <- c(2, 18, 25, 37, 45, 55, 70, 85, 100, 115, 150, 200)

# H2O(l) = H+ + OH-
AS01.H2O <- c(78.25, 79.34, 79.89, 80.90, 81.63, 82.59, 84.13, 85.78, 87.55, 89.42, 94.22, 102.21)
sout.H2O <- subcrt(c("H2O", "H+", "OH-"), c(-1, 1, 1), T = T)$out
# Tolerances set to lowest order of magnitute to pass
expect_true(maxdiff(sout.H2O$G/1000, AS01.H2O) < 0.01, info = info)

# AS01 Table 4.3 Reaction A1: H2(aq) + 0.5O2(aq) = H2O(l)
AS01.A1 <- c(-263.94, -263.45, -263.17, -262.62, -262.20, -261.63, -260.67, -259.60, -258.44, -257.18, -253.90, -248.44)
sout.A1 <- subcrt(c("H2", "O2", "H2O"), "aq", c(-1, -0.5, 1), T = T)$out
expect_true(maxdiff(sout.A1$G/1000, AS01.A1) < 0.01, info = info)

# AS01 Table 6.3 Reaction C7: 5S2O3-2 + H2O(l) + 4O2(aq) = 6SO4-2 + 2H+ + 4S(s)
AS01.C7 <- c(-1695.30, -1686.90, -1682.80, -1675.30, -1670.00, -1663.10, -1652.00, -1640.30, -1628.00, -1615.20, -1583.50, -1533.00)
s.C7 <- subcrt(c("S2O3-2", "H2O", "O2", "SO4-2", "H+", "S"), c("aq", "liq", "aq", "aq", "aq", "cr"), c(-5, -1, -4, 6, 2, 4), T = T)
sout.C7 <- s.C7$out
expect_true(maxdiff(sout.C7$G/1000, AS01.C7) < 0.06, info = info)
# We can also check that sulfur has expected polymorphic transitions
expect_equal(s.C7$polymorphs$sulfur, c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3), info = info)

info <- "Calculations using IAPWS-95 are possible"
oldwat <- water("IAPWS95")
# The tests pass if we get numeric values and no error
expect_silent(sb <- subcrt(c("H2O", "Na+"), T = c(-30, -20, 0, 10), P = 1)$out, info = info)
expect_true(all(sb$`Na+`$G < sb$water$G), info = info)
# Clean up
water(oldwat)

info <- "Phase stability limits give expected results"
expect_message(subcrt("gold", T = c(1300, 1350), P = 1, convert = FALSE), "subcrt: setting G to NA for gold\\(cr\\)", info = info)
# The reaction coefficients in the output should be unchanged 20171214
expect_equal(subcrt(c("bunsenite", "nickel", "oxygen"), c(-1, 1, 0.5))$reaction$coeff, c(-1, 1, 0.5), info = info) 
# Properties are NA only above (not at) the temperature limit for phase stability 20191111
sout <- subcrt("covellite", T = seq(780, 781, 0.5), P = 1, convert = FALSE)$out[[1]]
expect_equal(is.na(sout$G), c(FALSE, FALSE, TRUE), info = info)

# Use calories for comparisons with SUPCRT92
E.units("cal")

info <- "Calculations for K-feldspar are consistent with SUPCRT92"
# Use the superseded Helgeson et al., 1978 data
add.OBIGT("SUPCRT92", "K-feldspar")
T <- c(100, 100, 1000, 1000)
P <- c(5000, 50000, 5000, 50000)
SUPCRT_G <- c(-886628, -769531, -988590, -871493)
SUPCRT_H <- c(-932344, -815247, -865868, -748771)
SUPCRT_S <- c(62.6, 62.6, 150.6, 150.6)
SUPCRT_V <- c(108.9, 108.9, 108.9, 108.9)
SUPCRT_Cp <- c(56.7, 56.7, 80.3, 80.3)
CHNOSZ <- subcrt("K-feldspar", T = T, P = P)$out[[1]]
expect_equal(round(CHNOSZ$G), SUPCRT_G, info = info)
expect_equal(round(CHNOSZ$H), SUPCRT_H, info = info)
expect_equal(round(CHNOSZ$S, 1), SUPCRT_S, info = info)
expect_equal(round(CHNOSZ$V, 1), SUPCRT_V, info = info)
expect_equal(round(CHNOSZ$Cp, 1), SUPCRT_Cp, info = info)
OBIGT()

# Quartz tests re-activated after fixing dPdTtr to use parameters converted to Joules 20240211
info <- "Calculations for quartz are nearly consistent with SUPCRT92"
add.OBIGT("SUPCRT92")
# Using SUPCRT's equations, the alpha-beta transition occurs at
# 705 degC at 5000 bar and 1874 degC at 50000 bar,
# so here beta-quartz is stable only at T = 1000, P = 5000
T <- c(100, 100, 1000, 1000)
P <- c(5000, 50000, 5000, 50000)
SUPCRT_G <- c(-202778, -179719, -223906, -199129)
SUPCRT_H <- c(-214133, -191708, -199359, -177118)
SUPCRT_S <- c(12.3, 10.6, 31.8, 29.8)
SUPCRT_V <- c(22.5, 20.3, 23.7, 21.9)
SUPCRT_Cp <- c(12.3, 12.3, 16.9, 16.9)
CHNOSZ <- subcrt("quartz", T = T, P = P)$out[[1]]
# NOTE: It appears that, where alpha-quartz is stable above Ttr(Pr) but below Ttr(P),
# SUPCRT92 computes the heat capacity and its integrals using parameters of beta-quartz.
# (see e.g. the equation for CprdT under the (Cpreg .EQ. 2) case in the Cptrms subroutine of SUPCRT).
# ... is that incorrect?
# CHNOSZ's behavior is different - it only uses the lower-T polymorph below Ttr(P);
# so we get different values from SUPCRT at T = 1000, P = 50000 (4th condition) where alpha-quartz is stable
# (above Ttr@1 bar (575 degC), but below Ttr@50000 bar)
expect_equal(round(CHNOSZ$G)[-4], SUPCRT_G[-4], info = info)
expect_equal(round(CHNOSZ$H)[-4], SUPCRT_H[-4], info = info)
expect_equal(round(CHNOSZ$S, 1)[-4], SUPCRT_S[-4], info = info)
expect_equal(round(CHNOSZ$Cp, 1)[-4], SUPCRT_Cp[-4], info = info)
expect_equal(round(CHNOSZ$V, 1), SUPCRT_V, info = info)
OBIGT()

info <- "More calculations for quartz are nearly consistent with SUPCRT92"
add.OBIGT("SUPCRT92")
# Output from SUPCRT92 for reaction specified as "1 QUARTZ" run at 1 bar
# (SUPCRT shows polymorphic transition at 574.850 deg C, and does not give Cp values around the transition)
S92_1bar <- read.table(header = TRUE, text = "
    T       G       H    S       V
  573	-214507	-209517	24.7	23.3
  574	-214532	-209499	24.8	23.3
  575	-214557	-209192	25.1	23.7
  576	-214582	-209176	25.1	23.7
")
CHNOSZ_1bar <- subcrt("quartz", T = seq(573, 576), P = 1)$out[[1]]
expect_equal(round(CHNOSZ_1bar$G), S92_1bar$G, info = info)
expect_equal(round(CHNOSZ_1bar$H), S92_1bar$H, info = info)
expect_equal(round(CHNOSZ_1bar$S, 1), S92_1bar$S, info = info)
expect_equal(round(CHNOSZ_1bar$V, 1), S92_1bar$V, info = info)

# 5000 bar: SUPCRT shows polymorphic transition at 704.694 deg C
S92_5000bar <- read.table(header = TRUE, text = "
    T       G       H    S       V
  703	-215044	-204913	26.7	23.3
  704	-215071	-204895	26.7	23.3
  705	-215142	-204254	27.4	23.7
  706	-215170	-204238	27.5	23.7
")
CHNOSZ_5000bar <- subcrt("quartz", T = seq(703, 706), P = 5000)$out[[1]]
# NOTE: calculated values *below* the transition are different
expect_equal(CHNOSZ_5000bar$G, S92_5000bar$G, tolerance = 20, scale = 1, info = info)
expect_equal(CHNOSZ_5000bar$H, S92_5000bar$H, tolerance = 300, scale = 1, info = info)
expect_equal(CHNOSZ_5000bar$S, S92_5000bar$S, tolerance = 0.5, scale = 1, info = info)
expect_equal(CHNOSZ_5000bar$V, S92_5000bar$V, tolerance = 0.05, scale = 1, info = info)
OBIGT()

info <- "Duplicated species yield correct polymorphic transitions"
# If a mineral with polymorphic transitions is in both the basis and species lists,
# energy()'s call to subcrt() will have duplicated species.
# This wasn't working (produced NAs at low T) for a long time prior to 20171003.
s1 <- subcrt("chalcocite", T = c(100, 1000), P = 1000)
s2 <- subcrt(rep("chalcocite", 2), T = c(100, 1000), P = 1000)
expect_equal(s1$out[[1]]$logK, s2$out[[1]]$logK, info = info)
expect_equal(s1$out[[1]]$logK, s2$out[[2]]$logK, info = info)
## Another way to test it ...
basis(c("copper", "chalcocite"))
species("chalcocite")
a <- affinity(T = c(0, 1000, 2), P = 1)
expect_equal(as.numeric(a$values[[1]]), c(0, 0), info = info)

info <- "Reaction coefficients for repeated species are handled correctly"
# These were failing in version 1.1.3
s1 <- subcrt(c("quartz", "SiO2"), c(-1, 1))            
expect_equal(s1$reaction$coeff, c(-1, 1), info = info)
s2 <- subcrt(c("pyrrhotite", "pyrrhotite"), c(-1, 1))            
expect_equal(s2$reaction$coeff, c(-1, 1), info = info)
# These were failing in version 1.1.3-28
s3 <- subcrt(c("SiO2", "SiO2"), c(-1, 1))
expect_equal(s3$reaction$coeff, c(-1, 1), info = info)
s4 <- subcrt(c("H2O", "H2O", "H2O", "H2O", "H2O"), c(-2, 1, -3, 1, 3))
expect_equal(s4$reaction$coeff, c(-2, 1, -3, 1, 3), info = info)
# The reaction properties here should add up to zero
expect_equal(unique(s4$out$logK), 0, info = info)

info <- "Properties of HKF species below 0.35 g/cm3 are NA and give a warning"
wtext <- "below minimum density for applicability of revised HKF equations \\(2 T,P pairs\\)"
expect_warning(s1 <- subcrt(c("Na+", "quartz"), T = 450, P = c(400, 450, 500)), wtext, info = info) 
expect_equal(sum(is.na(s1$out$`Na+`$logK)), 2, info = info)
expect_equal(sum(is.na(s1$out$quartz$logK)), 0, info = info)
# Use exceed.rhomin to go below the minimum density
s2 <- subcrt(c("Na+", "quartz"), T = 450, P = c(400, 450, 500), exceed.rhomin = TRUE)
expect_equal(sum(is.na(s2$out$`Na+`$logK)), 0, info = info)

info <- "Combining minerals with polymorphic transitions and aqueous species with IS > 0 does not mangle output"
# s2 was giving quartz an extraneous loggam column and incorrect G and logK 20181107
add.OBIGT("SUPCRT92")
s1 <- subcrt(c("quartz", "K+"), T = 25, IS = 1)
s2 <- subcrt(c("K+", "quartz"), T = 25, IS = 1)
expect_true(identical(colnames(s1$out[[1]]), c("T", "P", "rho", "logK", "G", "H", "S", "V", "Cp", "polymorph")), info = info)
expect_true(identical(colnames(s2$out[[2]]), c("T", "P", "rho", "logK", "G", "H", "S", "V", "Cp", "polymorph")), info = info)
expect_true(identical(colnames(s1$out[[2]]), c("T", "P", "rho", "logK", "G", "H", "S", "V", "Cp", "loggam", "IS")), info = info)
expect_true(identical(colnames(s2$out[[1]]), c("T", "P", "rho", "logK", "G", "H", "S", "V", "Cp", "loggam", "IS")), info = info)
# Another one ... pyrrhotite was getting a loggam
expect_true(identical(colnames(subcrt(c("iron", "Na+", "Cl-", "OH-", "pyrrhotite"), T = 25, IS = 1)$out$pyrrhotite),
  c("T", "P", "rho", "logK", "G", "H", "S", "V", "Cp", "polymorph")), info = info)

info <- "Argument checking handles some types of invalid input"
expect_error(subcrt("H2O", -1, "liq", "xxx"), "invalid property name: xxx", info = info)
# Before version 1.1.3-63, having more than one invalid property gave a mangled error message
expect_error(subcrt("H2O", -1, "liq", c(1, 2)), "invalid property names: 1 2", info = info)

# Added on 20230620
info <- "Polymorphs are used by default"
sres_poly <- subcrt("pyrrhotite")
expect_equal(unique(sres_poly$out[[1]]$polymorph), c(1, 2, 3), info = info)
info <- "Polymorphs work for named species or numeric indices"
iPo <- info("pyrrhotite")
sres_poly1 <- subcrt(iPo)
expect_identical(sres_poly, sres_poly1, info = info)
info <- "Automatic identification of polymorphs can be turned off"
sres_nopoly <- subcrt("pyrrhotite", use.polymorphs = FALSE)
expect_null(sres_nopoly$out[[1]]$polymorph, info = info)
info <- "Warning is issued above beyond the transition temperature"
expect_warning(sres_nopoly <- subcrt("pyrrhotite", use.polymorphs = FALSE), "above T limit of 411 K", info = info)
info <- "Gibbs energy is silently extrapolated with exceed.Ttr = TRUE"
expect_silent(sres_nopoly_extrap <- subcrt("pyrrhotite", use.polymorphs = FALSE, exceed.Ttr = TRUE), info = info)
expect_false(anyNA(sres_nopoly_extrap$out[[1]]$G), info = info)

# Added on 20230621
info <- "Arguments 2 and 3 can't both be character"
expect_error(subcrt(c("hydrogen", "H2"), c("gas", "aq"), "G"), info = info)

# Added on 20230818
info <- "exceed.Ttr works for basis species in automatically balanced reactions"
basis(c("gypsum", "SO4-2", "H2O", "H+", "O2"))
# Defaults for subcrt() go above the temperature limit for gypsum, so use exceed.Ttr to calculate logK
automatic_reaction <- subcrt("Ca+2", 1, exceed.Ttr = TRUE)$out
expect_false(any(is.na(automatic_reaction$logK)), info = info)
# Check that logK is identical for the reaction entered manually
manual_reaction <- subcrt(c("gypsum", "Ca+2", "SO4-2", "H2O"), c(-1, 1, 1, 2), exceed.Ttr = TRUE)$out
expect_equal(automatic_reaction$logK, manual_reaction$logK, info = info)

# Added on 20231115
info <- "Cp equation limits give expected results"
expect_warning(sout <- subcrt("acanthite", T = 1000:1001, P = 1, convert = FALSE)$out[[1]], "above T limit of 1000 K", info = info)
expect_false(any(is.na(sout$G)), info = info)
info <- "exceed.Ttr doesn't interfere with polymorphic transitions"
# Stable polymorphs of pyrrhotite at default T,P conditions of subcrt()
polymorph <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3)
s1 <- subcrt("pyrrhotite")$out[[1]]
expect_equal(s1$polymorph, polymorph, info = info)
s2 <- subcrt("pyrrhotite", exceed.Ttr = TRUE)$out[[1]]
expect_equal(s2$polymorph, polymorph, info = info)

# Added on 20231204
info <- "Auto-balanced reactions apply ionic strength correction"
basis("CHNOS+")
sres <- subcrt("acetate", 1, IS = 1)
expect_length(sres$out$loggam, 15, info = info)

# Added on 20240206
info <- "High-temperature polymorph is not shown as stable below the transition temperature"
# This checks that the below-transition temperature code in subcrt() is working.
# If not, then cr2 is incorrectly identified as stable at 25 and 103 °C
# (that is, cr2 has a lower ΔG° than cr at those temperatures,
#  but should only be shown as stable above the transition temperature of 377 K)
T <- c(25, 50, 103, 104)
sout <- subcrt("carrollite", T = T, P = 1)$out[[1]]
expect_equal(sout$polymorph, c(1, 1, 1, 2), info = info)

info <- "Subzero degree C calculations are possible"
# Set default units
E.units("J")
# Start with H2O
s.H2O <- subcrt("H2O", T = c(-20.1, seq(-20, 0)), P = 1)$out$water
# We should get something at -20 deg C
expect_equal(s.H2O$G[2], convert(-56001, "J"), tolerance = 1, scale = 1, info = info)
# Following SUPCRT92, an input temperature of 0 is converted to 0.01
expect_equal(s.H2O$T[22], 0.01, info = info)
# The following test fails CRAN checks with Intel oneAPI 2023.x compilers
# (Expected TRUE, got FALSE) 20240211
if(!at_home()) exit_file("Skipping tests on CRAN")
# We shouldn't get anything at -20.1 deg C
expect_true(is.na(s.H2O$G[1]), info = info)

# References

# Amend, J. P. and Shock, E. L. (2001) 
#   Energetics of overall metabolic reactions of thermophilic and hyperthermophilic Archaea and Bacteria. 
#   FEMS Microbiol. Rev. 25, 175--243. https://doi.org/10.1016/S0168-6445(00)00062-0

# Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K. (1978) 
#   Summary and critique of the thermodynamic properties of rock-forming minerals. 
#   Am. J. Sci. 278-A, 1--229. http://www.worldcat.org/oclc/13594862
