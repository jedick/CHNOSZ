# Load default settings for CHNOSZ
reset()

info <- "Non-zero ionic strength transforms variables from activity to molality"
# What happens with activity coefficients when using subcrt() to calculate affinity,
# and in the rest of the main workflow of CHNOSZ?
# 20171025 first version
# 20181106 include non-zero activity coefficient of CO2(aq)

# The maximum absolute pairwise difference between x and y
maxdiff <- function(x, y) max(abs(y - x))

### First get the activity coefficients of H+ and HCO3-
## The long way...
wprop <- water(c("A_DH", "B_DH"), P = 1)
speciesprops <- subcrt(c("H+", "HCO3-", "CO2"), T = 25)$out
nonid <- nonideal(c("H+", "HCO3-", "CO2"), speciesprops, IS = 1, T = 298.15, P = 1, A_DH = wprop$A_DH, B_DH = wprop$B_DH)
# Compare with a precalculated value:
expect_true(maxdiff(nonid[[2]]$loggam, -0.1868168) < 1e-7, info = info)
## The short way...
out1 <- subcrt(c("H+", "HCO3-", "CO2"), T = 25, IS = 1)$out
loggam_HCO3 <- out1[[2]]$loggam
loggam_CO2 <- out1[[3]]$loggam
expect_equal(nonid[[2]]$loggam, loggam_HCO3, info = info)
expect_equal(nonid[[3]]$loggam, loggam_CO2, info = info)
## Take-home message -1: with default settings, the activity coefficient of H+ is always 1

### How do activity coefficient affect the value of G?
# Let's step back and look at the *standard Gibbs energy* at IS = 0
out0 <- subcrt(c("H+", "HCO3-", "CO2"), T = 25)$out
# The adjusted standard Gibbs energy is less than the standard Gibbs energy
# by an amount determined by the activity coefficient
expect_equal(out1[[2]]$G - out0[[2]]$G, -convert(loggam_HCO3, "G"), info = info)
expect_equal(out1[[3]]$G - out0[[3]]$G, -convert(loggam_CO2, "G"), info = info)
## Take-home message 0: setting IS in subcrt() gives adjusted standard Gibbs energy

# What is the equilibrium constant for the reaction CO2 + H2O = H+ + HCO3-?
# (this is the standard state property at IS = 0)
logK <- subcrt(c("CO2", "H2O", "H+", "HCO3-"), c(-1, -1, 1, 1), T = 25)$out$logK
# We get logK = -6.34468 (rounded)
expect_true(maxdiff(logK, -6.34468) < 1e-6, info = info)

### What is the affinity of the reaction at pH=7 and molalities of HCO3- and CO2 = 10^-3?
## Case 1: ionic strength = 0, so gamma = 0 and activity = molality
# First calculate it by hand from 2.303RTlog(K/Q)
# logQ = (logaH+ + logaHCO3-) - (logaH2O + logaCO2)
logQ0 <- (-7 + -3) - (0 + -3) # i.e. 7
# Convert logs to energy according to G = -2.303RTlogK
# (we negate that to get affinity)
A0manual <- -convert(logK - logQ0, "G")
# Now the run calculation with subcrt
A0subcrt <- subcrt(c("CO2", "H2O", "H+", "HCO3-"), c(-1, -1, 1, 1), T = 25, logact = c(-3, 0, -7, -3))$out$A
# We get the same affinity!
expect_equal(A0subcrt, A0manual, info = info)
## Case 2: ionic strength = 1, so activity = molality * gamma
logaHCO3 <- -3 + loggam_HCO3
logaCO2 <- -3 + loggam_CO2
logQ1 <- (-7 + logaHCO3) - (0 + logaCO2)
A1manual <- -convert(logK - logQ1, "G")
A1subcrt <- subcrt(c("CO2", "H2O", "H+", "HCO3-"), c(-1, -1, 1, 1), T = 25, logact = c(-3, 0, -7, -3), IS = 1)$out$A
expect_equal(A1subcrt, A1manual, info = info)
## Take-home message 1: using subcrt with IS not equal to zero, the "logact"
## argument is logmolal in affinity calculations for charged aqueous species

### Now calculate the affinities using affinity()
basis("CHNOS+")  # pH = 7, logaCO2 = -3
species(c("CO2", "HCO3-"))  # logactivities = -3
## Case 1: IS = 0
a0 <- affinity()
# That gives us values in log units; convert to energy
# (HCO3- is species #2)
A0affinity <- -convert(a0$values[[2]], "G")
expect_equal(A0affinity[[1]], A0subcrt, info = info)
## Case 2: IS = 1
a1 <- affinity(IS = 1)
A1affinity <- -convert(a1$values[[2]], "G")
expect_equal(A1affinity[[1]], A1subcrt, info = info)
## Take-home message 2: using affinity() with IS not equal to zero, the "logact"
## set by species() is logmolal in affinity calculations for charged aqueous species

### Now, swap HCO3- for CO2 in the basis
swap.basis("CO2", "HCO3-")
basis("HCO3-", -3)
a0 <- affinity()
a1 <- affinity(IS = 1)
# Look at HCO3 formation affinity:
# they're both zero at IS = 0 or 1
expect_equal(a0$values[[2]][1], 0, info = info)
expect_equal(a1$values[[2]][1], 0, info = info)
# Look at CO2 formation affinity:
ACO2_0affinity <- -convert(a0$values[[1]], "G")
ACO2_1affinity <- -convert(a1$values[[1]], "G")
# What's that?? we're looking at the reverse of the above
# reaction, i.e. H+ + HCO3- = CO2 + H2O
# so, logK = 6.345
logKrev <- -logK
logQrev0 <- -logQ0
logQrev1 <- -logQ1
ACO2_0manual <- -convert(logKrev - logQrev0, "G")
ACO2_1manual <- -convert(logKrev - logQrev1, "G")
expect_equal(ACO2_0manual, ACO2_0affinity[[1]], info = info)
expect_equal(ACO2_1manual, ACO2_1affinity[[1]], info = info)
## Take-home message 3: using affinity() with IS not equal to zero, the "logact"
## set by basis() is logmolal in affinity calculations for charged aqueous species

### Now look at equilibrate()
e0 <- equilibrate(a0)
e1 <- equilibrate(a1)
# Using the equilibrated values, calculate affinity of the reaction CO2 + H2O = H+ + HCO3-
# Case 1: IS = 0
logact_HCO3 <- e0$loga.equil[[2]]
logact_CO2 <- e0$loga.equil[[1]]
logQeq0 <- (-7 + logact_HCO3) - (logact_CO2 + 0)
Aeq0 <- -convert(logK - logQeq0, "G") # zero!
expect_equal(Aeq0[[1]], 0, info = info)
# Case 2: IS = 1
logact_HCO3 <- e1$loga.equil[[2]]
logact_CO2 <- e1$loga.equil[[1]]
# Here, loga.equil is the *molality*, so we must multiply by loggam
logQeq1 <- (-7 + logact_HCO3 + loggam_HCO3) - (logact_CO2 + loggam_CO2 + 0)
Aeq1 <- -convert(logK - logQeq1, "G") # zero!
expect_equal(Aeq1[[1]], 0, info = info)
## Take-home message 4: using affinity() with IS not equal to zero, the "loga.equil"
## returned by equilibrate() is logmolal for speciation calculations with charged aqueous species

# Finally, what is loga.balance?
a.balance <- 10^e1$loga.balance
m.total <- sum(10^unlist(e1$loga.equil))
expect_equal(a.balance, m.total, info = info)
## Take-home message 5: using affinity() with IS not equal to zero, the "loga.balance"
## used by equilibrate() is the logarithm of total molality of the balancing basis species
