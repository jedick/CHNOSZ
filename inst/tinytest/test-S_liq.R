# test-S_liq.R
# Fit JANAF data for S(l) for OBIGT
# 20260612 jmd

# Source data: https://janaf.nist.gov/tables/S-004.html

# Load data
Gf_298 <- 432     # J/mol
Hf_298 <- 1854    # J/mol
S_298  <- 36.825  # J/K/mol
Cp_298 <- 22.531  # J/K/mol

# Tables for Cp fits
T1_vals <- c(298.15, 300, 350, 388.360, 400, 432.020)
Cp1_vals <- c(22.531, 22.707, 27.434, 31.058, 32.162, 53.829)
T2_vals <- c(432.020, 450, 500, 600, 700, 800, 882.117)
Cp2_vals <- c(53.808, 43.046, 37.986, 34.308, 32.681, 31.699, 31.665)

# Entropy values for checking results
S1_vals <- c(36.825, 36.965, 40.821, 43.859, 44.793, 47.431)
S2_vals <- c(47.431, 49.308, 53.532, 60.078, 65.241, 69.530, 72.624)

# Fit Maier-Kelley Cp equation
# Cp = a + b*T + c*T^-2
# Pre-lambda
T1_2 <- T1_vals ^ -2
T1_0.5 <- T1_vals ^ -0.5
Cp1_lm <- lm(Cp1_vals ~ T1_vals + T1_2 + T1_0.5)
coeff1 <- signif(Cp1_lm$coefficients, 5)
# Post-lambda
T2_2 <- T2_vals ^ -2
T2_0.5 <- T2_vals ^ -0.5
Cp2_lm <- lm(Cp2_vals ~ T2_vals + T2_2 + T2_0.5)
coeff2 <- signif(Cp2_lm$coefficients, 5)

# Format today's date for ISO 8601, e.g. 2026-06-04
date <- format(Sys.time(), "%Y-%m-%d")
# Calculate molar volume (molecular weight / density)
V <- round(32.06 / 1.80, 2)

#############################################
### The following is adapted from FAQ.Rmd ###
#############################################

# The formula of the new mineral and literature reference
formula <- "S"
ref1 <- "JANAF98.2"
ref2 <- "OBIGT26.1"
# Thermodynamic properties of polymorph 1 at 25 °C (298.15 K)
G1 <- Gf_298
H1 <- Hf_298
S1 <- S_298
Cp1 <- Cp_298
# Heat capacity coefficients for polymorph 1
a1 <- coeff1[1]
b1 <- coeff1[2]
c1 <- coeff1[3]
d1 <- coeff1[4]
# Use same molar volume for both polymorphs
V1 <- V2 <- V
# Transition temperature
Ttr_val <- max(T1_vals)
# Transition enthalpy (J/mol)
DHtr <- 0
# Heat capacity coefficients for polymorph 2
a2 <- coeff2[1]
b2 <- coeff2[2]
c2 <- coeff2[3]
d2 <- coeff2[4]
# Maximum temperature of polymorph 2
T2 <- max(T2_vals)

# Use the temperature (Ttr) and enthalpy of transition (DHtr) to calculate the entropy of transition (DStr)
DGtr <- 0
TDStr <- DHtr - DGtr
DStr <- TDStr / Ttr_val

# The name of the new species
name <- "S_liq"

# Start new database entries
mod.OBIGT(name, abbrv = NA, formula = formula, state = "cr", ref1 = ref1, ref2 = ref2, date = date,
  model = "CGL", E_units = "J", G = 0, H = 0, S = 0, V = V1, Cp = Cp1,
  a = a1, b = b1, c = c1, d = d1, e = 0, f = 0, lambda = 0, T = Ttr_val)
# Use CGL_Ttr for the high-T polymorph to indicate that it undergoes a phase transition (i.e. vaporization) at T2
mod.OBIGT(name, abbrv = NA, formula = formula, state = "cr2", ref1 = ref1, ref2 = ref2, date = date,
  model = "CGL_Ttr", E_units = "J", G = 0, H = 0, S = 0, V = V2,
  a = a2, b = b2, c = c2, d = d2, e = 0, f = 0, lambda = 0, T = T2)

# Use K and J as the temperature and energy units for subcrt() (J is already the default)
T.units("K")
E.units("J")

# Calculate the entropy change of each polymorph from 298.15 to Ttr
DS1 <- subcrt(name, "cr", P = 1, T = Ttr_val, use.polymorphs = FALSE)$out[[1]]$S
DS2 <- subcrt(name, "cr2", P = 1, T = Ttr_val)$out[[1]]$S
# Obtain the difference of entropy between the polymorphs at 298.15 K
DS298 <- round(DS1 + DStr - DS2, 3)

# Put the values of S° at 298.15 into OBIGT, then calculate the changes
# of all thermodynamic properties of each polymorph between 298.15 K and Ttr
mod.OBIGT(name, state = "cr", S = S1)
mod.OBIGT(name, state = "cr2", S = S1 + DS298)
D1 <- subcrt(name, "cr", P = 1, T = Ttr_val, use.polymorphs = FALSE)$out[[1]]
D2 <- subcrt(name, "cr2", P = 1, T = Ttr_val)$out[[1]]

# Add up the contributions to ΔG°f and ΔH°f
DG298 <- round(D1$G + DGtr - D2$G, 2)
DH298 <- round(D1$H + DHtr - D2$H, 2)
# Make final changes to database
mod.OBIGT(name, state = "cr", G = G1, H = H1)
mod.OBIGT(name, state = "cr2", G = G1 + DG298, H = H1 + DH298)

if(FALSE){
  # Plot Cp and G of polymorphs
  T <- c(T1_vals, T2_vals)
  Tlow <- T[T < 600]
  sout <- subcrt(c("S_liq", "S_liq"), c("cr", "cr2"), T = Tlow, P = 1, use.polymorphs=FALSE)$out

  # Start plot with G
  par(mfrow = c(1, 2))
  plot(Tlow, sout[[2]]$G, type = "b", col = 2, ylab = "G")
  lines(Tlow, sout[[1]]$G, type = "b")
  legend("topright", c("cr", "cr2"), pch = 1, lty = 1, col = c(1, 2))
  plot(Tlow, sout[[2]]$Cp, type = "b", col = 2, ylab = "Cp", ylim = c(20, 80))
  lines(Tlow, sout[[1]]$Cp, type = "b")
  # Add experimental Cp
  points(T, c(Cp1_vals, Cp2_vals), pch = 19)
}

# Tests adapted from FAQ.Rmd
info <- "Entropy of transition is calculated correctly"
expect_equal(D2$S - D1$S, DStr, tolerance = 1e-3, scale = 1, info = info)
info <- "G, H, and S at 25 °C for each polymorph are self-consistent"
cr_parameters <- info(info(name, "cr"))
expect_true(abs(check.GHS(cr_parameters)) < 1, info = info)
cr2_parameters <- info(info(name, "cr2"))
expect_true(abs(check.GHS(cr2_parameters)) < 1, info = info)

# Tests added on 20260613

# Compute properties at T values in JANAF table
T <- c(T1_vals, T2_vals)
sout <- subcrt("S_liq", T = T, P = 1)$out[[1]]

info <- "Polymorphs are identified correctly"
expect_equal(sout$polymorph, c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2), info = info)

info <- "MAE of Cp is less than 1 J/K/mol"
Cp <- c(Cp1_vals, Cp2_vals)
expect_true(mean(abs(sout$Cp - Cp)) < 1, info = info)

info <- "MAE of S° is less than 1 J/K/mol"
S <- c(S1_vals, S2_vals)
expect_true(mean(abs(sout$S - S)) < 1, info = info)

info <- "Database values match those calculated here"
# Get the values we made
iS_liq <- info("S_liq", c("cr", "cr2"))
our_values <- thermo()$OBIGT[iS_liq, ]
# Get the default database values
reset()
db_values <- thermo()$OBIGT[iS_liq, ]
# Make the dates the same
our_values$date <- db_values$date
expect_equal(our_values, db_values, info = info)
