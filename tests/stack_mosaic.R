# Test variations on stacked mosaic diagrams culminating with stack_mosaic() wrapper function
# 20220722 v1 adapted from Mosaic Stacking 2 from multi-metal.Rmd

library(CHNOSZ)

# Use low resolution for routine package checking
res <- 20
## Increase resolution for nicer-looking PDF (the one kept in the test directory)
#res <- 200
# Define system
pH <- c(0, 14, res)
O2 <- c(-48, -33, res)
T <- 200
logmS <- -2
m_NaCl <- 0.1
logm_aq <- -6 # for both Fe- and Cu-bearing aq species
# Define basis species
S.aq <- c("H2S", "HS-", "HSO4-", "SO4-2")
# Define minerals
Fe.cr <- c("pyrite", "pyrrhotite", "magnetite", "hematite")
Fe.abbrv <- c("Py", "Po", "Mag", "Hem")
FeCu.cr <- c("chalcopyrite", "bornite")
Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
FeCu.abbrv <- c("Ccp", "Bn", "Cu", "Cpr", "Tnr", "Cct", "Cv")
# Define aqueous species
iFe.aq <- retrieve("Fe", c("S", "O", "H", "Cl"), "aq")
Fe.aq <- info(iFe.aq)$name
iCu.aq <- retrieve("Cu", c("S", "O", "H", "Cl"), "aq")
Cu.aq <- info(iCu.aq)$name
nacl <- NaCl(T = T, P = "Psat", m_tot = m_NaCl)

setup <- function() {
  reset()
  # Setup basis species
  basis(c("Cu+", "pyrite", "H2S", "oxygen", "H2O", "H+", "Cl-"))
  basis("H2S", logmS)
  basis("Cl-", log10(nacl$m_Cl))
}

ref_Cu_craq <- function() {
  # Load Cu-bearing minerals
  species(Cu.cr)
  # Add aqueous species 20210220
  species(iCu.aq, logm_aq, add = TRUE)

  mCu <- mosaic(list(S.aq), pH = pH, O2 = O2, T = T, IS = nacl$IS)
  diagram(mCu$A.species)
}

ref_FeCu_cr <- function() {
  # Load Fe-bearing minerals
  species(Fe.cr)
  mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T, IS = nacl$IS)
  dFe <- diagram(mFe$A.species, lwd = 0, names = FALSE, plot.it = FALSE)

  # Load Cu-bearing minerals
  species(c(FeCu.cr, Cu.cr))
  # Mosaic with all Fe species as basis species
  mFeCu <- mosaic(list(S.aq, Fe.cr), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dFe$predominant))

  diagram(mFeCu$A.species)
}

setup_FeCu <- function() {
  # Load Fe-bearing minerals
  species(Fe.cr)
  # Add aqueous species 20210220
  species(iFe.aq, logm_aq, add = TRUE)
  mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T, IS = nacl$IS)
  dFe <- diagram(mFe$A.species, lwd = 0, names = FALSE, plot.it = FALSE)

  # Load Cu-bearing minerals
  species(c(FeCu.cr, Cu.cr))
  # Add aqueous species 20210220
  species(iCu.aq, logm_aq, add = TRUE)
  
  # Return Fe diagram for stacked mosaic calculations
  dFe
}


test_FeCu_old <- function() {
  dFe <- setup_FeCu()

  # NOTE: limitation in mosaic() when using both solid and aqueous basis species (i.e. c(Fe.cr, Fe.aq)):
  #   Only the activity of the first-defined basis species is used for all basis species.
  #   The first-defined basis species is pyrite (logact = 0), but logact < 0 for the aq Fe species.
  # Workaround: Adjust standard Gibbs energies of aq Fe species to virtually change their activities
  DG_J <- convert(-logm_aq, "G", T = convert(T, "K"))
  # We should use calories here because the database values are in calories 20220604
  stopifnot(all(info(iFe.aq)$E_units == "cal"))
  DG_cal <- convert(DG_J, "cal")
  G.orig <- info(iFe.aq)$G
  G.new <- G.orig + DG_cal
  mod.OBIGT(iFe.aq, G = G.new)

  # Mosaic with all Fe species as basis species
  mFeCu <- mosaic(list(S.aq, c(Fe.cr, Fe.aq)), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dFe$predominant))

  diagram(mFeCu$A.species)
}

test_FeCu_new <- function() {
  dFe <- setup_FeCu()

  # Mosaic with all Fe species as basis species
  # Use loga_aq argument to control the activity of aqueous species in mosaic calculation 20220722
  # c(NA, logm_aq) means to use:
  #   basis()'s value for logact of aqueous S species
  #   logm_aq for logact of aqueous Fe species
  mFeCu <- mosaic(list(S.aq, c(Fe.cr, Fe.aq)), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dFe$predominant), loga_aq = c(NA, logm_aq))

  diagram(mFeCu$A.species)
}

test_FeCu_stack <- function() {
  # Aqueous S species
  bases <- S.aq
  # Fe-bearing minerals
  species1 <- Fe.cr
  # Cu-bearing and FeCu-bearing minerals
  species2 <- Cu.cr
  species12 <- FeCu.cr

  # Fe-bearing aqueous species
  species1 <- c(species1, Fe.aq)
  # Cu-bearing aqueous species
  species2 <- c(species2, Cu.aq)

  sm <- stack_mosaic(bases, species1, species2, species12, pH = pH, O2 = O2, T = T, IS = nacl$IS, loga_aq = logm_aq, plot.it = FALSE)
  diagram(sm[[2]])
}

## Setup plot
pdf("stack_mosaic.pdf", width = 8, height = 9)
mat <- matrix(c(1,6,2, 6,3,6, 4,6,5), nrow = 3, byrow = TRUE)
layout(mat, widths = c(2, 3, 2), heights = c(2, 3, 2, 3))
par(mar = c(4, 4, 3, 1))

## Run tests
# This first one is a simple mosaic diagram (no stacking)
setup()
ref_Cu_craq()
title("Cu\ncr+aq", font.main = 1)

setup()
ref_FeCu_cr()
title("Fe-Cu\ncr", font.main = 1)

T1result <- "(TODO: looks ok)"

setup()
d_FeCu_old <- test_FeCu_old()
title("Fe-Cu\ncr+aq", font.main = 1)

T2result <- "(TODO: looks ok)"
T3result <- "(TODO: looks different)"

setup()
d_FeCu_new <- test_FeCu_new()
title("Fe-Cu\ncr+aq", font.main = 1, xpd = NA)

if(identical(d_FeCu_old$predominant, d_FeCu_new$predominant)) T4result <- "(OK)" else T4result <- "(FAILED)"

setup()
d_FeCu_stack <- test_FeCu_stack()
title("Fe-Cu\ncr+aq", font.main = 1, xpd = NA)

if(identical(d_FeCu_new$predominant, d_FeCu_stack$predominant)) T5result <- "(OK)" else T5result <- "(FAILED)"

## Add arrows and text
plot.new()
text(0.5, 0.99, "mosaic()\nWITHOUT 'loga_aq'\n(pre-2.0.0)", font = 2, xpd = NA)

arrows(0.35, 0.9, 0.65, 0.9)
text(0.5, 0.93, paste0("Test 1: Cu mineral reactions unaffected\n", T1result))

arrows(0.1, 0.72, 0.25, 0.65)
text(0.05, 0.62, paste0("Test 2: aq reactions unaffected\n", T2result), xpd = NA)

arrows(0.9, 0.72, 0.75, 0.65)
text(0.9, 0.62, paste0("Test 3: cr reactions unaffected\n", T3result), xpd = NA)

abline(h = 0.25, lty = 2, col = 8)
text(0.4, 0.23, "mosaic()\nWITH 'loga_aq'", font = 2)
lines(c(0.5, 0.5), c(0.15, 0.25), lty = 2, col = 8)
text(0.6, 0.23, "stack_mosaic()", font = 2)

arrows(0.25, 0.4, 0.1, 0.3)
text(0.05, 0.35, paste0("Test 4: identical diagram\n", T4result), xpd = NA)

arrows(0.35, 0.1, 0.65, 0.1)
text(0.5, 0.05, paste0("Test 5: identical diagram\n", T5result))

dev.off()
