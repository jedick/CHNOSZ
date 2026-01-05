# CHNOSZ/demo/sum_S.R
# 20240604 Fe-S-O-C     log a (sum S) - pH      modified from mosaic.R
# 20241130 Fe-S-O + Au  logfO2 - log a (sum S)  added NaCl and solubility contours

library(CHNOSZ)

# Define conditions
res <- 500
loga_S <- c(-6, 0, res)
logf_O2 <- c(-45, -20, res)
pH <- 5
T <- 300
#P <- "Psat"
P <- 1000
#IS <- 0.8

# Calculate solution composition for 0.8 mol/kg NaCl
# Based on activity of Cl- from Fig. 18 of Skirrow and Walsh (2002)
m_NaCl = 0.8
NaCl <- NaCl(m_NaCl = m_NaCl, T = T, P = P, pH = pH)
IS <- NaCl$IS

# Setup chemical system
# Use oxygen instead of O2 to get the gas
basis(c("Fe", "SO4-2", "oxygen", "H+", "Cl-", "H2O"))
basis("pH", pH)
basis("Cl-", log10(NaCl$m_Clminus))
# Add minerals as formed species
species(c("pyrrhotite", "pyrite", "hematite", "magnetite"))

# List basis species to swap through
bases <- list( c("H2S", "HS-", "S3-", "SO2", "HSO4-", "SO4-2") )
# Calculate mosaic for Fe minerals
m <- mosaic(bases, "SO4-2" = loga_S, O2 = logf_O2, T = T, P = P, IS = IS)

# Loop over metals for solubility contours
for(metal in c("Fe", "Au")) {

  # Plot diagram for Fe minerals and show sulfur species
  d <- diagram(m$A.species, bold = TRUE, lwd = 2)
  diagram(m$A.bases[[1]], add = TRUE, col = 8, col.names = 8, lty = 4, dx = -2.5, dy = c(-4, 0, 0, 6), italic = TRUE)
  # Add water stability limit
  water.lines(d)

  if(metal == "Au") {
    # Put Au as the first basis species for solubility() calculation
    basis(c("Au", "SO4-2", "oxygen", "H+", "Cl-", "H2O"))
    basis("pH", pH)
    # This is the species whose solubility we want to calculate
    species("Au")
    # These are the aqueous species that form by dissolving Au
    iaq <- retrieve("Au", state = "aq", ligands = c("O", "H", "S", "Cl"))
    # Make orange contours for Au
    cols <- c(6, 7)
  } else {
    # For Fe, the basis species and formed species are already set up properly for solubility()
    # These are the aqueous species that form by dissolving Fe minerals
    iaq <- retrieve("Fe", state = "aq", ligands = c("O", "H", "S", "Cl"))
    # Make red contours for Fe
    cols <- c(2, 2)
  }
  # Calculate solubility and convert log molality to ppb
  s <- solubility(iaq, bases = bases, "SO4-2" = loga_S, O2 = logf_O2, T = T, P = P, IS = IS, in.terms.of = metal)
  sp <- convert(s, "ppb")
  # Plot twice to get deeper colors
  for(col in cols) {
    diagram(sp, contour.method = "flattest", levels = 10^(-3:5), add = TRUE, col = col, lty = 3, lwd = 1.8, cex = 1.2)
  }

  # Add legend and title
  T_P <- describe.property(c("T", "P"), c(T, P))
  P_Cl <- c(bquote(pH == .(pH)), bquote(NaCl == .(m_NaCl)~mol~kg^-1))
  legend <- lex(T_P, P_Cl)
  legend("topright", legend, bty = "n")
  main <- paste("Fe minerals with summed aq species: S (x-axis) and ppb", metal, "(contours)")
  title(main, font.main = 1)

}
