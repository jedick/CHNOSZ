# CHNOSZ/demo/contour.R
# Gold solubility contours on logfO2-pH diagram
# After Williams-Jones et al. (2009) doi:10.2113/gselements.5.5.281
#   and Ding et al. (2023) doi:10.1016/j.oregeorev.2022.105231
#   with optional data from Tagirov et al. (2024) doi:10.1016/j.gca.2024.08.022
# 20181107 initial version
# 20190415 cleanup for demo
# 20250215 add IGEM data and adjust conditions to reproduce plot from Ding et al.

library(CHNOSZ)

# Define plot resolution
res <- 600
# Define temperature (degrees C), pressure (bar)
T <- 250
P <- 300
# Define molality of NaCl and total S
m_NaCl <- 1
sum_S <- 0.01
# Define ranges of pH and logfO2
pH <- c(1, 10)
O2 <- c(-40, -30)
# Make smooth (TRUE) or sharp (FALSE) transitions between basis species 
blend <- TRUE

for(i in 1:2) {

  # Use the default database for the first diagram
  if(i == 1) reset()
  # Use IGEM data for the second diagram
  if(i == 2) add.OBIGT("IGEM")

  # Set up system
  basis(c("Au", "Cl-", "H2S", "H2O", "oxygen", "H+"))
  iaq <- retrieve("Au", c("H", "O", "S", "Cl"), state = "aq")
  species(iaq)
  basis("H2S", log10(sum_S))
  # Calculate solution composition for given NaCl molality
  NaCl <- NaCl(m_NaCl = m_NaCl, T = T, P = P)
  basis("Cl-", log10(NaCl$m_Clminus))
  # Calculate affinity with changing basis species across diagram
  bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
  m <- mosaic(bases, pH = c(pH, res), O2 = c(O2, res), T = T, P = P, IS = NaCl$IS, blend = blend)
  # Show predominance fields for S-species
  diagram(m$A.bases, col = 8, col.names = 8, lty = 3, italic = TRUE)
  # Show predominance fields for Au-species
  d <- diagram(m$A.species, col = 4, col.names = 4, lty = 2, bold = TRUE, add = TRUE)
  # Add dot-dash line for water stability limit
  water.lines(d, lty = 4)

  # Calculate and plot solubility of Au (use named 'bases' argument to trigger mosaic calculation)
  species("Au")
  s <- solubility(iaq, bases = bases, pH = c(pH, res), O2 = c(O2, res), T = T, P = P, IS = NaCl$IS, blend = blend)
  # Convert to ppb
  s <- convert(s, "ppb")
  diagram(s, levels = c(1, 10, 100, 1000), col = 1, add = TRUE)
  # Add legend and title
  dP <- describe.property(c("T", "P"), c(T, P))
  lexpr <- lex(dP, lNaCl(m_NaCl), lS(sum_S))
  legend("topright", lexpr, bty = "n")
  if(i == 1) title(main = ("Gold solubility (ppb), after Ding et al., 2023 (default data)"), font.main = 1)
  if(i == 2) title(main = ("Gold solubility (ppb), after Ding et al., 2023 (optional data)"), font.main = 1)

}
