# CHNOSZ/demo/contour.R
# Gold solubility contours on logfO2-pH diagram
# 20181107 initial version
# 20190415 cleanup for demo
# - After Williams-Jones et al. (2009) doi:10.2113/gselements.5.5.281
# 20250215 adjust conditions to reproduce plot from Ding et al.
# - After Ding et al. (2023) doi:10.1016/j.oregeorev.2022.105231
# 20260610 make plot with IGEM data (Au2S2-2 from Tagirov et al., 2025)
# - After Rubtsova et al. (2026) doi:10.1134/S1075701525600690

library(CHNOSZ)

# Define base plot resolution
res <- 200
# Make smooth (TRUE) or sharp (FALSE) transitions between basis species 
blend <- TRUE

for(author in c("Ding", "Rubtsova")) {

  if(author == "Ding") {
    # Define temperature (degrees C), pressure (bar)
    T <- 250
    P <- 300
    # Define molality of NaCl and total S
    m_NaCl <- 1
    Stot <- 0.01
    # Define ranges of pH and logfO2
    pH <- c(1, 11)
    O2 <- c(-40, -30)
    # Define solubility unit, contour levels and colors
    solubility_unit <- "ppb"
    levels <- c(0.01, 0.1, 1, 10, 100, 1000, 10000)
    title <- "Au solubility (ppb), after Ding et al. (2023)"
  }

  if(author == "Rubtsova") {
    # Define temperature (degrees C), pressure (bar)
    T <- 200
    P <- "Psat"
    # Define molality of NaCl and total S
    m_NaCl <- 0.1
    Stot <- 0.5
    # Define ranges of pH and logfO2
    pH <- c(1, 11)
    O2 <- c(-45, -30)
    # Define solubility unit, contour levels and colors
    solubility_unit <- "ppm"
    levels <- c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 2500, 10000)
    title <- "Au solubility (ppm), after Rubtsova et al. (2026)"
  }

  # Set up system
  basis(c("Au", "Cl-", "H2S", "H2O", "oxygen", "H+"))
  iaq <- retrieve("Au", c("H", "O", "S", "Cl"), state = "aq")
  basis("H2S", log10(Stot))
  # Calculate solution composition for given NaCl molality
  NaCl <- NaCl(m_NaCl = m_NaCl, T = T, P = P)
  basis("Cl-", log10(NaCl$m_Clminus))
  # Define S-bearing basis species to swap through for mosaic
  bases <- c("H2S", "HS-", "S3-", "SO2", "HSO4-", "SO4-2")

  # Load Au(cr) to calculate its solubility
  species("Au")
  # Use 'iaq' for dissolved Au species and 'bases' for mosaic of S-bearing basis species
  s <- solubility(iaq, bases = bases, pH = c(pH, res), O2 = c(O2, res), T = T, P = P, IS = NaCl$IS, blend = blend)
  # Convert to ppb or ppm
  s <- convert(s, solubility_unit)
  # Start plot
  d <- diagram(s, levels = levels)
  # Overlay colors - take logarithms to spread colors out evenly
  col <- hcl.colors(length(levels) - 1, "Blue-Yellow 3", rev = TRUE)
  image(d$vals$pH, d$vals$O2, log10(d$plotvals[[1]]), col = col, breaks = log10(levels), add = TRUE)
  # Re-add contour lines
  diagram(s, levels = levels, cex = 1.1, add = TRUE)

  # Calculate affinity with changing basis species across diagram
  species(iaq)
  m <- mosaic(bases, pH = c(pH, 2*res), O2 = c(O2, 2*res), T = T, P = P, IS = NaCl$IS, blend = blend)
  # Show predominance fields for S-bearing species
  diagram(m$A.bases, col = 8, col.names = 8, lty = 3, italic = TRUE, add = TRUE)
  # Show predominance fields for Au-bearing species
  # NOTE: balance = 1 is used for plotting predominance fields of polynuclear species (i.e. Au2S2-2)
  d <- diagram(m$A.species, col = 4, col.names = 4, lwd = 2, bold = TRUE, add = TRUE, balance = 1)
  # Add dot-dash line for water stability limit
  water.lines(d, lty = 4)

  # Add legend and title
  dP <- describe.property(c("T", "P"), c(T, P))
  lexpr <- c(dP, lNaCl(m_NaCl), lS(Stot))
  legend("topright", legend = lexpr, bty = "n")
  title(main = title, font.main = 1)

}
