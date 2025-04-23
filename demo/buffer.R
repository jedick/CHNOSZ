# CHNOSZ/demo/buffer.R
# Calculate buffered activities of basis species using two methods
# Ater Figure 6 of Schulte and Shock, 1995 (doi:10.1007/BF01581580)
library(CHNOSZ)

# Use Helgeson et al. (1978) minerals for closer reproduction 20201110
add.OBIGT("SUPCRT92")

b.species <- c("Fe", "CO2", "H2O", "N2", "H2", "H2S", "SiO2")
b.state <- c("cr", "gas", "liq", "gas", "gas", "aq", "aq")
b.logact <- c(0, 1, 0, 0, 0, 0, 0)
basis(b.species, b.state, b.logact)
xlim <- c(0, 350)
thermo.plot.new(xlim = xlim, ylim = c(-4, 4), xlab = axis.label("T"), ylab = axis.label("H2"))
# Method 1: in buffer(), assign name of buffer to basis species
bufferline <- function(buffer, ixlab) {
  basis("H2", buffer)
  a <- affinity(T = xlim, P = 300, return.buffer = TRUE, exceed.Ttr = TRUE)
  lines(a$vals[[1]], a$H2, col = 3, lwd = 2)
  text(a$vals[[1]][ixlab], a$H2[ixlab] + 0.2, buffer, font = 2)
}
bufferline("FeFeO", 40)
bufferline("QFM", 70)
bufferline("PPM", 204)
bufferline("HM", 102)
# Method 2: in diagram(), use the `type` argument
basis("H2", 0)
for(logact in c(-6, -10, -15)) {
  species(c("formaldehyde", "HCN"), logact)
  a <- affinity(T = xlim, P = 300)
  d <- diagram(a, type = "H2", lty = c(3, 2), add = TRUE)
  text(a$vals[[1]][13], mean(sapply(d$plotvals, c)[13, ]), logact)
}
# Add legends and title
legend("topright", legend = c("minerals", "formaldehyde", "HCN"),
  lty = c(1, 3, 2), lwd = c(2, 1, 1), col = c(3, 1, 1), bg = "white", cex = 0.9)
legend("bottomright", legend = c(describe.property("P", 300),
  describe.basis(c(2,4))), bg = "white", cex = 0.9)
title(main = paste("Mineral buffers and activities of aqueous species",
                 "(Schulte and Shock, 1995)", sep = "\n"), cex.main = 0.9)

# Reset OBIGT database
reset()
