# CHNOSZ/demo/NaCl.R
# NaCl dissocation logK f(T,P)
# after Shock et al., 1992, Fig. 1
#   Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) 
#   Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: 
#   Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 degrees C and 5 kbar. 
#   J. Chem. Soc. Faraday Trans. 88, 803-826. https://doi.org/10.1039/FT9928800803
# 20121111 first version

library(CHNOSZ)

## Uncomment these lines to make the plot with the g-function disabled
#mod.OBIGT("Cl-", z = 0)
#mod.OBIGT("Na+", z = 0)

# Read the data and start a new plot
expt <- read.csv(system.file("extdata/misc/SOJSH.csv", package = "CHNOSZ"), as.is = TRUE)
expt$pch[expt$pch == 0] <- 15
expt$pch[expt$pch == 1] <- 16
expt$pch[expt$pch == 2] <- 17
expt$pch[expt$pch == 5] <- 18
thermo.plot.new(xlim = c(0, 1000), ylim = c(-5.5, 1), xlab = axis.label("T"), ylab = axis.label("logK"))

# Use viridis palette
expt_Ps <- unique(expt$P)
cols <- hcl.colors(length(expt_Ps), palette = "Harmonic")
# Loop over pressures
for(i in 1:length(expt_Ps)) {
  this_expt <- subset(expt, P == expt_Ps[i])
  # Plot experimental logK
  points(this_expt$T, this_expt$logK, cex = 1.5, pch = this_expt$pch, col = cols[i])
}

# We'll be at 9 distinct pressure conditions, including Psat
# Psat is repeated to show "not considered" region
# (T >= 355 degC; Fig. 6 of Shock et al., 1992)
P <- c(list("Psat", "Psat"), as.list(seq(500, 4000, by = 500)))
cols <- c(cols[1], cols)
# The range of temperature for each pressure
T <- list()
T[[1]] <- seq(0, 354, 1)
T[[2]] <- seq(354, 370, 1)
T[[3]] <- seq(265, 465, 1)
T[[4]] <- seq(285, 760, 1)
T[[5]] <- seq(395, 920, 1)
T[[6]] <- T[[7]] <- T[[8]] <- T[[9]] <- T[[10]] <- seq(400, 1000, 1)

# Define reaction for calculating logK
species <- c("NaCl", "Na+", "Cl-")
coeffs <- c(-1, 1, 1)
logK <- numeric()
# Loop over pressures
for(i in 1:length(T)) {
  s <- suppressWarnings(subcrt(species, coeffs, T = T[[i]], P = P[[i]]))
  if(i == 2) lty <- 3 else lty <- 1
  lines(s$out$T, s$out$logK, lwd = 2, lty = lty, col = cols[i])
  # Keep the calculated values for each experimental condition (excluding Psat)
  iexpt <- which(P[[i]] == expt$P)
  Texpt <- expt$T[iexpt]
  if(i > 2) logK <- c(logK, splinefun(s$out$T, s$out$logK)(Texpt))
}

# Add title, labels, and legends
title(describe.reaction(s$reaction, states = 1))
text(150, -0.1, quote(italic(P)[sat]), cex = 1.2)
text(462, -4, "500 bar")
text(620, -4.3, "1000 bar")
text(796, -4.3, "1500 bar")
text(813, -1.3, "4000 bar")
legend("bottomleft", legend = unique(expt$source), pch = unique(expt$pch), pt.cex = c(1.2, 1.2, 1.4, 1.6), bty = "n")
#mtitle(c(describe.reaction(s$reaction), expression(italic(P)[sat]~"and 500-4000 bar")))
l1 <- quote("Revised HKF model with " * italic(g) * " function (Shock et al., 1992)")
l2 <- "Non-recommended region (Shock et al., 1992, Fig. 6)"
legend("topright", as.expression(c(l1, l2)), lty = c(1, 3), bty = "n")

# Test for average divergence (excluding Psat)
expt_test <- expt[!expt$P %in% "Psat", ]
stopifnot(mean(abs(logK - expt_test$logK)) < 0.09)
