# CHNOSZ/demo/dehydration.R
# Plot logK vs temperature for some dehydration reactions
# 20151126 first version
# 20250527 add state labels for amino acid and protein backbone

library(CHNOSZ)

# Define T values
T <- seq(0, 175)

# Define reactions
reactants <- c("[AABB]", "[AABB]", "malate-2", "goethite", "gypsum", "epsomite", "ethanol")
products <- c("[UPBB]", "[PBB]", "fumarate-2", "hematite", "anhydrite", "hexahydrite", "ethylene")
rstate <- c("aq", "cr", "aq", "cr", "cr", "cr", "aq")
pstate <- c("aq", "cr", "aq", "cr", "cr", "cr", "gas")
rcoeff <- c(-1, -1, -1, -2, -0.5, -1, -1)
pcoeff <- c(1, 1, 1, 1, 0.5, 1, 1)
# Assign position (index on x-axis) and rotation of the names
ilab <- c(140, 120, 60, 60, 20, 130, 120)
srt <- c(15, 25, 20, 12, 8, 15, 30)

# Setup plot
plot(range(T), c(-3, 1), type = "n", xlab = axis.label("T"), ylab = axis.label("logK"))
title(main = "Dehydration reactions")

for(i in 1:length(reactants)) {

  # Calculate standard thermodynamic properties of reaction
  # Suppress warnings about T limits of Cp equations for
  #   goethite, gypsum, epsomite, and hexahydrite 20250527
  s <- suppressWarnings(subcrt(c(reactants[i], products[i], "H2O"),
    c(rstate[i], pstate[i], "liq"), c(rcoeff[i], pcoeff[i], 1), T = T))

  # Plot line for logK
  lines(T, s$out$logK)

  # Add labels for reactant and product
  for(j in 1:2) {
    # Remove charge from names
    name <- gsub("-.*", "", s$reaction$name[j])
    # Add state for amino acid and protein backbones
    if(grepl("BB", name)) name <- paste(name, rstate[i])
    # Put reactant and product labels on opposite sides of line
    if(j == 1) dy <- 0.08 else dy <- -0.03
    if(j == 1) dx <- 0 else dx <- 5
    text(T[ilab[i]] + dx, s$out$logK[ilab[i]] + dy, name, adj = 1, srt = srt[i])
  }

}

# Plot dotted line at logK = 0
abline(h = 0, lty = 3)
