# CHNOSZ/demo/arsenic.R
# Eh-pH diagram for the system As-O-H-S,
# After Lu and Zhu, 2011 (doi:10.1007/s12665-010-0652-x)
# 20190415 extracted from go-IU.R; use retrieve()

library(CHNOSZ)

# Define temperature (degrees C), pressure (bar), grid resolution
res <- 400
T <- 25
P <- 1
# Change this to FALSE to make sharp transitions between the basis species,
# giving a diagram with straight lines around the AsS(OH)HS- wedge
blend <- TRUE

# Set basis species
basis(c("As", "H2O", "H2S", "H+", "e-"))
basis(c("H2S"), c(-3))
# Find and set formed species
iaq <- retrieve("As", c("S", "O", "H"), "aq")
icr <- retrieve("As", c("S", "O", "H"), "cr")
species(c(iaq, icr))
# Set activities of aqueous species
species(1:length(iaq), -5)

# The possible S-bearing basis species
bases <- c("H2S", "HS-", "S3-", "SO2", "HSO4-", "SO4-2")
# Calculate affinties of formation reactions using the speciated S basis species
m <- mosaic(bases, pH = c(0, 14, res), Eh = c(-0.8, 0.8, res), T = T, P = 1, blend = blend)
# Adjust name of realgar
m$A.species$species$name <- gsub(",alpha", "", m$A.species$species$name)
# Make the plot!
diagram(m$A.species)
# Add legend and title
dprop <- describe.property(c("T", "P"), c(T, P))
legend("bottomleft", legend = dprop, bty = "n")
t1 <- quote("As-O-H-S, "~list(Sigma*S == 10^-3~M, Sigma*As == 10^-5~M))
t2 <- "After Lu and Zhu, 2011 Fig. 2b"
mtitle(as.expression(c(t1, t2)), cex = 0.95)
