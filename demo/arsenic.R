# CHNOSZ/demo/arsenic.R

## Eh-pH diagram for the system As-O-H-S,
## After Lu and Zhu, 2011 (doi:10.1007/s12665-010-0652-x)
## 20190415 extracted from go-IU.R; use retrieve()

library(CHNOSZ)

# Define temperature (degrees C), pressure (bar), grid resolution
res <- 400
T <- 25
P <- 1
# Define As and total S activity
loga_As <- -5
logStot <- -3
# Change this to FALSE to make sharp transitions between the basis species,
# giving a diagram with straight lines around the AsS(OH)HS- wedge
blend <- TRUE

# Set basis species
basis(c("As", "H2O", "H2S", "H+", "e-"))
basis("H2S", logStot)
# Find and set formed species
iaq <- retrieve("As", c("S", "O", "H"), "aq")
icr <- retrieve("As", c("S", "O", "H"), "cr")
species(c(iaq, icr))
# Set activities of aqueous species
species(1:length(iaq), loga_As)

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
l_S <- bquote(log~italic(m)*S[tot] == .(logStot))
l_As <- bquote(log~italic(m)*As == .(loga_As))
legend("bottomleft", legend = c(dprop, l_S, l_As), bty = "n")
title("As-S-O-H, after Lu and Zhu (2011)", font.main = 1)

## logfO2-pH diagram for Fe-As-O-S minerals and aqueous Fe
## After Ding et al. (2023) (doi:10.1016/j.oregeorev.2022.105231)
## 20260612 first version

# Define conditions
T <- 250
P <- 300
Stot <- 0.01
Astot <- 0.005
m_NaCl <- 1
# Log activities of aqueous species
logact <- -5

# Set basis species
# NOTE: This must include the first species in S_species and As_species 
basis(c("Fe", "H2O", "H2S", "AsH3", "Cl-", "H+", "oxygen"))
basis("H2S", log10(Stot))
basis("AsH3", log10(Astot))
basis("Cl-", log10(m_NaCl))
# Set formed species
#iaq <- retrieve("Fe", c("As", "S", "O", "H"), "aq")
#icr <- retrieve("Fe", c("As", "S", "O", "H"), "cr")
#species(c(iaq, icr))
iaq <- retrieve("Fe", c("S", "O", "H", "Cl"), "aq")
species(iaq, logact)
cr <- c("pyrite", "pyrrhotite", "magnetite", "hematite", "arsenopyrite")
species(cr, add = TRUE)

# Speciate both S and As
S_species <- c("H2S", "HS-", "S3-", "SO2", "HSO4-", "SO4-2")
# Get only As species that have non-NA G at T and P
As_species <- names(retrieve("As", state = "aq", T = T, P = P))
bases <- list(S_species, As_species)
m <- mosaic(bases, pH = c(1, 11, res), O2 = c(-40, -30, res), T = T, P = P, IS = m_NaCl, blend = blend)
# Use abbreviations for minerals
names <- species()$name
abbrv <- info(info(cr))$abbrv
names[species()$state == "cr"] <- abbrv
# Make the base plot
diagram(m$A.species, names = names)
# Overlay S species
diagram(m$A.bases[[1]], add = TRUE, lty = 2, italic = TRUE, col = 8, col.names = 8)
# Overlay As species
dx <- rep(0, length(As_species))
dx[As_species == "As(OH)3"] <- 2.8
diagram(m$A.bases[[2]], add = TRUE, lty = 2, italic = TRUE, col = 2, col.names = 2, dx = dx)

# Add a legend and title
l_TP <- c(lT(T), lP(P))
l_S <- bquote(italic(m)*S[tot] == .(Stot))
l_As <- bquote(italic(m)*As[tot] == .(Astot))
l_NaCl <- bquote(italic(m)*NaCl == .(m_NaCl))
legend("bottomleft", legend = c(l_TP, l_S, l_As, l_NaCl), bg = "white")
title("Fe-As-S-O-H, after Ding et al. (2023)", font.main = 1)
