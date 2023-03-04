# CHNOSZ/demo/potassium.R
# Moved from berman.Rd 20200727
### Compare mineral stabilities predicted with the Berman and Helgeson datasets
### on a T - log(K+/H+) diagram, after Sverjensky et al., 1991
### (doi:10.1016/0016-7037(91)90157-Z)
library(CHNOSZ)

## Set up the system: basis species
basis(c("K+", "Al+3", "quartz", "H2O", "O2", "H+"))
# Use pH = 0 so that aK+ = aK+/aH+
basis("pH", 0)
# Load the species
species(c("K-feldspar", "muscovite", "kaolinite",
          "pyrophyllite", "andalusite"), "cr")

## Start with the data from Helgeson et al., 1978
add.OBIGT("SUPCRT92")
# Calculate affinities in aK+ - temperature space
# exceed.Tr: enable calculations above stated temperature limit of pyrophyllite
res <- 400
a <- suppressWarnings(affinity(`K+` = c(0, 5, res), T = c(200, 650, res), P = 1000, exceed.Ttr = TRUE))
# Make base plot with colors and no lines
diagram(a, xlab = ratlab("K+", molality = TRUE), lty = 0, fill = "terrain")
# Add the lines, extending into the low-density region (exceed.rhomin = TRUE)
a <- affinity(`K+` = c(0, 5, res), T = c(200, 650, res), P = 1000, 
              exceed.Ttr = TRUE, exceed.rhomin = TRUE)
diagram(a, add = TRUE, names = FALSE, col = 2, lwd = 1.5, lty = 2)
# The list of references:
ref1 <- thermo.refs(species()$ispecies)$key

## Now use the (default) data from Berman, 1988
# This resets the thermodynamic database
# without affecting the basis and species settings
OBIGT()
# Check that we have Berman's quartz
# and not coesite or some other phase of SiO2
iSiO2 <- rownames(basis()) == "SiO2"
stopifnot(info(basis()$ispecies[iSiO2])$name == "quartz")
# Berman's dataset doesn't have the upper temperature limits,
# so we don't need exceed.Ttr here
a <- affinity(`K+` = c(0, 5, res), T = c(200, 650, res), P = 1000, exceed.rhomin = TRUE)
diagram(a, add = TRUE, names = FALSE, col = 4, lwd = 1.5)
# The list of references:
ref2 <- thermo.refs(species()$ispecies)$key
ref2 <- paste(ref2, collapse = ", ")
# Add legend and title
legend("top", "low-density region", text.font = 3, bty = "n")
legend("topleft", describe.property(c("P", "IS"), c(1000, 1)), bty = "n", inset = c(0, 0.1))
legend("topleft", c("Helgeson et al., 1978", "Berman, 1988 and\nSverjensky et al., 1991"),
       lty = c(2, 1), lwd = 1.5, col = c(2, 4), bty = "n", inset = c(0.25, 0.08))
title(main = syslab(c("K2O", "Al2O3", "SiO2", "H2O", "HCl")), line = 1.8)
title(main = "After Sverjensky et al., 1991",
      line = 0.3, font.main = 1)

# Cleanup for next example
reset()
