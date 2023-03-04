# CHNOSZ/demo/saturation.R
## Make equilibrium activity diagrams including saturation limits
## and using activity ratios as variables
# first version (activity_ratios.R) 20170217
# keep one diagram and add saturation lines 20190127
library(CHNOSZ)

# The ratios are calculated with pH = 0 (activity of H+ = 1), so (activity of the ion) is equal to
# (activity of the ion) / [(activity of H+) ^ (charge of the ion)]

# NOTE: Bowers et al. use more complicated variables
# (involving the hydration numbers of H2O and the ion)
# with subsequently different axis ranges

## H2O-CaO-MgO-SiO2 at 300 degree C and 1000 bar
## Helgeson et al., 1969, p. 136 (http://www.worldcat.org/oclc/902423149)
## Bowers et al., 1984, p. 246 (http://www.worldcat.org/oclc/224591948)
par(cex = 1.4)
basis(c("SiO2", "Ca+2", "Mg+2", "carbon dioxide", "H2O", "O2", "H+"))
species(c("quartz", "talc", "chrysotile", "forsterite", "tremolite",
          "diopside", "wollastonite", "monticellite", "merwinite"))
# Calculate the chemical affinities of formation reactions
a <- affinity("Mg+2" = c(4, 10, 500), "Ca+2" = c(5, 15, 500), T = 300, P = 1000)
diagram(a, xlab = ratlab("Mg+2"), ylab = ratlab("Ca+2"), fill = "terrain", yline = 1.7)

# Add saturation limits for specified CO2 fugacity
basis("CO2", -1)
species(c("calcite", "dolomite", "magnesite", "brucite"))
# Use argument recall feature to rerun affinity over the same range of conditions
a <- affinity(a)
diagram(a, type = "saturation", add = TRUE, contour.method = c("edge", "edge", "flattest", "flattest"), lty = 2, cex = 1.4, col = "blue3")

# Add title and legend
title(main = syslab(c("H2O", "CO2", "CaO", "MgO", "SiO2")))
dprop <- describe.property(c("T", "P"), c(300, 1000))
dbasis <- describe.basis(4)
legend("bottomright", c(dprop, dbasis), bty = "n", cex = 0.9)
