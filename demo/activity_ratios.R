## Equilibrium activity diagrams for minerals using activity ratios as variables
## These are made with pH = 0 (activity of H+ = 1), so (activity of the ion) is equal to
## (activity of the ion) / [(activity of H+) ^ (charge of the ion)]

opar <- par(mfrow = c(2, 2))
res <- 200
fill <- "terrain"

data(thermo)
## get data for gibbsite from the SUPCRT92 database, without loading minerals that
## have been superseded by the Berman dataset (the default in CHNOSZ)
add.obigt("SUPCRT92", "gibbsite")
## or, get data for gibbsite and Al-bearing aqueous species from the SUPCRTBL database
#add.obigt("SUPCRTBL")

## K2O-Al2O3-SiO2-H2O, 25 degree C, 1 bar
## Steinmann et al., 1994 (http://ccm.geoscienceworld.org/content/42/2/197)
## Garrels and Christ, p. 361 (http://www.worldcat.org/oclc/517586)
## https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-75.html
basis(c("Al+3", "pseudo-H4SiO4", "K+", "H2O", "H+", "O2"))
species(c("gibbsite", "muscovite", "kaolinite", "pyrophyllite", "K-feldspar"))
a <- affinity(H4SiO4 = c(-6, -2, res), `K+` = c(-3, 6, res))
diagram(a, ylab = ratlab("K+"), fill = fill, yline = 1.7)
title(main = syslab(c("K2O", "Al2O3", "SiO2", "H2O")))
legend("bottomleft", describe.property(c("T", "P"), c(25, 1)), bty = "n")

## H2O-CaO-MgO-SiO2 at 300 degree C and 1000 bar
## Helgeson et al., 1969, p. 136 (http://www.worldcat.org/oclc/902423149)
## Bowers et al., 1984, p. 246 (http://www.worldcat.org/oclc/224591948)
basis(c("H2O", "Ca+2", "Mg+2", "SiO2", "O2", "H+"))
species(c("quartz", "talc", "chrysotile", "forsterite", "tremolite",
          "diopside", "wollastonite", "monticellite", "merwinite"))
# calculate the chemical affinities of formation reactions
a <- affinity("Mg+2" = c(4, 9, res), "Ca+2" = c(5, 14, res), T = 300, P = 1000)
diagram(a, xlab = ratlab("Mg+2"), ylab = ratlab("Ca+2"), fill = fill, yline = 1.7)
title(main = syslab(c("H2O", "CaO", "MgO", "SiO2")))
legend("bottomright", describe.property(c("T", "P"), c(300, 1000)), bty = "n")
# note: Bowers et al. use more complicated variables
# (involving the hydration numbers of H2O and the ion)
# with accordingly different axis ranges

## MgO-CaO-SiO2-H2O at 300 degree C and Psat
## Russell et al., 2010 (https://doi.org/10.1111/j.1472-4669.2010.00249.x)
basis(c("Mg+2", "Ca+2", "SiO2", "H2O", "O2", "H+"))
species(c("brucite", "chrysotile", "talc", "tremolite", "diopside", "akermanite"))
a <- affinity(SiO2 = c(-10, 0, res), `Ca+2` = c(0, 20, res), T = 300)
diagram(a, ylab = ratlab("Ca+2"), fill = fill, yline = 1.7)
title(main = syslab(c("MgO", "CaO", "SiO2", "H2O")))
legend("topright", describe.property(c("T", "P"), c(300, 85.84)), bty = "n")

## CaO-MgO-SiO2-H2O and
## CaO-Al2O3-MgO-SiO2-H2O at 300 degree C and 500 bar
## Bach and Klein, 2009 (https://doi.org/10.1016/j.lithos.2008.10.022)
basis(c("Ca+2", "Al+3", "Mg+2", "SiO2", "H2O", "O2", "H+"))
species(c("clinochlore", "clinozoisite", "prehnite", "grossular"))
a <- affinity(SiO2 = c(-5, 0, res), `Ca+2` = c(6, 11, res), T = 300, P = 500)
diagram(a, ylab = ratlab("Ca+2"), balance = "Al+3", fill = fill, yline = 1.7)
# (Hmmm... where is clinochlore? it doesn't appear on our diagram)
species(delete = TRUE)
species(c("brucite", "chrysotile", "talc", "tremolite", "diopside"))
a <- affinity(SiO2 = c(-5, 0, res), `Ca+2` = c(6, 11, res), T = 300, P = 500)
diagram(a, add = TRUE, col = "blue", col.names = "blue")
title(main = syslab(c("CaO", "Al2O3", "MgO", "SiO2", "H2O")))
legend("topright", describe.property(c("T", "P"), c(300, 500)), bty = "n")

par(opar)
