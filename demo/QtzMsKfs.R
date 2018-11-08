# CHNOSZ/demo/QtzMsKfs.R
# T - log(K+/H+) diagram, after Sverjensky et al., 1991
# (doi:10.1016/0016-7037(91)90157-Z)
# 20171009 diagram added to berman.Rd
# 20181108 moved to demo/QtzMsKfs.R; add molality calculations

# this demo compares diagrams made using the Berman and Helgeson datasets,
# and shows the use of nonideal calculations to set molalities in the basis species

## set up the system: basis species
basis(c("K+", "Al+3", "quartz", "H2O", "O2", "H+"))
# use pH = 0 so that aK+ = aK+/aH+
basis("pH", 0)
# load the species
species(c("K-feldspar", "muscovite", "kaolinite",
          "pyrophyllite", "andalusite"), "cr")
## the "b_gamma" equation gets closer to the published diagram than "B-dot"
thermo$opt$nonideal <<- "bgamma"

## start with the data from Helgeson et al., 1978
add.obigt("SUPCRT92")
# calculate affinities in aK+ - temperature space
# exceed.Tr: we go above stated temperature limit of pyrophyllite
# (this is above its stability field on the diagram, so pyrophyllite doesn't appear in this region,
# but its properties are needed needed to calculate relative stabilities of all minerals)
res <- 400
a <- affinity(`K+` = c(0, 5, res), T = c(200, 650, res), P = 1000, exceed.Ttr = TRUE)
# make base plot with colors and no lines
diagram(a, xlab = ratlab("K+", use.molality = TRUE), lty = 0, fill = "terrain")
# add the lines, extending into the low-density region (exceed.rhomin = TRUE)
a <- affinity(`K+` = c(0, 5, res), T = c(200, 650, res), P = 1000, exceed.Ttr = TRUE, exceed.rhomin = TRUE)
diagram(a, add = TRUE, names = NULL, col = "red", lty = 2, lwd = 1.5)
# calculate and plot the lines for 1 molal chloride
a <- affinity(`K+` = c(0, 5, res), T = c(200, 650, res), P = 1000, exceed.Ttr = TRUE, exceed.rhomin = TRUE, IS = 1)
diagram(a, add = TRUE, names = NULL, col = "red", lwd = 1.5)
# the list of references:
ref1 <- thermo.refs(species()$ispecies)$key

## now use the (default) data from Berman, 1988
# this resets the thermodynamic database
# without affecting the basis and species settings
data(OBIGT)
# we can check that we have Berman's quartz
# and not coesite or some other phase of SiO2
iSiO2 <- rownames(basis()) == "SiO2"
stopifnot(info(basis()$ispecies[iSiO2])$name == "quartz")
# Berman's dataset doesn't have the upper temperature limits, so we don't need exceed.Ttr here
a <- affinity(`K+` = c(0, 5, res), T = c(200, 650, res), P = 1000, exceed.rhomin = TRUE)
diagram(a, add = TRUE, names = NULL, col = "blue", lty = 2, lwd = 1.5)
a <- affinity(`K+` = c(0, 5, res), T = c(200, 650, res), P = 1000, exceed.rhomin = TRUE, IS = 1)
diagram(a, add = TRUE, names = NULL, col = "blue", lwd = 1.5)
# the list of references:
ref2 <- thermo.refs(species()$ispecies)$key
ref2 <- paste(ref2, collapse = ", ")

# add experimental points for 1000 bar (Table 1 of Sverjensky et al., 1991)
expt.T <- c(300, 400, 500, 550,  # KFs-Ms-Qtz
            400, 450, 500, 550,  # Ms-And-Qtz
            300, 350,            # Ms-P-Qtz
            300, 600)            # Kaol-Ms-Qtz, KFs-And-Qtz
expt.KH <- c(3.50, 2.75, 1.95, 1.40, 1.60, 1.57, 1.47, 1.38, 1.94, 1.80, 1.90, 0.63)
points(expt.KH, expt.T, pch = 19, cex = 1.2)
# add legend and title
legend("top", "low-density region", text.font = 3, bty = "n")
legend("top", describe.property(c(NA, NA, "P", "IS"), c(NA, NA, 1000, 1)), bty = "n")
legend("left", c(ref1, ref2, "ion molality", "ion activity", "experiments"),
       lty = c(1, 1, 1, 2, 0), lwd = 1.5, col = c(2, 4, 1, 1, 1), pch = c(NA, NA, NA, NA, 19), bty = "n")
title(main = syslab(c("K2O", "Al2O3", "SiO2", "H2O", "HCl")), line = 1.8)
title(main = "Helgeson and Berman minerals, after Sverjensky et al., 1991", line = 0.3, font.main = 1)
