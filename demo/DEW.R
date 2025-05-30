# CHNOSZ/demo/DEW.R
# Demo for the Deep Earth Water (DEW) model in CHNOSZ 20170927

library(CHNOSZ)

# Set up subplots
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(3.0, 3.5, 2.5, 1.0), mgp = c(1.7, 0.3, 0), las = 1, tcl = 0.3, xaxs = "i", yaxs = "i")

# Activate DEW model
oldwat <- water("DEW")

###########
#### Plot 1: quartz solubility at high pressure
## after Figure 7D of Sverjensky et al., 2014a [SHA14]
## (Geochim. Cosmochim. Acta, https://doi.org/10.1016/j.gca.2013.12.019)
###########

# Load SiO2 and Si2O4 data taken from DEW spreadsheet
iSi <- add.OBIGT("DEW", c("SiO2", "Si2O4"))
# Print the data references to confirm we got the right ones
thermo.refs(iSi)
# Set temperature ranges for different pressures
# data.frame is used to make P and T the same length
PT0.5 <- data.frame(P = 500, T = seq(200, 550, 10))
PT1.0 <- data.frame(P = 1000, T = seq(200, 700, 10))
PT2.0 <- data.frame(P = 2000, T = seq(200, 700, 10))
PT5.0 <- data.frame(P = 5000, T = seq(200, 850, 10))
PT10.0 <- data.frame(P = 10000, T = seq(200, 825, 10))
PT20.0 <- data.frame(P = 20000, T = seq(200, 800, 10))
PT <- rbind(PT0.5, PT1.0, PT2.0, PT5.0, PT10.0, PT20.0)
# Reaction 1: quartz = SiO2(aq) [equivalent to quartz + 3 H2O = Si(OH)4]
SiO2_logK <- suppressWarnings(subcrt(c("quartz", "SiO2"), c("cr", "aq"), c(-1, 1), P = PT$P, T = PT$T))$out$logK
# Reaction 2: 2 quartz = Si2O4(aq) [equivalent to 2 quartz + 3 H2O = Si2O(OH)6]
Si2O4_logK <- suppressWarnings(subcrt(c("quartz", "Si2O4"), c("cr", "aq"), c(-2, 1), P = PT$P, T = PT$T))$out$logK
# Plot the sum of molalities (== activities) for each pressure
plot(c(200, 1000), c(-2.5, 0.5), type = "n", xlab = axis.label("T"), ylab = "log molality")
for(P in unique(PT$P)) {
  icond <- PT$P == P
  SiO2_logm <- SiO2_logK[icond]
  Si2O4_logm <- Si2O4_logK[icond]
  logm <- log10(10^SiO2_logm + 10^Si2O4_logm)
  lines(PT$T[icond], logm)
  # add text label
  lastT <- tail(PT$T[icond], 1)
  Pkb <- paste(format(P/1000, nsmall = 1), "kb")
  text(lastT+25, tail(logm, 1), Pkb, adj = 0)
}
t1 <- quote("Solubility of"~alpha*"-quartz")
t2 <- "after Sverjensky et al., 2014a"
mtitle(as.expression(c(t1, t2)))
# TODO: lines are a little low at highest P and T ...

###########
#### Plot 2: correlations between non-solvation volume and HKF a1 parameter
## after Figures 12B and 12C of Sverjensky et al., 2014a [SHA14]
###########

# We can't use the DEW water model at 25 degC
reset()
# Load the fitted parameters for species as used by SHA14
# TODO: also use their Ca+2??
# NOTE: don't load NaCl, NH4+, or HS- here because the DEW spreadsheet lists a1 from the correlation
iDEW <- add.OBIGT("DEW", c("CO3-2", "BO2-", "MgCl+", "SiO2", "HCO3-", "Si2O4"))
# Override the check in subcrt() for the DEW water model 20220929
lapply(iDEW, mod.OBIGT, model = "HKF")
# Set up the plot
V0nlab <- expression(Delta * italic(V) * degree[n]~~(cm^3~mol^-1))
a1lab <- expression(italic(a)[1]%*%10~~(cal~mol~bar^-1))
plot(c(-25, 50), c(-4, 12), type = "n", xlab = V0nlab, ylab = a1lab)
# A function to get the HKF parameters, calculate nonsolvation volume, plot points, labels, error bars, and correlation lines
plotfun <- function(species, col, pch, cex, dy, error, xlim, corrfun) {
  # Get HKF a1 parameter
  a1 <- info(info(species), check.it = FALSE)$a1 * 10
  # Set omega to zero to calculate non-solvation volume 20220326
  mod.OBIGT(species, omega = 0)
  Vn <- do.call(rbind, subcrt(species, T = 25)$out)$V
  points(Vn, a1, col = col, pch = pch, cex = cex)
  for(i in 1:length(species)) text(Vn[i], a1[i] + dy, expr.species(species[i]))
  arrows(Vn, a1 - error, Vn, a1 + error, length = 0.03, angle = 90, code = 3, col = col)
  lines(xlim, corrfun(xlim), col = col)
}
# Monovalent ions: Na+, K+, Cl-, Br-
monofun <- function(Vn) 2.0754 + 0.10871 * Vn
# For easier reading, set y-offset to NA so the labels aren't plotted
plotfun(c("Na+", "K+", "Cl-", "Br-"), "red", 19, 1, NA, 0.5, c(-7, 35), monofun)
# Divalent ions: Mg+2, Ca+2, CO3-2, SO4-2
difun <- function(Vn) 3.5321 + 0.23911 * Vn
plotfun(c("Mg+2", "Ca+2", "CO3-2", "SO4-2"), "black", 15, 1, 1.2, 0.7, c(-20, 25), difun)
# Complexes and neutral molecules: BO2-, MgCl+, SiO2, NaCl, HCO3-, Si2O4, NH4+, HS-
compfun <- function(Vn) 1.5204 + 0.19421 * Vn
plotfun(c("MgCl+", "SiO2", "NaCl", "HCO3-", "Si2O4"), "blue1", 18, 1.5, 1, 0.5, c(-20, 50), compfun)
# For easier reading, put some labels below the points
plotfun(c("BO2-", "NH4+", "HS-"), "blue1", 18, 1.5, -1.2, 0.5, c(-20, 50), compfun)
# Include an empty subscript for better spacing between the lines
t1 <- quote("Correlations between non-solvation"[])
t2 <- quote("volume and HKF "*italic(a)[1]*" parameter")
mtitle(as.expression(c(t1, t2)))

# Clean up for next plot 
OBIGT()
water("DEW")

###########
#### Plot 3: logfO2-pH diagram for aqueous inorganic and organic carbon species at high pressure
## after Figure 1b of Sverjensky et al., 2014b [SSH14]
## (Nature Geoscience, https://doi.org/10.1038/NGEO2291)
###########

# Define system
basis("CHNOS+")
inorganic <- c("CO2", "HCO3-", "CO3-2", "CH4")
organic <- c("formic acid", "formate", "acetic acid", "acetate", "propanoic acid", "propanoate")
species(c(inorganic, organic), 0)
fill <- c(rep("aliceblue", length(inorganic)), rep("aquamarine", length(organic)))

# A function to make the diagrams
dfun <- function(T = 600, P = 50000, res = 300) {
  a <- affinity(pH = c(0, 10, res), O2 = c(-24, -12, res), T = T, P = P)
  # Set total C concentration to 0.03 molal
  # (from EQ3NR model for eclogite [Supporting Information of SSH14])
  e <- equilibrate(a, loga.balance = log10(0.03))
  diagram(e, fill = fill)
  dp <- describe.property(c("     T", "     P"), c(T, P), digits = 0)
  legend("bottomleft", legend = dp, bty = "n")
}

water("DEW")
## Not run: make diagram using CHNOSZ default database
## (not recommended for high P)
#dfun()
# Make diagram using CO2, HCO3-, CO3-2, CH4, and acetic acid data from DEW spreadsheet
# (the acetate field disappears if we also use the data for acetate from the spreadsheet 20200629)
#add.OBIGT("DEW", c("CO2", "HCO3-", "CO3-2", "CH4", "acetic acid"))
add.OBIGT("DEW")
dfun()
mtitle(c("Inorganic and organic species", "C[total] = 0.03 molal"))

###########
#### Plot 4: speciation of carbon as a function T, logfO2 and pH (added 20171008)
## after SSH14 Fig. 3
###########

# Conditions:
# T = 600, 700, 800, 900, 1000 degC
# P = 5.0GPa (50000 bar)
# fO2 = QFM - 2
# pH set by jadeite + kyanite + coesite
# output from EQ3NR calculations (SSH14 Supporting Information)
# dissolved carbon: 0.03, 0.2, 1, 4, 20 molal
# true ionic strength: 0.39, 0.57, 0.88, 1.45, 2.49
# pH: 3.80, 3.99, 4.14, 4.25, 4.33
## Activate DEW model
reset()
water("DEW")
# Add species data for DEW
inorganics <- c("CH4", "CO2", "HCO3-", "CO3-2")
organics <- c("formic acid", "formate", "acetic acid", "acetate", "propanoic acid", "propanoate")
add.OBIGT("DEW")
## Set basis species
basis(c("Fe", "SiO2", "CO3-2", "H2O", "oxygen", "H+"))
## Calculate logfO2 in QFM buffer
basis("O2", "QFM")
T <- seq(600, 1000, 100)
buf <- affinity(T = T, P = 50000, return.buffer = TRUE)
## Add species
species(c(inorganics, organics))
## Generate spline functions from IS, pH, and molC values at every 100 degC
IS <- c(0.39, 0.57, 0.88, 1.45, 2.49)
pH <- c(3.80, 3.99, 4.14, 4.25, 4.33)
molC <- c(0.03, 0.2, 1, 4, 20)
## Use extended Debye-Huckel equation with b_gamma set to zero
nonideal("bgamma0")
## Calculate affinities on the T-logfO2-pH-IS transect
a <- affinity(T = T, O2 = buf$O2 - 2, IS = IS, pH = pH, P = 50000)
## Calculate metastable equilibrium activities using the total
## Carbon molality as an approximation of total activity
e <- equilibrate(a, loga.balance = log10(molC))
## Make the diagram; don't plot names of low-abundance species
names <- c(inorganics, organics)
names[c(4, 5, 7, 9)] <- ""
col <- rep("black", length(names))
col[c(1, 3, 6, 8, 10)] <- c("red", "darkgreen", "purple", "orange", "navyblue")
diagram(e, alpha = "balance", names = names, col = col, ylim = c(0, 0.8), ylab = "carbon fraction", spline.method = "natural")

## Add legend and title
ltxt1 <- "P = 50000 bar"
ltxt2 <- substitute(logfO2=="QFM-2", list(logfO2 = axis.label("O2")))
pH <- seq(3.8, 4.3, length.out = length(T))
legend("left", legend = as.expression(c(ltxt1, ltxt2)), bty = "n")
t1 <- "Aqueous carbon speciation"
t2 <- "after Sverjensky et al., 2014b"
mtitle(c(t1, t2))

### Additional checks

# The maximum absolute pairwise difference between x and y
maxdiff <- function(x, y) max(abs(y - x))
## Check that we're within 0.1 of the QFM-2 values used by SSH14
stopifnot(maxdiff((buf$O2-2), c(-17.0, -14.5, -12.5, -10.8, -9.4)) < 0.1)

# Here are the logKs of aqueous species dissociation reactions at 600 degC and 50000 bar,
# values from EQ3NR output in Supporting Information of the paper (p. 103-109):
inorganic.logK <- c(24.4765, -9.0784, -5.3468, 0)
organic.logK <- c(1.7878, 2.5648, 15.3182, 16.9743, 30.4088, 28.9185)
# Calculate equilibrium constants of the reactions in CHNOSZ; use a negative sign to change from formation to dissociation
logK.calc <- -unlist(affinity(T = 600, P = 50000, property = "logK")$values)
logK.calc - c(inorganic.logK, organic.logK)
## Except for acetate, we're within 0.021 of the logK values used by SSH14
stopifnot(maxdiff(logK.calc[-8], c(inorganic.logK, organic.logK)[-8]) < 0.021)

## Check that we get similar activity coefficients
# Activity coefficients for monovalent species from EQ3NR output
loggamma <- c(-0.15, -0.18, -0.22, -0.26, -0.31)
# Activity coefficients calculated in CHNOSZ
sres <- subcrt("propanoate", T = seq(600, 1000, 100), P = 50000, IS = c(0.39, 0.57, 0.88, 1.45, 2.49))
stopifnot(maxdiff(sres$out[[1]]$loggam, loggamma) < 0.023)
# If m_star in nonideal() was zero, we could decrease the tolerance here
#stopifnot(maxdiff(sres$out[[1]]$loggam, loggamma) < 0.004)

# Reset OBIGT database
reset()

par(opar)
