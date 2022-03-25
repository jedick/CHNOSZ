# CHNOSZ/demo/AD.R
# Calculate Henry's constant using the Akinfiev-Diamond model 20190220
# Add volume and heat capacity 20220207
library(CHNOSZ)

# Start with default settings
reset()

########################
### HENRY'S CONSTANT ###
########################

# Function to plot natural logarithm of Henry's constant
lines.KH <- function(name = "CO2", T = 1:373, P = "Psat", HKF = FALSE, altH2S = FALSE) {
  # use AD or HKF model?
  if(HKF) add.OBIGT("inorganic_aq")
  if(!HKF) add.OBIGT("AD")
  # use alternative parameters for H2S? (AD03 Table 1)
  if(altH2S) mod.OBIGT("H2S", state = "aq", a = -11.2303, b = 12.6104, c = -0.2102)
  # get properties of aq - gas reaction
  sres <- suppressWarnings(subcrt(c(name, name), c("aq", "gas"), c(-1, 1), T = T, P = P))
  # calculate natural logarithm of Henry's constant in mole-fraction units
  ln_KH <- log(1000/18.0153) + log(10) * sres$out$logK
  # plot with units of reciprocal temperature (1000/K)
  TK <- convert(T, "K")
  if(HKF) lty <- 2 else lty <- 1
  if(HKF) col <- 2 else col <- 1
  if(altH2S) lty <- 3
  lines(1000/TK, ln_KH, lty = lty, col = col)
}

# Set up plot
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
par(mar = c(3.5, 3.5, 2.5, 1))
par(mgp = c(2.4, 1, 0))

ylab <- quote(ln~italic(K[H]))
xlab <- quote(1000 / list(italic(T), K))

# CO2 (Fig. 1a of AD03)
plot(0, 0, xlim = c(1, 4), ylim = c(4, 10), xlab = xlab, ylab = ylab)
lines.KH("CO2", 1:373, "Psat")
lines.KH("CO2", seq(100, 650, 10), 500)
lines.KH("CO2", 1:373, "Psat", HKF = TRUE)
lines.KH("CO2", seq(100, 650, 10), 500, HKF = TRUE)
dat <- read.csv(system.file("extdata/cpetc/AD03_Fig1a.csv", package = "CHNOSZ"))
points(dat$x, dat$y, pch = dat$pch)
text(3.5, 7.8, quote(italic(P)[sat]))
text(3.05, 9.2, "500 bar")
legend("bottom", c("Data (AD03, Fig. 1a)", "AD model", "Revised HKF model"),
       lty = c(0, 1, 2), pch = c(1, NA, NA), col = c(1, 1, 2), bty = "n")
title(main = syslab(c("CO2", "H2O"), dash = " - "))

# H2 (Fig. 1b of AD03)
plot(0, 0, xlim = c(1, 4), ylim = c(8, 12), xlab = xlab, ylab = ylab)
lines.KH("H2", 1:373, "Psat")
lines.KH("H2", seq(100, 650, 10), 1000)
lines.KH("H2", 1:373, "Psat", HKF = TRUE)
lines.KH("H2", seq(100, 650, 10), 1000, HKF = TRUE)
text(3.4, 11.4, quote(italic(P)[sat]))
text(1.5, 11, "1000 bar")
dat <- read.csv(system.file("extdata/cpetc/AD03_Fig1b.csv", package = "CHNOSZ"))
points(dat$x, dat$y, pch = dat$pch)
legend("bottomright", c("Data (AD03, Fig. 1b)", "AD model", "Revised HKF model"),
       lty = c(0, 1, 2), pch = c(1, NA, NA), col = c(1, 1, 2), bty = "n")
title(main = syslab(c("H2", "H2O"), dash = " - "))

# H2S (Fig. 1c of AD03)
plot(0, 0, xlim = c(1, 4), ylim = c(4, 9), xlab = xlab, ylab = ylab)
lines.KH("H2S", 1:373, "Psat")
lines.KH("H2S", seq(100, 650, 10), 1000)
lines.KH("H2S", 1:373, "Psat", altH2S = TRUE)
lines.KH("H2S", seq(100, 650, 10), 1000, altH2S = TRUE)
lines.KH("H2S", 1:373, "Psat", HKF = TRUE)
lines.KH("H2S", seq(100, 650, 10), 1000, HKF = TRUE)
dat <- read.csv(system.file("extdata/cpetc/AD03_Fig1c.csv", package = "CHNOSZ"))
points(dat$x, dat$y, pch = dat$pch)
text(3.4, 6.9, quote(italic(P)[sat]))
text(3.1, 8.6, "1000 bar")
legend("bottom", c("Data (AD03, Fig. 1c)", "AD model", "AD model (alt. H2S)", "Revised HKF model"),
       lty = c(0, 1, 3, 2), pch = c(1, NA, NA, NA), col = c(1, 1, 1, 2), bty = "n")
title(main = syslab(c("H2S", "H2O"), dash = " - "))

# CH4 (Fig. 1d of AD03)
plot(0, 0, xlim = c(1.5, 4), ylim = c(8, 12), xlab = xlab, ylab = ylab)
lines.KH("CH4", 1:350, "Psat")
lines.KH("CH4", 1:350, "Psat", HKF = TRUE)
dat <- read.csv(system.file("extdata/cpetc/AD03_Fig1d.csv", package = "CHNOSZ"))
points(dat$x, dat$y, pch = dat$pch)
text(3.4, 11, quote(italic(P)[sat]))
legend("bottomright", c("Data (AD03, Fig. 1d)", "AD model", "Revised HKF model"),
       lty = c(0, 1, 2), pch = c(1, NA, NA), col = c(1, 1, 2), bty = "n")
title(main = syslab(c("CH4", "H2O"), dash = " - "))

##############
### VOLUME ###
##############

# Function to plot volume
lines.V <- function(species = "CO2", T = seq(300, 440, 1), P = 280, HKF = FALSE) {
  # Use AD or HKF model?
  if(HKF) add.OBIGT("inorganic_aq")
  if(!HKF) add.OBIGT("AD")
  # Calculate V; use exceed.rhomin to allow HKF calculations in low-density region
  V <- subcrt(species, T = T, P = P, exceed.rhomin = HKF)$out[[1]]$V
  if(HKF) lty <- 2 else lty <- 1
  if(HKF) col <- 2 else col <- "gray20"
  lines(T, V, lty = lty, col = col, lwd = 1.5)
}

# Read file with V data for four species
file <- system.file("extdata/cpetc/HWM96_V.csv", package = "CHNOSZ")
dat <- read.csv(file)
# Use data near 280 bar
dat <- dat[abs(dat$P - 28) < 0.1, ]
# Setup plot
par(mfrow = c(2, 2))
par(mar = c(4, 4, 3, 1))
par(cex = 1.2)
par(mgp = c(2.5, 1, 0))
# Loop over species
for(species in c("CH4", "CO2", "H2S", "NH3")) {
  thisdat <- dat[dat$species == species, ]
  T <- convert(thisdat$T, "C")
  if(species %in% c("CH4", "CO2")) ylim <- c(0, 2200)
  if(species %in% c("H2S", "NH3")) ylim <- c(0, 1200)
  plot(T, thisdat$V, xlim = c(300, 440), ylim = ylim, xlab = axis.label("T"), ylab = axis.label("V"))
  lines.V(species)
  lines.V(species, HKF = TRUE)
  legend("topleft", legend = expr.species(species, use.state = TRUE), bty = "n")
}
par(mfrow = c(1, 1))
plot.window(c(0, 1), c(0, 1))
box(col = "transparent")
legend <- c("Hn\u011bdkovsk\u00fd et al. (1996)", "AD model", "Revised HKF model")
legend(0.26, 0.54, legend, pch = c(1, NA, NA), lty = c(NA, 1, 2), col = c(1, "gray20", 2), lwd = 1.5, inset = c(-10, -10), bty = "n", text.font = 2)
title(main = "Aqueous nonelectrolytes (280 bar)")

#####################
### HEAT CAPACITY ###
#####################

# Function to plot heat capacity
lines.Cp <- function(species = "CO2", T = seq(300, 440, 1), P = 280, HKF = FALSE) {
  # Use AD or HKF model?
  if(HKF) add.OBIGT("inorganic_aq")
  if(!HKF) add.OBIGT("AD")
  # Calculate Cp; use exceed.rhomin to allow HKF calculations in low-density region
  Cp <- subcrt(species, T = T, P = P, exceed.rhomin = HKF)$out[[1]]$Cp
  if(HKF) lty <- 2 else lty <- 1
  if(HKF) col <- 2 else col <- "gray20"
  lines(T, Cp, lty = lty, col = col, lwd = 1.5)
}

# Read file with Cp data for four species
file <- system.file("extdata/cpetc/HW97_Cp.csv", package = "CHNOSZ")
dat <- read.csv(file)
# Setup plot
par(mfrow = c(2, 2))
par(mar = c(4, 4, 3, 1))
par(cex = 1.2)
par(mgp = c(2.5, 1, 0))
# Loop over species
for(species in c("CH4", "CO2", "H2S", "NH3")) {
  thisdat <- dat[dat$species == species, ]
  T <- convert(thisdat$T, "C")
  if(species %in% c("CH4", "CO2")) ylim <- c(-40000, 40000)
  if(species %in% c("H2S", "NH3")) ylim <- c(-20000, 20000)
  plot(T, thisdat$Cp, xlim = c(300, 440), ylim = ylim, xlab = axis.label("T"), ylab = axis.label("Cp"))
  lines.Cp(species)
  lines.Cp(species, HKF = TRUE)
  legend("topleft", legend = expr.species(species, use.state = TRUE), bty = "n")
}
par(mfrow = c(1, 1))
plot.window(c(0, 1), c(0, 1))
box(col = "transparent")
legend <- c("Hn\u011bdkovsk\u00fd and Wood (1997)", "AD model", "Revised HKF model")
legend(0.26, 0.54, legend, pch = c(1, NA, NA), lty = c(NA, 1, 2), col = c(1, "gray20", 2), lwd = 1.5, inset = c(-10, -10), bty = "n", text.font = 2)
title(main = "Aqueous nonelectrolytes (280 bar)")

par(opar)
