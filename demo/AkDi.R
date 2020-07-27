# CHNOSZ/demo/AkDi.R
# calculations using the Akinfiev-Diamond model 20190220
# after Fig. 1 of Akinfiev and Diamond, 2003
library(CHNOSZ)

# function to plot natural logarithm of Henry's constant
lines.KH <- function(name = "CO2", T = 1:373, P = "Psat", HKF = FALSE, altH2S = FALSE) {
  # use AkDi or HKF model?
  if(!HKF) add.OBIGT("AkDi")
  # use alternative parameters for H2S? (AD03 Table 1)
  if(altH2S) mod.OBIGT("H2S", state="aq", a=-11.2303, b=12.6104, c=-0.2102)
  # get properties of aq - gas reaction
  sres <- suppressWarnings(subcrt(c(name, name), c("aq", "gas"), c(-1, 1), T = T, P = P))
  # calculate natural logarithm of Henry's constant in mole-fraction units
  ln_KH <- log(1000/18.0153) + log(10) * sres$out$logK
  # plot with units of reciprocal temperature (1000/K)
  TK <- convert(T, "K")
  lty <- 1
  if(altH2S) lty <- 2
  if(HKF) lty <- 3
  if(HKF) col <- "red" else col <- "black"
  lines(1000/TK, ln_KH, lty = lty, col = col)
  reset()
}

# set up plot
par(mfrow=c(2, 2))
ylab <- quote(ln~italic(K[H]))
xlab <- quote(1000 / list(italic(T), K))

# CO2 (Fig. 1a of AD03)
plot(0, 0, xlim=c(1, 4), ylim=c(4, 10), xlab=xlab, ylab=ylab)
lines.KH("CO2", 1:373, "Psat")
lines.KH("CO2", seq(100, 650, 10), 500)
lines.KH("CO2", 1:373, "Psat", HKF = TRUE)
lines.KH("CO2", seq(100, 650, 10), 500, HKF = TRUE)
dat <- read.csv(system.file("extdata/cpetc/AD03_Fig1a.csv", package="CHNOSZ"))
points(dat$x, dat$y, pch=dat$pch)
text(3.5, 7.8, quote(italic(P)[sat]))
text(3.05, 9.2, "500 bar")
legend("bottom", c("Data (AD03, Fig. 1a)", "AkDi model", "HKF model"), lty=c(0, 1, 3), pch=c(1, NA, NA), col=c(1, 1, 2), bty="n")
title(main=syslab(c("CO2", "H2O"), dash = " - "))

# H2 (Fig. 1b of AD03)
plot(0, 0, xlim=c(1, 4), ylim=c(8, 12), xlab=xlab, ylab=ylab)
lines.KH("H2", 1:373, "Psat")
lines.KH("H2", seq(100, 650, 10), 1000)
lines.KH("H2", 1:373, "Psat", HKF = TRUE)
lines.KH("H2", seq(100, 650, 10), 1000, HKF = TRUE)
text(3.4, 11.4, quote(italic(P)[sat]))
text(1.5, 11, "1000 bar")
dat <- read.csv(system.file("extdata/cpetc/AD03_Fig1b.csv", package="CHNOSZ"))
points(dat$x, dat$y, pch=dat$pch)
legend("bottomright", c("Data (AD03, Fig. 1b)", "AkDi model", "HKF model"), lty=c(0, 1, 3), pch=c(1, NA, NA), col=c(1, 1, 2), bty="n")
title(main=syslab(c("H2", "H2O"), dash = " - "))

# H2S (Fig. 1c of AD03)
plot(0, 0, xlim=c(1, 4), ylim=c(4, 9), xlab=xlab, ylab=ylab)
lines.KH("H2S", 1:373, "Psat")
lines.KH("H2S", seq(100, 650, 10), 1000)
lines.KH("H2S", 1:373, "Psat", altH2S = TRUE)
lines.KH("H2S", seq(100, 650, 10), 1000, altH2S = TRUE)
lines.KH("H2S", 1:373, "Psat", HKF = TRUE)
lines.KH("H2S", seq(100, 650, 10), 1000, HKF = TRUE)
dat <- read.csv(system.file("extdata/cpetc/AD03_Fig1c.csv", package="CHNOSZ"))
points(dat$x, dat$y, pch=dat$pch)
text(3.4, 6.9, quote(italic(P)[sat]))
text(3.1, 8.6, "1000 bar")
legend("bottom", c("Data (AD03, Fig. 1c)", "AkDi model", "AkDi model (alt. H2S)", "HKF model"), lty=c(0, 1, 2, 3), pch=c(1, NA, NA, NA), col=c(1, 1, 1, 2), bty="n")
title(main=syslab(c("H2S", "H2O"), dash = " - "))

# CH4 (Fig. 1d of AD03)
plot(0, 0, xlim=c(1.5, 4), ylim=c(8, 12), xlab=xlab, ylab=ylab)
lines.KH("CH4", 1:350, "Psat")
lines.KH("CH4", 1:350, "Psat", HKF = TRUE)
dat <- read.csv(system.file("extdata/cpetc/AD03_Fig1d.csv", package="CHNOSZ"))
points(dat$x, dat$y, pch=dat$pch)
text(3.4, 11, quote(italic(P)[sat]))
legend("bottomright", c("Data (AD03, Fig. 1d)", "AkDi model", "HKF model"), lty=c(0, 1, 3), pch=c(1, NA, NA), col=c(1, 1, 2), bty="n")
title(main=syslab(c("CH4", "H2O"), dash = " - "))

