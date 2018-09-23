# CHNOSZ/demos/adenine.R
# plot thermodynamic data and model fits for aqueous adenine 20170725

# notable functions in this demo:
# EOSregress() - to regress HKF coefficients from Cp data
# add.obigt() - to add thermodynamic data that are not in the default database

# LCT17 = Lowe et al., 2017 (J. Chem. Thermodyn., doi:10.1016/j.jct.2017.04.005)
# LH06 = LaRowe and Helgeson, 2006 (Geochim. Cosmochim. Acta, doi:10.1016/j.gca.2006.04.010)
# HKF = Helgeson-Kirkham-Flowers equations (e.g. Am. J. Sci., doi:10.2475/ajs.281.10.1249)

# adenine Cp0 and V0 from LCT17 Table 4
AdH <- data.frame(
  T = seq(288.15, 363.15, 5),
  V = c(89.1, 89.9, 90.6, 91.3, 92, 92.7, 93.1, 93.6, 94.1, 94.9, 95.4, 95.9, 96.3, 96.9, 97.1, 97.8),
  V_SD = c(1.1, 1.3, 1.1, 1, 1.1, 1, 0.9, 0.9, 0.8, 0.6, 0.7, 0.7, 0.7, 0.4, 1.1, 0.7),
  Cp = c(207, 212, 216, 220, 224, 227, 230, 234, 237, 241, 243, 245, 248, 250, 252, 255),
  Cp_SD = c(5, 7, 8, 7, 8, 7, 6, 7, 6, 5, 6, 6, 5, 4, 7, 5)
)
# functions to calculate V and Cp using density model (LCT17 Eq. 28)
Vfun <- function(v1, v2, k, T) {
  # gas constant (cm3 bar K-1 mol-1)
  R <- 83.144598
  # isothermal compressibility (bar-1)
  beta <- water("beta", TK)$beta
  v1 + v2 / (T - 228) - k * R * beta
}
Cpfun <- function(c1, c2, k, T) {
  # gas constant (J K-1 mol-1)
  R <- 8.3144598
  # isobaric temperature derivative of expansivity (K-2)
  daldT <- water("daldT", TK)$daldT
  c1 + c2 / (T - 228) ^ 2 - k * R * T * daldT
}
# set up units (used for axis labels and HKF calculations)
E.units("J")
T.units("K")
# temperature and pressure points for calculations
TK <- seq(275, 425)
P <- water("Psat", TK)$Psat
# set up plots
layout(matrix(1:3), heights=c(1, 8, 8))
# title at top
par(mar=c(0, 0, 0, 0), cex=1)
plot.new()
text(0.5, 0.5, "Heat capacity and volume of aqueous adenine\n(Lowe et al., 2017)", font=2)
# settings for plot areas
par(mar = c(4, 4, 0.5, 1), mgp = c(2.4, 0.5, 0))
# location of x-axis tick marks
xaxp <- c(275, 425, 6)

### Cp plot (LCT17 Figures 4 and 12) ###
plot(AdH$T, AdH$Cp, type = "p", xlim = range(TK), ylim = c(150, 350),
     xlab = axis.label("T"), ylab = axis.label("Cp0"), 
     pch = 5, tcl = 0.3, xaxs = "i", yaxs = "i", las = 1, xaxp = xaxp
)
# error bars (arrows trick from https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars)
arrows(AdH$T, AdH$Cp - AdH$Cp_SD, AdH$T, AdH$Cp + AdH$Cp_SD, length = 0.05, angle = 90, code = 3)
# get LH06 predictions using HKF model;
# this version of adenine parameters has been superseded by LCT17,
# so we add them by hand
mod.obigt("adenine-old", formula="C5H5N5", a1=21.5046, a2=8.501, a3=-2.6632, a4=-5.3561, c1=87.88, c2=-15.87, omega=0.065)
LH06 <- subcrt("adenine-old", T = TK)$out$adenine
lines(TK, LH06$Cp, lty = 3)
# density model (parameters from LCT17 Table 11)
lines(TK, Cpfun(160.4, -653239, -7930.3, TK), lty = 2)

# regress HKF parameters
# specify the terms in the HKF equations
var <- c("invTTheta2", "TXBorn")
# build a data frame with T, P, and Cp columns
Cpdat <- data.frame(T = AdH[, "T"], P = 1, Cp = AdH[, "Cp"])
# convert Cp data from J to cal
Cpdat$Cp <- convert(Cpdat$Cp, "cal")
# regress HKF parameters from Cp data
HKFcoeffs <- EOSregress(Cpdat, var)$coefficients
# get predictions from the fitted model
Cpfit <- EOScalc(HKFcoeffs, TK, P)
# plot the fitted model
lines(TK, convert(Cpfit, "J"), lwd = 2, col = "green3")
# format coefficients for legend; use scientific notation for c2 and omega
coeffs <- format(signif(HKFcoeffs, 4), scientific = TRUE)
# keep c1 out of scientific notation
coeffs[1] <- signif(HKFcoeffs[[1]], 4)
ipos <- which(coeffs >= 0)
coeffs[ipos] <- paste("+", coeffs[ipos], sep = "")
fun.lab <- as.expression(lapply(1:length(coeffs), function(x) {
  EOSlab(names(HKFcoeffs)[x], coeffs[x])
}))
# add legend: regressed HKF coefficients
legend("topleft", legend = fun.lab, pt.cex = 0.1, box.col = "green3")
# add legend: lines
legend("bottomright", lty = c(3, 2, 1), lwd = c(1, 1, 2), col = c("black", "black", "green3"), bty = "n",
  legend = c("HKF model (LaRowe and Helgeson, 2006)",
  "density model (Lowe et al., 2017)", "HKF model (fit using CHNOSZ)")
)

### V plot (LCT17 Figures 3 and 11) ###
plot(AdH$T, AdH$V, type = "p", xlim = range(TK), ylim = c(85, 105),
     xlab = axis.label("T"), ylab = axis.label("V0"), 
     pch = 5, tcl = 0.3, xaxs = "i", yaxs = "i", las = 1, xaxp = xaxp
)
axis(3, labels = FALSE, tcl = 0.3, xaxp = xaxp)
axis(4, labels = FALSE, tcl = 0.3)
arrows(AdH$T, AdH$V - AdH$V_SD, AdH$T, AdH$V + AdH$V_SD, length = 0.05, angle = 90, code = 3)
# HKF model with coefficients from LH06
lines(TK, LH06$V, lty = 3)
# density model with coefficients from LCT17
lines(TK, Vfun(73.9, -917.0, -7930.3, TK), lty = 2)
# HKF heat capacity coefficients from LCT17
LCT17 <- subcrt("adenine", T = TK)$out$adenine
lines(TK, LCT17$V, lwd = 2, col = "royalblue")
legend("bottomright", lty = c(3, 2, 1), lwd = c(1, 1, 2), col = c("black", "black", "royalblue"), bty = "n",
  legend=c("HKF model (LaRowe and Helgeson, 2006)",
  "density model (Lowe et al., 2017)", "HKF model (fit by Lowe et al., 2017 using CHNOSZ)")
)
# reset database and computational settings
data(thermo)
