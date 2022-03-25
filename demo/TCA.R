# TCA.R 20171010
# Reproduce Fig. 6 in Canovas and Shock, 2016:
# Plots of the standard partial molal Gibbs energy of reaction for each step in
# the citric acid cycle for temperatures to 500 degrees C and pressures to 5 kbar.
library(CHNOSZ)

# These plots use calories 20220325
E.units("cal")

# species in reactions
NADox <- "NAD(ox)-"; NADred <- "NAD(red)-2"
ADP <- "ADP-3"; ATP <- "ATP-4"
species <- list(
  c("oxaloacetate-2", "pyruvate", "H2O", NADox, "citrate-3", NADred, "CO2", "H+"),
  c("citrate-3", "cis-aconitate-3", "H2O"),
  c("cis-aconitate-3", "H2O", "isocitrate-3"),
  c("isocitrate-3", NADox, "a-ketoglutarate-2", NADred, "CO2"),
  c("a-ketoglutarate-2", ADP, "HPO4-2", NADox, "succinate-2", ATP, NADred, "CO2"),
  c("succinate-2", "fumarate-2", "H2"),
  c("fumarate-2", "H2O", "malate-2"),
  c("malate-2", NADox, "oxaloacetate-2", NADred, "H+"),
  c("pyruvate", NADox, ADP, "HPO4-2", "H2O", "CO2", NADred, "H+", ATP, "H2")
)
# reaction coefficients
coeffs <- list(
  c(-1, -1, -1, -1, 1, 1, 1, 1),
  c(-1, 1, 1),
  c(-1, -1, 1),
  c(-1, -1, 1, 1, 1),
  c(-1, -1, -1, -1, 1, 1, 1, 1),
  c(-1, 1, 1),
  c(-1, -1, 1),
  c(-1, -1, 1, 1, 1),
  c(-1, -4, -1, -1, -2, 3, 4, 2, 1, 1)
)
# species names
oxal <- quote(Oxaloacetate^-2)
pyr <- quote(Pyruvate^-"")
h2o <- quote(H[2]*O)
nox <- quote(NAD[ox]^-"")
cit <- quote(Citrate^-3)
nred <- quote(NAD[red]^-2)
co2 <- quote(CO[2*(italic(aq))])
hplus <- quote(H^+"")
iso <- quote(Isocitrate^-3)
aco <- quote(italic(cis)*"-Aconitate"^-3)
ket <- quote(alpha*"-Ketoglutarate"^-2)
adp <- quote(ADP^-3)
hpo4 <- quote(HPO[4]^-2)
suc <- quote(Succinate^-2)
atp <- quote(ATP^-4)
fum <- quote(Fumarate^-2)
h2 <- quote(H[2*(italic(aq))])
mal <- quote(Malate^-2)
# the reaction double arrow
eq <- "\u21cc"
sublist <- list(oxal=oxal, pyr=pyr, h2o=h2o, nox=nox, cit=cit, nred=nred,
                co2=co2, hplus=hplus, aco=aco, iso=iso, ket=ket, adp=adp,
                hpo4=hpo4, suc=suc, atp=atp, fum=fum, h2=h2, mal=mal, eq=eq)
# reaction titles
rtitle <- list(
  c(substitute("        "*oxal + pyr + h2o + nox ~eq~ "", sublist), substitute(cit + nred + co2 + hplus, sublist)),
  substitute(cit ~eq~ aco + h2o, sublist),
  substitute(aco + h2o ~eq~ iso*"   ", sublist),
  c(substitute(iso + nox ~eq~ "   ", sublist), substitute(ket + nred + co2*"   ", sublist)),
  c(substitute(ket + adp + hpo4 + nox ~eq~ "", sublist), substitute("      "*suc + atp + nred + co2, sublist)),
  c(substitute(suc ~eq~ "", sublist), substitute(fum + h2, sublist)),
  substitute(fum + h2o ~eq~ mal, sublist),
  c(substitute(mal + nox ~eq~ "                          ", sublist), substitute(oxal + nred + hplus * "             ", sublist)),
  c(substitute(pyr + 4*nox + adp + hpo4 + 2*h2o ~eq~ "                    ", sublist),
    substitute(3*co2 + 4*nred + 2*hplus + atp + h2 * "           ", sublist))
)
# set up plot
opar <- par(no.readonly = TRUE)
par(mfrow=c(3, 3))
ylims <- list(
  c(-10, 45), c(1, 6),   c(-2.5, 7.5),
  c(-35, 5),  c(-9, 5),  c(5, 28),
  c(-1.5, 6),   c(14, 18), c(20, 80)
)
# loop over reactions
for(i in seq_along(species)) {
  thermo.plot.new(xlim=c(0, 500), ylim=ylims[[i]], xlab=axis.label("T"),
                  ylab=axis.label("DrG0", prefix="k"), mar=c(3.0, 3.5, 3.5, 2.0))
  # loop over isobars
  for(P in seq(500, 5000, 500)) {
    T <- seq(0, 500, 10)
    if(P==500) T <- seq(0, 350, 10)
    if(P==5000) T <- seq(100, 500, 10)
    # calculate and plot standard Gibbs energy
    sout <- subcrt(species[[i]], coeffs[[i]], T=T, P=P)$out
    lines(T, sout$G/1000)
  }
  if(is.list(rtitle[[i]])) mtitle(as.expression(rtitle[[i]]), spacing = 1.6, cex=0.8)
  else mtitle(as.expression(rtitle[[i]]), line=0.4, cex=0.8)
}
# make an overall title
par(xpd=NA)
text(-70, 284, "Citric Acid Cycle, after Canovas and Shock, 2016", font=2, cex=1.5)
par(xpd=FALSE)
par(opar)

# Reset the units
reset()
