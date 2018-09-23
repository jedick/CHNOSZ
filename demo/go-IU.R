# CHNOSZ/demo/go-IU.R  20171018
# diagrams using data from the SUPCRTBL compilation
# (BL = Bloomington campus of Indiana University)

## set up plotting area
par(mfrow=c(2, 2))

## start with default database
data(thermo)

###########
### plot 1: boehmite - kaolinite equilibrium
###########
## experimental data from Table 1 of Hemley et al., 1980
# doi:10.2113/gsecongeo.75.2.210
xT <- c(200, 200, 200, 200, 250, 250, 300, 300, 300, 300)
xlogaSiO2 <- -c(2.54, 2.59, 2.65, 2.77, 2.21, 2.32, 1.90, 1.95, 1.94, 1.90)
## set up basis species so that axis.label shows activity of SiO2
basis(c("Al2O3","SiO2", "H2O", "O2"))
T <- 125:350
thermo.plot.new(xlim=range(T), ylim=c(-3.5, -1.5), xlab = axis.label("T"), ylab=axis.label("SiO2"))
points(xT, xlogaSiO2)
basis(delete=TRUE)
## first calculation: after SUPCRT92 (SiO2(aq) from SHS89)
add.obigt("SUPCRT92") # gets kaolinite and boehmite from HDNB78
r1 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T=T, P=1000, exceed.Ttr = TRUE) 
# we need exceed.Ttr = TRUE because the T limit for boehmite is 500 K (Helgeson et al., 1978)
## second calculation: kaolinite from Berman, 1988
data(thermo)
add.obigt("SUPCRT92", "boehmite") # gets only boehmite from HDNB78
r2 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T=T, P=1000, exceed.Ttr = TRUE) 
## third calculation: boehmite from Hemingway et al., 1991
add.obigt("SUPCRTBL", "boehmite")
r3 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T=T, P=1000) 
## fourth calculation: SiO2 from Apps and Spycher, 2004
add.obigt("SUPCRTBL", "SiO2")
r4 <- subcrt(c("boehmite", "H2O", "SiO2", "kaolinite"), c(-1, -0.5, -1, 0.5), T=T, P=1000) 
## log activity of SiO2 is -ve logK
lines(T, -r1$out$logK)
lines(T, -r2$out$logK, lty=2)
lines(T, -r3$out$logK, lty=2, col="red")
lines(T, -r4$out$logK, col="red")
## add points calculated using the SUPCRTBL package
points(seq(125, 350, 25), -c(3.489, 3.217, 2.967, 2.734, 2.517, 2.314, 2.124, 1.946, 1.781, 1.628), pch=4, col="red")
## add labels, legend, and title
text(182.5, -3.17, "SUPCRT92", srt=48, cex=0.7, font=2)
text(154, -3.04, "circa SUPCRTBL", srt=45, cex=0.7, font=2, col="red")
legend("topleft", lty=c(1, 2, 2, 1, 0), pch=c(NA, NA, NA, NA, 4),
       col=c("black", "black", "red", "red", "red"), bty="n", cex=0.9,
       legend=c("Kln,Bhm:HDNB78; SiO2(aq):SHS89", "Kln:Ber88 (SUPCRTBL uses HP11)",
                "+ Bhm:HRA91", "+ SiO2(aq):AS04", "SUPCRTBL"))
legend("bottomright", pch=1, legend="Hemley et al., 1980", bty="n", cex=0.9)
mtitle(c("Kaolinite - Boehmite", "After Zhu and Lu, 2009 Fig. A1"), cex=0.95)
# Zhu and Lu, 2009: doi:10.1016/j.gca.2009.03.015
# Helgeson et al., 1978 (HDNB78): http://www.worldcat.org/oclc/13594862
# Shock et al., 1989 (SHS89): doi:10.1016/0016-7037(89)90341-4
# Berman, 1988 (Ber88): doi:10.1093/petrology/29.2.445
# Holland and Powell, 2011 (HP11): 10.1111/j.1525-1314.2010.00923.x
# Hemingway et al., 1991 (HRA91): http://pubs.er.usgs.gov/publication/70016664
# Apps and Spycher, 2004 (AS04): Bechtel SAIC Company, LLC ANL-NBS-HS-000043 REV 00 (DOC.20041118.0004)

###########
### plot 2: dawsonite solubility
###########
## experimental data from Benezeth et al., 2007 Table 5
# doi:10.1016/j.gca.2007.07.003
# (averages for each temperature in a single run)
T <- c(100.1, 100.1, 150.1, 100.1, 150.1, 99.8, 99.8, 200.7, 99.8, 50.1, 75.1, 100.3, 150.1)
logK <- -c(14.825, 14.735, 13.625, 14.79, 13.665, 14.725, 14.1775, 12.74, 14.4925, 16.8625, 15.61, 14.51, 13.455)
plot(T, logK, xlim=c(25, 250), ylim=c(-18, -10), xlab=axis.label("T"), ylab=axis.label("logK"))
# this gets us dawsonite and Al(OH)4-
add.obigt("SUPCRTBL")
T <- 0:250
# calculation 1: dawsonite with non-zero Cp
species <- c("dawsonite", "H2O", "Al(OH)4-", "HCO3-", "Na+", "H+")
coeffs <- c(-1, -2, 1, 1, 1, 1)
Daw1 <- subcrt(species, coeffs, T=T)
# calculation 2: dawsonite with 0 Cp
mod.obigt("dawsonite", Cp=0)
Daw2 <- subcrt(species, coeffs, T=T)
## plot the calculated logKs
lines(T, Daw1$out$logK, col="red")
lines(T, Daw2$out$logK, col="red", lty=2)
## add points calculated using the SUPCRTBL package
points(seq(25, 250, 25), c(-17.829, -16.523, -15.402, -14.425, -13.568, -12.815, -12.154, -11.581, -11.094, -10.699), pch=4, col="red")
## add labels, legend, and title
text(182.5, -3.17, "SUPCRT92\n(CHNOSZ default)", srt=43, cex=0.7, font=2)
text(145, -3.1, "SUPCRTBL", srt=41.5, cex=0.7, font=2, col="red")
legend("topleft", lty=c(1, 2, 0), pch=c(NA, NA, 4), col="red", bty="n", cex=0.9,
       legend=c("Daw Cp != 0", "Daw Cp = 0", "SUPCRTBL"))
legend("bottomright", pch=1, legend="Ben\u00e9z\u00e9th et al., 2007", bty="n", cex=0.9)
mtitle(c("Dawsonite - aqueous species", "After Zimmer et al., 2016 Fig. 2"), cex=0.95)
# doi:10.1016/j.cageo.2016.02.013

###########
### plot 3: Eh-pH diagram for As-O-H-S
###########
add.obigt("SUPCRTBL")
#basis(c("Fe", "As", "H2O", "H2S", "H+", "e-"))
#basis(c("Fe", "H2S"), c(-6, -3))
basis(c("As", "H2O", "H2S", "H+", "e-"))
basis(c("H2S"), c(-3))
As_aq <- c("H3AsO4", "H2AsO4-", "HAsO4-2", "AsO4-3", "H3AsO3", "H2AsO3-", "HAsO3-2", "AsO3-3")
AsS_aq <- c("AsS(OH)HS-", "As3S4(HS)2-")
As_cr <- "As"
AsS_cr <- c("realgar,alpha", "realgar,beta", "orpiment", "orpiment,amorphous")
FeAs_cr <- c("arsenopyrite", "scorodite", "ferric arsenate,amorphous")
#species(c(As_aq, AsS_aq, As_cr, AsS_cr, FeAs_cr))
species(c(As_aq, AsS_aq, As_cr, AsS_cr))
species(c(As_aq, AsS_aq), -5)
## a simple diagram, but using only H2S at all pH
#a <- affinity(pH=c(0, 14), Eh=c(-1, 1.5))
#diagram(a)
# the S basis species depends on pH
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
# calculate affinties of formation reactions using the speciated S basis species
res <- 300
# we "blend" the transitions with pH, unlike LZ11's diagram where
# it appears the S-basis species are switched in an on/off fashion
m <- mosaic(bases, pH=c(0, 14, res), Eh=c(-0.8, 0.8, res), blend=TRUE)
# adjust colors and names
fill <- rev(heat.colors(nrow(species())))
fill[11:15] <- "darkgrey"
m$A.species$species$name <- gsub(",alpha", "", m$A.species$species$name)
diagram(m$A.species, fill=fill)
dprop <- describe.property(c("T", "P"), c(25, 1))
legend("bottomleft", legend=dprop, bty="n")
t1 <- quote("As-O-H-S, "~list(Sigma*S==10^-3~M, Sigma*As==10^-5~M))
t2 <- "After Lu and Zhu, 2011 Fig. 2b"
# doi:10.1007/s12665-010-0652-x
mtitle(as.expression(c(t1, t2)), cex=0.95)

###########
### plot 4: aqueous Al species
###########
add.obigt("SUPCRTBL")
basis(c("Al+3", "F-", "H+", "O2", "H2O"))
AlOH <- c("Al(OH)4-", "Al(OH)3", "Al(OH)2+", "AlOH+2")
Al <- "Al+3"
AlF <- c("AlF+2", "AlF2+", "AlF3", "AlF4-")
AlOHF <- c("Al(OH)2F2-", "Al(OH)2F", "AlOHF2")
species(c(AlOH, Al, AlF, AlOHF), "aq")
res <- 300
a <- affinity(pH=c(0.5, 6.5, res), `F-`=c(-2, -9, res), T=200)
diagram(a, fill=rev(cm.colors(nrow(species()))))
dprop <- describe.property(c("T", "P"), c(200, "Psat"))
legend("topright", legend=dprop, bty="n")
mtitle(c("Aqueous aluminum species",
         "After Tagirov and Schott, 2001 Fig. 4d"), cex=0.95)
# doi:10.1016/S0016-7037(01)00705-0

###########
### clean up: restore thermodynamic database to default
###########
data(thermo)
