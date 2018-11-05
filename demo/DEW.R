# demo for the Deep Earth Water (DEW) model in CHNOSZ 20170927

# set up subplots
par(mfrow = c(2, 2), mar=c(3.0, 3.5, 2.5, 1.0), mgp=c(1.7, 0.3, 0), las=1, tcl=0.3, xaxs="i", yaxs="i")

# activate DEW model
oldwat <- water("DEW")

###########
#### plot 1: quartz solubility at high pressure
## after Figure 7D of Sverjensky et al., 2014a [SHA14]
## (Geochim. Cosmochim. Acta, https://doi.org/10.1016/j.gca.2013.12.019)
###########

# load SiO2 and Si2O4 data taken from DEW spreadsheet
iSi <- add.obigt("DEW_aq", c("SiO2", "Si2O4"))
# print the data references to confirm we got the right ones
thermo.refs(iSi)
# set temperature ranges for different pressures
# data.frame is used to make P and T the same length
PT0.5 <- data.frame(P=500, T=seq(200, 550, 10))
PT1.0 <- data.frame(P=1000, T=seq(200, 700, 10))
PT2.0 <- data.frame(P=2000, T=seq(200, 700, 10))
PT5.0 <- data.frame(P=5000, T=seq(200, 850, 10))
PT10.0 <- data.frame(P=10000, T=seq(200, 825, 10))
PT20.0 <- data.frame(P=20000, T=seq(200, 800, 10))
PT <- rbind(PT0.5, PT1.0, PT2.0, PT5.0, PT10.0, PT20.0)
# reaction 1: quartz = SiO2(aq) [equivalent to quartz + 3 H2O = Si(OH)4]
SiO2_logK <- subcrt(c("quartz", "SiO2"), c("cr", "aq"), c(-1, 1), P=PT$P, T=PT$T)$out$logK
# reaction 2: 2 quartz = Si2O4(aq) [equivalent to 2 quartz + 3 H2O = Si2O(OH)6]
Si2O4_logK <- subcrt(c("quartz", "Si2O4"), c("cr", "aq"), c(-2, 1), P=PT$P, T=PT$T)$out$logK
# plot the sum of molalities (== activities) for each pressure
plot(c(200, 1000), c(-2.5, 0.5), type="n", xlab=axis.label("T"), ylab="log molality")
for(P in unique(PT$P)) {
  icond <- PT$P == P
  SiO2_logm <- SiO2_logK[icond]
  Si2O4_logm <- Si2O4_logK[icond]
  logm <- log10(10^SiO2_logm + 10^Si2O4_logm)
  lines(PT$T[icond], logm)
  # add text label
  lastT <- tail(PT$T[icond], 1)
  Pkb <- paste(format(P/1000, nsmall=1), "kb")
  text(lastT+25, tail(logm, 1), Pkb, adj=0)
}
t1 <- quote("Solubility of"~alpha*"-quartz")
t2 <- "after Sverjensky et al., 2014a"
mtitle(as.expression(c(t1, t2)))
# TODO: lines are a little low at highest P and T ...

###########
#### plot 2: correlations between non-solvation volume and HKF a1 parameter
## after Figures 12B and 12C of Sverjensky et al., 2014a [SHA14]
###########

# load the fitted parameters for species as used by SHA14
# TODO: also use their Ca+2??
# NOTE: don't load NaCl, NH4+, or HS- here because the DEW spreadsheet lists a1 from the correlation
add.obigt("DEW", c("CO3-2", "BO2-", "MgCl+", "SiO2", "HCO3-", "Si2O4"))
# set up the plot
V0nlab <- expression(Delta * italic(V) * degree[n]~~(cm^3~mol^-1))
a1lab <- expression(italic(a)[1]%*%10~~(cal~mol~bar^-1))
plot(c(-25, 50), c(-4, 12), type="n", xlab=V0nlab, ylab=a1lab)
# a function to get the HKF parameters, calculate nonsolvation volume, plot points, labels, error bars, and correlation lines
plotfun <- function(species, col, pch, cex, dy, error, xlim, corrfun) {
  # get HKF parameters
  par <- info(info(species))
  a1 <- par$a1 * 10
  # get the nonsolvation volume
  Vn <- unlist(hkf("V", par, contrib="n")$aq)
  points(Vn, a1, col=col, pch=pch, cex=cex)
  for(i in 1:length(species)) text(Vn[i], a1[i]+dy, expr.species(species[i]))
  arrows(Vn, a1 - error, Vn, a1 + error, length = 0.03, angle = 90, code = 3, col=col)
  lines(xlim, corrfun(xlim), col=col)
}
# monovalent ions: Na+, K+, Cl-, Br-
monofun <- function(Vn) 2.0754 + 0.10871 * Vn
# for easier reading, set y-offset to NA so the labels aren't plotted
plotfun(c("Na+", "K+", "Cl-", "Br-"), "red", 19, 1, NA, 0.5, c(-7, 35), monofun)
# divalent ions: Mg+2, Ca+2, CO3-2, SO4-2
difun <- function(Vn) 3.5321 + 0.23911 * Vn
plotfun(c("Mg+2", "Ca+2", "CO3-2", "SO4-2"), "black", 15, 1, 1.2, 0.7, c(-20, 25), difun)
# complexes and neutral molecules: BO2-, MgCl+, SiO2, NaCl, HCO3-, Si2O4, NH4+, HS-
compfun <- function(Vn) 1.5204 + 0.19421 * Vn
plotfun(c("MgCl+", "SiO2", "NaCl", "HCO3-", "Si2O4"), "blue1", 18, 1.5, 1, 0.5, c(-20, 50), compfun)
# for easier reading, put some labels below the points
plotfun(c("BO2-", "NH4+", "HS-"), "blue1", 18, 1.5, -1.2, 0.5, c(-20, 50), compfun)
# include an empty subscript for better spacing between the lines
t1 <- quote("Correlations between non-solvation"[])
t2 <- quote("volume and HKF "*italic(a)[1]*" parameter")
mtitle(as.expression(c(t1, t2)))

###########
#### plot 3: logfO2-pH diagram for aqueous inorganic and organic carbon species at high pressure
## after Figure 1b of Sverjensky et al., 2014b [SSH14]
## (Nature Geoscience, https://doi.org/10.1038/NGEO2291)
###########

# define system with loga.species = 0
basis("CHNOS+")
species(c("CO2", "HCO3-", "CO3-2", "acetic acid", "acetate", "CH4"))
species(1:6, 0)

# a function to make the diagrams
dfun <- function(T = 600, P = 50000, res=300) {
  a <- affinity(pH = c(0, 10, res), O2 = c(-24, -12, res), T = T, P = P)
  diagram(a, limit.water = FALSE, fill=tail(topo.colors(7), -1))
  dp <- describe.property(c("     T", "     P"), c(T, P), digits=0)
  legend("bottomleft", legend=dp, bty="n")
}

data(OBIGT)
## (not run) make diagram using CHNOSZ default database
#dfun()
#t1 <- quote("CHNOSZ default database"[])
#t2 <- quote("(not recommended for high"~italic(P)*")")
#mtitle(as.expression(c(t1, t2)))
# make diagram using CO2, HCO3-, CO3-2, and methane data from DEW spreadsheet
add.obigt("DEW_aq", c("CO2", "HCO3-", "CO3-2", "methane"))
dfun()
CO2quote <- quote(list(CO[2], HCO[3]^"-", CO[3]^"-2"))
DEWexpr <- substitute("DEW data for"~x, list(x=CO2quote))
mtitle(as.expression(c(DEWexpr, "and methane")))

###########
#### plot 4: speciation of carbon as a function T, logfO2 and pH (added 20171008)
## after SSH14 Fig. 3
###########

# conditions:
# T = 600, 700, 800, 900, 1000 degC
# P = 5.0GPa (50000 bar)
# fO2 = QFM - 2
# pH set by jadeite + kyanite + coesite
# output from EQ3NR calculations (SSH14 Supporting Information)
# dissolved carbon: 0.03, 0.2, 1, 4, 20 molal
# true ionic strength: 0.39, 0.57, 0.88, 1.45, 2.49
# pH: 3.80, 3.99, 4.14, 4.25, 4.33
## activate DEW model
data(thermo)
water("DEW")
# add species data for DEW
inorganics <- c("methane", "CO2", "HCO3-", "CO3-2")
organics <- c("formic acid", "formate", "acetic acid", "acetate", "propanoic acid", "propanoate")
# skip updating acetate because the new data from the DEW spreadsheet give different logK
add.obigt("DEW", c(inorganics, organics[-4]))
## set basis species
basis(c("Fe", "SiO2", "CO3-2", "H2O", "oxygen", "H+"))
## calculate logfO2 in QFM buffer
basis("O2", "QFM")
T <- seq(600, 1000, 100)
buf <- affinity(T=T, P=50000, return.buffer=TRUE)
## add species
species(c(inorganics, organics))
## generate spline functions from IS, pH, and molC values at every 100 degC
IS <- c(0.39, 0.57, 0.88, 1.45, 2.49)
pH <- c(3.80, 3.99, 4.14, 4.25, 4.33)
molC <- c(0.03, 0.2, 1, 4, 20)
## use extended Debye-Huckel equation with b_gamma set to zero
nonideal("bgamma0")
## calculate affinities on the T-logfO2-pH-IS transect
a <- affinity(T = T, O2 = buf$O2 - 2, IS = IS, pH = pH, P = 50000)
## calculate metastable equilibrium activities using the total
## carbon molality as an approximation of total activity
e <- equilibrate(a, loga.balance = log10(molC))
## make the diagram; don't plot names of low-abundance species
names <- c(inorganics, organics)
names[c(4, 5, 7, 9)] <- ""
col <- rep("black", length(names))
col[c(1, 3, 6, 8, 10)] <- c("red", "darkgreen", "purple", "orange", "navyblue")
if(packageVersion("CHNOSZ") > "1.1.3") {
  diagram(e, alpha = "balance", names = names, col = col, ylim = c(0, 0.8), ylab="carbon fraction", spline.method="natural")
} else {
  diagram(e, alpha = "balance", names = names, col = col, ylim = c(0, 0.8), ylab="carbon fraction")
}


## add legend and title
ltxt1 <- "P = 50000 bar"
ltxt2 <- substitute(logfO2=="QFM-2", list(logfO2 = axis.label("O2")))
pH <- seq(3.8, 4.3, length.out = length(T))
legend("left", legend = as.expression(c(ltxt1, ltxt2)), bty = "n")
t1 <- "Aqueous carbon speciation"
t2 <- "after Sverjensky et al., 2014b"
mtitle(c(t1, t2))

### additional checks

## check that we're within 0.1 of the QFM-2 values used by SSH14
stopifnot(maxdiff((buf$O2-2), c(-17.0, -14.5, -12.5, -10.8, -9.4)) < 0.1)

# Here are the logKs of aqueous species dissociation reactions at 600 degC and 50000 bar,
# values from EQ3NR output in Supporting Information of the paper (p. 103-109):
inorganic.logK <- c(24.4765, -9.0784, -5.3468, 0)
organic.logK <- c(1.7878, 2.5648, 15.3182, 16.9743, 30.4088, 28.9185)
# calculate equilibrium constants of the reactions in CHNOSZ; use a negative sign to change from formation to dissociation
logK.calc <- -unlist(affinity(T = 600, P = 50000, property = "logK")$values)
logK.calc - c(inorganic.logK, organic.logK)
## check that we're within 0.021 of the logK values used by SSH14
stopifnot(maxdiff(logK.calc, c(inorganic.logK, organic.logK)) < 0.021)

## check that we get similar activity coefficients
# activity coefficients for monovalent species from EQ3NR output
loggamma <- c(-0.15, -0.18, -0.22, -0.26, -0.31)
# activity coefficients calculated in CHNOSZ
sres <- subcrt("propanoate", T = seq(600, 1000, 100), P = 50000, IS = c(0.39, 0.57, 0.88, 1.45, 2.49))
stopifnot(maxdiff(sres$out[[1]]$loggam, loggamma) < 0.023)
# if m_star in nonideal() was zero, we could decrease the tolerance here
#stopifnot(maxdiff(sres$out[[1]]$loggam, loggamma) < 0.004)

###########
### all done!
# reset the database and previous water computational option
data(OBIGT)
water(oldwat)
###########
