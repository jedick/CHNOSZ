## affinities of metabolic reactions
## after Amend and Shock, 2001, Fig. 7
##  Amend, J. P. and Shock, E. L. (2001) Energetics of overall metabolic reactions of thermophilic and hyperthermophilic Archaea and Bacteria. \emph{FEMS Microbiol. Rev.} \bold{25}, 175--243. \url{https://doi.org/10.1016/S0168-6445(00)00062-0}
# use aq state for all basis species (including O2)
basis(c("CO2", "H2", "NH3", "O2", "H2S", "H+"), "aq")
# we're going to make H2O
species("H2O")
# a function to create the plots
doplot <- function(T) {
  res <- 20
  # calculate affinity/2.303RT as a function of loga(H2) and loga(O2)
  a <- affinity(H2=c(-10, 0, res), O2=c(-10, 0, res), T=T)
  T.K <- convert(T, "K")                   # temperature in Kelvin
  acal <- convert(a$values[[1]], "G", T.K) # affinity (cal/mol)
  akJ <- convert(acal, "J")/1000           # affinity (kJ/mol)
  # now contour the values
  xyvals <- seq(-10, 0, length.out=res)
  contour(x=xyvals, y=xyvals, z=t(akJ), levels=seq(-150, -250, -20),
    labcex=1, xlab=axis.label("H2"), ylab=axis.label("O2"))
  # show the temperature
  legend("topleft", bg="white", cex=1,
    legend=describe.property("T", T, digits=0, ret.val=TRUE) )
}
# plot layout with space for title at top
layout(matrix(c(1, 1, 2, 3, 4, 5), ncol=2, byrow=TRUE), heights=c(1, 4, 4))
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
# we use subcrt() to generate a reaction for titling the plot
rxnexpr <- describe.reaction(subcrt("H2O", 1)$reaction, states="all")
# also in the title is the property with its units
E.units("J")
Gexpr <- axis.label("DGr", prefix="k")[[2]]
text(0.5, 0.6, substitute(paste(G~~"for"~~r), list(G=Gexpr, r=rxnexpr)), cex=2)
text(0.5, 0.2, "after Amend and Shock, 2001 Figure 7", cex=2)
# now make the plots
par(mar=c(3, 3, 0.5, 0.5), cex=1.3, mgp=c(2, 1, 0))
sapply(c(25, 55, 100, 150), doplot)
# affinity() can handle the three dimensions simultaneously
print(affinity(H2=c(-10, 0, 3), O2=c(-10, 0, 3), T=c(25, 150, 4))$values)
# this is so the plots in the next examples show up OK
E.units("cal")
layout(matrix(1))
par(opar)

## amino acid synthesis at low and high temperatures
## after Amend and Shock, 1998
##  Amend, J. P. and Shock, E. L. (1998) Energetics of amino acid synthesis in hydrothermal ecosystems. \emph{Science} \bold{281}, 1659--1662. \url{https://doi.org/10.1126/science.281.5383.1659}
# select the basis species and species of interest
# and set their activities, first for the 18 degree C case
basis(c("H2O", "CO2", "NH4+", "H2", "H+", "H2S"),
  log10(c(1, 1e-4, 5e-8, 2e-9, 5e-9, 1e-15)))
species(sort(aminoacids("Z")),
  log10(c(3.9, 0.7, 1.1, 3.3, 0.5, 3.8, 1.0, 5.8, 1.2, 0.7,
  0.8, 1.0, 2.8, 0.5, 0.5, 4.6, 5.8, 0.6, 0.9, 2.8)/1e9))
T <- 18
TK <- convert(T, "K")
# calculate A/2.303RT (dimensionless), convert to G of reaction (cal/mol)
a <- affinity(T=T)
G.18.cal <- convert(unlist(a$values), "G", T=TK)
# covvert to kJ/mol
G.18.kJ <- convert(G.18.cal, "J")/1000
# the 100 degree C case
basis(c("H2O", "CO2", "NH4+", "H2", "H+", "H2S"),
  log10(c(1, 2.2e-3, 2.9e-6, 3.4e-4, 1.9e-6, 1.6e-3)))
species(1:20, log10(c(2.8e-9, 5.0e-10, 7.9e-10, 2.4e-9, 3.6e-10,
  2.7e-9, 7.2e-10, 4.2e-9, 8.6e-10, 5.0e-10, 5.7e-10, 7.2e-10, 2.0e-9,
  3.6e-10,3.6e-10, 3.3e-9, 4.2e-9, 4.3e-10, 6.5e-10, 2.0e-9)))
T <- 100
TK <- convert(T, "K")
a <- affinity(T=T)
G.100.cal <- convert(unlist(a$values), "G", T=TK)
G.100.kJ <- convert(G.100.cal, "J")/1000
# the average oxidation states of carbon
Z.C <- ZC(thermo()$obigt$formula[thermo()$species$ispecies])
# put everything together a la Table 3 in the paper
print(out <- data.frame(G.18=G.18.kJ, G.100=G.100.kJ, Z.C=Z.C))
# make a plot; set units to get correct label
E.units("J")
plot(out$Z.C, out$G.18, pch=20, xlim=c(-1.1, 1.1), ylim=c(-200, 500), 
  xlab=axis.label("ZC"), ylab=axis.label("DGr", prefix="k"))
points(out$Z.C, out$G.100, col="red", pch=20)
legend("topleft", pch=c(20, 20), col=c("black", "red"),
  legend=describe.property(c("T", "T"), c(18, 100)))
title(main="Amino acid synthesis, after Amend and Shock, 1998")
# 9 amino acids have negative delta Gr under hydrothermal conditions
# (cf. AS98 with 11; we are using more recent thermodynamic data)
stopifnot(sum(out$G.100 < 0)==9)
# reset units and species to run next examples
E.units("cal")
species(delete=TRUE)
