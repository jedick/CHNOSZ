## Demo for revisit(): CV of equilibrium activities of proteins in
## subcellular compartment of S. cerevisiae

# get the proteins in the requested location
loc <- "cell.periphery"
y <- yeastgfp(loc)
# get the amino acid compositions of the proteins
aa <- yeast.aa(y$protein)
# don't use those with NA abundance or sequence
ina <- is.na(y$abundance) | is.na(aa$chains)
aa <- aa[!ina, ]
# add these proteins to CHNOSZ's inventory
ip <- add.protein(aa)
# set up a the chemical system
basis("CHNOS+")
# calculate affinities of formation in logfO2 space
a <- affinity(O2=c(-85, -60), iprotein=ip)
# show the equilibrium activities
opar <- par(mfrow=c(2, 2))
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
diagram(e, names=NULL)
# make a title
expr <- as.expression(substitute(x~loc~"proteins in"~
  italic("S. cerevisiae"), list(x=length(ip), loc=loc)))
mtitle(c("Equilibrium activities of", expr))
# show the coefficient of variation
revisit(e, "CV")
mtitle(c("CV of equilibrium activities of", expr))
# calculate affinities in logfO2-logaH2O space
a <- affinity(O2=c(-85, -65), H2O=c(-5, 5), iprotein=ip)
# show the predominances
diagram(a, normalize=TRUE, fill="heat")
# calculate the equilibrium activities
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
# show the coefficient of variation
r <- revisit(e, "CV")
mtitle(c("CV of equilibrium activities of", expr))
par(opar)
