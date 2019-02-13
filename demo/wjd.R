## carbon-containing compounds in prebiological atmospheres
## Dayhoff et al., 1964 (https://doi.org/10.1126/science.146.3650.1461)
#pdf("dayhoff.pdf", width=6, height=6)
# read formulas and Gibbs energies
file <- system.file("extdata/adds/DLEN67.csv", package="CHNOSZ")
dlen <- read.csv(file, as.is=TRUE, row.names=1)
# turn formulas into a stoichiometric matrix
A <- i2A(dlen$formula)
# assemble Gibbs energies/RT at 500 K
T <- 500  # K
R <- 1.9872  # gas constant, cal K^-1 mol^-1
G0.RT <- 1000 * dlen$G500 / R / T
# a function to minimize Gibbs energy for system with 
# given mole fraction of carbon (xC)
min.atmos <- function(xC) {
  # the bulk composition C:H:N:O
  B <- c(xC, 100-40-xC, xC, 40)
  # guess the initial composition
  Y <- guess(A, B)
  w <- wjd(A=A, G0.RT=G0.RT, Y=Y, P=1, imax=90, Gfrac=1e-14)
  if(!is.near.equil(w)) cat(paste("not near equilibrium for xC=", xC, "\n"))
  return(w)
}
# vary carbon content
xCs <- seq(8, 47, 1)
Xs <- sapply(xCs, function(xC) min.atmos(xC)$X)
# normalize the mole numbers to mole fractions
Xs <- t(t(Xs)/colSums(Xs))
plot(-10, 0, xlim=c(0, 55), ylim=c(-25, 1), xlab="mole percent C", ylab="log10 mole fraction")
for(i in 1:nrow(Xs)) lines(xCs, log10(Xs[i, ]))
text(48, log10(Xs[, length(xCs)]), dlen$formula, adj=0)
text(35, log10(Xs[, 27]) + 0.5, dlen$formula, adj=0)
text(7, log10(Xs[, 1]), dlen$formula, adj=1)
title(main="Prebiological atmospheres (Dayhoff et al., 1964)")
#dev.off()

## run.wjd with proteins: cell periphery of yeast
# get the proteins in the requested location
y <- yeastgfp("cell.periphery")
# get the amino acid compositions of the proteins
aa <- yeast.aa(y$protein)
# don't use those with NA abundance or sequence
ina <- is.na(y$abundance) | is.na(aa$chains)
aa <- aa[!ina, ]
# let's try normalizing the proteins to single residues
# columns 6:25 are the actual amino acid counts
aa.625 <- aa[, 6:25]
aa[, 6:25] <- aa.625 / rowSums(aa.625)
# add proteins to thermo$protein
add.protein(aa)
# add proteins to thermo$obigt
iobigt <- info(paste(aa$protein, aa$organism, sep="_"))
# use equal initial abundances, with total equal to yeastGFP abundances
Y <- rep(mean(y$abundance[!ina]), length(y$abundance[!ina]))
# run the Gibbs energy minimization
w <- run.wjd(iobigt, Y=Y, imax=100)
# make a log-log plot
plot(log10(y$abundance[!ina]), log10(w$X), xlim=c(1.5, 5), ylim=c(1.5, 5),
  xlab="log10(abundance) reported in YeastGFP study",
  ylab="log10(abundance) calculated using Gibbs energy minimization")
# get the element potentials (tolerating "close enough" to equilibrium)
emu <- equil.potentials(w, tol=1e7)
# then the logarithms of activities of the basis species
basis("CHNOS")
bl <- basis.logact(emu)
# make a title and legend
title(main="Relative abundances of proteins: yeast cell periphery")
basis(names(bl), bl)
legend("topleft", describe.basis(digits=2))
