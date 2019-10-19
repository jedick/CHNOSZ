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
