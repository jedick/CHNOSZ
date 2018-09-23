\encoding{UTF-8}
\name{swap.basis}
\alias{basis.elements}
\alias{element.mu}
\alias{basis.logact}
\alias{swap.basis}
\alias{ibasis}
\title{Swap Basis Species}
\description{
  Swap the basis species defining a chemical system. One basis species is replaced by a new one with a different chemical formula.
}

\usage{
  swap.basis(species, species2, T = 25)
  basis.elements(basis = get("thermo")$basis)
  element.mu(basis = get("thermo")$basis, T = 25)
  basis.logact(emu, basis = get("thermo")$basis, T = 25)
  ibasis(species)
}

\arguments{
  \item{basis}{dataframe, a basis definition}
  \item{T}{numeric, temperature in Kelvin}
  \item{emu}{numeric, chemical potentials of elements}
  \item{species}{character, names or formulas of species, or numeric, indices of species}
  \item{species2}{character or numeric, a species to swap in to the basis definition}
}

\details{
\code{swap.basis} allows to change the basis definition by swapping out a basis species for a new one.
Specify the names or formulas of the old and replacement basis species in the first argument.
When the basis definition is changed, any species of interest that were present are deleted, unless the new basis definition has exactly the same elements as before.
In that case, the species are kept; also, the activities of the new basis species are set in order to maintain the chemical potentials of the elements at \code{T} \degC and 1 bar.

The other functions are supporting functions: \code{basis.elements} returns the stoichiometric matrix for the current basis definition, \code{element.mu} calculates the chemical potentials of elements corresponding to the activities of the basis species, \code{basis.logact} does the inverse operation, and \code{ibasis} returns the index in the basis set for a given species index (in \code{\link{thermo}$obigt}), name or formula.
}

\seealso{ \code{\link{basis}}, and \code{\link{mosaic}} }

\examples{
\dontshow{data(thermo)}
## swapping basis species
# start with a preset basis definition
b1 <- basis("CHNOS+")
# swap H2(aq) for O2(gas)
(b2 <- swap.basis("O2", "H2"))
# the logarithm of activity calculated for H2
# is equal to the one calculated from the equilibrium constant
# for H2O = H2 + 0.5O2
logK <- subcrt(c("oxygen","H2","H2O"), c(-0.5,-1,1), T=25)$out$logK
# the equilibrium value of logaH2 
# (for logaH2O = 0 and logfO2 = -80)
(logaH2 <- -logK + 40)
stopifnot(all.equal(logaH2, b2$logact[5]))
# put O2 back in
b3 <- swap.basis("H2", "oxygen")
# we have returned to starting point
stopifnot(all.equal(b1$logact, b3$logact))

## demonstrating the interconvertibility between 
## chemical potentials of elements and logarithms 
## of activities of basis species at high temperature
#basis("CHNOS+")
#bl1 <- basis()$logact
#emu <- element.mu(T=100)
#bl2 <- basis.logact(emu, T=100)
## note that basis.logact produces a named array
#stopifnot(all.equal(bl1, as.numeric(bl2)))

## swapping basis species while species are defined
## and using numeric species indices
basis("MgCHNOPS+") 
# load some Mg-ATP species
species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"))
# swap in CO2(g) for CO2(aq)
iCO2g <- info("CO2", "gas")
swap.basis("CO2", iCO2g)
a1 <- affinity()
# swap in CH4(g) for CO2(g)
iCH4g <- info("CH4", "gas")
swap.basis(iCO2g, iCH4g)
a2 <- affinity()
# the equilibrium fugacity of CH4 is *very* low
# swap in CO2(aq) for CH4(g)
iCO2a <- info("CO2", "aq")
swap.basis(iCH4g, iCO2a)
a3 <- affinity()
# swapping the basis species didn't affect the affinities
# of the formation reactions of the species, since
# the chemical potentials of the elements were unchanged
stopifnot(all.equal(a1$values, a2$values))
stopifnot(all.equal(a1$values, a3$values))
}

\concept{Extended workflow}
