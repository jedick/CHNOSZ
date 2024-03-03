\encoding{UTF-8}
\name{swap.basis}
\alias{basis.elements}
\alias{element.mu}
\alias{basis.logact}
\alias{swap.basis}
\alias{ibasis}
\title{Swap basis species}
\description{
  Swap the basis species defining a chemical system. One basis species is replaced by a new one with a different chemical formula.
}

\usage{
  swap.basis(species, species2, T = 25)
  basis.elements(basis = thermo()$basis)
  element.mu(basis = thermo()$basis, T = 25)
  basis.logact(emu, basis = thermo()$basis, T = 25)
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

The other functions are supporting functions: \code{basis.elements} returns the stoichiometric matrix for the current basis definition, \code{element.mu} calculates the chemical potentials of elements corresponding to the activities of the basis species, \code{basis.logact} does the inverse operation, and \code{ibasis} returns the index in the basis set for a given species index (in \code{\link{thermo}$OBIGT}), name or formula.
}

\seealso{ \code{\link{basis}}, and \code{\link{mosaic}} }

\examples{
\dontshow{reset()}
## Swapping basis species
# Start with a preset basis definition
b1 <- basis("CHNOS+")
# Swap H2(aq) for O2(gas)
b2 <- swap.basis("O2", "H2")
# Put oxygen back in
b3 <- swap.basis("H2", "oxygen")

# Interconversion of chemical potentials of elements and
# logarithms of activities of basis species at high temperature
basis("CHNOS+")
bl1 <- basis()$logact
emu <- element.mu(T = 100)
bl2 <- basis.logact(emu, T = 100)
# There's no difference
round(bl2 - bl1, 10)
}

\concept{Extended workflow}
