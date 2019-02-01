\encoding{UTF-8}
\name{mosaic}
\alias{mosaic}
\title{Chemical Affinities with Changing Basis Species}
\description{
Calculate chemical affinities of formation reactions of species using basis species that change with the conditions.
}

\usage{
  mosaic(bases, bases2 = NULL, blend = FALSE, mixing = FALSE, ...)
}

\arguments{
  \item{bases}{character, basis species to be changed in the calculation, or list, containing vectors for each group of changing basis species}
  \item{bases2}{character, second set of changing basis species}
  \item{blend}{logical, use relative abundances of basis species?}
  \item{mixing}{logical, include a term for the Gibbs energy of mixing?}
  \item{...}{additional arguments to be passed to \code{\link{affinity}}}
}

\details{

\code{mosaic} can be used to calculate the reaction affinities when the basis species listed in \code{bases} change in relative abundance over the range of conditions, due to e.g. ionization, complexation or redox reactions.
Chemical activity or predominance diagrams constructed by assembling sub-diagrams corresponding to the predominant basis species have sometimes been described as \dQuote{mosaic diagrams} in the literature.
This is a way to \dQuote{speciate the basis species}.
For example, the speciation of sulfur (\samp{SO4-2}, \samp{HSO4-}, \samp{HS-} and \samp{H2S}) as a function of Eh and pH affects the formation affinities, and therefore relative stabilities of iron oxide and sulfide minerals.

The function calculates the affinities using each basis species listed in \code{bases} in turn, changing them via \code{\link{swap.basis}}.
The first species listed in \code{bases} should be in the current basis definition.
The arguments in \code{...} are passed to \code{affinity} to specify the conditions.
If \code{blend} is FALSE (the default), the function returns the affinities calculated using the single predominant basis species in \code{bases} at each condition.
If \code{blend} is TRUE, the function combines the affinities of the formation reactions in proportion to the relative abundances of the basis species at each condition.
Additionally, if \code{mixing} is TRUE, a term is included to account for the Gibbs energy of mixing.
See the second example in \code{\link{solubility}} for a numerical test of the calculations using \code{blend} and \code{mixing}.

The basis species listed in \code{bases} should all be related to the first basis species there (i.e. all share the same element).
A second, independent set of basis species can be provided in \code{bases2} (for example \samp{CO3-2}, \samp{HCO3-}, \samp{CO2}, if the first set of basis species are the sulfur-bearing ones listed above).
The function then works recursively, by calling itself instead of \code{affinity}, so that the inner loop changes the basis species in \code{bases2}.
In this way, all possible combinations of the two sets of basis species are used in the calculation.

A more flexible method of specifying multiple sets of basis species is now available.
Instead of using \code{bases} and \code{bases2}, supply a list for just the \code{bases} argument.
The list should contain any number of vectors specifying the groups of basis species.
All combinations of basis species in these groups are used for the calculations.
This overcomes the prior limitation of only having two changing groups of basis species.

}

\value{
A list containing \code{A.species} (affinities of formation of the species with changing basis species) and \code{A.bases} (affinities of formation of the basis species in terms of the first basis species), each having same structure as the list returned by \code{\link{affinity}}.
If \code{bases2} is provided, the list also contains \code{A.bases2} (affinities of formation of the second set of basis species, in terms of the first one in that set).
}

\seealso{
\code{demo("mosaic")}, extending the example below by addition of carbonate species in \code{bases2}, with \code{blend} set to TRUE, and using thermodynamic data from Garrels and Christ, 1965.
}

\examples{
\dontshow{data(thermo)}# Fe-minerals and aqueous species in Fe-S-O-H system
# speciate SO4-2, HSO4-, HS-, H2S as a function of Eh and pH
# after Garrels and Christ, 1965 Figure 7.20
pH <- c(0, 14, 200)
Eh <- c(-1, 1, 200)
T <- 25
basis(c("FeO", "SO4-2", "H2O", "H+", "e-"))
basis("SO4-2", -6)
species(c("Fe+2", "Fe+3"), -6)
species(c("pyrrhotite", "pyrite", "hematite", "magnetite"))
# the basis species we'll swap through
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
# calculate affinities using the predominant basis species
m1 <- mosaic(bases, pH=pH, Eh=Eh, T=T)
# make a diagram and add water stability lines
d <- diagram(m1$A.species, lwd=2)
water.lines(d, col="seagreen", lwd=1.5)
# show lines for Fe(aq) = 10^-4 M
species(c("Fe+2", "Fe+3"), -4)
m2 <- mosaic(bases, pH=pH, Eh=Eh, T=T)
diagram(m2$A.species, add=TRUE, names=NULL)
title(main=paste("Iron oxides and sulfides in water, log(total S) = -6",
  "After Garrels and Christ, 1965", sep="\n"))
# we could overlay the basis species predominance fields
#diagram(m1$A.bases, add=TRUE, col="blue", col.names="blue", lty=3)
}

\references{
  Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{http://www.worldcat.org/oclc/517586}
}

\concept{Extended workflow}
