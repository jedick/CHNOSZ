\encoding{UTF-8}
\name{mosaic}
\alias{mosaic}
\title{Chemical Affinities with Changing Basis Species}
\description{
Calculate chemical affinities of formation reactions of species using basis species that change with the conditions.
}

\usage{
  mosaic(bases, bases2 = NULL, blend = TRUE, ...)
}

\arguments{
  \item{bases}{character, basis species to be changed in the calculation, or list, containing vectors for each group of changing basis species}
  \item{bases2}{character, second set of changing basis species}
  \item{blend}{logical, use relative abundances of basis species?}
  \item{...}{additional arguments to be passed to \code{\link{affinity}}}
}

\details{

\code{mosaic} can be used to calculate the affinities of formation of species when the relative abundances of the basis species listed in \code{bases} change over the range of conditions, due to e.g. ionization, complexation or redox reactions.
This is a way to \dQuote{speciate the basis species}.
For example, the speciation of sulfur (\samp{SO4-2}, \samp{HSO4-}, \samp{HS-} and \samp{H2S}) as a function of Eh and pH affects the formation affinities, and therefore relative stabilities of iron oxide and sulfide minerals.
Chemical activity diagrams constructed by assembling sub-diagrams corresponding to the predominant basis species can described as \dQuote{mosaic diagrams}.

The function calculates the affinities using all combination of basis species given as vector arguments to \code{bases} and \code{bases2}.
The first species listed in each group should be in the current basis definition, and all the basis species in each group should be related to the first basis species there (i.e. all share the same element).
A second, independent set of basis species can be provided in \code{bases2} (for example \samp{CO3-2}, \samp{HCO3-}, \samp{CO2}, if the first set of basis species are the sulfur-bearing ones listed above).
The arguments in \code{...} are passed to \code{affinity} to specify the variable conditions, such as temperature, pressure, and activities of other basis species.

If \code{blend} is TRUE (the default), the relative abundances of the basis species in each group are calculated using \code{\link{equilibrate}}, with the total activity taken from the corresponding basis species in the incoming \code{\link{basis}} definition.
Then, the function calculates overall affinities of the formation reactions of each species by combining reactions written using individual basis species in proportion to the relative abundances of the basis species.
If \code{blend} is FALSE, the function returns the affinities calculated using the single predominant basis species in \code{bases} at each condition.

A more flexible method of specifying multiple sets of basis species is now available.
Instead of using \code{bases} and \code{bases2}, supply a list for just the \code{bases} argument.
The list should contain any number of vectors specifying the groups of basis species.
All combinations of basis species in these groups are used for the calculations.
This overcomes the prior limitation of only having two changing groups of basis species.

}

\value{
A list containing \code{A.species} (affinities of formation of the species with changing basis species) and \code{A.bases} (affinities of formation of the basis species in terms of the first basis species), each having same structure as the list returned by \code{\link{affinity}}.
If \code{blend} is TRUE, the output also contains \code{E.bases} (the output of \code{\link{equilibrate}} for each group of basis species)
If \code{bases2} is provided, the list also contains \code{A.bases2} (affinities of formation of the second group of basis species).
}

\seealso{
\code{demo("mosaic")}, extending the example below by addition of carbonate species in \code{bases2}, and using thermodynamic data from Garrels and Christ, 1965.
The help page of \code{\link{solubility}} has an example combining \code{mosaic} with solubility calculations.
}

\examples{
\dontshow{reset()}# Fe-minerals and aqueous species in Fe-S-O-H system
# speciate SO4-2, HSO4-, HS-, H2S as a function of Eh and pH
# after Garrels and Christ, 1965 Figure 7.20
pH <- c(0, 14, 250)
Eh <- c(-1, 1, 250)
T <- 25
basis(c("FeO", "SO4-2", "H2O", "H+", "e-"))
basis("SO4-2", -6)
species(c("Fe+2", "Fe+3"), -6)
species(c("pyrrhotite", "pyrite", "hematite", "magnetite"))
# the basis species we'll swap through
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
# calculate affinities using the relative abundances of the basis species
# NOTE: set blend = FALSE for sharp transitions between the basis species
# (looks more like the diagram in GC65)
m1 <- mosaic(bases, pH = pH, Eh = Eh, T = T)
# make a diagram and add water stability lines
d <- diagram(m1$A.species, lwd = 2)
water.lines(d, col = "seagreen", lwd = 1.5)
# show lines for Fe(aq) = 10^-4 M
species(c("Fe+2", "Fe+3"), -4)
m2 <- mosaic(bases, pH = pH, Eh = Eh, T = T)
diagram(m2$A.species, add = TRUE, names = NULL)
title(main=paste("Iron oxides and sulfides in water, log(total S) = -6",
  "After Garrels and Christ, 1965", sep="\n"))
legend("bottomleft", c("log(act_Fe) = -4", "log(act_Fe) = -6"), lwd = c(2, 1), bty = "n")
# we could overlay the basis species predominance fields
#diagram(m1$A.bases, add=TRUE, col="blue", col.names="blue", lty=3)
}

\references{
  Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{http://www.worldcat.org/oclc/517586}
}

\concept{Extended workflow}
