\encoding{UTF-8}
\name{mosaic}
\alias{mosaic}
\title{Chemical Affinities with Changing Basis Species}
\description{
Calculate chemical affinities of formation reactions of species using basis species that change with the conditions.
}

\usage{
  mosaic(bases, blend = TRUE, stable = list(), loga_aq = NULL, ...)
}

\arguments{
  \item{bases}{character, basis species to be changed in the calculation, or list, each element of which defines an independent group of changing basis species}
  \item{blend}{logical, use relative abundances of basis species?}
  \item{stable}{list, previously determined stable species}
  \item{loga_aq}{numeric, activities of aqueous species (overrides current values in \code{\link{basis}})}
  \item{...}{additional arguments to be passed to \code{\link{affinity}}}
}

\details{

\code{mosaic} calculates the affinities of formation of species when the relative abundances of the basis species listed in \code{bases} change over the range of conditions, due to e.g. ionization, complexation or redox reactions.
This is a way to \dQuote{speciate the basis species}.
For example, the speciation of sulfur (\samp{SO4-2}, \samp{HSO4-}, \samp{HS-} and \samp{H2S}) as a function of Eh and pH affects the formation affinities, and therefore relative stabilities of iron oxide and sulfide minerals.
Chemical activity diagrams constructed by assembling sub-diagrams corresponding to the predominant (i.e. most stable) basis species can described as \dQuote{mosaic diagrams}.

The function calculates the affinities using all combination of basis species given as a vector argument to \code{bases}.
Or, a list can be provided contain any number of vectors that specify different groups of basis species.
All combinations of basis species in these groups are used for the calculations.

The first species listed in each group should be in the current basis definition, and all the basis species in each group should be related to the first basis species there (i.e. all share the same element).
The arguments in \code{...} are passed to \code{affinity} to specify the variable conditions, such as temperature, pressure, and activities of other basis species.

\code{blend} can be a single logical value or multiple values, in order to control the calculations for individual groups of basis species.
If \code{blend} is TRUE (the default), the relative abundances of the basis species in each group are calculated using \code{\link{equilibrate}}, with the total activity taken from the corresponding basis species in the incoming \code{\link{basis}} definition.
Then, the function calculates overall affinities of the formation reactions of each species by combining reactions written using individual basis species in proportion to the relative abundances of the basis species.

If \code{blend} is FALSE, the function returns the affinities calculated using the single predominant basis species in \code{bases} at each condition.
This is appropriate when minerals, rather than aqueous species, are used as the changing basis species.
Note, however, that \code{mosaic} is not internally recursive: the stabilities of one group of basis species (e.g. minerals) are not affected by changes in another group (e.g. aqueous species).

(\dQuote{Stacked mosaic diagrams}) are useful for making diagrams for multi-metal systems.
By using the stable minerals in one calculation as the new basis species in a subsequent calculation, a series of stacked \code{mosaic} diagrams with increasing complexity can be made.
Specifically, this is done by supplying previously calculated stabilities (from the \code{predominant} element of the output of \code{\link{diagram}}) as an element of the list in the \code{stable} argument whose position corresponds to the appropriate group of basis species.
Note that a value in any position of the \code{stable} list forces \code{blend = FALSE} for the corresponding group of basis species, so there is no need to explicity change the \code{blend} argument.

The activities of mosaiced basis species in each group are taken from the current \code{\link{basis}} definition.
Generally it makes sense to set the activity of minerals to 1 (logact = 0) and the activity of aqueous species to some smaller value.
For mosaic stacking calculations where the mosaiced basis species include both minerals and aqueous species, the \code{loga_aq} argument specifies the activity of aqueous species to be used \emph{in each group}.
That is, there should be one value of \code{loga_aq} for each group of basis species; use NA to indicate that the activity comes from the current \code{\link{basis}} definition.
See the Mosaic Stacking 2 section of the vignette \viglink{multi-metal} for an example.
}

\value{
A list containing \code{A.species} (affinities of formation of the species with changing basis species) and \code{A.bases} (affinities of formation of the basis species in terms of the first basis species), each having same structure as the list returned by \code{\link{affinity}}.
If \code{bases} is a list, \code{A.bases} is also a list, each element of which corresponds to a group of basis species in \code{bases}.
If \code{blend} is TRUE, the output also contains \code{E.bases} (the output of \code{\link{equilibrate}} for each group of basis species)
}

\seealso{
\code{demo("mosaic")}, which extends the example below with carbonate species in order to plot a siderite field.
To calculate mineral solubilities with mosaic calculations that account for ligand speciation, use \code{bases} as the first argument to \code{\link{solubility}}.
\code{\link{stack_mosaic}} implements calculations for bimetallic systems.
}

\examples{
\dontshow{reset()}# Fe-minerals and aqueous species in Fe-S-O-H system
# Speciate SO4-2, HSO4-, HS-, H2S as a function of Eh and pH
# After Garrels and Christ, 1965 Figure 7.20
pH <- c(0, 14)
Eh <- c(-1, 1)
T <- 25
basis(c("FeO", "SO4-2", "H2O", "H+", "e-"))
basis("SO4-2", -6)
species(c("Fe+2", "Fe+3"), -6)
species(c("pyrrhotite", "pyrite", "hematite", "magnetite"), add = TRUE)
# The basis species we'll swap through
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
# Calculate affinities using the relative abundances of the basis species
# NOTE: set blend = FALSE for sharp transitions between the basis species
# (looks more like the diagram in GC65)
m1 <- mosaic(bases, pH = pH, Eh = Eh, T = T)
# Make a diagram and add water stability lines
d <- diagram(m1$A.species, lwd = 2)
water.lines(d, col = "seagreen", lwd = 1.5)
# Show lines for Fe(aq) = 10^-4 M
species(c("Fe+2", "Fe+3"), -4)
m2 <- mosaic(bases, pH = pH, Eh = Eh, T = T)
diagram(m2$A.species, add = TRUE, names = FALSE)
title(main=paste("Iron oxides and sulfides in water, log(total S) = -6",
  "After Garrels and Christ, 1965", sep="\n"))
legend("bottomleft", c("log(act_Fe) = -4", "log(act_Fe) = -6"), lwd = c(2, 1), bty = "n")
# We could overlay the predominance fields for the basis species
#diagram(m1$A.bases, add = TRUE, col = "blue", col.names = "blue", lty = 3)
}

\references{
  Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{https://www.worldcat.org/oclc/517586}
}

\concept{Extended workflow}
