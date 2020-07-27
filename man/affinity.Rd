\encoding{UTF-8}
\name{affinity}
\alias{affinity}
\title{Chemical Affinities of Formation Reactions}
\description{
Calculate the chemical affinities of formation reactions of species.
}

\usage{
  affinity(..., property = NULL, sout = NULL, exceed.Ttr = FALSE,
    exceed.rhomin = FALSE, return.buffer = FALSE, return.sout = FALSE,
    balance = "PBB", iprotein = NULL, loga.protein = -3, transect = NULL)
}

\arguments{
  \item{...}{numeric, zero or more named arguments, used to identify the variables of interest in the calculations. For argument recall, pass the output from a previous calculation of \code{affinity} as an unnamed first argument.}
  \item{property}{character, the property to be calculated. Default is \samp{A}, for chemical affinity of formation reactions of species of interest}
  \item{sout}{list, output from \code{\link{subcrt}}}
  \item{exceed.Ttr}{logical, allow \code{\link{subcrt}} to compute properties for phases beyond their transition temperature?}
  \item{exceed.rhomin}{logical, allow \code{\link{subcrt}} to compute properties of species in the HKF model below 0.35 g cm\S{-3}?}
  \item{return.buffer}{logical. If \code{TRUE}, and a \code{\link{buffer}} has been associated with one or more basis species in the system, return the values of the activities of the basis species calculated using the buffer. Default is \code{FALSE}.}
  \item{return.sout}{logical, return only the values calculated with \code{\link{subcrt}}?}
  \item{balance}{character. This argument is used to identify a conserved basis species (or \samp{PBB}) in a chemical activity buffer. Default is \samp{PBB}.}
  \item{iprotein}{numeric, indices of proteins in \code{\link{thermo}$protein} for which to calculate properties}
  \item{loga.protein}{numeric, logarithms of activities of proteins identified in \code{iprotein}}
  \item{transect}{logical, force a transect calculation, even for three or fewer values of the variables?}
}

\details{
\code{affinity} calculates the chemical affinities of reactions to form the species of interest from the basis species.
The equation used to calculate chemical affinity (\emph{A}), written for base-10 (decimal) logarithms, is \eqn{A/(2.303RT)}=\eqn{\log (K/Q)}{log(K/Q)}, where \eqn{K} is the equilibrium constant of the reaction, \eqn{Q} is the activity product of the species in the reaction, and 2.303 is the conversion factor from natural to decimal logarithms.
The calculation of chemical affinities relies on the current definitions of the \code{\link{basis}} species and \code{\link{species}} of interest.
Calculations are possible at single values of temperature, pressure, ionic strength and chemical activities of the basis species, or as a function of one or more of these variables. 

The argument \code{property} can be changed to calculate other thermodynamic properties of formation reactions. 
Valid properties are \samp{A} or NULL for chemical affinity, \samp{logK} or \samp{logQ} for logarithm of equilibrium constant and reaction activity product, or any of the properties available in \code{\link{subcrt}} except for \samp{rho}.
The properties returned are those of the formation reactions of the species of interest from the basis species.
It is also possible to calculate the properties of the species of interest themselves (not their formation reactions) by setting the \code{property} to \samp{G.species}, \samp{Cp.species}, etc.
Except for \samp{A}, the properties of proteins or their reactions calculated in this manner are restricted to nonionized proteins.

Zero, one, or more leading arguments to the function identify which of the chemical activities of basis species, temperature, pressure and/or ionic strength to vary.
The names of each of these arguments may be the formula of any of the basis species of the system, or \samp{T}, \samp{P}, \samp{pe}, \samp{pH}, \samp{Eh}, or \samp{IS} (but names may not be repeated).
The names of charged basis species such as \samp{K+} and \samp{SO4-2} should be quoted when used as arguments.
The value of each argument is of the form \code{c(min, max)} or \code{c(min, max, res)} where \code{min} and \code{max} refer to the minimimum and maximum values of variable identified by the name of the argument, and \code{res} is the resolution, or number of points along which to do the calculations; \code{res} is assigned a default value of 128 if it is missing.
For any arguments that refer to basis species, the numerical values are the logarithms of activity (or fugacity for gases) of that basis species.

If \samp{T}, \samp{P}, and/or \samp{IS} are not among the \code{var}s, their constant values can be supplied in \code{T}, \code{P}, or \code{IS} (in mol kg\eqn{^{-1}}{^-1}).
The units of \samp{T} and \samp{P} are those set by \code{\link{T.units}} and \code{\link{P.units}} (on program start-up these are \degC and bar, respectively). 
\code{sout}, if provided, replaces the call to \code{\link{subcrt}}, which can greatly speed up the calculations if this intermediate result is stored by other functions.
\code{exceed.Ttr} is passed to \code{\link{subcrt}} so that the properties of mineral phases beyond their transition temperatures can optionally be calculated.

If one or more buffers are assigned to the definition of \code{\link{basis}} species, the logarithms of activities of these basis species are taken from the buffer (see \code{\link{buffer}}).

The \code{iprotein} and \code{loga.protein} arguments can be used to compute the chemical affinities of formation reactions of proteins that are not in the current \code{\link{species}} definition.
\code{iprotein} contains the indices (rownumbers) of desired proteins in \code{\link{thermo}$protein}.
This uses some optimizations to calculate the properties of many proteins in a fraction of the time it would take to calculate them individually.

When the length(s) of the variables is(are) greater than 3, the function enters the \samp{transect} mode of operation.
In this mode of operation, instead of performing the calculations on an \eqn{n}{n}-dimensional grid, the affinities are calculated on a transect of changing T, P, and/or chemical activity of basis species.

Argument recall is invoked by passing a previous result of \code{affinity} as the first argument.
The function then calls itself using the settings from the previous calculation, with additions or modifications indicated by the remaining arguments in the current function call.
}

\value{
A list, elements of which are \code{fun} the name of the function (\samp{affinity}), \code{args} all of the arguments except for \samp{sout} (these are used for argument recall), \code{sout} output from \code{\link{subcrt}}, \code{property} name of the calculated property (\samp{A} for chemical affinity), \code{basis} and \code{species} definition of basis species and species of interest in effect at runtime, \code{T} and \code{P} temperature and pressure, in the system units of Kelvin and bar, set to \code{numeric()} (length=0) if either one is a variable, \code{vars} the names of the variables, \code{vals} the values of the variables (a list, one element for each variable), \code{values} the result of the calculation (a list, one element for each species, with names taken from the species index in \code{\link{thermo}$OBIGT}).
The elements of the lists in \code{vals} and \code{values} are arrays of \eqn{n}{n} dimensions, where \eqn{n}{n} is the number of variables.
The values of chemical affinity of formation reactions of the species are returned in dimensionless units (for use with decimal logarithms, i.e., \emph{A}/\eqn{2.303RT}).

Names other than \samp{T} or \samp{P} in \code{vars} generally refer to basis species, and the corresponding \code{vals} are the logarithms of activity or fugacity.
However, if one or more of \code{pe}, \code{Eh} or \code{pH} is among the variables of interest, \code{vals} holds the values of the those variables as indicated.
}

\seealso{
\code{\link{ionize.aa}}, activated if proteins are among the species of interest and \samp{H+} is in the basis.
\code{\link{equilibrate}} for using the results of \code{affinity} to calculate equilibrium activities of species, and \code{\link{diagram}} to plot the results.
\code{\link{demo}("saturation")} for an example using the argument recall feature.
}

\examples{
\dontshow{reset()}## set up a system and calculate
## chemical affinities of formation reactions
basis(c("SiO2", "MgO", "H2O", "O2"), c(-5, -5, 0, 999))
species(c("quartz","enstatite","forsterite"))
# chemical affinities (A/2.303RT) at 25 deg C and 1 bar
affinity()
# at higher temperature and pressure
affinity(T=500, P=2000)
# at 25 temperatures and pressures,
# some are in the low-density region so we suppress warnings
suppressWarnings(affinity(T=c(500, 1000, 5), P=c(1000, 5000, 5)))
# equilibrium constants of formation reactions
affinity(property="logK")
# standard molal Gibbs energies of species,
# user units (default: cal/mol)
affinity(property="G.species")
# standard molal Gibbs energies of reactions
affinity(property="G")
# a T,P-transect
# (fluid pressure from Helgeson et al., 2009 Fig. 7)
affinity(T=c(25, 110, 115, 215), P=c(11, 335, 500, 1450))
}

\references{
Helgeson, H. C., Richard, L, McKenzie, W. F., Norton, D. L. and Schmitt, A. (2009) A chemical and thermodynamic model of oil generation in hydrocarbon source rocks. \emph{Geochim. Cosmochim. Acta} \bold{73}, 594--695. \url{https://doi.org/10.1016/j.gca.2008.03.004}
}

\concept{Main workflow}
