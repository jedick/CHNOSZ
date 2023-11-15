\encoding{UTF-8}
\name{subcrt}
\alias{subcrt}
\title{Properties of Species and Reactions}
\description{
  Calculate the standard molal thermodynamic properties of one or more species or a reaction between species as a function of temperature and pressure.
}

\usage{
  subcrt(species, coeff = 1, state = NULL,
    property = c("logK","G","H","S","V","Cp"),
    T = seq(273.15,623.15,25), P = "Psat", grid = NULL, 
    convert = TRUE, exceed.Ttr = FALSE, exceed.rhomin = FALSE,
    logact = NULL, autobalance = TRUE, use.polymorphs = TRUE, IS = 0)
}

\arguments{
  \item{species}{character, name or formula of species, or numeric, rownumber of species in \code{thermo()$OBIGT}}
  \item{coeff}{numeric, reaction coefficients on species}
  \item{state}{character, state(s) of species}
  \item{property}{character, property(s) to calculate}
  \item{T}{numeric, temperature(s) of the calculation}
  \item{P}{numeric, pressure(s) of the calculation, or character, \samp{Psat}}
  \item{grid}{character, type of \code{P}\eqn{\times}{x}\code{T} grid to produce (NULL, the default, means no gridding)}
  \item{exceed.Ttr}{logical, calculate Gibbs energies of mineral phases and other species beyond their transition temperatures?}
  \item{exceed.rhomin}{logical, return properties of species in the HKF model below 0.35 g cm\S{-3}?}
  \item{logact}{numeric, logarithms of activities of species in reaction}
  \item{convert}{logical, are units of T, P, and energy settable by the user (default) (see \code{\link{T.units}})?}
  \item{autobalance}{logical, attempt to automatically balance reaction with basis species?}
  \item{use.polymorphs}{logical, automatically identify available polymorphs in OBIGT and use the stable one at each value of \code{T}?}
  \item{IS}{numeric, ionic strength(s) at which to calculate adjusted molal properties, mol kg\eqn{^{-1}}{^-1}}
}

\details{
\code{subcrt} calculates the standard molal thermodynamic properties of species and reactions as a function of temperature and pressure.
For each of the \code{species} (a formula or name), optionally identified in a given \code{state}, the standard molal thermodynamic properties and equations-of-state parameters are retrieved via \code{\link{info}} (except for \H2O liquid).
The standard molal properties of the species are computed using the thermodynamic \code{model} given for each species (see \code{\link{thermo}}).
This function also calculates the thermodynamic properties of \emph{reactions} by summing those of the respective species.
This functionality is inspired the \acronym{SUPCRT92} package (Johnson et al., 1992).

\code{T} and \code{P} denote the temperature and pressure for the calculations.
The only valid non-numeric value is \samp{Psat} for \code{P}, which is the default (see \code{\link{water}}).
For calculations below 273.16 K, \code{P} set to 1, as \Psat is not defined at subzero (\degC) temperatures.
At temperatures above the critical point of water, \code{P} must be set to a numeric value; unless \code{exceed.rhomin} is TRUE, this should correspond to a fluid density \ge 0.35 g cm\S{-3}.

Argument \code{\link{grid}} if present can be one of \code{T} or \code{P} to perform the computation of a \code{T}\eqn{\times}{x}\code{P} or \code{P}\eqn{\times}{x}\code{T} grid.
The \code{property}s to be calculated can be one or more of those shown below:

 \tabular{lll}{
         \code{rho}     \tab Density of water                      \tab g cm\eqn{^{-3}}{^-3} \cr
         \code{logK}    \tab Logarithm of equilibrium constant     \tab dimensionless \cr
         \code{G}       \tab Gibbs energy                          \tab (cal | J) mol\eqn{^{-1}}{^-1} \cr
         \code{H}       \tab Enthalpy                              \tab (cal | J) mol\eqn{^{-1}}{^-1} \cr
         \code{S}       \tab Entropy                               \tab (cal | J) K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \cr
         \code{V}       \tab Volume                                \tab cm\eqn{^3} mol\eqn{^{-1}}{^-1} \cr
         \code{Cp}      \tab Heat capacity                         \tab (cal | J) K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \cr
         \code{E}       \tab Exapansibility                        \tab cm\eqn{^3} K\eqn{^{-1}}{^-1} \cr
         \code{kT}      \tab Isothermal compressibility            \tab cm\eqn{^3} bar\eqn{^{-1}}{^-1} \cr
 }

If \code{convert} is \code{TRUE} (the default), the input values of \code{T} and \code{P} are interpreted to have the units given by \code{\link{T.units}} and \code{\link{P.units}} (default: \degC and bar), and the output values of \code{G}, \code{H}, \code{S} and \code{Cp} are based on the units given in \code{\link{E.units}} (default: Joules).
If \code{convert} is \code{FALSE}, the user units (\code{\link{T.units}}, \code{\link{P.units}}, and \code{\link{E.units}}) are ignored, and \code{T} and \code{P} are taken to be in Kelvin and bar, and the returned values of \code{G}, \code{H}, \code{S} and \code{Cp} are in Joules.

A chemical reaction is defined if \code{coeff} is given.
In this mode the standard molal properties of species are summed according to the stoichiometric \code{coeff}icients, where negative values denote reactants.
An unbalanced reaction is signalled if the amount of any element on the reactant and product sides differs by more than 1e-7; in this case, \code{subcrt} prints the missing composition needed to balance the reaction and produces a warning but computes a result anyway.
Alternatively, if \code{autobalance} is \code{TRUE}, the \code{\link{basis}} species of a system were previously defined, and all elements in the reaction are represented by the basis species, an unbalanced reaction given in the arguments to \code{subcrt} will be balanced automatically.
The auto balancing doesn't change the reaction coefficients of any species in the reaction that are not among the basis species.

If \code{logact} is provided, the chemical affinities of reactions are calculated.
\code{logact} indicates the logarithms of activities (fugacities for gases) of species in the reaction; if there are fewer values of \code{logact} than number of species those values are repeated as necessary.
If the reaction was unbalanced to start, the logarithms of activities of any basis species added to the reaction are taken from the current definition of the \code{\link{basis}} species.
Columns appended to the output are \code{logQ} for the log10 of the activity product of the reaction, and \code{A} for the chemical affinity, in the units set by \code{\link{E.units}}.
Note that \code{\link{affinity}} provides related functionality but is geared toward the properties of formation reactions of species from the basis species and can be performed in more dimensions.
Calculations of chemical affinity in \code{subcrt} can be performed for any reaction of interest; however, they are currently limited to constant values of the logarithms of activities of species in the reactions, and hence of \code{logQ}, across the computational range.

If \code{IS} is set to a single value other than zero, \code{\link{nonideal}} is used to calculate the adjusted properties (\code{G}, \code{H}, \code{S} and \code{Cp}) of charged aqueous species at the given ionic strength.
To perform calculations at a single \code{P} and \code{T} and for multiple values of ionic strength, supply these values in \code{IS}.
Calculations can also be performed on a \code{P}-\code{IS}, \code{T}-\code{IS} or \code{P,T}-\code{IS} grid.
Values of \code{logK} of reactions calculated for \code{IS} not equal to zero are consistent with the adjusted Gibbs energies of the charged aqueous species.

If \code{thermo()$opt$varP} is \code{TRUE}, standard Gibbs energies of gases will be converted from a standard state at 1 bar (as used in SUPCRT) to a variable pressure standard state (see chapter 12 in Anderson and Crerar, 1993).
This is useful for constructing e.g. boiling curves for organic compounds.
}

\section{Warning}{
Although \acronym{SUPCRT92} prohibits calculations above 350 \degC at \Psat (\dQuote{beyond range of applicability of aqueous species equations}), CHNOSZ does not impose this limitation, and allows calculations up to the critical temperature (373.917 \degC) at \Psat.
Interpret calculations between 350 \degC and the critical temperature at \Psat at your own risk.
The discontinuity in the value of \logK at \Psat that is apparent in \code{demo("NaCl")} demonstrates one unexpected result.

\code{NA}s are produced for calculations at \samp{Psat} when the temperature exceeds the critical temperature of \H2O.
In addition, properties of species using the revised HKF equations are set to \code{NA} wherever the density of \H2O < 0.35 g/cm\S{3} (threshold just above the critical isochore; Johnson et al., 1992).
Both of these situations produce warnings, which are stored in the \samp{warnings} element of the return value.

\code{NA}s are also output if the T, P conditions are otherwise beyond the capabilities of the water equations of state derived from SUPCRT92 (H2O92D.f), but the messages about this are produced by \code{\link{water.SUPCRT92}} rather than \code{subcrt}.
}

\section{Limitations}{
Note that \code{E} and \code{kT} can only be calculated for aqueous species and only if the option (\code{\link{thermo}()$opt$water}) for calculations of properties using \code{water} is set to \code{IAPWS}.
On the other hand, if the \code{water} option is \samp{SUPCRT} (the default), \code{E} and \code{kT} can be calculated for water but not for aqueous species.
(This is not an inherent limitation in either formulation, but it is just a matter of implementation.)
}

\section{Note on polymorphic transitions}{
Minerals with polymorphic transitions (denoted in OBIGT by having states \samp{cr} (lowest-\T phase), \samp{cr2}, etc.) can be specified by name with \samp{cr} for the \code{state} or by using a numeric species index for the lowest-\T polymorph.
If \code{use.polymorphs} is TRUE, \code{subcrt} uses the transition temperatures calculated from those at P = 1 bar given in OBIGT together with functions of the entropies and volumes of transitions (see \code{\link{dPdTtr}}) to determine the stable polymorph at each grid point and uses the properties of that polymorph in the output.
A \code{polymorph} column is added to the output to indicate the stable polymorph at each \T-\P condition.
If \code{exceed.Ttr} is \code{FALSE} (the default), output values of Gibbs energy are assigned NA beyond the transition temperature of the highest-temperature polymorph.
Set \code{exceed.Ttr} to \code{TRUE} to identify the stable polymorphs by comparing their extrapolated Gibbs energies instead of the tabulated transition temperatures.
This is generally not advised, as extrapolated Gibbs energies may not reliably determine the stable polymorph at extreme temperatures.
}

\value{
A list of length two or three.
If the properties of a reaction were calculated, the first element of the list (named \samp{reaction}) contains a dataframe with the reaction parameters; the second element, named \samp{out}, is a dataframe containing the calculated properties.
Otherwise, the properties of species (not reactions) are returned: the first element, named \samp{species}, contains a dataframe with the species identification; the second element, named \samp{out}, is itself a list, each element of which is a dataframe of properties for a given species.
If minerals with polymorphic transitions are present, a third element (a dataframe) in the list indicates for all such minerals the stable phase at each grid point.
}

\seealso{
\code{\link{info}} can be used to find species in the thermodynamic database.
\code{\link{makeup}} is used by \code{subcrt} for parsing formulas to check mass balance of reactions.
\code{demo("ORP")} and \code{\link{nonideal}} for examples using the \code{IS} argument.
}

\examples{
\dontshow{reset()}
## Properties of species
subcrt("water")
# Change temperature
subcrt("water", T = seq(0, 100, 20))
# Change temperature and pressure 	
T <- seq(500, 1000, 100)
P <- seq(5000, 10000, 1000)
subcrt("water", T = T, P = P)
# Temperature-pressure grid
subcrt("water", T = c(500, 1000), P = c(5000, 10000), grid = "P")

## Properties of reactions
subcrt(c("glucose", "ethanol", "CO2"), c(-1, 2, 2), T = 25)
# Use CO2(gas) (or just change "CO2" to "carbon dioxide")
subcrt(c("glucose", "ethanol", "CO2"), c(-1, 2, 2), c("aq", "aq", "gas"), T = 25)

## Automatically balance reactions
# First define the basis species
basis(c("CO2", "H2O", "NH3", "H2S", "O2"))
# Auto-balance adds the required amount of H2O and O2
subcrt(c("ethanol", "glucose"), c(-3, 1), T = 37)
# An example with H+
basis(c("H2O", "H2S", "O2", "H+"))
subcrt(c("HS-", "SO4-2"), c(-1, 1), T = 100)

## Mineral polymorphs
# Properties of the stable polymorph at each temperature
subcrt("pyrrhotite")
# Reactions automatically use stable polymorph
subcrt(c("pyrite", "pyrrhotite", "H2O", "H2S", "O2"), c(-1, 1, -1, 1, 0.5))
# Extrapolated properties of the lowest-T polymorph (metastable at higher temperatures)
subcrt(c("pyrrhotite"), use.polymorphs = FALSE, exceed.Ttr = TRUE)

## Messages about problems with the calculation
# Above the T, P limits for the H2O equations of state
subcrt("alanine", T = c(2250, 2251), P = c(30000, 30001), grid = "T")
# Psat is not defined above the critical point
# (suppressWarnings is used so that checks of examples don't raise warnings)
suppressWarnings(subcrt("alanine", T = seq(0, 5000, by = 1000)))

## Minerals with polymorphic transitions
# Compare calculated values of heat capacity of iron with
# values from Robie and Hemingway, 1995
T.units("K")
# We set pressure here otherwise subcrt uses Psat (saturation 
# vapor pressure of H2O above 100 degrees C) which can't be 
# calculated above the critical point of H2O (~647 K)
s <- subcrt("Fe", T=seq(300, 1800, 20), P=1)
plot(s$out[[1]]$T, s$out[[1]]$Cp, type="l",
  xlab=axis.label("T"), ylab=axis.label("Cp"))
# Add points from RH95
RH95 <- read.csv(system.file("extdata/cpetc/RH95.csv", package="CHNOSZ"))
points(RH95[,1], RH95[,2])
title(main=paste("Heat capacity of Fe(cr)\n",
  "(points - Robie and Hemingway, 1995)"))
# Reset the units to default values
T.units("C")

## Subzero (degrees C) calculations
# Uncomment the following to try IAPWS95 instead of SUPCRT92
#water("IAPWS95")
# The limit for H2O92D.f (from SUPCRT92) is currently -20 deg C
# but we go to -30 knowing properties will become NA
sb <- subcrt(c("H2O", "Na+"), T = seq(-30, 10), P = 1)$out
# Start plot with extra room on right
opar <- par(mar=c(5, 4, 4, 4))
# Plot Delta G
plot(sb$water$T, sb$water$G, ylim = c(-264000, -234000),
  xlab = axis.label("T"), ylab = axis.label("DG0"))
points(sb$`Na+`$T, sb$`Na+`$G, pch = 2)
# Add Cp
# change y-axis
par("usr" = c(par("usr")[1:2], -400, 100))
axis(4)
mtext(axis.label("Cp0"), side = 4, line = 3)
points(sb$water$T, sb$water$Cp, pch = 16)
points(sb$`Na+`$T, sb$`Na+`$Cp, pch = 17)
legend("topleft", c("H2O Cp", "H2O G", "Na+ Cp", "Na+ G"), pch = c(16, 1, 17, 2))
H2O <- expr.species("H2O")
Na <- expr.species("Na+")
degC <- expr.units("T")
title(main = substitute(H2O~and~Na~to~-20~degC, list(H2O = H2O, Na = Na, degC = degC)))
par(opar)

## Calculations using a variable-pressure standard state
thermo("opt$varP" = TRUE)
# Calculate the boiling point of octane at 2 and 20 bar
# We need exceed.Ttr = TRUE because the liquid is metastable
# at high temperatures (also, the gas is metastable at low
# temperatures, but that doesn't produce NA in the output)
sout2 <- subcrt(rep("octane", 2), c("liq", "gas"),
  c(-1, 1), T = seq(-50, 300, 0.1), P = 2, exceed.Ttr = TRUE)$out
sout20 <- subcrt(rep("octane", 2), c("liq", "gas"),
  c(-1, 1), T = seq(-50, 300, 0.1), P = 20, exceed.Ttr = TRUE)$out
# Find T with the Gibbs energy of reaction that is closest to zero
Tvap2 <- sout2$T[which.min(abs(sout2$G))]
Tvap20 <- sout20$T[which.min(abs(sout20$G))]
# Compare these with experimental values (Fig. 1 of Helgeson et al., 1998)
Tvap2.exp <- 156
Tvap20.exp <- 276
# Reset varP to FALSE (the default)
thermo("opt$varP" = FALSE)
}

\references{
  Anderson, G. M. and Crerar, D. A. (1993) \emph{Thermodynamics in Geochemistry: The Equilibrium Model}, Oxford University Press. \url{https://www.worldcat.org/oclc/803272549}

  Johnson, J. W., Oelkers, E. H. and Helgeson, H. C. (1992) SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\degC. \emph{Comp. Geosci.} \bold{18}, 899--947. \doi{10.1016/0098-3004(92)90029-Q}

  Helgeson, H. C., Owens, C. E., Knox, A. M. and Richard, L. (1998) Calculation of the standard molal thermodynamic properties of crystalline, liquid, and gas organic molecules at high temperatures and pressures. \emph{Geochim. Cosmochim. Acta} \bold{62}, 985--1081. \doi{10.1016/S0016-7037(97)00219-6}

  LaRowe, D. E. and Helgeson, H. C. (2007) Quantifying the energetics of metabolic reactions in diverse biogeochemical systems: electron flow and ATP synthesis. \emph{Geobiology} \bold{5}, 153--168. \doi{10.1111/j.1472-4669.2007.00099.x}

  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \degC and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \doi{10.1039/FT9928800803}
}

\concept{Main workflow}
