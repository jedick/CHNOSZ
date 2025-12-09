\encoding{UTF-8}
\name{IAPWS95}
\alias{IAPWS95}
\title{Properties of water from IAPWS-95}
\description{
Calculate thermodynamic properties of water following the IAPWS-95 formulation.
}

\usage{
  IAPWS95(property, T = 298.15, rho = 1000)
}

\arguments{
  \item{property}{character, name(s) of property(s) to calculate}
  \item{T}{numeric, temperature (K)}
  \item{rho}{numeric, density (kg m\eqn{^{-3}}{^-3})}
}

\details{

\code{IAPWS95} provides an implementation of the IAPWS-95 formulation for properties (including pressure) calculated as a function of temperature and density.

The \code{IAPWS95} function returns values of thermodynamic properties in specific units (per gram).
The IAPWS-95 formulation follows the triple point convention used in engineering (values of internal energy and entropy are taken to be zero at the triple point).

For \code{IAPWS95} the upper temperature limit of validity is 1000 \degC, but extrapolation to much higher temperatures is possible (Wagner and Pruss, 2002).
Valid pressures are from the greater of zero bar or the melting pressure at temperature to 10000 bar (with the provision for extrapolation to more extreme conditions).
The function does not check these limits and will attempt calculations for any range of input parameters, but may return \code{NA} for properties that fail to be calculated at given temperatures and pressures and/or produce warnings or even errors when problems are encountered.
}

\value{
  A data frame the number of rows of which corresponds to the number of input temperature, pressure and/or density values.
}

\seealso{
\code{\link{util.water}} for properties along the saturation curve (\code{WP02.auxiliary}) and calculation of density from pressure and temperature (\code{rho.IAPWS92}).
\code{\link{water.IAPWS95}} is a wrapper around \code{IAPWS95} and the utility functions, which converts the specific units to molar quantities, and is used in higher-level functions (\code{\link{water}}).
}

\examples{
# Calculate pressure for given temperature and density
IAPWS95("P", T = 500, rho = 838.0235)
}

\references{
Wagner W, Pruss A. 2002. The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. \emph{J Phys Chem Ref Data} \bold{31}: 387--535. \doi{10.1063/1.1461829}
}

\concept{Water properties}
