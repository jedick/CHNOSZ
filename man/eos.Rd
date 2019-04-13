\encoding{UTF-8}
\name{eos}
\alias{eos}
\alias{hkf}
\alias{cgl}
\alias{AkDi}
\title{Equations of State}
\description{
Calculate thermodynamic properties using the revised Helgeson-Kirkham-Flowers (HKF) or Akinfiev-Diamond (AkDi) equations of state for aqueous species, or using a generic heat capacity equation for crystalline, gas, and liquid species.
}

\usage{
  cgl(property = NULL, parameters = NULL, T = 298.15, P = 1)
  hkf(property = NULL, parameters = NULL, T = 298.15, P = 1,
    contrib = c("n", "s", "o"), H2O.props = "rho")
  AkDi(property = NULL, parameters = NULL, T = 298.15, P = 1, isPsat = TRUE)
}

\arguments{
  \item{property}{character, name(s) of properties to calculate}
  \item{parameters}{dataframe, species parameters as one or more rows from \code{thermo()$obigt}}
  \item{T}{numeric, temperature(s) at which to calculate properties (K)}
  \item{P}{numeric, pressure(s) at which to calculate properties (bar)}
  \item{contrib}{character, which contributions to consider in the revised HKF equations equations of state: (\code{n})onsolvation, (\code{s})olvation (the \eqn{\omega}{omega} terms), or (o)rigination contributions (i.e., the property itself at 25 \degC and 1 bar). Default is \code{c("n","s","o")}, for all contributions}
  \item{H2O.props}{character, properties to calculate for water}
  \item{isPsat}{logical, is this a calculation along the liquid-vapor saturation curve (Psat)?}
}

\details{
The equations of state permit the calculation of the standard molal properties of species as a function of temperature and pressure.
The \code{property} argument is required and refers to one or more of \samp{G}, \samp{H}, \samp{S}, \samp{Cp} and \samp{V}, and for aqueous species only, \samp{kT} and \samp{E}.
The units of these properties are the first ones shown in the description for \code{\link{subcrt}}.
The names of the properties are matched without regard to case. 

\code{hkf} implements the revised HKF equations of state (Helgeson et al., 1981; Tanger and Helgeson, 1988; Shock and Helgeson, 1988).
The equations-of-state parameters are \code{a1}, \code{a2}, \code{a3}, \code{a4}, \code{c1}, \code{c2}, \code{omega} and \code{Z}; the units of these parameters are as indicated for \code{\link{thermo}$obigt}, without the order of magnitude multipliers.
Note that the equation-of-state parameter \code{Z} (appearing in the \eqn{g}{g}-function for the temperature derivatives of the omega parameter; Shock et al., 1992) is taken from \code{thermo()$obigt} and \emph{not} from the \code{\link{makeup}} of the species.
\code{H2O.props} is an optional argument that lists the properties of water that should be returned; it is used by \code{\link{subcrt}} so that the time-consuming water calculations are only performed once.

The temperature and pressure derivatives of the \code{omega} parameter for charged species (where \code{Z != 0}, but not for the aqueous proton, H+) are calculated using the \eqn{g}{g}- and \eqn{f}{f}-functions described by Shock et al., 1992 and Johnson et al., 1992.
If the IAPWS-95 or DEW equations are activated (see \code{\link{water}}), only the \eqn{g}{g}-function (applicable to \samp{G}), but not its derivatives (needed for \samp{H}, \samp{S}, \samp{Cp}, and \samp{V}), is calculated.

The parameters in the \code{cgl} equations of state for crystalline, gas and liquid species (except liquid water) include \code{V}, \code{a}, \code{b}, \code{c}, \code{d}, \code{e}, \code{f} and \code{lambda}.
The terms denoted by \code{a}, \code{b} and \code{c} correspond to the Maier-Kelley equation for heat capacity (Maier and Kelley, 1932); the additional terms are useful for representing heat capacities of minerals (Robie and Hemingway, 1995) and gaseous or liquid organic species (Helgeson et al., 1998).
The standard molal volumes (\samp{V}) of species in these calculations are taken to be independent of temperature and pressure.

For both \code{hkf} and \code{cgl}, if at least one equations-of-state parameter for a species is provided, any NA values of the other parameters are reset to zero.
If all equations-of-state parameters are NA, but values of \samp{Cp} and/or \samp{V} are available, those values are used in the integration of \samp{G}, \samp{H} and \samp{S} as a function of temperature. 

\code{AkDi} provides the Akinfiev-Diamond model for aqueous species (Akinfiev and Diamond, 2003).
To run this code, the database must also include the corresponding gasesous species (with the same name or chemical formula).
Currently, only the standard chemical potential (Gibbs energy) is calculated.

}

\section{Warning}{
The range of applicability of the revised HKF equations of state for aqueous species corresponds to the stability region of liquid water or the supercritical fluid with density greater than 0.35 g/cm3, and between 0 to 1000 \degC and 1 to 5000 bar (Tanger and Helgeson, 1988; Shock and Helgeson, 1988).
The \code{hkf} function does not check these limits and will compute properties as long as the requisite electrostatic properties of water are available.
There are conceptually no temperature limits (other than 0 Kelvin) for the validity of the \code{cgl} equations of state.
However, the actual working upper temperature limits correspond to the temperatures of phase transitions of minerals or to those temperatures beyond which extrapolations from experimental data become highly uncertain.
These temperature limits are stored in the thermodynamic database for some minerals, but \code{cgl} ignores them; however, \code{\link{subcrt}} warns if they are exceeded.
}

\value{
A list of length equal to the number of species (i.e., number rows of \code{parameters}).
Each element of the list contains a dataframe, each column of which corresponds to one of the specified properties; the number of rows is equal to the number of pressure-temperature points.
Furthermore, in \code{hkf}, the output is a list consisting of the above-described object (named \samp{aq}) and a data frame of the calculated properties of water (named \samp{H2O}).
}

\seealso{
\code{\link{info}} for retrieving equations of state parameters from the thermodynamic database, \code{\link{water}} for equations of state of water, \code{\link{subcrt}} for interactive use of these equations.
}

\examples{
\dontshow{reset()}
## aqueous species
CH4aq <- info(info("methane", "aq"))	
hkf(property = "Cp", parameters = CH4aq)
# the non-solvation heat capacity
hkf(property = "Cp", parameters = CH4aq, contrib = "n")	
# at different temperature and pressure
hkf(property = "Cp", parameters = CH4aq, T = c(373.15,473.15), P = 1000)

## crystalline, gas, liquid species
CH4gas <- info(info("methane", "gas"))	
cgl(property = "Cp", parameters = CH4gas)
# melting and vaporization of octane
C8H18par <- info(info(rep("octane", 3), c("cr", "liq", "gas")))
myT <- seq(200, 420, 10)
DG0f <- cgl(property = "G", parameters = C8H18par, T = myT, P = 1)
cbind(T = myT, which.pmax(DG0f, pmin = TRUE))  # 1 = cr, 2 = liq, 3 = gas
# compare that result with the tabulated transition temperatures
print(C8H18par)
}

\references{

  Akinfiev, N. N. and Diamond, L. W. (2003) Thermodynamic description of aqueous nonelectrolytes at infinite dilution over a wide range of state parameters. \emph{Geochim. Cosmochim. Acta} \bold{67}, 613--629. \url{https://doi.org/10.1016/S0016-7037(02)01141-9}

  Helgeson, H. C., Kirkham, D. H. and Flowers, G. C. (1981) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures. IV. Calculation of activity coefficients, osmotic coefficients, and apparent molal and standard and relative partial molal properties to 600\degC and 5 Kb. \emph{Am. J. Sci.} \bold{281}, 1249--1516. \url{https://doi.org/10.2475/ajs.281.10.1249}

  Helgeson, H. C., Owens, C. E., Knox, A. M. and Richard, L. (1998) Calculation of the standard molal thermodynamic properties of crystalline, liquid, and gas organic molecules at high temperatures and pressures. \emph{Geochim. Cosmochim. Acta} \bold{62}, 985--1081. \url{https://doi.org/10.1016/S0016-7037(97)00219-6}

  Maier, C. G. and Kelley, K. K. (1932) An equation for the representation of high-temperature heat content data. \emph{J. Am. Chem. Soc.} \bold{54}, 3243--3246. \url{https://doi.org/10.1021/ja01347a029}

  Robie, R. A. and Hemingway, B. S. (1995) \emph{Thermodynamic Properties of Minerals and Related Substances at 298.15 K and 1 Bar (\eqn{10^5} Pascals) Pressure and at Higher Temperatures}. U. S. Geol. Surv., Bull. 2131, 461 p. \url{http://www.worldcat.org/oclc/32590140}

  Shock, E. L. and Helgeson, H. C. (1988) Calculation of the thermodynamic and transport properties of aqueous species at high pressures and temperatures: Correlation algorithms for ionic species and equation of state predictions to 5 kb and 1000\degC. \emph{Geochim. Cosmochim. Acta} \bold{52}, 2009--2036. \url{https://doi.org/10.1016/0016-7037(88)90181-0}
  
  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \degC and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{https://doi.org/10.1039/FT9928800803}

  Tanger, J. C. IV and Helgeson, H. C. (1988) Calculation of the thermodynamic and transport properties of aqueous species at high pressures and temperatures: Revised equations of state for the standard partial molal properties of ions and electrolytes. \emph{Am. J. Sci.} \bold{288}, 19--98. \url{https://doi.org/10.2475/ajs.288.1.19}

}

\concept{Thermodynamic calculations}
