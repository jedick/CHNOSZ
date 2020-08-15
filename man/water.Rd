\encoding{UTF-8}
\name{water}
\alias{water}
\alias{water.SUPCRT92}
\alias{water.IAPWS95}
\alias{water.DEW}
\title{Properties of Water}
\description{
Calculate thermodynamic and electrostatic properties of water.
}

\usage{
  water(property = NULL, T = 298.15, P = "Psat", P1 = TRUE)
  water.SUPCRT92(property=NULL, T = 298.15, P = 1, P1 = TRUE)
  water.IAPWS95(property=NULL, T = 298.15, P = 1)
  water.DEW(property=NULL, T = 373.15, P = 1000)
}

\arguments{
  \item{property}{character, computational setting or property(s) to calculate}
  \item{T}{numeric, temperature (K)}
  \item{P}{numeric, pressure (bar), or \samp{Psat} for vapor-liquid saturation}
  \item{P1}{logical, output pressure of 1 bar below 100 \degC instead of calculated values of \samp{Psat}?}
}

\details{

These functions compute the thermodynamic (Gibbs energy and it derivatives) and electrostatic (dielectric constant and its derivatives) properties of liquid or supercritical \H2O as a function of temperature and pressure using equations of state taken from the literature.
The high-level function \code{water} performs different computations, depending on the setting of \code{\link{thermo}$opt$water}:

\describe{

  \item{\samp{SUPCRT92} (default) or \samp{SUPCRT}}{Thermodynamic and electrostatic properties are calculated using a FORTRAN subroutine taken from the \acronym{SUPCRT92} software package (Johnson et al., 1992). See more information below.}

  \item{\samp{IAPWS95} or \samp{IAPWS}}{Thermodynamic properties are calculated using an implementation in \R code of the \acronym{IAPWS-95} formulation (Wagner and Pruss, 2002), and electrostatic properties are calculated using the equations of Archer and Wang, 1990. See \code{\link{IAPWS95}} and more information below.}

  \item{\samp{DEW}}{Thermodynamic and electrostatic properties are calculated using the Deep Earth Water (\acronym{DEW}) model (Sverjensky et al., 2014). The defaults for \code{T} and \code{P} reflect the minimum values for applicability of the model; calculations at lower \code{T} and/or \code{P} points fall back to using \samp{SUPCRT92}. See \code{\link{DEW}}.}

}

Calling the function with no arguments returns the current computational setting.
Use e.g. \code{water("DEW")} to make the setting; the previous setting (at the time of the function call) is returned invisibly.
Subsequent calculations with \code{water}, or other functions such as \code{subcrt} and \code{affinity}, will use that setting.

The allowed \code{property}s for \code{water} are one or more of those given below, depending on the computational setting; availability is shown by an asterisk.
Note that some of the properties that can actually be calculated using the different formulations are not implemented here.
Except for \code{rho}, the units are those used by Johnson and Norton, 1991.

  \tabular{llllll}{
     Property \tab Description \tab Units \tab IAPWS95 \tab SUPCRT92 \tab DEW \cr
     \code{A} \tab Helmholtz energy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \tab NA \cr
     \code{G} \tab Gibbs energy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \tab * \cr
     \code{S} \tab Entropy \tab cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \tab * \tab * \tab NA \cr
     \code{U} \tab Internal energy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \tab NA \cr
     \code{H} \tab Enthalpy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \tab NA \cr
     \code{Cv} \tab Isochoric heat capacity \tab cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \tab * \tab * \tab NA \cr
     \code{Cp} \tab Isobaric heat capacity \tab cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \tab * \tab * \tab NA \cr
     \code{Speed} \tab Speed of sound \tab cm s\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \code{alpha} \tab Coefficient of isobaric expansivity \tab K\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \code{beta} \tab Coefficient of isothermal compressibility \tab bar\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \code{epsilon} \tab Dielectric constant \tab dimensionless \tab NA \tab * \tab * \cr
     \code{visc} \tab Dynamic viscosity \tab g cm\eqn{^{-1}}{^-1} s\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \code{tcond} \tab Thermal conductivity \tab cal cm\eqn{^{-1}}{^-1} s\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \code{tdiff} \tab Thermal diffusivity \tab cm\eqn{^2} s\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \code{Prndtl} \tab Prandtl number \tab dimensionless \tab NA \tab * \tab NA \cr
     \code{visck} \tab Kinematic viscosity \tab cm\eqn{^2} s\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \code{albe} \tab Isochoric expansivity \tab bar K\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \tab -compressibility \tab \tab \tab \tab \cr
     \code{ZBorn} \tab Z Born function \tab dimensionless \tab NA \tab * \tab NA \cr
     \code{YBorn} \tab Y Born function \tab K\eqn{^{-1}}{^-1} \tab * \tab * \tab NA \cr
     \code{QBorn} \tab Q Born function \tab bar\eqn{^{-1}}{^-1} \tab * \tab * \tab * \cr
     \code{daldT} \tab Isobaric temperature derivative \tab K\eqn{^{-2}}{^-2} \tab NA \tab * \tab NA \cr
     \tab of expansibility \tab \tab \tab \tab \cr
     \code{XBorn} \tab X Born function \tab K\eqn{^{-2}}{^-2} \tab * \tab * \tab NA \cr
     \code{NBorn} \tab N Born function \tab bar\eqn{^{-2}}{^-2} \tab * \tab NA \tab NA \cr
     \code{UBorn} \tab U Born function \tab bar\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} \tab * \tab NA \tab NA \cr
     \code{V} \tab Volume \tab cm\eqn{^3} mol\eqn{^{-1}}{^-1} \tab * \tab * \tab * \cr
     \code{rho} \tab Density \tab kg m\eqn{^3} \tab * \tab * \tab * \cr     
     \code{Psat} \tab Saturation vapor pressure \tab bar \tab * \tab * \tab NA \cr
     \code{E} \tab Isobaric expansivity \tab cm\eqn{^3} K\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \code{kT} \tab Isothermal compressibility \tab cm\eqn{^3} bar\eqn{^{-1}}{^-1} \tab NA \tab * \tab NA \cr
     \code{de.dT} \tab Temperature derivative \tab K\eqn{^{-1}}{^-1} \tab * \tab NA \tab NA \cr
     \tab of dielectric constant \tab \tab \tab \tab \cr
     \code{de.dP} \tab Pressure derivative \tab bar\eqn{^{-1}}{^-1} \tab * \tab NA \tab NA \cr
     \tab of dielectric constant \tab \tab \tab \tab \cr
     \code{P} \tab Pressure \tab bar \tab * \tab NA \tab NA \cr
     \code{A_DH} \tab A Debye-Huckel parameter \tab kg\eqn{^{0.5}}{0.5} mol\eqn{^{-0.5}}{-0.5} \tab * \tab * \tab * \cr
     \code{B_DH} \tab B Debye-Huckel parameter \tab kg\eqn{^{0.5}}{0.5} mol\eqn{^{-0.5}}{-0.5} cm\eqn{^{-1}}{-1} \tab * \tab * \tab * \cr
  }

Call \code{water.SUPCRT92}, \code{water.IAPWS95}, or \code{water.DEW} with no arguments to list the available properties.

\code{water.SUPCRT92} interfaces to the FORTRAN subroutine taken from the \acronym{SUPCRT92} package (H2O92D.F) for calculating properties of water.
These calculations are based on data and equations of Levelt-Sengers et al., 1983, Haar et al., 1984, and Johnson and Norton, 1991, among others (see Johnson et al., 1992).
A value of \code{P} set to \samp{Psat} refers to one bar below 100 \degC, otherwise to the vapor-liquid saturation pressure at temperatures below the critical point (\samp{Psat} is not available at temperatures above the critical point).
\code{water.SUPCRT92} provides a limited interface to the FORTRAN subroutine; some functions provided there are not made available here (e.g., using variable density instead of pressure, or calculating the properties of steam).

The stated temperature limits of validity of calculations in \code{water.SUPCRT92} are from the greater of 0 \degC or the melting temperature at pressure, to 2250 \degC (Johnson et al., 1992).
Valid pressures are from the greater of zero bar or the melting pressure at temperature to 30000 bar.
The present functions do not check these limits and will attempt calculations for any range of input parameters, but may return \code{NA} for properties that fail to be calculated at given temperatures and pressures and/or produce warnings or even errors when problems are encountered.

Starting with version 0.9-9.4, a check for minimum pressure (in valTP function in H2O92D.f) has been bypassed so that properties of H2O can be calculated using \code{water.SUPCRT92} at temperatures below the 0.01 \degC triple point.
A primary check is still enforced (Tbtm), giving a minimum valid temperature of 253.15 K.

\code{water.IAPWS95} is a wrapper around \code{\link{IAPWS95}}, \code{\link{rho.IAPWS95}} and \code{\link{water.AW90}}.
\code{water.IAPWS95} provides for calculations at specific temperature and pressure; density, needed for \code{IAPWS95}, is inverted from pressure using \code{rho.IAPWS95}.
The function also contains routines for calculating the Born functions as numerical derivatives of the static dielectric constant (from \code{water.AW90}).
For compatibility with geochemical modeling conventions, the values of Gibbs energy, enthalpy and entropy output by \code{IAPWS95} are converted by \code{water.IAPWS95} to the triple point reference state adopted in \code{SUPCRT92} (Johnson and Norton, 1991; Helgeson and Kirkham, 1974).
\code{water.IAPWS95} also accepts setting \code{P} to \samp{Psat}, with the saturation pressure calculated from \code{\link{WP02.auxiliary}}; by default the returned properties are for the liquid, but this can be changed to the vapor in \code{\link{thermo}$opt$IAPWS.sat}.

\code{A_DH} and \code{B_DH} are solvent parameters in the \dQuote{B-dot} (extended Debye-Huckel) equation (Helgeson, 1969; Manning, 2013).
}


\value{
A data frame, the number of rows of which corresponds to the number of input temperature-pressure pairs.
}

\seealso{
For calculating properties of reactions, \code{\link{subcrt}} coordinates the calculation of properties among \code{water} and \code{\link{hkf}} and \code{\link{cgl}} for other species.
}

\examples{\dontshow{reset()}
## calculations along saturation curve
T <- seq(273.15, 623.15, 25)
# liquid density, from SUPCRT92
water("rho", T=T, P="Psat")
# values of the saturation pressure, Gibbs energy	
water(c("Psat", "G"), T=T, P="Psat") 
# derivatives of the dielectric constant (Born functions)
water(c("QBorn", "YBorn", "XBorn"), T=T, P="Psat")
# now at constant pressure
water(c("QBorn", "YBorn", "XBorn"), T=T, P=2000)

## comparing the formulations
T <- convert(c(25, 100, 200, 300), "K")
# IAPWS-95
oldwat <- water("IAPWS95")
water(water.IAPWS95(), T=T)
# Deep Earth Water (DEW)
water("DEW")
water(water.DEW(), T=T, P=1000)
# SUPCRT92 (the default)
water(oldwat)
water(water.SUPCRT92(), T=T)

## calculating Q Born function
# after Table 22 of Johnson and Norton, 1991
T <- rep(c(375, 400, 425, 450, 475), each=5)
P <- rep(c(250, 300, 350, 400, 450), 5)
w <- water("QBorn", T=convert(T, "K"), P=P)
# the rest is to make a neat table
w <- as.data.frame(matrix(w[[1]], nrow=5))
colnames(w) <- T[1:5*5]
rownames(w) <- P[1:5]
print(w)
}

\references{
Archer, D. G. and Wang, P. M. (1990) The dielectric constant of water and Debye-Hückel limiting law slopes. \emph{J. Phys. Chem. Ref. Data} \bold{19}, 371--411. \url{https://doi.org/10.1063/1.555853}

Haar, L., Gallagher, J. S. and Kell, G. S. (1984) \emph{NBS/NRC Steam Tables}. Hemisphere, Washington, D. C., 320 p. \url{https://www.worldcat.org/oclc/301304139}

Helgeson, H. C. and Kirkham, D. H. (1974) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures. I. Summary of the thermodynamic/electrostatic properties of the solvent. \emph{Am. J. Sci.} \bold{274}, 1089--1098. \url{https://doi.org/10.2475/ajs.274.10.1089}

Helgeson, H. C. (1969) Thermodynamics of hydrothermal systems at elevated temperatures and pressures. \emph{Am. J. Sci.} \bold{267}, 729--804. \url{https://doi.org/10.2475/ajs.267.7.729}

Johnson, J. W. and Norton, D. (1991) Critical phenomena in hydrothermal systems: state, thermodynamic, electrostatic, and transport properties of H\eqn{_2}{2}O in the critical region. \emph{Am. J. Sci.} \bold{291}, 541--648. \url{https://doi.org/10.2475/ajs.291.6.541}

Johnson, J. W., Oelkers, E. H. and Helgeson, H. C. (1992) SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\degC. \emph{Comp. Geosci.} \bold{18}, 899--947. \url{https://doi.org/10.1016/0098-3004(92)90029-Q}

Levelt-Sengers, J. M. H., Kamgarparsi, B., Balfour, F. W. and Sengers, J. V. (1983) Thermodynamic properties of steam in the critical region. \emph{J. Phys. Chem. Ref. Data} \bold{12}, 1--28. \url{https://doi.org/10.1063/1.555676}

Manning, C. E. (2013) Thermodynamic modeling of fluid-rock interaction at mid-crustal to upper-mantle conditions. \emph{Rev. Mineral. Geochem.} \bold{76}, 135--164. \url{https://doi.org/10.2138/rmg.2013.76.5}

Sverjensky, D. A., Harrison, B. and Azzolini, D. (2014) Water in the deep Earth: The dielectric constant and the solubilities of quartz and corundum to 60 kb and 1,200 \degC. \emph{Geochim. Cosmochim. Acta} \bold{129}, 125--145. \url{https://doi.org/10.1016/j.gca.2013.12.019}

Wagner, W. and Pruss, A. (2002) The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. \emph{J. Phys. Chem. Ref. Data} \bold{31}, 387--535. \url{https://doi.org/10.1063/1.1461829}
}

\concept{Water properties}
