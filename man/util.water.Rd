\encoding{UTF-8}
\name{util.water}
\alias{util.water}
\alias{WP02.auxiliary}
\alias{rho.IAPWS95}
\alias{water.AW90}
\title{Functions for Properties of Water and Steam}
\description{
  Utility functions for properties of water and steam.
}

\usage{
  WP02.auxiliary(property, T = 298.15)
  rho.IAPWS95(T = 298.15, P = 1, state="", trace=0)
  water.AW90(T = 298.15, rho = 1000, P = 0.1)
}

\arguments{
  \item{property}{character, property to calculate}
  \item{T}{numeric, temperature (K)}
  \item{P}{numeric, pressure (units of bar, except MPa for \code{water.AW90})}
  \item{state}{character, state or phase of \H2O}
  \item{trace}{integer number}
  \item{rho}{numeric, density (kg m\eqn{^{-3}}{^-3})}
}

\details{

Auxiliary equations to the IAPWS-95 formulation (Wagner and Pruß, 2002) are provided in \code{WP02.auxiliary}.
The \code{property} for this function can be one of \samp{P.sigma} (saturation vapor pressure in MPa), \samp{dP.sigma.dT} (derivative of saturation vapor pressure with respect to temperature), or \samp{rho.liquid} or \samp{rho.vapor} (density of liquid or vapor in kg m\eqn{^{-3}}{^-3}).

\code{rho.IAPWS95} implements a root-finding technique (using \code{\link{uniroot}}) to determine the values of density for the stable phase of \H2O at the given temperature and pressure.
The \code{state} option is used internally in order to determine the stable phase at conditions close to saturation (0.9999*\Psat <= P <= 1.00005*\Psat, where \Psat is the saturation pressure calculated by \code{\link{WP02.auxiliary}}).
Alternatively, the user can specify a \code{state} of \samp{liquid} or \samp{vapor} to force the calculation of density for the corresponding phase, even if it is metastable (e.g. superheated water, supercooled steam; this option has no effect in the supercritical region).
The \code{state} is set in calls by \code{\link{water.IAPWS95}} to the value in \code{thermo$opt$IAPWS.sat} (default \samp{liquid}) so that higher-level functions (\code{\link{water}}, \code{\link{subcrt}}) take properties for that state along the saturation curve.
Diagnostic messages are printed if \code{trace} is positive (it is also included in the call to \code{uniroot}).

\code{water.AW90} provides values of the static dielectric constant (\code{epsilon}) calculated using equations given by Archer and Wang, 1990.
}

\examples{\dontshow{data(thermo)}
# calculate density of stable phase at 500 K, 500 bar
rho <- rho.IAPWS95(T=500, P=500)
# calculate pressure (= 50 MPa) at this density
IAPWS95("P", T=500, rho=rho)
# calculate dielectric constant 	
water.AW90(T=500, rho=rho, P=50)

# density along saturation curve
T <- seq(273.15, 623.15, 25)
WP02.auxiliary(T=T)  # liquid from WP02
WP02.auxiliary("rho.vapor", T)  # vapor from WP02	

# WP02.auxiliary gives a close estimate of saturation pressure...
T <- 445:455
P.sigma <- WP02.auxiliary("P.sigma", T)
# ... but alternates between being just on the liquid or vapor side
# (low rho: steam; high rho: water)
rho.IAPWS95(T, convert(P.sigma, "bar"))
# thermo$opt$IAPWS.sat allows for choosing liquid or vapor or ""
thermo$opt$IAPWS.sat <<- ""
# shows artifactual vapor-liquid transition
water.IAPWS95("V", T, "Psat")
# the calculated Psat, while not exact, should be close enough for most
# geochemical calculations to select liquid or vapor
oldwat <- water("IAPWS95")
thermo$opt$IAPWS.sat <<- "vapor"
V.vapor <- subcrt("water", T=convert(445:455, "C"))$out[[1]]$V
thermo$opt$IAPWS.sat <<- "liquid" # the default
V.liquid <- subcrt("water", T=convert(445:455, "C"))$out[[1]]$V
stopifnot(all(V.vapor > V.liquid))
water(oldwat)
}

\references{
Archer, D. G. and Wang, P. M. (1990) The dielectric constant of water and Debye-Hückel limiting law slopes. \emph{J. Phys. Chem. Ref. Data} \bold{19}, 371--411. \url{https://srd.nist.gov/JPCRD/jpcrd383.pdf}

Wagner, W. and Pruß, A. (2002) The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. \emph{J. Phys. Chem. Ref. Data} \bold{31}, 387--535. \url{https://doi.org/10.1063/1.1461829}
}

\concept{Water properties}
