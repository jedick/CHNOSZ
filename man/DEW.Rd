\encoding{UTF-8}
\name{DEW}
\alias{DEW}
\alias{calculateDensity}
\alias{calculateGibbsOfWater}
\alias{calculateEpsilon}
\alias{calculateQ}
\title{Deep Earth Water (DEW) Model}
\description{
Calculate thermodynamic properties of water using the Deep Earth Water (DEW) model.
}

\usage{
  calculateDensity(pressure, temperature, error = 0.01)
  calculateGibbsOfWater(pressure, temperature)
  calculateEpsilon(density, temperature)
  calculateQ(density, temperature)
}

\arguments{
  \item{pressure}{numeric, pressure (bar)}
  \item{temperature}{numeric, temperature (\degC)}
  \item{error}{numeric, residual error for bisection calculation}
  \item{density}{numeric, density (g/cm^3)}
}

\details{
The Deep Earth Water (DEW) model, described by Sverjensky et al., 2014, extends the applicability of the revised HKF equations of state to 60 kbar.
This implementation of DEW is based on the VBA macro code in the May, 2017 version of the DEW spreadsheet downloaded from \url{http://www.dewcommunity.org/}.
The spreadsheet provides multiple options for some calculations; here the default equations for density of water (Zhang and Duan, 2005), dielectric constant (Sverjensky et al., 2014) and Gibbs energy of water (integral of volume, equation created by Brandon Harrison) are used.

Comments in the original code indicate that \code{calculateGibbsOfWater} is valid for 100 \le \T \le 1000 \degC and \P \ge 1000 bar.
Likewise, the power function fit of the dielectric constant (epsilon) is valid for 100 \le \T \le 1200 \degC and \P \ge 1000 bar (Sverjensky et al., 2014).
}

\value{
The calculated values of density, Gibbs energy, and the Q Born coefficient have units of g/cm^3, cal/mol, and bar^-1 (epsilon is dimensionless).
}

\seealso{
\code{\link{water.DEW}}; use \code{water("DEW")} to activate these equations for the main functions in CHNOSZ.
}

\examples{
pressure <- c(1000, 60000)
temperature <- c(100, 1000)
calculateGibbsOfWater(pressure, temperature)
(density <- calculateDensity(pressure, temperature))
calculateEpsilon(density, temperature)
calculateQ(density, temperature)
}

\references{
Sverjensky, D. A., Harrison, B. and Azzolini, D. (2014) Water in the deep Earth: The dielectric constant and the solubilities of quartz and corundum to 60 kb and 1,200 \degC. \emph{Geochim. Cosmochim. Acta} \bold{129}, 125--145. \doi{10.1016/j.gca.2013.12.019}

Zhang, Z. and Duan, Z. (2005) Prediction of the \emph{PVT} properties of water over wide range of temperatures and pressures from molecular dynamics simulation. \emph{Phys. Earth Planet. Inter.} \bold{149}, 335--354. \doi{10.1016/j.pepi.2004.11.003}
}

\concept{Water properties}
