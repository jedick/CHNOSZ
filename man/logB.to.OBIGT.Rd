\encoding{UTF-8}
\name{logB.to.OBIGT}
\alias{logB.to.OBIGT}
\title{Fit Thermodynamic Parameters to Formation Constants (\logB)}
\description{
Fit thermodynamic parameters to experimental formation constants for an aqueous species and add the parameters to OBIGT.
}

\usage{
  logB.to.OBIGT(logB, species, coeffs, T, P, optimize.omega = TRUE, tolerance = 0.05)
}

\arguments{
  \item{logB}{Values of \logB}
  \item{species}{Species in the formation reaction (the formed species must be last)}
  \item{coeffs}{Coefficients of the formation reaction}
  \item{T}{Temperature (units of \code{\link{T.units}})}
  \item{P}{Temperature (\samp{Psat} or numerical values with units of \code{\link{P.units}})}
  \item{optimize.omega}{logical, optimize the omega parameter (relevant for charged species)?}
  \item{tolerance}{Tolerance for checking equivalence of input and calculated \logB values}
}

\details{

This function can be used to update the \code{\link{OBIGT}} thermodynamic database with parameters fit to formation constants of aqueous species as a function of temperature.
The \code{logB} argument should have decimal logarithm of equilibrium constants (\logB) for a formation reaction.
The formation reaction is defined by \code{species} and \code{coeffs}.
All species in the formation reaction must be available in OBIGT except for the \emph{last} species, which is the newly formed species.

The function works as follows.
First, the properties of the known species are combined with \logB to calculate the standard Gibbs energy (G[T]) of the formed species at each value of \code{T} and \code{P}.
Then, the values of G[T] are used to retrieve the thermodynamic parameters \code{G}, \code{S}, \code{c1}, \code{c2}, and \code{omega}.
The first two of these are standard-state values at 25 \degC and 1 bar, and the last three are parameters in the revised HKF equations (e.g. Eq. B25 of Shock et al., 1992).
Volumetric parameters (\a1 to \a4) in the HKF equations currently aren't included, so the resulting fits are valid only at the pressure corresponding to the given values of \logB.
The fitted parameters for the formed species are then added to OBIGT.
\code{\link{all.equal}} is used to test for approximate equivalence of the input and calculated equilibrium constants; if the mean absolute difference exceeds \code{tolerance}, an error occurs.

The number of parameters used in the fit is not more than the number of data values (\emph{n}), and the parameters are added in the order stated above.
If \emph{n} is less than 5, then the values of the unused parameters are set to 0 for addition to OBIGT.
For instance, a single value of \logB would be fit only with \code{G}, implying that computed values of G[T] have no temperature dependence.

For uncharged species, the value of \omega is a constant, but for charged species, it is a function of \T and \P as described by the \dQuote{g function} (Shock et al., 1992).
By default an optimization step is used to refine the value of \code{omega} at 25 \degC and 1 bar for charged species.
In case the data are not sufficient to constrain \code{omega} (see first example), this step can be skipped by setting \code{optimize.omega} to FALSE; then, it is assumed that \omega is constant during the fitting.
In OBIGT, a constant omega for a charged species can be enforced by setting the \code{z} parameter to 0, but that is not done by this function.

}

\seealso{
\code{logB.to.OBIGT} calls \code{\link{mod.OBIGT}} with \code{zap = TRUE} to clear the parameters of a formed species if it already exists in the OBIGT database.
If preexisting parameters (e.g. volumetric HKF coefficients) weren't cleared, they would interfere with the fitting done here, which uses only selected parameters.
}

\value{
  The species index of the new species in OBIGT.
}

\examples{
## CoHS+ from Migdisov et al. (2011)
logB <- c(6.24, 6.02, 5.84, 5.97, 6.52)
T <- c(120, 150, 200, 250, 300)
P <- "Psat"
species <- c("Co+2", "HS-", "CoHS+")
coeffs <- c(-1, -1, 1)
opar <- par(mfrow = c(2, 1))
for(optimize.omega in c(TRUE, FALSE)) { 
  # Fit the parameters with or without variable omega
  inew <- logB.to.OBIGT(logB, species, coeffs, T, P, optimize.omega = optimize.omega)
  # Print the new database entry
  info(inew)
  # Plot experimental logB
  plot(T, logB, "n", c(100, 320), c(5.8, 6.8), xlab = axis.label("T"), ylab = quote(log~beta))
  points(T, logB, pch = 19, cex = 2)
  # Plot calculated values
  Tfit <- seq(100, 320, 10)
  sres <- subcrt(species, coeffs, T = Tfit)
  lines(sres$out$T, sres$out$logK, col = 4)
  title(describe.reaction(sres$reaction))
  legend <- c("Migdisov et al. (2011)", paste0("logB.to.OBIGT(optimize.omega =",optimize.omega,")"))
  legend("top", legend, pch = c(19, NA), lty = c(0, 1), col = c(1, 4), pt.cex = 2, bg = "#FFFFFFB0")
}
par(opar)

## ZnCl+ from Mei et al. (2015)
# These are values for 5000 bar calculated from the modified Ryzhenko-Bryzgalin (RB) model
logB <- c(-1.93, -1.16, -0.38, 0.45, 1.15, 1.76, 2.30, 2.80, 3.26, 3.70, 4.12, 4.53, 4.92)
species <- c("Zn+2", "Cl-", "ZnCl+")
coeffs <- c(-1, -1, 1)
T <- c(25, 60, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600)
P <- 5000
# Note: ZnCl+ is in the default OBIGT database;
# logB.to.OBIGT() "zaps" the previous parameters before adding the fitted ones
inew <- logB.to.OBIGT(logB, species, coeffs, T, P)
# Plot RB and logB.to.OBIGT values
plot(T, logB, xlab = axis.label("T"), ylab = axis.label("logB"), pch = 19, cex = 2)
Tfit <- seq(25, 600, 25)
sres <- subcrt(species, coeffs, T = Tfit, P = P)
lines(sres$out$T, sres$out$logK, col = 4)
title(describe.reaction(sres$reaction), line = 3)
title("5000 bar", font.main = 1, line = 1)
legend <- c("Mei et al. (2015)", "logB.to.OBIGT()")
legend("topleft", legend, pch = c(19, NA), lty = c(0, 1), col = c(1, 4), pt.cex = 2)
}

\references{
Migdisov, Art. A., Zezin, D. and Williams-Jones, A. E. (2011) An experimental study of Cobalt (II) complexation in Cl\S{-} and H\s{2}S-bearing hydrothermal solutions. \emph{Geochim. Cosmochim. Acta} \bold{75}, 4065--4079. \doi{10.1016/j.gca.2011.05.003}

Mei, Y., Sherman, D. M., Liu, W., Etschmann, B., Testemale, D. and Brugger, J. (2015) Zinc complexation in chloride-rich hydrothermal fluids (25--600 \degC): A thermodynamic model derived from \emph{ab initio} molecular dynamics. \emph{Geochim. Cosmochim. Acta} \bold{150}, 264--284. \doi{10.1016/j.gca.2014.09.023}

Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \degC and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \doi{10.1039/FT9928800803}
}
