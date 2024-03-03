\encoding{UTF-8}
\name{logB.to.OBIGT}
\alias{logB.to.OBIGT}
\title{Fit thermodynamic parameters to formation constants (\logB)}
\description{
Fit thermodynamic parameters to experimental formation constants for an aqueous species and add the parameters to OBIGT.
}

\usage{
  logB.to.OBIGT(logB, species, coeffs, T, P, npar = 3,
    optimize.omega = FALSE, tolerance = 0.05)
}

\arguments{
  \item{logB}{Values of \logB}
  \item{species}{Species in the formation reaction (the formed species must be last)}
  \item{coeffs}{Coefficients of the formation reaction}
  \item{T}{Temperature (units of \code{\link{T.units}})}
  \item{P}{Temperature (\samp{Psat} or numerical values with units of \code{\link{P.units}})}
  \item{npar}{numeric, number of parameters to use in fitting}
  \item{optimize.omega}{logical, optimize the omega parameter (relevant for charged species)?}
  \item{tolerance}{Tolerance for checking equivalence of input and calculated \logB values}
}

\details{

This function updates the \code{\link{OBIGT}} thermodynamic database with parameters fit to formation constants of aqueous species as a function of temperature.
The \code{logB} argument should have decimal logarithm of formation constants for an aqueous complex (\logB).
The formation reaction is defined by \code{species} and \code{coeffs}.
All species in the formation reaction must be present in OBIGT except for the \emph{last} species, which is the newly formed species.

The function works as follows.
First, the properties of the known species are combined with \logB to calculate the standard Gibbs energy (G[T]) of the formed species at each value of \code{T} and \code{P}.
Then, \code{\link{lm}} is used to fit one or more of the thermodynamic parameters \code{G}, \code{S}, \code{c1}, \code{c2}, and \code{omega} to the values of G[T].
The first two of these parameters are standard-state values at 25 \degC and 1 bar, and the last three are parameters in the revised HKF equations (e.g. Eq. B25 of Shock et al., 1992).
The fitted parameters for the formed species are then added to OBIGT.
Finally, \code{\link{all.equal}} is used to test for approximate equivalence of the input values of \logB and calculated equilibrium constants; if the mean absolute difference exceeds \code{tolerance}, an error occurs.

To avoid overfitting, only the first three parameters (\code{G}, \code{S}, and \code{c1}) are used by default.
More parameters (up to 5) or fewer (down to 1) can be selected by changing \code{npar}.
Volumetric parameters (\a1 to \a4) in the HKF equations currently aren't included, so the resulting fits are valid only at the input pressure values.

Independent of \code{npar}, the number of parameters used in the fit is not more than the number of data values (\emph{n}).
If \emph{n} is less than 5, then the values of the unused parameters are set to 0 for addition to OBIGT.
For instance, a single value of \logB would be fit only with \code{G}, implying that computed values of G[T] have no temperature dependence.

The value of \omega is a constant in the revised HKF equations for uncharged species, but for charged species, it is a function of \T and \P as described by the \dQuote{g function} (Shock et al., 1992).
An optimization step is available to refine the value of \code{omega} at 25 \degC and 1 bar for charged species taking its temperature dependence into account for the fitting.
However, representative datasets are not well represented by variable \code{omega} (see first example), so this step is skipped by default.
Furthermore, \code{logB.to.OBIGT} by default also sets the \code{z} parameter in OBIGT to 0; this enforces a constant \omega for charged species in calculations with \code{\link{subcrt}} (note that this is a parameter for the HKF equations and does not affect reaction balancing).
Set \code{optimize.omega} to TRUE to override the defaults and activate variable \omega for charged species; this takes effect only if \code{npar == 5} and \emph{n} > 5.

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
for(o.o in c(TRUE, FALSE)) { 
  # Fit the parameters with or without variable omega
  inew <- logB.to.OBIGT(logB, species, coeffs, T, P, npar = 5, optimize.omega = o.o)
  # Print the new database entry
  info(inew)
  # Plot experimental logB
  plot(T, logB, "n", c(100, 320), c(5.8, 6.8),
    xlab = axis.label("T"), ylab = quote(log~beta))
  points(T, logB, pch = 19, cex = 2)
  # Plot calculated values
  Tfit <- seq(100, 320, 10)
  sres <- subcrt(species, coeffs, T = Tfit)
  lines(sres$out$T, sres$out$logK, col = 4)
  title(describe.reaction(sres$reaction))
  legend <- c("Migdisov et al. (2011)",paste0("logB.to.OBIGT(optimize.omega = ",o.o,")"))
  legend("top", legend, pch = c(19, NA), lty = c(0, 1), col = c(1, 4),
    pt.cex = 2, bg = "#FFFFFFB0")
}
par(opar)
# NB. Optimizing omega leads to unphysical oscillations in the logK (first plot)

## ZnCl+ from Mei et al. (2015)
# Values for 5000 bar calculated with the modified Ryzhenko-Bryzgalin (RB) model
logB <- c(-1.93,-1.16,-0.38,0.45,1.15,1.76,2.30,2.80,3.26,3.70,4.12,4.53,4.92)
species <- c("Zn+2", "Cl-", "ZnCl+")
coeffs <- c(-1, -1, 1)
T <- c(25, 60, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600)
P <- 5000
# Note: ZnCl+ is in the default OBIGT database;
# logB.to.OBIGT() "zaps" the previous parameters before adding the fitted ones
inew <- logB.to.OBIGT(logB, species, coeffs, T, P, npar = 5)
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
