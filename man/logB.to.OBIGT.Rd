\encoding{UTF-8}
\name{logB.to.OBIGT}
\alias{logB.to.OBIGT}
\title{Fit Thermodynamic Parameters to Formation Constants (\logB)}
\description{
Fit thermodynamic parameters to experimental formation constants for an aqueous species and add the parameters to OBIGT.
}

\usage{
  logB.to.OBIGT(logB, species, coeff, T, P, tolerance = 0.05)
}

\arguments{
  \item{logB}{Values of \logB}
  \item{species}{Species in the formation reaction (the formed species must be last)}
  \item{coeff}{Coefficients of the formation reaction}
  \item{T}{Temperature (units of \code{\link{T.units}})}
  \item{P}{Temperature (\samp{Psat} or numerical values with units of \code{\link{P.units}})}
  \item{tolerance}{Tolerance for checking equivalence of input and calculated \logB values}
}

\details{

\code{logB.to.OBIGT} starts with the decimal logarithms of equilibrium constants (\logB) for a formation reaction.
The formation reaction is defined by \code{species} and \code{coeff}.
All species in the formation reaction must be available in OBIGT except for the \emph{last} species, which must be the newly formed species.
The properties of the known species are combined with \logB to calculate the standard Gibbs energy (G[T]) of the formed species as a function of \T and \P.

The values of G[T] are fit using the thermodynamic parameters \code{G}, \code{S}, and \code{Cp} (i.e., standard-state values at 25 °C and 1 bar).
The equation is this: G[T] = G[Tr] + -S[Tr] * (T - Tr) + Cp[Tr] (T - Tr - T * ln(T / Tr)), where G, S, and Cp are standard-state Gibbs energy, entropy, and isobaric heat capacity, T is the temperature in Kelvin, Tr is the reference temperature (298.15 K), and ln stands for the natural logarithm.
Note that the equation assumes a constant heat capacity and no volume change, and it may be unsuitable for extrapolation beyond the original \T and \P range.

The fitted \code{G}, \code{S}, and \code{Cp} for the formed species are added to OBIGT.
\code{\link{all.equal}} is used to test the near equivalence of the input and calculated equilibrium constants; if the mean absolute difference exceeds \code{tolerance}, an error occurs.

}

\seealso{
\code{logB.to.OBIGT} calls \code{\link{mod.OBIGT}} with \code{zap = TRUE} to clear the parameters of a formed species if it already exists in the OBIGT database.
If preexisting parameters (e.g. HKF coefficients) weren't cleared, they would interfere with the fitting done here, which uses only three parameters (\code{G}, \code{S}, and \code{Cp}).
}

\examples{
## CoHS+ from Migdisov et al. (2011)
logB <- c(6.24, 6.02, 5.84, 5.97, 6.52)
species <- c("Co+2", "HS-", "CoHS+")
coeff <- c(-1, -1, 1)
T <- c(120, 150, 200, 250, 300)
P <- "Psat"
inew <- logB.to.OBIGT(logB, species, coeff, T, P)
# Print the new database entry
# (Note: H is calculated from G and S)
info(inew)
# Make a plot of experimental logB
xlim <- c()
plot(T, logB, "n", c(100, 320), c(5.8, 6.8), xlab = axis.label("T"), ylab = quote(log~beta))
points(T, logB, pch = 19, cex = 2)
# Plot fitted values
Tfit <- seq(100, 320, 10)
sres <- subcrt(species, coeff, T = Tfit)
lines(sres$out$T, sres$out$logK, col = 4)
title(describe.reaction(sres$reaction))
legend <- c("Migdisov et al. (2011)", "logB.to.OBIGT()")
legend("topleft", legend, pch = c(19, NA), lty = c(0, 1), col = c(1, 4), pt.cex = 2)

## ZnCl+ from Mei et al. (2015)
# These are values for 5000 bar calculated from the modified Ryzhenko-Bryzgalin model
logB <- c(-1.93, -1.16, -0.38, 0.45, 1.15, 1.76, 2.30, 2.80, 3.26, 3.70, 4.12, 4.53, 4.92)
species <- c("Zn+2", "Cl-", "ZnCl+")
coeff <- c(-1, -1, 1)
T <- c(25, 60, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600)
P <- 5000
# Note: ZnCl+ is in the default OBIGT database;
# logB.to.OBIGT() "zaps" the previous parameters before adding the fitted ones
inew <- logB.to.OBIGT(logB, species, coeff, T, P)
# Plot MRB and logB.to.OBIGT fitted values
plot(T, logB, xlab = axis.label("T"), ylab = axis.label("logB"), pch = 19, cex = 2)
Tfit <- seq(25, 600, 25)
sres <- subcrt(species, coeff, T = Tfit, P = P)
lines(sres$out$T, sres$out$logK, col = 4)
title(describe.reaction(sres$reaction), line = 3)
title("5000 bar", font.main = 1, line = 1)
legend <- c("Mei et al. (2015)", "logB.to.OBIGT()")
legend("topleft", legend, pch = c(19, NA), lty = c(0, 1), col = c(1, 4), pt.cex = 2)
}

\references{
Migdisov, Art. A., Zezin, D. and Williams-Jones, A. E. (2011) An experimental study of Cobalt (II) complexation in Cl\S{-} and H\s{2}S-bearing hydrothermal solutions. \emph{Geochim. Cosmochim. Acta} \bold{75}, 4065--4079. \doi{10.1016/j.gca.2011.05.003}

Mei, Y., Sherman, D. M., Liu, W., Etschmann, B., Testemale, D. and Brugger, J. (2015) Zinc complexation in chloride-rich hydrothermal fluids (25--600 \degC): A thermodynamic model derived from \emph{ab initio} molecular dynamics. \emph{Geochim. Cosmochim. Acta} \bold{150}, 264--284. \doi{10.1016/j.gca.2014.09.023}
}
