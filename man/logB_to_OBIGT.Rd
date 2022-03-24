\encoding{UTF-8}
\name{logB_to_OBIGT}
\alias{logB_to_OBIGT}
\title{Fit Thermodynamic Parameters to Formation Constants (\logB)}
\description{
Fit thermodynamic parameters to experimental formation constants for an aqueous species and add the parameters to OBIGT.
}

\usage{
  logB_to_OBIGT(logB, species, coeff, T, P, tolerance = 0.05)
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

\code{logB_to_OBIGT} starts with the decimal logarithms of equilibrium constants (\logB) for a formation reaction.
The formation reaction is defined by \code{species} and \code{coeff}; all species in the formation reaction must be available in OBIGT except for the last species, which is the newly formed species.
The properties of the known species are combined with \logB to calculate the standard Gibbs energy (G[T]) of the formed species as a function of \T and \P.

The values of G[T] are fit using the thermodynamic parameters \emph{G}, \emph{S}, and \emph{Cp} (i.e., standard-state values at 25 °C and 1 bar).
The equation is this: G[T] = G[Tr] + -S[Tr] * (T - Tr) + Cp[Tr] (T - Tr - T * log(T / Tr)), where G, S, and Cp are standard-state Gibbs energy, entropy, and isobaric heat capacity, T is the temperature in Kelvin, Tr is the reference temperature (298.15 K), and log stands for the natural logarithm.
Note that the equation assumes a constant heat capacity and no volume change, and it may be unsuitable for extrapolation beyond the original \T and \P range.

The fitted \emph{G}, \emph{S}, and \emph{Cp} for the formed species are added to OBIGT.
\code{\link{all.equal}} is used to test the near equivalence of the input and calculated equilibrium constants; if the differences exceed \code{tolerance}, the code generates an error.

}

\seealso{
See \code{\link{add.OBIGT}} for more information about modifying the database.
}


\examples{
# Formation constants of CoHS+ from Migdisov et al. (2011)
logB <- c(6.24, 6.02, 5.84, 5.97, 6.52)
species <- c("Co+2", "HS-", "CoHS+")
coeff <- c(-1, -1, 1)
T <- c(120, 150, 200, 250, 300)
P <- "Psat"
inew <- logB_to_OBIGT(logB, species, coeff, T, P)
# Print the new database entry
# (Note: H is calculated from G and S)
info(inew)
# Make a plot of experimental logB
xlim <- c()
plot(T, logB, "n", c(100, 320), c(5.8, 6.8), xlab = axis.label("T"))
points(T, logB, pch = 19, cex = 2)
# Add fitted values
Tfit <- seq(100, 320, 10)
sres <- subcrt(species, coeff, T = Tfit)
lines(sres$out$T, sres$out$logK, col = 4)
title(describe.reaction(sres$reaction))
legend <- c("Migdisov et al. (2011)", "logB_to_OBIGT()")
legend("topleft", legend, pch = c(19, NA), lty = c(0, 1), col = c(1, 4), pt.cex = 2)
}

\references{
Migdisov, Art. A., Zezin, D. and Williams-Jones, A. E. (2011) An experimental study of Cobalt (II) complexation in Cl\S{-} and H\s{2}S-bearing hydrothermal solutions. \emph{Geochim. Cosmochim. Acta} \bold{75}, 4065--4079. \doi{10.1016/j.gca.2011.05.003}
}
