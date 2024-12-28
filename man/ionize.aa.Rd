\encoding{UTF-8}
\name{ionize.aa}
\alias{ionize.aa}
\title{Properties of ionization of proteins}
\description{
  Calculate the charges of proteins and contributions of ionization to the thermodynamic properties of proteins.
}

\usage{
  ionize.aa(aa, property = "Z", T = 25, P = "Psat", pH = 7,
    ret.val = NULL, suppress.Cys = FALSE)
}

\arguments{
  \item{aa}{data frame, amino acid composition in the format of \code{thermo()$protein}}
  \item{property}{character, property to calculate}
  \item{T}{numeric, temperature in \degC}
  \item{P}{numeric, pressure in bar, or \samp{Psat} for vapor pressure of \H2O above 100 \degC}
  \item{pH}{numeric, pH}
  \item{ret.val}{character, return the indicated value from intermediate calculations}
  \item{suppress.Cys}{logical, suppress (ignore) the ionization of the cysteine groups?}
}

\details{
  The properties of ionization of proteins calculated by this function take account of the standard molal thermodynamic properties of ionizable amino acid sidechain groups and the terminal groups in proteins ([AABB]) and their equations of state parameters taken from Dick et al., 2006. The values of the ionization constants (pK) are calculated as a function of temperature, and the charges and the ionization contributions of other thermodynamic properties to the proteins are calculated additively, without consideration of electrostatic interactions, so they are best applied to the unfolded protein reference state.

  For each amino acid composition in \code{aa}, the additive value of the \code{property} is calculated as a function of \code{T}, \code{P} and \code{pH}. \code{property} can be NULL to denote net charge, or if not NULL is one of the properties available in \code{\link{subcrt}}, or is \samp{A} to calculate the dimensionless chemical affinity (A/2.303RT) of the ionization reaction for the protein. If \code{ret.val} is one of \samp{pK}, \samp{alpha}, or \samp{aavals} it indicates to return the value of the ionization constant, degree of formation, or the values of the \code{property} for each ionizable group rather than taking their sums for the amino acid compositions in \code{aa}.
}

\value{
  The function returns a matrix (possibly with only one row or column) with number of rows corresponding to the longest of \code{T}, \code{P} or \code{pH} (values of any of these with shorter length are recycled) and a column for each of the amino acid compositions in \code{aa}.
}

\seealso{ \code{\link{pinfo}}, \code{\link{affinity}} }

\examples{\dontshow{reset()}
## Heat capacity of LYSC_CHICK as a function of T
pH <- c(5, 9, 3)
T <- seq(0, 100)
# Cp of non-ionized protein
Cp.nonion <- subcrt("LYSC_CHICK", T = T)$out[[1]]$Cp
plot(T, Cp.nonion, xlab = axis.label("T"), type = "l",
  ylab = axis.label("Cp"), ylim = c(20000, 35000))
# Cp of ionization and ionized protein
aa <- pinfo(pinfo("LYSC_CHICK"))
for(pH in c(5, 9, 3)) {
  Cp.ionized <- Cp.nonion + ionize.aa(aa, "Cp", T = T, pH = pH)[, 1]
  lines(T, Cp.ionized, lty = 2)
  text(80, Cp.ionized[70], paste("pH =", pH) )
}
# Makhatadze and Privalov's group contributions
T <- c(5, 25, 50, 75, 100, 125)
points(T, MP90.cp("LYSC_CHICK", T))
# Privalov and Makhatadze's experimental values
e <- read.csv(system.file("extdata/misc/PM90.csv", package = "CHNOSZ"))
points(e$T, e$LYSC_CHICK, pch = 16)
legend("bottomright", pch = c(16, 1, NA, NA), lty = c(NA, NA, 1, 2),
  legend = c("PM90 experiment", "MP90 groups", 
  "DLH06 groups no ion", "DLH06 groups ionized"))
title("Heat capacity of unfolded LYSC_CHICK")
}

\references{
  Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \doi{10.5194/bg-3-311-2006}

  Makhatadze, G. I. and Privalov, P. L. (1990) Heat capacity of proteins. 1. Partial molar heat capacity of individual amino acid residues in aqueous solution: Hydration effect. \emph{J. Mol. Biol.} \bold{213}, 375--384. \doi{10.1016/S0022-2836(05)80197-4}

  Privalov, P. L. and Makhatadze, G. I. (1990) Heat capacity of proteins. II. Partial molar heat capacity of the unfolded polypeptide chain of proteins: Protein unfolding effects. \emph{J. Mol. Biol.} \bold{213}, 385--391. \doi{10.1016/S0022-2836(05)80198-6}
}

\concept{Protein properties}
