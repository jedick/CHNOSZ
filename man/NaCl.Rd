\encoding{UTF-8}
\name{NaCl}
\alias{NaCl}
\title{Simple NaCl-Water Solution}
\description{
Calculate speciation and ionic strength of aqueous solutions with a given molality of NaCl.
}

\usage{
  NaCl(m_tot = 1, T = 25, P = "Psat", pH = NA, attenuate = FALSE)
}

\arguments{
  \item{m_tot}{numeric, total molality of NaCl (single value)}
  \item{T}{numeric, temperature in \degC}
  \item{P}{numeric, pressure in bar}
  \item{pH}{numeric, pH}
  \item{attenuate}{logical, halve changes of variables in each step?}
}

\details{
Thermodynamic models for metal solubility and speciation involving chloride complexes are commonly specified in terms of amount of NaCl rather than activity (or molality) of Cl\S{-} as an independent variable.
This function calculates distribution of species and ionic strength in a simple aqueous solution given a total amount (\code{m_tot}, in mol/kg) of NaCl.
The aqueous Cl-bearing species considered in the system are Cl\S{-}, NaCl, and optionally HCl.
Na\S{+} is present as a basis species, but the formation of Na-bearing species such as NaOH is not considered.
The activity coefficients of charged species are calculated using the extended Debye-Hückel equation (see \code{\link{nonideal}}) via the \code{IS} argument of \code{\link{affinity}}.
The function first sets the molality of Na\S{+} and ionic strength equal to \code{m_tot}, then calculates the distribution of Cl-bearing species.
Based on mass balance of Na atoms, the molality of NaCl is then used to recalculate the molality of Na\S{+}, followed by ionic strength.
To find a solution, the function iterates until the change of molality of Na\S{+} and ionic strength are both less than \code{m_tot} / 100.

At very high NaCl concentrations, which are beyond the applicability limits of the extended Debye-Hückel model and therefore not recommended for normal use, the iterations tend to oscillate without converging.
Setting \code{attenuate} to TRUE, which halves the amount of change in each step, may help with convergence.
If a solution is not found after 100 iterations, the function stops with an error.

If \code{pH} is NA (the default), then HCl is not included in the calculation and its molality in the output is also assigned NA.
Note that only a single value is accepted for \code{m_tot}, but the other numeric arguments can have length > 1, allowing multiple combinations \code{T}, \code{P}, and \code{pH} in a single function call.
However, due to limitations in \code{\link{affinity}}, only one of \code{T} and \code{P} can have length > 1.
}

\section{Warning}{
It is important to keep in mind the ionic strength limits of the Debye-Hückel equation, but this function doesn't enforce them.
Furthermore, metal-ligand complexing is not calculated by this function, so metal solubility and speciation calculations will be accurate only for relatively insoluble metals in NaCl-dominated solutions.
}

\value{
A list with components \samp{IS} (ionic strength calculated from molalities of Na\S{+} and Cl\S{-}), \samp{m_Cl}, \samp{m_Cl}, \samp{m_NaCl}, and \samp{m_HCl} (molalities of Na\S{+}, Cl\S{-}, NaCl, and HCl).
}

\seealso{
This function is used in a few demos (\code{demo("contour")}, \code{demo("gold")}, \code{demo("minsol")}, \code{demo("sphalerite")}).
\code{demo("yttrium")} uses the \code{pH} argument.
}

\examples{\dontshow{reset()}
# Ionic strength calculated with HCh version 3.7 (Shvarov and Bastrakov, 1999)
# at 1000 bar, 100, 200, and 300 degress C, and 1 to 6 molal NaCl
m.HCh <- 1:6
IS.HCh <- list(`100` = c(0.992, 1.969, 2.926, 3.858, 4.758, 5.619),
               `300` = c(0.807, 1.499, 2.136, 2.739, 3.317, 3.875),
               `500` = c(0.311, 0.590, 0.861, 1.125, 1.385, 1.642))
# Total molality in the calculation with NaCl()
m_tot <- seq(1, 6, 0.5)
N <- length(m_tot)
# Where we'll put the calculated values
IS.calc <- data.frame(`100` = numeric(N), `300` = numeric(N), `500` = numeric(N))
# NaCl() is *not* vectorized over m_tot, so we use a loop here
for(i in 1:length(m_tot)) {
  NaCl.out <- NaCl(m_tot[i], c(100, 300, 500), P = 1000)
  IS.calc[i, ] <- NaCl.out$IS
}
# Plot ionic strength from HCh and NaCl() as points and lines
col <- c("black", "red", "orange")
plot(c(1,6), c(0,6), xlab = "NaCl (mol/kg)", ylab = axis.label("IS"), type = "n")
for(i in 1:3) {
  # NOTE: the differences are probably mostly due to different models
  # for the properties of NaCl(aq) (HCh: B.Ryhzenko model;
  # CHONSZ: revised HKF with parameters from Shock et al., 1997)
  points(m.HCh, IS.HCh[[i]], col = col[i])
  lines(m_tot, IS.calc[, i], col = col[i])
}
# Add legend and title
dprop <- describe.property(rep("T", 3), c(100, 300, 500))
legend("topleft", dprop, lty = 1, pch = 1, col = col)
title(main="H2O + NaCl; HCh (points) and 'NaCl()' (lines)")
}

\references{
Shvarov, Y. and Bastrakov, E. (1999) HCh: A software package for geochemical equilibrium modelling. User's Guide. \emph{Australian Geological Survey Organisation} \bold{1999/25}. \url{http://pid.geoscience.gov.au/dataset/ga/25473}
}

\concept{Extended workflow}
