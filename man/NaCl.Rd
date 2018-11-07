\encoding{UTF-8}
\name{NaCl}
\alias{NaCl}
\title{Simple NaCl-Water Solution}
\description{
Calculate speciation and ionic strength of aqueous solutions with a given molality of NaCl.
}

\usage{
  NaCl(T = seq(100, 500, 100), P = 1000, m_tot = 2)
}

\arguments{
  \item{T}{numeric, temperature in \degC}
  \item{P}{numeric, pressure in bar (single value)}
  \item{m_tot}{numeric, total molality of NaCl (single value)}
}

\details{
This function calculates speciation (ion activities) and ionic strength in aqueous solutions given a total amount (\code{m_tot}, in mol/kg) of NaCl.
The function is written for quick calculations along a temperature range (\code{T}) at constant pressure (\code{P}).
The only reaction considered is Na\S{+} + Cl\S{-} = NaCl(aq).
The algorithm starts by calculating the equilibrium constant (\emph{K}) of the reaction and assuming complete dissociation of NaCl(aq).
This also permits calculating the ionic strength from the molalities of Na\S{+} and Cl\S{-}.
Then, \code{\link{uniroot}} is used to find the equilibrium molality of Cl\S{-}; that is, where the affinity of the reaction (log(\emph{K}/\emph{Q})) becomes zero.
The activity quotient (\emph{Q}) is evaluated taking account of activity coefficients of Na\S{+}, Cl\S{-}, and NaCl(aq) calculated for the nominal ionic strength (see \code{\link{nonideal}}).
The calculated molality of Cl\S{-} yields a new estimate of the ionic strength of the system.
The calculations are iterated until the deviation in ionic strength at all temperatures is less than 0.01.
}

\section{Warning}{
This function provides only a first-order estimate of the solution composition, and is intended for solubility calculations of relatively insoluble metals in NaCl-dominated solutions.
The formation of other species such as HCl or NaOH is not accounted for.
}

\value{
A list with components \samp{IS} (\dQuote{true} ionic strength from concentrations of unpaired ions), \samp{m_Cl} (molality of Cl\S{-}), \samp{gam_Na}, and \samp{gam_Cl} (activity coefficients of Na\S{+} and Cl\S{-}).
}

\seealso{
\code{demo("gold")} for an application of this function.
}

\examples{\dontshow{data(thermo)}
# ionic strength of solution and activity coefficient of Cl-
# from HCh version 3.7 (Shvarov and Bastrakov, 1999) at 1000 bar,
# 100, 200, and 300 degress C, and 1 to 6 molal NaCl
m.HCh <- 1:6
IS.HCh <- list(`100`=c(0.992, 1.969, 2.926, 3.858, 4.758, 5.619),
               `300`=c(0.807, 1.499, 2.136, 2.739, 3.317, 3.875),
               `500`=c(0.311, 0.590, 0.861, 1.125, 1.385, 1.642))
gam_Cl.HCh <- list(`100`=c(0.565, 0.545, 0.551, 0.567, 0.589, 0.615),
                  `300`=c(0.366, 0.307, 0.275, 0.254, 0.238, 0.224),
                  `500`=c(0.19, 0.137, 0.111, 0.096, 0.085, 0.077))
# total molality in the calculation with NaCl()
m_tot <- seq(1, 6, 0.5)
N <- length(m_tot)
# where we'll put the calculated values
IS.calc <- data.frame(`100`=numeric(N), `300`=numeric(N), `500`=numeric(N))
gam_Cl.calc <- data.frame(`100`=numeric(N), `300`=numeric(N), `500`=numeric(N))
# NaCl() is *not* vectorized over m_tot, so we use a loop here
for(i in 1:length(m_tot)) {
  NaCl.out <- NaCl(c(100, 300, 500), P=1000, m_tot=m_tot[i])
  IS.calc[i, ] <- NaCl.out$IS
  gam_Cl.calc[i, ] <- NaCl.out$gam_Cl
}
# plot ionic strength from HCh and NaCl() as points and lines
opar <- par(mfrow=c(2, 1))
col <- c("black", "red", "orange")
plot(c(1,6), c(0,6), xlab="NaCl (mol/kg)", ylab=axis.label("IS"), type="n")
for(i in 1:3) {
  # NOTE: the differences are probably mostly due to different models
  # for the properties of NaCl(aq) (HCh: B.Ryhzenko model;
  # CHONSZ: revised HKF with parameters from Shock et al., 1997)
  points(m.HCh, IS.HCh[[i]], col=col[i])
  lines(m_tot, IS.calc[, i], col=col[i])
}
# add 1:1 line, legend, and title
abline(0, 1, lty=3)
dprop <- describe.property(rep("T", 3), c(100, 300, 500))
legend("topleft", dprop, lty=1, pch=1, col=col)
title(main="H2O + NaCl; HCh (points) and 'NaCl()' (lines)")
plot(c(1,6), c(0,0.8), xlab="NaCl (mol/kg)", ylab=expression(gamma[Cl^"-"]), type="n")
# plot activity coefficient (gamma)
for(i in 1:3) {
  points(m.HCh, gam_Cl.HCh[[i]], col=col[i])
  lines(m_tot, gam_Cl.calc[, i], col=col[i])
}
# we should be fairly close
stopifnot(maxdiff(unlist(gam_Cl.calc[seq(1,11,2), ]), unlist(gam_Cl.HCh)) < 0.033)
par(opar)
}

\references{
Shvarov, Y. and Bastrakov, E. (1999) HCh: A software package for geochemical equilibrium modelling. User's Guide. \emph{Australian Geological Survey Organisation} \bold{1999/25}. \url{http://pid.geoscience.gov.au/dataset/ga/25473}
}

\concept{Extended workflow}
