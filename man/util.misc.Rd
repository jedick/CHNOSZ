\encoding{UTF-8}
\name{util.misc}
\alias{util.misc}
\alias{dPdTtr}
\alias{Ttr}
\alias{GHS_Tr}
\alias{unitize}
\title{Functions for Miscellaneous Tasks}
\description{
Calculate \eqn{dP/dT}{dP/dT} and temperature of phase transitions; scale logarithms of activity to a desired total activity.
}

\usage{
  dPdTtr(ispecies, ispecies2 = NULL)
  Ttr(ispecies, ispecies2 = NULL, P = 1, dPdT = NULL)
  GHS_Tr(ispecies, Htr)
  unitize(logact = NULL, length = NULL, logact.tot = 0)
}

\arguments{
  \item{ispecies}{numeric, species index of a mineral phase}
  \item{ispecies2}{numeric, species index of next mineral phase (the default is ispecies + 1)}
  \item{P}{numeric, pressure (bar)}
  \item{dPdT}{numeric, values of (\eqn{dP/dT}{dP/dT}) of phase transitions (\code{Ttr})}
  \item{Htr}{numeric, enthalpy(ies) of transition (cal/mol)}
  \item{logact}{numeric, logarithms of activity}
  \item{length}{numeric, numbers of residues}
  \item{logact.tot}{numeric, logarithm of total activity}
}

\details{
\code{dPdTtr} returns values of \eqn{(dP/dT)_{Ttr}}{(dP/dT)Ttr}, where \eqn{Ttr}{Ttr} represents the transition temperature, of the phase transition at the high-\eqn{T}{T} stability limit of the \code{ispecies} in \code{thermo$obigt} (other than checking that the names match, the function does not check that the species in fact represent different phases of the same mineral).
\code{dPdTtr} takes account of the Clapeyron equation, \eqn{(dP/dT)_{Ttr}}{(dP/dT)Ttr}=\eqn{{\Delta}S/{\Delta}V}{DS/DV}, where \eqn{{\Delta}S}{DS} and \eqn{{\Delta}V}{DV} represent the changes in entropy and volume of phase transition, and are calculated using \code{subcrt} at Ttr from the standard molal entropies and volumes of the two phases involved.
Using values of \code{dPdT} calculated using \code{dPdTtr} or supplied in the arguments, \code{Ttr} returns as a function of \code{P} values of the upper transition temperature of the mineral phase represented by \code{ispecies}.

\code{GHS_Tr} can be used to calculate values of G, H, and S at Tr for the cr2, cr3, and cr4 phases in the database.
It combines the given \code{Htr} (enthalpies of transition) with the database values of GHS @ Tr only for the phase that is stable at 298.15 K (cr) and the transition temperatures and Cp coefficients for higher-temperature phases, to calculate the GHS @ Tr (i.e. low-temperature metastable conditions) of the phases that are stable at higher temperatures.

\code{unitize} scales the logarithms of activities given in \code{logact} so that the logarithm of total activity of residues is equal to zero (i.e. total activity of residues is one), or to some other value set in \code{logact.tot}.
\code{length} indicates the number of residues in each species.
If \code{logact} is NULL, the function takes the logarithms of activities from the current species definition.
If any of those species are proteins, the function gets their lengths using \code{protein.length}.
}

\examples{\dontshow{data(thermo)}
# we need the Helgeson et al., 1978 minerals for this example
add.obigt("SUPCRT92")
# that replaces the existing enstatite with the first phase;
# the other phases are appended to the end of thermo$obigt
i1 <- info("enstatite")
i2 <- info("enstatite", "cr2")
i3 <- info("enstatite", "cr3")
# (dP/dT) of transitions
dPdTtr(i1, i2)  # first transition
dPdTtr(i2, i3)  # second transition
# temperature of transitions (Ttr) as a function of P
Ttr(i1, i2, P=c(1,10,100,1000))
Ttr(i2, i3, P=c(1,10,100,1000))
# restore default database
data(OBIGT)

# calculate the GHS at Tr for the high-temperature phases of iron
# using transition enthalpies from the SUPCRT92 database (sprons92.dat)
Htr <- c(326.0, 215.0, 165.0)
iiron <- info("iron")
GHS_Tr(iiron, Htr)
# the results calculated above are stored in the database ...
info(1:3 + iiron)[, c("G", "H", "S")]
# ... meaning that we can recalculate the transition enthalpies using subcrt()
sapply(info(0:2 + iiron)$T, function(T) {
  # a very small T increment around the transition temperature
  T <- convert(c(T-0.01, T), "C")
  # use suppressMessages to make the output less crowded
  substuff <- suppressMessages(subcrt("iron", T=T, P=1))
  diff(substuff$out$iron$H)
})

## scale logarithms of activity
# suppose we have two proteins whose lengths are 100 and 
# 200; what are the logarithms of activity of the proteins 
# that are equal to each other and that give a total 
# activity of residues equal to unity?
logact <- c(-3,-3)  # could be any two equal numbers
length <- c(100,200)
logact.tot <- 0
loga <- unitize(logact,length,logact.tot)
# the proteins have equal activity
stopifnot(identical(loga[1],loga[2]))
# the sum of activity of the residues is unity
stopifnot(isTRUE(all.equal(sum(10^loga * length),1)))
## now, what if the activity of protein 2 is ten
## times that of protein 1?
logact <- c(-3,-2)
loga <- unitize(logact,length,logact.tot)
# the proteins have unequal activity
stopifnot(isTRUE(all.equal(loga[2]-loga[1],1)))
# but the activities of residues still add up to one
stopifnot(isTRUE(all.equal(sum(10^loga * length),1)))
}

\concept{Thermodynamic calculations}
