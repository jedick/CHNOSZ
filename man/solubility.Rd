\encoding{UTF-8}
\name{solubility}
\alias{solubility}
\title{Equilibrium chemical activities of species}
\description{
Calculate chemical activities of aqueous species in equilibrium with a mineral or gas.
}

\usage{
  solubility(iaq, ..., in.terms.of = NULL, dissociate = FALSE, find.IS = FALSE)
}

\arguments{
  \item{iaq}{character (names) or numeric (species indices) of aqueous species}
  \item{...}{arguments for \code{\link{affinity}} or \code{\link{mosaic}} (i.e. plotting variables)}
  \item{in.terms.of}{character, express the total solubility in terms of moles of this species}
  \item{dissociate}{logical, does the mineral undergo a dissociation reaction?}
  \item{find.IS}{logical, find the equilibrium ionic strength by iteration?}
}

\details{
\code{solubility} calculates the activities of aqueous species in equilibrium with one or more minerals or gases.
The minerals or gases should be loaded as the formed \code{\link{species}}, and the aqueous species (including ions and/or neutral complexes) that can be produced by dissolution should be listed in the \code{iaq} argument.
The definitions of plotting variables should be provided in \code{...}, which are passed as arguments to \code{\link{affinity}}, or to \code{\link{mosaic}} if the first one is named \code{bases}.

It must be possible to obtain a valid set of basis species by substituting each of the minerals or gases in the first position of the current \code{\link{basis}} defintion, and all of the aqueous species must include that basis species in their formation reactions.
(This essentially means that all minerals, gases and aqueous species must share a common element, which is what the reactions are balanced on.)

For a single mineral or gas, the output of \code{solubility} can be used by \code{\link{diagram}} with \code{type = "loga.balance"} (the default when \code{type} is \samp{auto}) to plot the sum of activities of aqueous species.
This value corresponds to the solubility of the mineral or gas; \code{in.terms.of} can be used to express this value in terms of another species or element.
For example, for dissolution of gaseous S\s{2}, \code{in.terms.of = "S"} gives the total amount of S in solution, which is twice the amount of S\s{2} dissolved.
Likewise, the solubility of corundum (Al\s{2}O\s{3}) can be expressed in terms of the moles of Al\S{+3} in solution (see the vignette \viglink{anintro}).


For multiple minerals, the function calculates the solubilities for each of the minerals separately; these are stored in the \code{loga.equil} element of the output.
The overall \emph{minimum} solubility among all the minerals at each point is stored in \code{loga.balance}.
This corresponds to the total activity of dissolved species in equilibrium with the most stable mineral.

Note that \code{\link{diagram}} by default plots \code{loga.balance}.
Use \code{type = "loga.equil"} in the function call to plot the solubilities for the individual minerals or gases, or the activities of aqueous species when there is one mineral or gas being dissolved.
}

\section{Backward Compatibility}{
For compatibility with previous versions of the function, the \code{iaq} argument can be the output of \code{\link{affinity}} or \code{\link{mosaic}} for aqueous species.
The examples for ionic strength and dissociation reactions were designed for this calling style.

In this case the (single) mineral or gas being dissolved is taken from the current \code{\link{basis}} species.
Usually, the basis species should be set up so that the first basis species represents the substance being dissolved (a mineral such as CaCO\s3 or gas such as CO\s2).
This is treated as the conserved basis species, so it must be present in all of the formation reactions of the aqueous species.

The \code{\link{species}} should be defined to represent one set of aqueous species (including ions and/or neutral complexes) formed in solution, all involving the conserved basis species.
For a dissociation reaction, the second basis species should be used to represent the counterion (cation or anion).
Other variables (pH, ionic strength, activities of other basis species) should be defined in the call to \code{\link{affinity}} to make \code{iaq}.
}

\section{Dissociation Reactions}{
The function perfoms some additional steps to calculate the solubility of something that dissociates (not just dissolves).
For example, the dissolution of calcite (CaCO\s3), involves the release of both calcium ions and different forms of carbonate in solution, depending on the pH.
The equilibrium calculation must take account of the \emph{total} activity of the shared ion (Ca\S{+2}), which is unknown at the start of the calculation.
The solution is found by recalculating the affinities, essentially working backward from the assumption that the dissociation didn't occur.
The resulting activities correspond to equilibrium considering the system-wide activity of Ca\S{+2}.

A \emph{not recommended} alternative is to set \code{dissociate} to a numeric value corresponding to the number of dissociated species (i.e. 2 for a 1:1 electrolyte).
This setting indicates to calculate activities on a per-reaction basis, where each reaction has its own (independent) activity of Ca\S{+2}.
That does not give a complete equilibrium in the system, but may be required to reproduce some published diagrams (see comment in the calcite example of \code{demo("solubility")}).
}

\section{Ionic Strength}{
Set \code{find.IS} to TRUE to determine the final ionic strength due to dissolution of a substance in pure water.
This works by calculating the ionic strength from the amounts of aqueous species formed, then re-running \code{affinity} with the calculated \code{IS} value.
Note that for dissociation reactions, the ionic strength is calculated from both the ions present in the species definition and the counter ion, which should be the second basis species.
The calculation is iterated until the ionic strength deviation at every point is lower than a preset tolerance (1e-4).
Alternatively, speciation of counterions (e.g. ionized forms of carbonate or sulfate) can also be performed by using the \code{\link{mosaic}} function instead of \code{affinity}; this is used in the second example below.
}

\section{Warning}{
This function has not been tested for systems that may form dimers or higher-order complexes (such as Au\s{2}S\s{2}\S{2-}).
Except for relatively simple systems, even after careful refinement, the results from CHNOSZ, which considers chemical activities as the independent variables, will not match the results from speciation-solubility (or Gibbs energy minimization) codes, where the system is defined by its bulk composition.
}

\seealso{
\code{\link{retrieve}} provides a way to list all of the aqueous species in the database that have the specified elements.

\code{demo("solubility")} shows solubilities of CO\s{2} and calcite calculated as a function of pH and T; note that for calcite, the \code{dissociate} argument is set to TRUE.
\code{demo("gold")} shows solubility calculations for Au in aqueous solutions with hydroxide, chloride, and hydrosulfide complexes.

Solubility calculations for multiple minerals are used for generating isosolubility (aka. equisolubility) lines in \code{demo("Pourbaix")} and \code{demo("minsol")}.
The latter demo combines the calculation of solubilities with a \code{\link{mosaic}} calculation to account for the speciation of aqueous sulfur.

Whereas \code{solubility} yields a stable equilibrium condition (the affinities of formation reactions of aqueous species are zero), \code{\link{equilibrate}} generates metastable equilibrium (the affinities of formation reactions are equal to each other, but not necessarily zero).
}

\examples{\dontshow{reset()}

## EXAMPLE 1

# Calculate solubility of a single substance:
# Gaseous SO2 with a given fugacity
# Define basis species (any S-bearing basis species should be first)
basis(c("sulfur", "oxygen", "H2O", "H+"))
basis("pH", 6)
# Load the substances (minerals or gases) to be dissolved
species("sulfur dioxide", -20)
# List the formed aqueous species
# We can use retrieve() to identify all the possible aqueous species
iaq <- retrieve("S", c("O", "H"), "aq")
# Place arguments for affinity() after the first argument of solubility()
s1 <- solubility(iaq, O2 = c(-56, -46), T = 125, in.terms.of = "S")

# Calculate overall solubility for multiple substances:
# Gaseous S2 and SO2 with a given fugacity
basis(c("sulfur", "oxygen", "H2O", "H+"))
basis("pH", 6)
species(c("S2", "sulfur dioxide"), -20)
s2 <- solubility(iaq, O2 = c(-56, -46), T = 125, in.terms.of = "S")

# Make expressions for legends
S_ <- expr.species("SO2", "gas", -20, TRUE)
pH_ <- quote(pH == 6)
T_ <- lT(125)
lexpr <- lex(S_, pH_, T_)

# Make diagrams from the results of solubility calculations
layout(matrix(c(1, 3, 2, 3), nrow = 2))
# Logarithm of activity of aqueous species in equilibrium with SO2(gas)
diagram(s1, type = "loga.equil", ylim = c(-15, 0))
diagram(s1, col = 3, lwd = 3, add = TRUE)
legend("topright", legend = lexpr, bty = "n")
# Logarithm of concentration (parts per million) of aqueous species
sppm <- convert(s1, "logppm")
diagram(sppm, type = "loga.equil", ylim = c(-10, 5))
diagram(sppm, col = 3, lwd = 3, add = TRUE)
legend("topright", legend = lexpr, bty = "n")
par(xpd = NA)
text(-58, 6.5, paste("Solubility of gaseous SO2 (green line) is",
  "sum of concentrations of aqueous species"), cex = 1.2, font = 2)
par(xpd = FALSE)

# Show overall (minimum) solubility of multiple gases
diagram(s2, col = 4, lwd = 3)
# Show solubilities of individual gases
names <- info(species()$ispecies)$formula
diagram(s2, type = "loga.equil", names = names, add = TRUE)
title("Minimum solubility (blue line) corresponds to the most stable gas")
layout(matrix(1))

## EXAMPLE 2

## Two ways to calculate pH-dependent solubility of calcite
## with ionic strength determination
## Method 1: CO2 and carbonate species as formed species
basis(c("CO2", "Ca+2", "H2O", "O2", "H+"))
species("calcite")
iaq <- info(c("CO2", "HCO3-", "CO3-2"))
# Ionic strength calculations don't converge below around pH = 3
sa0 <- solubility(iaq, pH = c(4, 14), dissociate = TRUE)
saI <- solubility(iaq, pH = c(4, 14), dissociate = TRUE, find.IS = TRUE)
## Method 2: CO2 and carbonate species as basis species
basis(c("Ca+2", "CO2", "H2O", "O2", "H+"))
species("calcite")
iaq <- info("Ca+2")
bases <- c("CO2", "HCO3-", "CO3-2")
sm0 <- solubility(iaq, bases = bases, pH = c(4, 14), dissociate = TRUE)
smI <- solubility(iaq, bases = bases, pH = c(4, 14), dissociate = TRUE, find.IS = TRUE)
## Plot the results
plot(0, 0, xlab="pH", ylab="solubility, log mol", xlim = c(4, 14), ylim = c(-5, 2))
# Method 1 with/without ionic strength
lines(saI$vals[[1]], saI$loga.balance, lwd = 5, col = "lightblue")
lines(sa0$vals[[1]], sa0$loga.balance, lwd = 5, col = "pink")
# Method 2 with/without ionic strength
lines(smI$vals[[1]], smI$loga.balance, lty = 2)
lines(sm0$vals[[1]], sm0$loga.balance, lty = 2)
legend("topright", c("I = 0", "I = calculated", "mosaic method"),
       col = c("pink", "lightblue", "black"), lwd = c(5, 5, 1), lty = c(1, 1, 2))
title(main = "Solubility of calcite: Ionic strength and mosaic method")
# The two methods give nearly equivalent results
stopifnot(all.equal(sa0$loga.balance, sm0$loga.balance))
stopifnot(all.equal(saI$loga.balance, smI$loga.balance, tolerance = 0.003))
## NOTE: the second method (using mosaic) is slower, but is
## more flexible; e.g. complexes with Ca+2 could be included
}

\concept{Main workflow}
