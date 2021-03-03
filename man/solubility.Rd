\encoding{UTF-8}
\name{solubility}
\alias{solubility}
\alias{solubilities}
\title{Equilibrium Chemical Activities of Species}
\description{
Calculate chemical activities of aqueous species in equilibrium with a mineral or gas.
}

\usage{
  solubility(aout, dissociation = NULL, find.IS = FALSE, in.terms.of = NULL,
    codeanal = FALSE)
  solubilities(a_cr, i_aq, in.terms.of = NULL, dissociation = NULL)
}

\arguments{
  \item{aout}{list, output from \code{\link{affinity}}}
  \item{dissociation}{logical, does the mineral undergo a dissociation reaction?}
  \item{find.IS}{logical, find the equilibrium ionic strength by iteration?}
  \item{in.terms.of}{character, express the total solubility in terms of moles of this species}
  \item{codeanal}{logical, print diagnostic messages and values of internal variables for code analysis}
  \item{a_cr}{list, output from \code{\link{affinity}} or \code{\link{mosaic}} for one or more minerals}
  \item{i_aq}{character or numeric, aqueous species that can be formed by dissolving the minerals}
}

\details{
\code{solubility} performs a simple task: given the values of \code{\link{affinity}} of formation reactions of species at specified activity, it works backward to find the activities of species that make the affinities zero.
This corresponds to complete equilibrium with all of the basis species.
Usually, the basis species should be set up so that the first basis species represents the substance being dissolved (a mineral such as CaCO\s3 or gas such as CO\s2).
Internally, this is treated as the conserved basis species, so it must be present in all of the formation reactions of the species.

To generate \code{aout}, the \code{\link{species}} should be defined to represent one set of aqueous species (including ions and/or neutral complexes) formed in solution, all involving the conserved basis species.
For a dissociation reaction, the second basis species should be used to represent the counterion (cation or anion).
Other variables (pH, ionic strength, activities of other basis species) should be defined in the call to \code{\link{affinity}} to make \code{aout}.

A more complex calculation involves setting \code{find.IS} to TRUE in order to determine the final ionic strength of dissolving a substance in pure water.
This works by calculating the ionic strength from the equilibrium solubility calculation, then re-running \code{affinity} with those values.
Note that for dissociation reactions, the ionic strength is calculated from both the ions present in the species definition and the counter ion, which should be the second basis species.
The calculation is iterated until the ionic strength deviation at every point is lower than a preset tolerance (1e-4).
Alternatively, speciation of counterions (e.g. ionized forms of carbonate or sulfate) can also be accomplished by using the \code{\link{mosaic}} function instead of \code{affinity}; this is used in the second example below.

The output of \code{solubility} has the same format as that of \code{equilibrate}, and can be used by \code{\link{diagram}} with \code{type = "loga.balance"} to plot the solubilities, or with \code{type = NULL} to plot the activities of species. 
The value of \code{loga.balance} reflects the activity (or molality) of the conserved basis species, i.e. the thing being dissolved.
Use \code{in.terms.of} to express this value in terms of another species.
For example, the solubility of corundum (Al\s{2}O\s{3}) can be expressed in terms of the moles of Al\S{+3} in solution (see the vignette \viglink{anintro}).

Whereas \code{solubility} is limited to calculations for a single mineral at a time, \code{solubilities} combines the results of solubility calculations for multiple minerals.
\code{a_cr} should contain the affinities of only minerals with any set of basis species, as long as the first one contains the main element.
\code{i_aq} gives the names or species indices of all the aqueous species that all have that element in combination with the other basis species (possibly generated using \code{\link{retrieve}}).
The function calculates the solubilities for each of the minerals present in code{a_cr}.
Then, it takes the \emph{minimum} solubility among all the minerals at each grid location.
This corresponds to the maximum activity of dissolved species.
At this point the solution is saturated (not oversaturated) with respect to the most stable mineral.

The result of \code{solubilities} is stored in the \code{loga.balance} component of the returned list object.
Only this component is used by \code{\link{diagram}} for plotting; the other components, while present, refer to the solubility calculations for only the last mineral.
This function is used for generating isosolubility (aka. equisolubility) lines in \code{demo("Pourbaix")} and \code{demo{"zinc"}}.
The latter demo combines the calculation of solubilities with a \code{\link{mosaic}} calculation to account for the speciation of aqueous sulfur.
}

\section{Dissociation Reactions}{
The function perfoms some additional steps to calculate the solubility of something that dissociates (not just dissolves).
For example, the dissolution of calcite (CaCO\s3), involves the release of both calcium ions and different forms of carbonate in solution, depending on the pH.
The equilibrium calculation must take account of the \emph{total} activity of the shared ion (Ca\S{+2}), which is unknown at the start of the calculation.
The solution is accomplished by recalculating the affinities, essentially working backward from the assumption that the dissociation didn't occur.
The resulting activities correspond to equilibrium considering the system-wide activity of Ca\S{+2}.

The function attempts to automatically detect whether dissociation reactions are involved.
If the formation reactions of all species includes \emph{both} the first and second basis species, the \code{dissociation} flag is set to TRUE.
An example reaction of this type can be found in demo(solubility): CaCO3 (first basis species) = Ca+2 (second basis species) + CO3-2 (first species).
Note that if the conserved basis species is not the first basis species, then the automatic detection of \code{dissociation} will always return FALSE.
Therefore, a reaction corresponding to Au (fourth basis species) + ... = ... gives \code{dissociation} = FALSE (see \code{demo(gold)}).
This algorithm for determining whether dissociation occurs is not perfect, so \code{dissociation} can be explicitly set in the arguments.
A \emph{not recommended} alternative is to set \code{dissociation} to a numeric value corresponding to the stoichiometry of released species (i.e. 2 for a 1:1 electrolyte).
This setting indicates to calculate activities on a per-reaction basis, where each reaction has its own (independent) activity of Ca\S{+2}.
That does not give a complete equilibrium in the system, but may be required to reproduce some published diagrams.
}

\section{Warning}{
This function has not been tested for systems that may form dimers or higher-order complexes (such as Au\s{2}S\s{2}\S{2-}).
Except for relatively simple systems, even after careful refinement, the results from CHNOSZ, which considers chemical activities as the independent variables, will not match the results from speciation-solubility (or Gibbs energy minimization) codes, where the system is defined by its bulk composition.
}

\seealso{
\code{demo("solubility")} shows \T-pH diagrams related to the CO\s{2} and calcite example here.
\code{demo("gold")} shows solubility calculations for Au in aqueous solutions with hydroxide, chloride, and hydrosulfide complexes.
In the latter demo, the conserved basis species is not the first one (this might have been done to prevent the incorrect automatic detection of dissociation reactions).

\code{\link{retrieve}} provides a way to list all of the aqueous species in the database that have the specified elements.
In contrast to \code{solubility}, \code{\link{equilibrate}} calculates equilibrium chemical activities of species given a constant value of \code{loga.balance} (the logarithm of total activity of the conserved basis species).
}

\examples{\dontshow{reset()}
## Solubility of CO2 and calcite as a function of pH
opar <- par(mfrow = c(1, 2))
## Set pH range and resolution, constant temperature and ionic strength
pH <- c(0, 14)
res <- 100
T <- 25
IS <- 0
## Start with CO2 (not a dissociation reaction)
basis(c("carbon dioxide", "H2O", "O2", "H+"))
# This is close to atmospheric PCO2
basis("CO2", -3.5)
species(c("CO2", "HCO3-", "CO3-2"))
a <- affinity(pH = c(pH, res), T = T, IS = IS)
s <- solubility(a)
# First plot total activity line
diagram(s, ylim = c(-10, 4), type = "loga.balance", lwd = 4, col = "green2")
# Add activities of species
diagram(s, ylim=c(-10, 4), add = TRUE, dy = 1)
# Add legend
lexpr <- as.expression(c("total", expr.species("CO2", state = "aq"),
  expr.species("HCO3-"), expr.species("CO3-2")))
legend("topleft", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
  col = c("green2", rep("black", 3)), legend = lexpr)
title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
  list(what = expr.species("CO2"), T = T)), line = 1.5)
mtext("cf. Fig. 4.5 of Stumm and Morgan, 1996")
## Now do calcite (a dissociation reaction)
basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
species(c("CO2", "HCO3-", "CO3-2"))
a <- affinity(pH = c(pH, res), T = T, IS = IS)
s <- solubility(a)
diagram(s, ylim = c(-10, 4), type = "loga.balance", lwd = 4, col = "green2")
diagram(s, add = TRUE, dy = 1)
legend("topright", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
  col = c("green2", rep("black", 3)), legend = lexpr)
title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
  list(what = "calcite", T = T)))
mtext("cf. Fig. 4A of Manning et al., 2013")
par(opar)

## Two ways to calculate pH-dependent solubility of calcite
## with ionic strength determination
## Method 1: CO2 and carbonate species as formed species
basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
species(c("CO2", "HCO3-", "CO3-2"))
# Ionic strength calculations don't converge below around pH = 3
a <- affinity(pH = c(4, 14))
sa0 <- solubility(a)
saI <- solubility(a, find.IS = TRUE)
## Method 2: CO2 and carbonate species as basis species
basis(c("calcite", "CO2", "H2O", "O2", "H+"))
species(c("Ca+2"))
m <- mosaic(c("CO2", "HCO3-", "CO3-2"), pH = c(4, 14))
sm0 <- solubility(m)
smI <- solubility(m, find.IS = TRUE)
## Plot the results
plot(0, 0, xlab="pH", ylab="solubility, log mol", xlim = c(4, 14), ylim = c(-5, 2))
# Method 1 with/without ionic strength
lines(a$vals[[1]], saI$loga.balance, lwd=5, col="lightblue")
lines(a$vals[[1]], sa0$loga.balance, lwd=5, col="pink")
# Method 2 with/without ionic strength
lines(a$vals[[1]], smI$loga.balance, lty=2)
lines(a$vals[[1]], sm0$loga.balance, lty=2)
legend("topright", c("I = 0", "I = calculated", "mosaic method"),
       col = c("pink", "lightblue", "black"), lwd = c(5, 5, 1), lty = c(1, 1, 2))
title(main = "Solubility of calcite: Ionic strength and mosaic method")
# The two methods give nearly equivalent results
stopifnot(all.equal(sa0$loga.balance, sm0$loga.balance))
stopifnot(all.equal(saI$loga.balance, smI$loga.balance, tolerance = 0.003))
## NOTE: the second method (using mosaic) takes longer, but is
## more flexible; e.g. complexes with Ca+2 could be included
}

\references{
Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \url{https://doi.org/10.2138/rmg.2013.75.5}

Stumm, W. and Morgan, J. J. (1996) \emph{Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters}, John Wiley & Sons, New York, 1040 p. \url{https://www.worldcat.org/oclc/31754493}
}

\concept{Main workflow}
