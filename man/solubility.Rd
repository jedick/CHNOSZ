\encoding{UTF-8}
\name{solubility}
\alias{solubility}
\title{Equilibrium Chemical Activities of Species}
\description{
Calculate chemical activities of species in equilibrium with a soluble basis species.
}

\usage{
  solubility(aout, dissociation = NULL, find.IS = FALSE)
}

\arguments{
  \item{aout}{list, output from \code{\link{affinity}}}
  \item{dissociation}{logical, does the mineral undergo a dissociation reaction?}
  \item{find.IS}{logical, find the equilibrium ionic strength by iteration?}
}

\details{
This function performs a simple task: from the values of \code{\link{affinity}} of formation reactions of species at given activity, it works backward to find the activities of species that make the affinities zero.
This corresponds to complete equilibrium with all of the basis species.
Usually, the basis species should be set up so that the first basis species represents the substance being dissolved (a mineral such as CaCO\s3 or gas such as CO\s2).
Internally, this is treated as the conserved basis species, so it must be present in all of the formation reactions of the species.
It is also possible to set the conserved basis species as other than the first one (see \code{demo(gold)}), but this implies that dissociation reactions are not occurring (see below).

The \code{\link{species}} should be defined to represent one set of ions (anions or cations or their complexes) formed in solution, all involving the conserved basis species.
For a dissociation reaction, the second basis species should be used to represent the counterion (cation or anion).

The function perfoms some additional steps to calculate the solubility of something that dissociates (not just dissolves).
For example, the dissolution of calcite (CaCO\s3), involves the release of both calcium ions and different forms of carbonate in solution, depending on the pH.
The equilibrium calculation must take account of the \emph{total} activity of the shared ion (Ca\S{+2}), which is unknown at the start of the calculation.
The solution is accomplished by recalculating the affinities, essentially working backward from the assumption that the dissociation didn't occur.
The resulting activities correspond to equilibrium considering the system-wide activity of Ca\S{+2}.

The function attempts to automatically detect whether dissociation reactions are involved by looking at the first species.
If the formation reaction of that species includes \emph{both} the first and second basis species, the \code{dissociation} flag is set to TRUE.
An example reaction of this type can be found in demo(solubility): CaCO3 (first basis species) = Ca+2 (second basis species) + CO3-2 (first species).
Note that if the conserved basis species is not the first basis species, then the automatic detection of \code{dissociation} will always return FALSE.
Therefore, a reaction corresponding to Au (fourth basis species) + ... = ... gives \code{dissociation} = FALSE (see \code{demo(gold)}).
This algorithm for determining whether dissociation occurs is prone to error, so \code{dissociation} can be explicitly set in the arguments.
A \emph{not recommended} alternative is to set \code{dissociation} to a numeric value corresponding to the stoichiometry of released species (i.e. 2 for a 1:1 electrolyte).
This setting indicates to calculate activities on a per-reaction basis, where each reaction has its own (independent) activity of Ca\S{+2}.
That does not give a complete equilibrium in the system, but may be required to reproduce some published diagrams.

Note that other variables (pH, ionic strength, activities of other basis species) should be defined in the preceding call to \code{\link{affinity}}.
However, for dissolving a substance in pure water, \code{find.IS} can be set to TRUE to determine the final ionic strength.
This works by calculating the ionic strength from the equilibrium solubility calculation, then re-running \code{affinity} with those values.
Note that for dissociation reactions, the ionic strength is calculated from both the ions present in the species definition and the counter ion, which should be the second basis species.
TODO: Speciation of counterions (e.g. ionized forms of carbonate or sulfate) can be accomplished by using the \code{\link{mosaic}} function instead of \code{affinity}.
The calculation is iterated until the ionic strength deviation at every point is lower than a preset tolerance (1e-4).

The output of \code{solubility} has the same format as that of \code{equilibrate}, and can be used by \code{\link{diagram}} with \code{type = "loga.balance"} to plot the solubilities, or with \code{type = NULL} to plot the activities of species. 
}

\section{Warning}{
This function has not been tested for systems that may form dimers or higher-order complexes (such as Au\s{2}S\s{2}\S{2-}).
Except for relatively simple systems, even after careful refinement, the results from CHNOSZ, which considers chemical activities as the independent variables, will not match the results from speciation-solubility (or Gibbs energy minimization) codes, where the system is defined by its bulk composition.
}

\seealso{
\code{demo("solubility")} adds \T-pH diagrams to the CO\s{2} and calcite example here.
\code{demo("gold")} shows solubility calculations for Au in aqueous solutions with hydroxide, chloride, and hydrosulfide complexes.
\code{\link{equilibrate}} calculates equilibrium chemical activities of species given a constant value of \code{loga.balance} (the logarithm of total activity of the conserved basis species).
}

\examples{\dontshow{data(thermo)}
## solubility of CO2 and calcite as a function of pH
opar <- par(mfrow = c(1, 2))
## set pH range and resolution, constant temperature and ionic strength
pH <- c(0, 14)
res <- 100
T <- 25
IS <- 0
## start with CO2 (not a dissociation reaction)
basis(c("carbon dioxide", "H2O", "O2", "H+"))
# ca. atmospheric PCO2
basis("CO2", -3.5)
species(c("CO2", "HCO3-", "CO3-2"))
a <- affinity(pH = c(pH, res), T = T, IS = IS)
s <- solubility(a)
# first plot total activity line
diagram(s, ylim = c(-10, 4), type = "loga.balance", lwd = 4, col = "green2")
# add activities of species
diagram(s, ylim=c(-10, 4), add = TRUE, dy = 1)
# add legend
lexpr <- as.expression(c("total", expr.species("CO2", state = "aq"),
  expr.species("HCO3-"), expr.species("CO3-2")))
legend("topleft", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
  col = c("green2", rep("black", 3)), legend = lexpr)
title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
  list(what = expr.species("CO2"), T = T)), line = 1.5)
mtext("cf. Fig. 4.5 of Stumm and Morgan, 1996")
## now do calcite (a dissociation reaction)
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

## two ways to calculate pH-dependent solubility of calcite
## with ionic strength determination
## method 1: CO2 and carbonate species as formed species
basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
species(c("CO2", "HCO3-", "CO3-2"))
# ionic strength calculations don't converge below around pH=3
a <- affinity(pH = c(3, 14))
sa0 <- solubility(a)
saI <- solubility(a, find.IS = TRUE)
## method 2: CO2 and carbonate species as basis species
basis(c("calcite", "CO2", "H2O", "O2", "H+"))
species(c("Ca+2"))
m <- mosaic(c("CO2", "HCO3-", "CO3-2"), pH = c(3, 14), blend = TRUE)
sm0 <- solubility(m)
smI <- solubility(m, find.IS = TRUE)
## plot the results
plot(0, 0, xlab="pH", ylab="solubility, log mol", xlim = c(3, 14), ylim = c(-5, 2))
# method 1 with/without ionic strength
lines(a$vals[[1]], saI$loga.balance, lwd=5, col="lightblue")
lines(a$vals[[1]], sa0$loga.balance, lwd=5, col="pink")
# method 2 with/without ionic strength
lines(a$vals[[1]], smI$loga.balance, lty=2)
lines(a$vals[[1]], sm0$loga.balance, lty=2)
legend("topright", c("I = 0", "I = calculated", "mosaic method"),
       col = c("pink", "lightblue", "black"), lwd = c(5, 5, 1), lty = c(1, 1, 2))
title(main = "Solubility of calcite: Ionic strength and mosaic method")
# the two methods give nearly equivalent results
stopifnot(all.equal(sa0$loga.balance, sm0$loga.balance))
stopifnot(all.equal(saI$loga.balance, smI$loga.balance, tolerance = 0.003))
## NOTE: the second method (using mosaic) takes longer, but is
## more flexible; e.g. complexes with Ca+2 could be included
}

\references{
Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \url{https://doi.org/10.2138/rmg.2013.75.5}

Stumm, W. and Morgan, J. J. (1996) \emph{Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters}, John Wiley & Sons, New York, 1040 p. \url{http://www.worldcat.org/oclc/31754493}
}

\concept{Main workflow}
