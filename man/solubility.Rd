\encoding{UTF-8}
\name{solubility}
\alias{solubility}
\title{Equilibrium Chemical Activities of Species}
\description{
Calculate chemical activities of species in equilibrium with a soluble basis species.
}

\usage{
  solubility(aout, balance = NULL, split = FALSE)
}

\arguments{
  \item{aout}{list, output from \code{\link{affinity}}}
  \item{balance}{character, basis species to conserve in reactions}
  \item{split}{logical, does the mineral undergo a dissociation reaction?}
}

\details{
This function performs a simple task: from the values of \code{\link{affinity}} of formation reactions of species at given activity, it works backward to find the activities of species that make the affinities zero.
This corresponds to complete equilibrium with all of the basis species.
For solubility calculations, the basis species should be set up so that formation reactions of species are balanced on the thing being dissolved (a mineral such as CaCO\s3 or gas such as CO\s2).
This is usually identified as the first basis species, but if there is any ambiguity, use \code{balance} to indicate the conserved basis species.

\code{split} should be set when calculating the solubility of something that dissociates (not just dissolves).
For example, to calculate the solubility of calcite (CaCO\s3), each form of carbonate is represented by a different reaction, but all of these reactions also release calcium ions.
The equilibrium calculation must take account of the \emph{total} activity of the shared ion (Ca\S{+2}); naturally, that value was unknown for the calculation of \code{\link{affinity}}.
The solution is accomplished by setting \code{split} to TRUE to recalculate the affinities (working backward, as if the split didn't occur).
Then, the resulting activities correspond to equilibrium considering the system-wide activity of Ca\S{+2}.
A \emph{not recommended} alternative is to set \code{split} to a numeric value (probably 2) to calculate activities on a per-reaction basis, where each reaction has its own activity of Ca\S{+2}.
That does not give a complete equilibrium, but may be required to reproduce some published diagrams.

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
# solubility of CO2 and calcite as a function of pH
opar <- par(mfrow = c(1, 2))

# set pH range and resolution, constant temperature and ionic strength
pH <- c(0, 14)
res <- 100
T <- 25
IS <- 0

# start with CO2
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

# now do calcite
basis(c("calcite", "Ca+2", "H2O", "O2", "H+"))
species(c("CO2", "HCO3-", "CO3-2"))
a <- affinity(pH = c(pH, res), T = T, IS = IS)
s <- solubility(a, split = TRUE)
diagram(s, ylim = c(-10, 4), type = "loga.balance", lwd = 4, col = "green2")
diagram(s, add = TRUE, dy = 1)
legend("topright", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
  col = c("green2", rep("black", 3)), legend = lexpr)
title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
  list(what = "calcite", T = T)))
mtext("cf. Fig. 4A of Manning et al., 2013")

par(opar)
}

\references{
Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \url{https://doi.org/10.2138/rmg.2013.75.5}

Stumm, W. and Morgan, J. J. (1996) \emph{Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters}, John Wiley & Sons, New York, 1040 p. \url{http://www.worldcat.org/oclc/31754493}
}

\concept{Main workflow}
