\encoding{UTF-8}
\name{solubility}
\alias{solubility}
\title{Equilibrium Chemical Activities of Species}
\description{
Calculate chemical activities of species in equilibrium with a soluble basis species.
}

\usage{
  solubility(eout, exp = 1)
}

\arguments{
  \item{eout}{list, output from \code{\link{equilibrate}}}
  \item{exp}{numeric, exponent characterizing the stoichiometry of the dissociation reaction}
}

\details{
Use this function to calculate solubilities of minerals (such as CaCO\s3) or gases (such as CO\s2).
Start by using \code{\link{equilibrate}} to calculate equilibrium chemical activities of species given a constant value of \code{loga.balance} (the logarithm of total activity of the balanced basis species).
Note that this produces affinities of formation reactions of species that are equal to each other, but are generally not equal to zero.
\code{solubility} adjusts \code{loga.balance} such that the affinities of the formation reaction become zero.
This corresponds to \dQuote{true} equilibrium for a solution in contact with the balanced basis species - i.e. the solubility of that species.

Normally, the balance is automatically identified as the first basis species that is present in all of the species.
If that is not adequate, it can be explicitly set via the \samp{balance} setting in \code{equilibrate}.

The value of \code{exp} should be changed when calculating solubility of species that dissociate (not just dissolve).
For example, set \code{exp} to 2 for calculating the solubility of calcite (CaCO\s3).

The output of \code{solubility} has the same format as that of \code{equilibrate}, and can be used by \code{\link{diagram}} with \code{type = "loga.balance"}. 
}

\section{Warning}{
This function has not been tested for systems that may form dimers or higher-order complexes (such as Au\s{2}S\s{2}\S{2-}).
The lower figures in \code{demo("gold")} are incomplete, as they do not account for other possible reactions not involving Au, particularly the decrease of Cl\S{-} concentration owing to the rising stability of the NaCl\s{(aq)} complex at high temperature.
Except for relatively simple systems, even after careful refinement, the results from CHNOSZ, which considers chemical activities as the independent variables, will not match the results from speciation-solubility (or Gibbs energy minimization) codes, where the system is defined by its bulk composition.
}

\seealso{
\code{demo("solubility")} adds \T-pH diagrams to the CO\s{2} and calcite example here.
\code{demo("gold")} shows solubility calculations for Au in aqueous solutions with hydroxide, chloride, and hydrosulfide complexes.
}

\examples{\dontshow{data(thermo)}
# solubility of CO2 and calcite as a function of pH
par(mfrow = c(1, 2))

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
e <- equilibrate(a)
s <- solubility(e)
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
e <- equilibrate(a)
s <- solubility(e, exp = 2)
diagram(s, ylim = c(-10, 4), type = "loga.balance", lwd = 4, col = "green2")
diagram(s, add = TRUE, dy = 1)
legend("topright", lty = c(1, 1:3), lwd = c(4, 2, 2, 2),
  col = c("green2", rep("black", 3)), legend = lexpr)
title(main = substitute("Solubility of"~what~"at"~T~degree*"C",
  list(what = "calcite", T = T)))
mtext("cf. Fig. 4A of Manning et al., 2013")
}

\references{
Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \url{https://doi.org/10.2138/rmg.2013.75.5}

Stumm, W. and Morgan, J. J. (1996) \emph{Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters}, John Wiley & Sons, New York, 1040 p. \url{http://www.worldcat.org/oclc/31754493}
}

\concept{Main workflow}
