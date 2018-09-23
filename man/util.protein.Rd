\encoding{UTF-8}
\name{util.protein}
\alias{util.protein}
\alias{group.formulas}
\alias{MP90.cp}
\title{Functions for Proteins (Other Calculations)}

\description{
  Return chemical formulas of groups in proteins, and calculate heat capacity using an additivity model from the literature.
}

\usage{
  MP90.cp(protein, T)
  group.formulas()
}

\arguments{
  \item{protein}{proteins specified in any format usable by \code{\link{pinfo}}}
  \item{T}{numeric, temperature in \degC}
}

\details{
\code{group.formulas} returns the chemical formulas of each of the 20 common amino acid residues in proteins, as well as the terminal -H and -H (treated as the [H2O] group).

\code{MP90.cp} takes \code{protein} (name of protein) and \code{T} (one or more temperatures in \degC and returns the additive heat capacity (J mol \eqn{^{-1}}{^-1}) of the unfolded protein using values of heat capacities of the residues taken from Makhatadze and Privalov, 1990.
Those authors provided values of heat capacity at six points between 5 and 125 \degC; this function interpolates (using \code{\link{splinefun}}) values at other temperatures.
}

\seealso{
\code{\link{ionize.aa}} for an example that compares \code{MP90.cp} with heat capacities calculated in CHNOSZ at different temperatures and pHs.
}

\references{
  Makhatadze, G. I. and Privalov, P. L. (1990) Heat capacity of proteins. 1. Partial molar heat capacity of individual amino acid residues in aqueous solution: Hydration effect \emph{J. Mol. Biol.} \bold{213}, 375--384. \url{https://doi.org/10.1016/S0022-2836(05)80197-4}
}

\concept{Protein properties}
