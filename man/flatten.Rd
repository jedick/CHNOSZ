\encoding{UTF-8}
\name{flatten}
\alias{flatten}
\alias{duplex}
\title{Combine Diagrams}
\description{
  Combine two diagrams for different systems by flattening them or using a secondary balancing constraint.
}

\usage{
  flatten(d1, d2)
  duplex(d1, d2, balance = NULL)
}

\arguments{
  \item{d1}{list, output of \code{\link{diagram}} for first system}
  \item{d2}{list, output of \code{diagram} for second system}
  \item{balance}{character or numeric, specification of secondary balancing coefficients}
}

\details{

These functions both make a new \code{\link{affinity}} object from two diagrams.
\code{flatten} identifies the intersection of predominance fields for all possible combinations of species (without interaction between the systems), while \code{duplex} creates a new diagram by comparing the affinities of reactions between species in both systems.
Both functions mask the non-predominant species by assigning them -Inf values of affinity, so the result can be used to make a new diagram that shows the combined system.

\code{flatten} makes a simple overlay of two diagrams using new "species" generated from all combinations of those in \code{d1} and \code{d2}.
The new names are formed from the \code{names} used in the source diagrams; for example if "Cp" and "Py" are predominant minerals at the same position in diagrams 1 and 2, the field for the flattened diagram will be labeled "Cp+Py".
The resulting affinities are simply the sum of affinities of the two species; they are assigned values of -Inf wherever one of the species is not predominant in either of the diagrams.

Diagrams for different systems likely use different \emph{primary} balancing coefficients, such as balancing on different metals.
\code{duplex} uses \emph{secondary} balancing coefficients, specified acording to \code{balance} (see \code{\link{equilibrate}} for a description of this argument), to determine the reactions between the species in the two systems.
The affinities of these reactions are then used \emph{only} to identify the predominant species at each grid point.
The \emph{returned} value of affinity are carried forward from those used to make the source diagrams (\samp{plotvals} in \code{d1} and \code{d2}), and therefore reflect the primary balancing coefficients.
The returned values are assigned -Inf wherever that species is determined to not predominate according to the secondary balancing.

Because \code{flatten} yields finite values of affinity for only a single species at any grid point, the final diagram can be made with any setting of \code{balance}.
However, for \code{duplex}, \code{balance} in the final diagram should be set to \samp{1} in order to preserve the primary balancing coefficients.

}

\value{
A list object with the same structure as the output from \code{\link{affinity}}, so it can be used as input to \code{diagram}.
}

\seealso{
A longer example is in the vignette \viglink{multi-metal}.
}

\examples{\dontshow{opar <- par(no.readonly = TRUE)}
par(mfrow = c(2, 2))
# Define basis species with Fe and Cu
basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
xlab <- ratlab("Fe+2", "Cu+")
# Calculate diagram for only Fe-bearing minerals
species(c("pyrite", "pyrrhotite", "magnetite", "hematite"))
aFe <- affinity("Fe+2" = c(0, 12), O2 = c(-40, -16), T = 400, P = 2000)
dFe <- diagram(aFe, xlab = xlab, main = "Fe-S-O-H")
# Calculate diagram for only Cu-bearing minerals
species(c("covellite", "chalcocite", "tenorite", "cuprite"))
aCu <- affinity(aFe)  # argument recall
dCu <- diagram(aCu, xlab = xlab, main = "Cu-S-O-H")
### flatten() diagram
ac <- flatten(dFe, dCu)
diagram(ac, xlab = xlab, main = "Cu-Fe-S-O-H with flatten()")
### duplex() diagram
ad <- duplex(dFe, dCu)
diagram(ad, xlab = xlab, balance = 1, main = "Cu-Fe-S-O-H with duplex()")
db <- describe.basis(ibasis = 3)
leg <- lex(lTP(400, 2000), db)
legend("bottomleft", legend = leg, bty = "n")
\dontshow{par(opar)}}

\concept{Extended workflow}
