\encoding{UTF-8}
\name{duplex}
\alias{duplex}
\title{Combine Diagrams for Two Metals}
\description{
  Combine diagrams from two systems using a secondary balancing constraint.
}

\usage{
  duplex(d1, d2, balance = NULL)
}

\arguments{
  \item{d1}{list, output of \code{\link{diagram}} for first system}
  \item{d2}{list, output of \code{diagram} for second system}
  \item{balance}{character or numeric, specification of secondary balancing coefficients}
}

\details{

This function makes a new \code{\link{affinity}} object from two diagrams.
Each of the diagrams likely uses different \emph{primary} balancing coefficients (e.g., balancing on different metals).
This function uses \emph{secondary} balancing coefficients to combine the diagrams.

See \code{\link{equilibrate}} for a description of the \code{balance} argument.

}

\value{
A list object with the same structure as the output from \code{\link{affinity}}, so it can be used as input to \code{diagram}.
}

\seealso{
A longer example is in the vignette \viglink{multi-metal}.
}

\examples{
def.par <- par(no.readonly = TRUE)
layout(matrix(c(1, 1, 2, 2, 0, 3, 3, 0), nrow = 4))
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
# Combine the diagrams
ad <- duplex(dFe, dCu)
diagram(ad, xlab = xlab, balance = 1, main = "Cu-Fe-S-O-H")
db <- describe.basis(ibasis = 3)
leg <- lex(lTP(400, 2000), db)
legend("bottomleft", legend = leg, bty = "n")
par(def.par)
}

\concept{Extended workflow}
