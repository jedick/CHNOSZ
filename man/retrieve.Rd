\encoding{UTF-8}
\name{retrieve}
\alias{retrieve}
\title{Retrieve Species by Element}
\description{
Retrieve species in the database containing one or more chemical elements.
}

\usage{
  retrieve(elements)
}

\arguments{
  \item{elements}{character, one or more chemical elements}
}

\details{
This function retrieves the species in the thermodynamic database (see \code{\link{thermo}}) that have all of the elements specified in \code{elements}.
The return value is a named numeric vector giving the species index (i.e. rownumber(s) of \code{thermo()$obigt}) with names corresponding to the chemical formulas of the species.

The first time the function is run, it uses \code{\link{i2A}} to build the stoichiometric matrix for the current database.
Following runs use the previously calculated stoichiometric matrix, unless a change to the database is detected, which triggers a recalculation.
}

\seealso{
\code{\link{info}}
}

\examples{
# species index of Ti-bearing minerals
retrieve("Ti")

# thermodynamic data for Au-Cl complexes
info(retrieve(c("Au", "Cl")))
}

\concept{Extended workflow}
