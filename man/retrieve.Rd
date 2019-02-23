\encoding{UTF-8}
\name{retrieve}
\alias{retrieve}
\title{Retrieve Species by Element}
\description{
Retrieve species in the database containing one or more chemical elements.
}

\usage{
  retrieve(..., include.electron = FALSE, include.groups = FALSE, state = NULL)
}

\arguments{
  \item{...}{list, one or more arguments, each of which is a character vector with the names of one or more chemical elements}
  \item{include.electron}{logical, include the electron in the result for chemical systems?}
  \item{include.groups}{logical, include groups in the result?}
  \item{state}{character, filter the result on these state(s).}
}

\details{
This function retrieves the species in the thermodynamic database (see \code{\link{thermo}}) that have all of the elements specified in the arguments.
A single argument is interpreted as a combination of one or more elements that must be present in each species.
The return value is a named numeric vector giving the species index (i.e. rownumber(s) of \code{thermo()$obigt}) with names corresponding to the chemical formulas of the species.
If the argument list is empty, then the function returns an empty (length 0) numeric value.

If more than one argument is present, all of the species identified by each argument are combined, then any species containing any other elements are excluded.
This can be used to retrieve all of the species in the database within a given chemical system.
When searching for charged species in a chemical system (using the element named \samp{Z}), the electron is excluded unless \code{include.electron} is TRUE (note that the electron has a chemical formula of \samp{(Z-1)}).

Groups used in group-additivity calculations, which have names with square brackets (e.g. [-CH2-]), are excluded unless \code{include.groups} is TRUE.
Results can be filtered on physical state by setting the \code{state} argument.

The first time the function is run, it uses \code{\link{i2A}} to build the stoichiometric matrix for the current database.
Following runs use the previously calculated stoichiometric matrix, unless a change to the database is detected, which triggers a recalculation.
}

\seealso{
\code{\link{info}}
}

\examples{
# species index of Ti-bearing minerals
retrieve("Ti")
# thermodynamic data for those minerals
info(retrieve("Ti"))

# thermodynamic data for Au-Cl complexes
info(retrieve(c("Au", "Cl")))

# all species that have Au
retrieve("Au")
# all species that have both Au and Cl
retrieve(c("Au", "Cl"))
# all species that have Au and/or Cl, and no other elements
retrieve("Au", "Cl")
# all species that have Au and/or Cl and/or Z, and no other elements
retrieve("Au", "Cl", "Z")
# include the electron with these species
retrieve("Au", "Cl", "Z", include.electron = TRUE)

# minerals in the system SiO2-MgO-CaO-CO2
retrieve("Si", "Mg", "Ca", "C", "O", state="cr")
}

\concept{Extended workflow}
