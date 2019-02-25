\encoding{UTF-8}
\name{retrieve}
\alias{retrieve}
\title{Retrieve Species by Element}
\description{
Retrieve species in the database containing one or more chemical elements.
}

\usage{
  retrieve(..., state = NULL, add.charge = TRUE,
    hide.electron = TRUE, hide.proton = TRUE, hide.groups = TRUE)
}

\arguments{
  \item{...}{list, one or more arguments, each of which is a character vector with the names of one or more chemical elements}
  \item{state}{character, filter the result on these state(s).}
  \item{add.charge}{logical, add charge to the system?}
  \item{hide.electron}{logical, exclude the electron from the result?}
  \item{hide.proton}{logical, exclude the proton from the result?}
  \item{hide.groups}{logical, exclude groups from the result?}
}

\details{
This function retrieves the species in the thermodynamic database (see \code{\link{thermo}}) that have all of the elements specified in the arguments.
A single argument is interpreted as a combination of one or more elements that must be present in each species.
Results can be filtered on physical state by setting the \code{state} argument.

If more than one argument is present, all of the species identified by each argument are combined, and all species containing any other elements are excluded.
This can be used to retrieve all of the species in the database within a given chemical system.
A chemical system includes charged species if \code{add.charge} is TRUE (the default) or the user supplies the \dQuote{element} of charge (\samp{Z}).

The electron and proton are excluded from the result unless \code{hide.electron} and/or \code{hide.proton} is FALSE.
Groups used in group-additivity calculations, which have names with square brackets (e.g. [-CH2-]), are also excluded unless \code{hide.groups} is FALSE.

The return value is a named numeric vector giving the species index (i.e. rownumber(s) of \code{thermo()$obigt}) with names corresponding to the chemical formulas of the species.
However, if electron is included in the result, the result uses its name (\samp{e-}) instead of its chemical formula (\samp{(Z-1)}).
If the argument list is empty, then the function returns an empty (length 0) numeric value.
A special argument value \samp{all} can be used to retrieve all species in the thermodynamic database, including filtering on state and hiding of the proton, electron, and/or groups.

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

# all species that have Au
retrieve("Au")
# all species that have both Au and Cl
retrieve(c("Au", "Cl"))
# all species that have Au and/or Cl, including charged species, but no other elements
retrieve("Au", "Cl")
# include the electron with these species
retrieve("Au", "Cl", hide.electron = FALSE)
# all uncharged species that have Au and/or Cl
retrieve("Au", "Cl", add.charge = FALSE)

# minerals in the system SiO2-MgO-CaO-CO2
retrieve("Si", "Mg", "Ca", "C", "O", state="cr")
}

\concept{Extended workflow}
