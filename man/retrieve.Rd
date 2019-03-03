\encoding{UTF-8}
\name{retrieve}
\alias{retrieve}
\title{Retrieve Species by Element}
\description{
Retrieve species in the database containing one or more chemical elements.
}

\usage{
  retrieve(..., state = NULL, add.charge = TRUE, hide.groups = TRUE,
    req1 = FALSE)
}

\arguments{
  \item{...}{list, one or more arguments, each of which is a character vector with the names of one or more chemical elements}
  \item{state}{character, filter the result on these state(s).}
  \item{add.charge}{logical, add charge to the system?}
  \item{hide.groups}{logical, exclude groups from the result?}
  \item{req1}{logical, include only species that contain the first element in the system?}
}

\details{
This function retrieves the species in the thermodynamic database (see \code{\link{thermo}}) that have all of the elements specified in the arguments.
A single argument is interpreted as a combination of one or more elements that must be present in each species.
Results can be filtered on physical state by setting the \code{state} argument.

If more than one argument is present, all of the species identified by each argument are combined, and all species containing any other elements are excluded.
This can be used to retrieve all of the species in the database within a given chemical system.
A chemical system includes charged species if \code{add.charge} is TRUE (the default) or the user supplies the \dQuote{element} of charge (\samp{Z}).
If \code{req1} is TRUE, the result corresponds to the intersection of all of the species in the system with those identified by the first argument (i.e. those bearing the first element).

Groups used in group-additivity calculations, which have names with square brackets (e.g. [-CH2-]), are excluded unless \code{hide.groups} is FALSE.

The return value is a named numeric vector giving the species index (i.e. rownumber(s) of \code{thermo()$obigt}) with names corresponding to the chemical formulas of the species.
However, if the electron is in the result, its name (\samp{e-}) is used instead of its chemical formula (\samp{(Z-1)}).
If the argument list is empty, then the function returns an empty (length 0) numeric value.
A special argument value \samp{all} can be used to retrieve all species in the thermodynamic database, including filtering on state and hiding of the groups.

To speed up operation, the function uses a precalculated stoichiometric matrix for the default database, which is loaded with the package (see \code{\link{thermo}}).
If the function detects a change to any chemical formulas in database, it triggers a recalculation using \code{\link{i2A}}.
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
# all uncharged species that have Au and/or Cl
retrieve("Au", "Cl", add.charge = FALSE)

# minerals in the system SiO2-MgO-CaO-CO2
retrieve("Si", "Mg", "Ca", "C", "O", state="cr")

# make an Eh-pH diagram for aqueous S-bearing species
basis("CHNOSe")
ispecies <- retrieve("S", "O", "H", req1 = TRUE, state = "aq")
species(ispecies)
a <- affinity(pH = c(0, 14), Eh = c(-1, 1))
diagram(a, fill = "terrain")
}

\concept{Extended workflow}
