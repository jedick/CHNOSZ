\encoding{UTF-8}
\name{retrieve}
\alias{retrieve}
\title{Retrieve species by element}
\description{
Retrieve species in the database containing one or more chemical elements.
}

\usage{
  retrieve(elements = NULL, ligands = NULL, state = NULL,
    T = NULL, P = "Psat", add.charge = TRUE, hide.groups = TRUE)
}

\arguments{
  \item{elements}{character, combination of elements, or list, elements in a chemical system}
  \item{ligands}{character, elements present in any ligands}
  \item{state}{character, filter the result on these state(s).}
  \item{T}{numeric, temperature where \DG0 of species must be not NA}
  \item{P}{numeric, pressure where \DG0 of species must be not NA}
  \item{add.charge}{logical, add charge to the system?}
  \item{hide.groups}{logical, exclude groups from the result?}
}

\details{
This function retrieves the species in the thermodynamic database (see \code{\link{thermo}}) that have the indicated \code{elements}.
A character value of \code{elements} is interpreted as a combination of one or more elements that must be present in each species.
A list value of \code{elements} is used to specify a chemical system -- the species must contain one or more of the indicated elements, but no other elements.
\code{ligands}, if present, gives the elements that may be present in any ligands; this can be used to retrieve all species in a system bearing the \code{elements} (usually a single metal).

The result includes charged species if \code{add.charge} is TRUE (the default) or the user supplies the \dQuote{element} of charge (\samp{Z}).
Results can be filtered on physical state by setting the \code{state} argument.
Groups used in group-additivity calculations, which have names with square brackets (e.g. [-CH2-]), are excluded unless \code{hide.groups} is FALSE.
A special argument value \samp{all} can be used to retrieve all species in the thermodynamic database, including filtering on state and hiding of the groups.

The return value is a named integer vector giving the species index (i.e. rownumber(s) of \code{thermo()$OBIGT}) with names corresponding to the chemical formulas of the species.
If the electron is in the result, its name (\samp{e-}) is used instead of its chemical formula (\samp{(Z-1)}).
An empty (length 0) integer value is returned if no \code{elements} are specified or no species are retrieved.

Set \code{T} (and optionally \code{P}) to require that species have non-NA values of \DG0 at this temperature and pressure.
These values are passed to \code{\link{subcrt}} and therefore have the units set in \code{\link{T.units}} and \code{\link{P.units}}.

To speed up operation, the function uses a precalculated stoichiometric matrix for the default database, which is loaded with the package (see \code{\link{thermo}}).
If the function detects a change to any chemical formulas in the database, it updates the stoichiometric matrix using \code{\link{i2A}}.
}

\seealso{
\code{\link{info}} for basic database searches; \viglink{anintro} for a diagram made with \code{retrieve}d species in the Mn-O-H system.
}

\examples{
# Species index of Ti-bearing minerals
retrieve("Ti")
# Thermodynamic data for those minerals
info(retrieve("Ti"))

# All species that have Au
retrieve("Au")
# All species that have both Au and Cl
retrieve(c("Au", "Cl"))
# Au-Cl system: species that have Au and/or Cl,
# including charged species, but no other elements
retrieve(list("Au", "Cl"))
# All Au-bearing species in the Au-Cl system
retrieve("Au", "Cl")
# All uncharged Au-bearing species in the Au-Cl system
retrieve("Au", "Cl", add.charge = FALSE)

# Minerals in the system SiO2-MgO-CaO-CO2
retrieve(list("Si", "Mg", "Ca", "C", "O"), state = "cr")
\dontshow{par(thermo()$opar)}
}

\concept{Extended workflow}
