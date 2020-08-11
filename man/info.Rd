\encoding{UTF-8}
\name{info}
\alias{info}
\title{Search the Thermodynamic Database}
\description{
  Search for species by name or formula, retrieve their thermodynamic properties and parameters, and add proteins to the thermodynamic database.
}

\usage{
  info(species = NULL, state = NULL, check.it=TRUE)
}

\arguments{
  \item{species}{character, names or formulas of species, or (for \code{info} only) numeric with same meaning as \code{ispecies}}
  \item{state}{character, physical states of the species}
  \item{check.it}{logical, check GHS and EOS parameters for self-consistency?}
}

\details{
\code{info} is the primary function used for querying the thermodynamic database (\code{\link{thermo}$OBIGT}).
It is often called recursively; first with a character value (or values) for \code{species} indicating the name(s) or formula(s) of the species of interest.
The result of this call is a numeric value, which can be provided as an argument in a second call to \code{info} in order to retrieve a data frame of the thermodynamic properties of the species.

The text of \code{species} is searched in the names, chemical formulas, and abbreviations (in the \samp{abbrv} column) in the thermodynamic database.
If the text of the \code{species} is matched, the index of that species is returned.
If there are multiple matches for the \code{species}, and \code{state} is NULL, the index of first match is returned.
The order of entries in the database is grouped by states in the order \samp{aq}, \samp{cr}, \samp{gas}, \samp{liq}.
Therefore, for substances represented by both aqueous and gaseous species, the index of the aqueous species is returned, unless \code{state} is set to \samp{gas}.
Note that names (not formulas) of inorganic species, such as \samp{oxygen} and \samp{methane}, are used only for the gas.

Names of species including an underscore character are indicative of proteins, e.g. \samp{LYSC_CHICK}.
If the name of a protein is provided to \code{info} and the composition of the protein can be found using \code{\link{pinfo}}, the thermodyamic properties and parameters of the nonionized protein (calculated using amino acid group additivity) are added to the thermodynamic database.
Included in the return value, as for other species, is the index of the protein in the thermodynamic database or \code{NA} if the protein is not found. Names of proteins and other species can be mixed.

If no exact matches are found, \code{info} searches the database for similar names or formulas using \code{\link{agrep}}.
If any of these are found, the results are summarized on the screen, but the function always returns NA in this case.

With a numeric argument, the rows of \code{thermo()$OBIGT} indicated by \code{ispecies} are returned, after removing any order-of-magnitude scaling factors.
If these species are all aqueous or are all not aqueous, the compounded column names used in \code{thermo()$OBIGT} are replaced with names appropriate for the corresponding equations of state.
A missing value of one of the standard molal Gibbs energy (\code{G}) or enthalpy (\code{H}) of formation from the elements or entropy (\code{S}) is calculated from the other two, if available.
If \code{check.it} is TRUE, several checks of self consistency among the thermodynamic properties and parameters are performed using \code{\link{checkGHS}} and \code{\link{checkEOS}}.
}


\seealso{ \code{\link{retrieve}} for searching species by element; \code{\link{check.OBIGT}} for checking self-consistency of each species. }

\examples{
\dontshow{reset()}
## summary of available data
info()

## species information
# search for something named (or whose formula is) "Fe"
si <- info("Fe")
# use the number to get the full entry
info(si)
# show data for the higher-temperature phases
info(si:(si+3))

## dealing with states
# standard order of precedence for names:
# aq > cr > gas > liq
info(c("ethanol", "adenosine"))  # aq, aq
# state argument overrides the default
info(c("ethanol", "adenosine"), state = c("gas", "cr"))
# exceptions: gases have precedence for names of methane and inorganic gases
info(c("methane", "oxygen")) # gas, gas

# formulas default to aqueous species, if available
i1 <- info(c("CH4", "CO2", "CS2", "MgO"))
info(i1)$state  # aq, aq, gas, cr
# state argument overrides the default
i2 <- info(c("CH4", "CO2", "MgO"), "gas") 
info(i2)$state  # gas, gas, NA

## partial name or formula searches
info("ATP")
info("thiol")
info("MgC")
# add an extra character to refine a search
# or to search using terms that have exact matches
info("MgC ")
info("acetate ")
info(" H2O")
}

\concept{Main workflow}
