\encoding{UTF-8}
\name{species}
\alias{species}
\title{Species of Interest}
\description{
  Define the species of interest in a system; modify their physical states and logarithms of activities.
}

\usage{
  species(species = NULL, state = NULL, delete = FALSE, add = FALSE,
    index.return = FALSE)
}

\arguments{
  \item{species}{character, names or formulas of species to add to the species definition; numeric, rownumbers of species to modify or delete}
  \item{state}{character, physical states; numeric, logarithms of activities or fugacities}
  \item{delete}{logical, delete the species identified by numeric values of \code{species} (or all species if that argument is missing)?}  
  \item{add}{logical, delete a previous species definition instead of adding to it?}  
  \item{index.return}{logical, return the affected rownumbers of \code{thermo()$species} instead of its contents?}
}

\details{
After defining the \code{\link{basis}} species of your system you can use \code{species} to identify the species of interest.
A species is uniquely identified by a combination of a \code{name} and \code{state}, which are columns of the thermodynamic database in \code{\link{thermo}()$OBIGT}.
For each match of \code{species} to the name, formula, or abbreviation of a species, and of \code{state} to the state (\samp{aq}, \samp{cr}, \samp{gas}, \samp{liq}), the species is added to the current species definition in \code{thermo()$species}.  

The \code{state} argument can be omitted, in which case the first matching species in any state is added (in many cases, this means the aqueous species).
If there are multiple matches for a species, the one that is in the state given by \code{thermo()$opt$state} is chosen, otherwise the matching (or \eqn{n}{n}'th matching duplicate) species is used.
Note that the \code{state}s of species representing phases of minerals that undergo phase transitions are coded as \samp{cr} (lowest-T phase), \samp{cr2}, \samp{cr3}, \code{...} (phases with increasing temperature).
If \code{state} is \samp{cr} when one of these minerals is matched, all the phase species are added.

The data frame in \code{thermo()$species} holds the species names and indices as well as the stoichiometric reaction coefficients for the formation reaction from the basis species and the logarithms of activities or fugacities that are used by \code{\link{affinity}}.
The default values for logarithms of activities are -3 for aqueous (\samp{aq}) species and 0 for others. 

To modify the logarithms of activities of species (logarithms of fugacities for gases) provide one or more numeric values of \code{species} referring to the rownumbers of the species dataframe, or \code{species} NULL, to modify all currently defined species.
The values in \code{state}, if numeric, are interpreted as the logarithms of activities, or if character are interpreted as new states for the species.
If \code{species} is numeric and \code{delete} is \code{TRUE}, these species are deleted from the dataframe; if the only argument is \code{delete} and it is \code{TRUE}, all the species are removed.

By default, when identifying new species, any previous species definition is removed.
Set \code{add} to TRUE to add species to an existing definition.
}

\value{
With no arguments or when adding species, \code{species} returns the value of \code{thermo()$species}, unless \code{index.return} is TRUE, when the function returns the rownumbers of \code{thermo()$species} having the new species.
With \samp{delete=TRUE}, the value is the definition that existed prior the deletion; with \samp{delete=TRUE} and \samp{species} not NULL, the number of species remaining after the selected ones have been deleted, or \code{NULL} if no species remain.
}

\seealso{
  Use \code{\link{info}} to search the thermodynamic database without adding species to the system. \code{\link{basis}} is a prerequisite for \code{\link{species}}.
}

\examples{\dontshow{reset()}
# Set up the basis species
basis("CHNOS")
# Define some aqueous species
species(c("CO2", "NH3"))
# Add some gases
species(c("CO2", "NH3"), "gas", add = TRUE)
# Delete the aqueous species
species(1:2, delete = TRUE)
# Modify the "logact" value
# (log10 of activity for aqueous species;
#  log10 of fugacity for gases)
species(1:2, c(-2, -5))
# Change the state to aqueous
species(1:2, "aq")
# Load a new species definition (deletes the old one first)
species(c("glycine", "alanine"))
# Delete all the species
species(delete = TRUE)

# changing the elements in the basis definition
# causes species to be deleted
basis(c("CaO", "CO2", "H2O", "SiO2", "MgO", "O2"))
species(c("dolomite", "quartz", "calcite", "forsterite"))
basis(c("CO2", "H2O", "O2"))
species()  # NULL
}

\concept{Main workflow}
