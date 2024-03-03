\encoding{UTF-8}
\name{util.data}
\alias{util.data}
\alias{thermo.refs}
\alias{check.EOS}
\alias{check.GHS}
\alias{check.OBIGT}
\alias{dumpdata}
\alias{RH2OBIGT}
\title{Functions for checking thermodynamic data}
\description{
Show table of references in a web browser or get individual references for species.
Check self consistency of individual entries in database.
}

\usage{
  thermo.refs(key = NULL, keep.duplicates = FALSE)
  check.EOS(eos, model, prop, return.difference = TRUE)
  check.GHS(ghs, return.difference = TRUE)
  check.OBIGT()
  dumpdata(file)
  RH2OBIGT(compound = NULL, state = "cr", 
    file = system.file("extdata/adds/RH98_Table15.csv", package = "CHNOSZ"))
}

\arguments{
  \item{key}{character, numeric, or list; bibliographic reference key(s)}
  \item{keep.duplicates}{logical, keep duplicated references?}
  \item{eos}{dataframe, equations-of-state parameters in the format of \code{thermo()$OBIGT}}
  \item{model}{character, thermodynamic model (see \code{\link{thermo}})}
  \item{prop}{character, property of interest (\samp{Cp} or \samp{V})}
  \item{return.difference}{logical, return the difference between database and calculated values?}
  \item{ghs}{dataframe, containing G, H and S, in the format of \code{thermo()$OBIGT}}
  \item{file}{character, path to a file}
  \item{compound}{character, name of compound(s) in group additivity calculation}
  \item{state}{character, physical state of species}
}

\details{

\code{thermo.refs} with default arguments uses \code{\link{browseURL}} to display the sources of thermodynamic data in \code{thermo()$refs}, with the URLs in that table showing as hyperlinks in the browser.
Otherwise, if \code{key} is character, the citation information for those reference \code{key}s (including URLs) are returned.
If \code{key} is numeric, the values refer to the species in those rows of \code{thermo()$OBIGT}, and the citation information for each listed reference (\code{thermo()$OBIGT$ref1}, \code{thermo()$OBIGT$ref2}) is returned.
If \code{key} is a list, it is interpreted as the result of a call to \code{subcrt}, and the citation information for each species involved in the calculation is returned.
Only unique references are returned, unless \code{keep.duplicates} is TRUE.
In that case, a single reference for each species is returned, ignoring anything in \code{thermo()$OBIGT$ref2}.

\code{check.EOS} calculates heat capacity (\code{prop = "Cp"}) or volume (\code{prop = "V"}) from equation-of-state parameters at 25 \degC and 1 bar.
\code{check.GHS} calculates G (standard molal Gibbs energy of formation from the elements) from H (standard molal enthalpy of formation) and S (standard molal entropy) at 25 \degC and 1 bar.
The calculated values of Cp, V, or G are then compared with the given values (i.e., database values).
If \code{return.difference} is TRUE (the default), the difference between the database and calculated values is returned.
If \code{return.difference} is FALSE, the difference is compared with a tolerance setting (see below).
If the absolute value of the difference exceeds the tolerance, the function prints a message and returns the calculated value (not the difference) of the property.
If the absolute value of the difference is less than the tolerance, the function returns NA with no message.

For \code{check.EOS}, the thermodynamic parameters should be provided in \code{eos}, which is a data frame with column names in the same format as \code{\link{thermo}$OBIGT}.
For \code{check.GHS}, the data frame should include G, H, S, and the chemical formula of the species.
The default tolerances are 1 J/K.mol or 1 cal/K.mol for Cp (depending on the \code{E_units} for the species), 1 cm3/mol for V, and 100 cal/mol for G.
These can be changed by setting \code{thermo()$opt$Cp.tol}, \code{thermo()$opt$V.tol}, and \code{thermo()$opt$G.tol}.

\code{check.OBIGT} is a function to check self-consistency of each entry in the thermodynamic database, using \code{check.EOS} and \code{check.GHS}.
The output is a table listing only species that exceed at least one of the tolerance limits, giving the species index (rownumber in `thermo()$OBIGT`), species name and state, and DCp, DV and DG, for the calculated differences (only those above the tolerances are given).
Values of DCp and DG are given in the units present in the data files.
This function is used to generate the file found at \code{extdata/thermo/OBIGT_check.csv}.

\code{dumpdata} returns all of the available data, from both the default and optional data files, or writes it to a file if \code{file} is not NULL.
The format is the same as \code{\link{thermo}$OBIGT}, except for a single prepended column named \samp{source}, giving the source of the data (\samp{OBIGT} refers to the default database, and \samp{DEW}, \samp{SLOP98}, and \samp{SUPCRT92} are the optional data files).

  \code{RH2OBIGT} implements a group additivity algorithm for standard molal thermodynamic properties and equations of state parameters of crystalline and liquid organic molecules from Richard and Helgeson, 1998. The names of the \code{compound}s and their physical \code{state} are searched for in the indicated \code{file}, that also contains chemical formulas and group stoichiometries; the names of the groups are stored in the column names of this file, and must be present in \code{\link{thermo}$OBIGT}. The default \code{file} (\code{\link{extdata}/thermo/RH98_Table15.csv}) includes data taken from Table 15 of Richard and Helgeson, 1998 for high molecular weight compounds in \samp{cr}ystalline and \samp{liq}uid states. An error is produced if any of the \code{compound}-\code{state} combinations is not found in the \code{file}, if any of the group names for a given \code{compound}-\code{state} combination is not found in \code{thermo()$OBIGT}, or if the chemical formula calculated from group additivity (with the aid of \code{\link{i2A}} and \code{\link{as.chemical.formula}}) is not identical to that listed in the \code{file}.

}

\value{
The values returned (\code{\link{invisible}}-y) by \code{mod.OBIGT} are the rownumbers of the affected species.
}

\seealso{ \code{\link{thermo}}, \code{\link{add.OBIGT}}, \code{\link{mod.buffer}} }

\examples{\dontshow{reset()}
# Citation information for Helgeson et al., 1998
thermo.refs("HOK+98")
# Two references for alanine
thermo.refs(info("alanine"))
# Three references for species in the reaction
s <- subcrt(c("O2","O2"),c("gas","aq"),c(-1,1))
thermo.refs(s) 
\dontrun{
## Marked dontrun because it opens a browser
# Show the contents of thermo()$refs
thermo.refs()
} 

## Calculate thermodynamic properties of organic compounds
## using group additivity, after Richard and Helgeson, 1998
RH2OBIGT()
}

\references{
Richard, L. and Helgeson, H. C. (1998) Calculation of the thermodynamic properties at elevated temperatures and pressures of saturated and aromatic high molecular weight solid and liquid hydrocarbons in kerogen, bitumen, petroleum, and other organic matter of biogeochemical interest. \emph{Geochim. Cosmochim. Acta} \bold{62}, 3591--3636. \doi{10.1016/S0016-7037(97)00345-1}
}

\concept{Thermodynamic data}
