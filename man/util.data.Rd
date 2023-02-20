\encoding{UTF-8}
\name{util.data}
\alias{util.data}
\alias{thermo.refs}
\alias{checkEOS}
\alias{checkGHS}
\alias{check.OBIGT}
\alias{dumpdata}
\alias{RH2OBIGT}
\title{Functions for Checking Thermodynamic Data}
\description{
Show table of references in a web browser or get individual references for species.
Check self consistency of individual entries in database.
}

\usage{
  thermo.refs(key = NULL, keep.duplicates = FALSE)
  checkEOS(eos, model, prop, ret.diff = FALSE)
  checkGHS(ghs, ret.diff = FALSE)
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
  \item{ret.diff}{logical, return the difference between calculated and tabulated values?}
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

\code{checkEOS} compares heat capacity and volume calculated from equation-of-state parameters with reference (tabulated) values at 25 \degC and 1 bar and prints a message and returns the calculated value if tolerance is exceeded.
The thermodynamic parameters should be provided in \code{eos}, which is a data frame with columns (and column names) in the same format as \code{\link{thermo}$OBIGT}.
The \code{prop}erty can be one of \samp{Cp} or \samp{V}.
The default tolerances, given in \code{thermo()$opt$Cp.tol} and \code{thermo()$opt$V.tol}, are 1 J/K.mol or 1 cal/K.mol for Cp and 1 cm3/mol for V.
If \code{ret.diff} is TRUE, the differences are returned irrespective of their values, and no messages are printed.

\code{checkGHS} compares G (standard molal Gibbs energy of formation from the elements) calculated from H (standard molal enthalpy of formation) and S (standard molal entropy) with reference (tabulated) values of G at 25 \degC and 1 bar.
A message is printed and the calculated difference is returned if it exceeds the value given in \code{thermo()$opt$G.tol}, which has a default value of 100 cal/mol.
The calculation requires that G, H and S, and the chemical formula of the species all be present.

\code{check.OBIGT} is a function to check self-consistency of each entry in the thermodynamic database, using \code{checkEOS} and \code{checkGHS}.
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
# citation information for Helgeson et al., 1998
thermo.refs("HOK+98")
# two references for alanine
thermo.refs(info("alanine"))
# three references for species in the reaction
s <- subcrt(c("O2","O2"),c("gas","aq"),c(-1,1))
thermo.refs(s) 
\dontrun{
## marked dontrun because it opens a browser
# show the contents of thermo()$refs
thermo.refs()
} 

## calculate thermodynamic properties of organic compounds
## using group additivity, after Richard and Helgeson, 1998
RH2OBIGT()
}

\references{
Richard, L. and Helgeson, H. C. (1998) Calculation of the thermodynamic properties at elevated temperatures and pressures of saturated and aromatic high molecular weight solid and liquid hydrocarbons in kerogen, bitumen, petroleum, and other organic matter of biogeochemical interest. \emph{Geochim. Cosmochim. Acta} \bold{62}, 3591--3636. \doi{10.1016/S0016-7037(97)00345-1}
}

\concept{Thermodynamic data}
