\encoding{UTF-8}
\name{add.obigt}
\alias{add.obigt}
\alias{mod.obigt}
\alias{today}
\title{Functions to Work with the Thermodynamic Database}
\description{
Add or modify species in the thermodynamic database.
}

\usage{
  add.obigt(file, species = NULL, force = TRUE, E.units = "cal")
  mod.obigt(...)
  today()
}

\arguments{
  \item{file}{character, path to a file}
  \item{species}{character, names of species to load from file}
  \item{force}{logical, force replacement of already existing species?}
  \item{E.units}{character, units of energy, \samp{cal} or \samp{J}}
  \item{...}{character or numeric, properties of species to modify in the thermodynamic database}
}

\details{
\code{\link{add.obigt}} is used to update the thermodynamic database (\code{\link{thermo}$obigt}) in the running session.
The format (column names) of the specified file must be the same as the \code{extdata/OBIGT/*.csv} files provided with CHNOSZ.

\code{file} is first matched against the names of files in the \code{extdata/OBIGT} directory packaged with CHNOSZ.
In this case, the filename suffixes are removed, so \samp{DEW_aq}, \samp{organic_aq}, and \samp{organic_cr} are valid names.
The function also accepts single matches with the state suffix dropped (\samp{DEW} but not \samp{organic}).
If there are no matches to a system file, then \code{file} is interpreted as the path a user-supplied file.

If \code{species} is NULL (default), all species listed in the file are used.
If \code{species} is given and matches the name(s) of species in the file, only those species are added to the database.

By default, species in the file replace any existing species having the same combination of name and state.
Set \code{force} to FALSE to avoid replacing species that are present in (\code{thermo$obigt}).

Given the default setting of \code{E.units}, the function does not perform any unit conversions.
If \code{E.units} is set to \samp{J}, then the thermodynamic parameters are converted from units of Joules to calories, as used in the CHNOSZ database. 

When adding species, there is no attempt made to keep the order of physical states in the database (aq-cr-liq-gas); the function simply adds new rows to the end of \code{thermo}$obigt.
As a result, retrieving the properties of an added aqueous species using \code{\link{info}} requires an explicit \code{state="aq"} argument to that function if a species with the same name is present in one of the cr, liq or gas states.

\code{mod.obigt} changes one or more of the properties of species or adds species to the thermodynamic database.
These changes are lost if you reload the database by calling \code{\link{data}(thermo)} or if you quit the \R session without saving it.
The name of the species to add or change must be supplied as the first argument of \code{...} or as a named argument (named \samp{name}).
When adding new species, a chemical formula should be included along with the values of any of the thermodynamic properties.
The formula is taken from the \samp{formula} argument, or if that is missing, is taken to be the same as the \samp{name} of the species.
An error results if the formula is not valid (i.e. can not be parsed by\code{\link{makeup}}).
Additional arguments refer to the name of the property(s) to be updated and are matched to any part of compound column names in \code{\link{thermo}$obigt}, such as \samp{z} or \samp{T} in \samp{z.T}.
Unless \samp{state} is specified as one of the properties, its value is taken from \code{thermo$opt$state}.
When adding species, properties that are not specified become NA (except for \samp{state}).
The values provided should be in the units specifed in the documentation for the \code{thermo} data object, including any order-of-magnitude scaling factors.

\code{today} returns the current date in the format adopted for \code{thermo$obigt} (inherited from SUPCRT-format data files) e.g. \samp{13.May.12} for May 13, 2012.
}

\value{
The values returned (\code{\link{invisible}}-y) are the indices of the added and/or modified species.
}

\seealso{ \code{\link{thermo}}, \code{\link{util.data}}, \code{\link{mod.buffer}} }

\examples{\dontshow{data(thermo)}
## modify an existing species (example only)
ialanine <- mod.obigt("alanine", state="cr", G=0, H=0, S=0)
# we have made the values of G, H, and S inconsistent
# with the elemental composition of alanine, so the following 
# now produces a message about that
info(ialanine)
## add a species
iCl2O <- mod.obigt("Cl2O", G=20970)
info(iCl2O)
# add a species with a name that is different from the formula
mod.obigt("buckminsterfullerene", formula="C60", state="cr", date=today())
# retrieve the species data (thermodynamic properties in this toy example are NA)
info(info("C60"))
# reset database
data(OBIGT)

# using add.obigt():
# compare stepwise stability constants of cadmium chloride complexes
# using data from Sverjensky et al., 1997 and Bazarkina et al., 2010
Cdspecies <- c("Cd+2", "CdCl+", "CdCl2", "CdCl3-", "CdCl4-2")
P <- c(1, seq(25, 1000, 25))
SSH97 <- lapply(1:4, function(i) {
  subcrt(c(Cdspecies[i], "Cl-", Cdspecies[i+1]),
    c(-1, -1, 1), T=25, P=P)$out$logK
})
file <- system.file("extdata/thermo/BZA10.csv", package="CHNOSZ")
add.obigt(file)
BZA10 <- lapply(1:4, function(i) {
  subcrt(c(Cdspecies[i], "Cl-", Cdspecies[i+1]),
    c(-1, -1, 1), T=25, P=P)$out$logK
})
# reset default database
data(OBIGT)
matplot(P, do.call(cbind, SSH97), type="l")
matplot(P, do.call(cbind, BZA10), type="l", add=TRUE, lwd=2)
legend("topleft", legend=c("", "", "Sverjensky et al., 1997",
  "Bazarkina et al., 2010"), lwd=c(0, 0, 1, 2), bty="n")
# make reaction labels
y <- c(1.8, 0.2, -0.5, -1)
invisible(lapply(1:4, function(i) {
  text(800, y[i], describe.reaction(subcrt(c(Cdspecies[i], "Cl-",
    Cdspecies[i+1]), c(-1, -1, 1), T=25, P=1)$reaction))
}))
}

\references{
Bazarkina, E. F., Zotov, A. V., and Akinfiev, N. N. (2010)  Pressure-dependent stability of cadmium chloride complexes: Potentiometric measurements at 1-1000 bar and 25°C. \emph{Geology of Ore Deposits} \bold{52}, 167--178. \url{https://doi.org/10.1134/S1075701510020054}

Sverjensky, D. A., Shock, E. L., and Helgeson, H. C. (1997) Prediction of the thermodynamic properties of aqueous metal complexes to 1000 °C and 5 kbar. \emph{Geochim. Cosmochim. Acta} \bold{61}, 1359--1412. \url{https://doi.org/10.1016/S0016-7037(97)00009-4}
}

\concept{Thermodynamic data}
