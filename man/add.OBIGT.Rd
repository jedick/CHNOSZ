\encoding{UTF-8}
\name{add.OBIGT}
\alias{add.OBIGT}
\alias{mod.OBIGT}
\alias{today}
\title{Functions to Work with the Thermodynamic Database}
\description{
Add or modify species in the thermodynamic database.
}

\usage{
  add.OBIGT(file, species = NULL, force = TRUE)
  mod.OBIGT(...)
  today()
}

\arguments{
  \item{file}{character, path to a file}
  \item{species}{character, names of species to load from file}
  \item{force}{logical, force replacement of already existing species?}
  \item{...}{character or numeric, properties of species to modify in the thermodynamic database}
}

\details{
\code{\link{add.OBIGT}} is used to update the thermodynamic database (\code{\link{thermo}$OBIGT}) in the running session.
The format (column names) of the specified file must be the same as the \code{extdata/OBIGT/*.csv} files provided with CHNOSZ.
Howvever, for backwards compatibility, the \code{E_units} column of the input can be missing; then, it is assigned a value of \samp{cal}.

\code{file} is first matched against the names of files in the \code{extdata/OBIGT} directory packaged with CHNOSZ.
In this case, the file suffixes are removed, so \samp{DEW}, \samp{organic_aq}, and \samp{organic_cr} are valid names.
If there are no matches to a system file, then \code{file} is interpreted as the path a user-supplied file.

If \code{species} is NULL (default), all species listed in the file are used.
If \code{species} is given and matches the name(s) of species in the file, only those species are added to the database.

By default, species in the file replace any existing species having the same combination of name and state.
Set \code{force} to FALSE to avoid replacing species that are present in (\code{thermo()$OBIGT}).

When adding (not replacing) species, there is no attempt made to keep the order of physical states in the database (aq-cr-liq-gas); the function simply adds new rows to the end of \code{thermo}$OBIGT.
As a result, retrieving the properties of an added aqueous species using \code{\link{info}} requires an explicit \code{state="aq"} argument to that function if a species with the same name is present in one of the cr, liq or gas states.

\code{mod.OBIGT} changes one or more of the properties of species or adds species to the thermodynamic database.
These changes are lost if you reload the database by calling \code{\link{data}(thermo)} or if you quit the \R session without saving it.
The name of the species to add or change must be supplied as the first argument of \code{...} or as a named argument (named \samp{name}).
When adding new species, a chemical formula should be included along with the values of any of the thermodynamic properties.
The formula is taken from the \samp{formula} argument, or if that is missing, is taken to be the same as the \samp{name} of the species.
An error results if the formula is not valid (i.e. can not be parsed by\code{\link{makeup}}).
Additional arguments refer to the name of the property(s) to be updated and are matched to any part of compound column names in \code{\link{thermo}$OBIGT}, such as \samp{z} or \samp{T} in \samp{z.T}.
Unless \samp{state} is specified as one of the properties, its value is taken from \code{thermo()$opt$state}.
When adding species, properties that are not specified become NA, except for \samp{state}, which takes a default value from \code{thermo()$opt$state}, and \samp{z.T}, which for aqueous species is set to the charge calculated from the chemical formula (otherwise, NA charge for newly added species would trigger the \code{\link{AkDi}} model).
The values provided should have energy units correponding to the current setting (\code{\link{E.units}}), but this can be overridden by giving a value for \samp{E_units} in the new data.
The values provided should also include any order-of-magnitude scaling factors (see \code{\link{thermo}}).

\code{today} returns the current date in the format adopted for \code{thermo()$OBIGT} (inherited from SUPCRT-format data files) e.g. \samp{13.May.12} for May 13, 2012.
}

\value{
The values returned (\code{\link{invisible}}-y) are the indices of the added and/or modified species.
}

\seealso{ \code{\link{thermo}}, \code{\link{util.data}}, \code{\link{mod.buffer}} }

\examples{\dontshow{reset()}
## modify an existing species (example only)
ialanine <- mod.OBIGT("alanine", state="cr", G=0, H=0, S=0)
# we have made the values of G, H, and S inconsistent
# with the elemental composition of alanine, so the following 
# now produces a message about that
info(ialanine)
## add a species
iCl2O <- mod.OBIGT("Cl2O", G=20970)
info(iCl2O)
# add a species with a name that is different from the formula
mod.OBIGT("buckminsterfullerene", formula="C60", state="cr", date=today())
# retrieve the species data (thermodynamic properties in this toy example are NA)
info(info("C60"))
# reset database
OBIGT()

# using add.OBIGT():
# compare stepwise stability constants of cadmium chloride complexes
# using data from Sverjensky et al., 1997 and Bazarkina et al., 2010
Cdspecies <- c("Cd+2", "CdCl+", "CdCl2", "CdCl3-", "CdCl4-2")
P <- c(1, seq(25, 1000, 25))
SSH97 <- lapply(1:4, function(i) {
  subcrt(c(Cdspecies[i], "Cl-", Cdspecies[i+1]),
    c(-1, -1, 1), T=25, P=P)$out$logK
})
file <- system.file("extdata/adds/BZA10.csv", package="CHNOSZ")
add.OBIGT(file)
BZA10 <- lapply(1:4, function(i) {
  subcrt(c(Cdspecies[i], "Cl-", Cdspecies[i+1]),
    c(-1, -1, 1), T=25, P=P)$out$logK
})
# reset default database
OBIGT()
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

# another use of add.OBIGT()
# compare Delta G of AABB = UPBB + H2O
# (Figure 9 of Kitadai, 2014)
E.units("J")
# default database has values from Kitadai, 2014
Kit14 <- subcrt(c("[AABB]", "[UPBB]", "H2O"), c(-1, 1, 1), T = seq(0, 300, 10))
# optional file OldAA has superseded values of [UPBB] from Dick et al., 2006
add.OBIGT("OldAA")
DLH06 <- subcrt(c("[AABB]", "[UPBB]", "H2O"), c(-1, 1, 1), T = seq(0, 300, 10))
xlab <- axis.label("T"); ylab <- axis.label("DG", prefix="k")
plot(Kit14$out$T, Kit14$out$G/1000, type = "l", ylim = c(10, 35),
     xlab = xlab, ylab = ylab)
lines(DLH06$out$T, DLH06$out$G/1000, lty = 2)
legend("topleft", c("Dick et al., 2006", "Kitadai, 2014"), lty = c(2, 1))
title(main = "AABB = UPBB + H2O; after Figure 9 of Kitadai, 2014")
# reset database *and* settings (units)
reset()

# Another use of add.OBIGT(): calculate Delta G of
# H4SiO4 = SiO2 + 2H2O using different data for SiO2.
# first, get H4SiO4 from Stefansson, 2001
add.OBIGT("AS04", "H4SiO4")
T <- seq(0, 350, 10)
s1 <- subcrt(c("H4SiO4", "SiO2", "H2O"), c(-1, 1, 2), T = T)
# now, get SiO2 from Apps and Spycher, 2004
add.OBIGT("AS04", "SiO2")
s2 <- subcrt(c("H4SiO4", "SiO2", "H2O"), c(-1, 1, 2), T = T)
# plot logK from the first and second calculations
plot(T, s1$out$G, type = "l", xlab = axis.label("T"),
  ylab = axis.label("DG"), ylim = c(-100, 600))
lines(T, s2$out$G, lty = 2)
# add title and legend
title(main = describe.reaction(s1$reaction))
stxt <- lapply(c("H4SiO4", "SiO2", "SiO2"), expr.species)
legend("top", legend = as.expression(stxt), bty = "n")
legend("topright", c("Stef\u00e1nsson, 2001", "Shock et al., 1989",
  "Apps and Spycher, 2004"), lty = c(0, 1, 2), bty = "n")
reset()
}

\references{
Apps, J. and Spycher, N. (2004) \emph{Data qualification for thermodynamic data used to support THC calculations}. DOC.20041118.0004 ANL-NBS-HS-000043 REV 00. Bechtel SAIC Company, LLC.

Bazarkina, E. F., Zotov, A. V., and Akinfiev, N. N. (2010)  Pressure-dependent stability of cadmium chloride complexes: Potentiometric measurements at 1-1000 bar and 25°C. \emph{Geology of Ore Deposits} \bold{52}, 167--178. \url{https://doi.org/10.1134/S1075701510020054}

Kitadai, N. (2014)  Thermodynamic prediction of glycine polymerization as a function of temperature and pH consistent with experimentally obtained results. \emph{J. Mol. Evol.} \bold{78}, 171--187. \url{https://doi.org/10.1007/s00239-014-9616-1}

Shock, E. L., Helgeson, H. C. and Sverjensky, D. A. (1989) Calculation of the thermodynamic and transport properties of aqueous species at high pressures and temperatures: Standard partial molal properties of inorganic neutral species. \emph{Geochim. Cosmochim. Acta} \bold{53}, 2157--2183. \url{https://doi.org/10.1016/0016-7037(89)90341-4}

Stef{\aacute}nsson, A. (2001) Dissolution of primary minerals of basalt in natural waters. I. Calculation of mineral solubilities from 0\degC to 350\degC. \emph{Chem. Geol.} \bold{172}, 225--250. \url{https://doi.org/10.1016/S0009-2541(00)00263-1}

Sverjensky, D. A., Shock, E. L., and Helgeson, H. C. (1997) Prediction of the thermodynamic properties of aqueous metal complexes to 1000 °C and 5 kbar. \emph{Geochim. Cosmochim. Acta} \bold{61}, 1359--1412. \url{https://doi.org/10.1016/S0016-7037(97)00009-4}
}

\concept{Thermodynamic data}
