\encoding{UTF-8}
\name{add.OBIGT}
\alias{add.OBIGT}
\alias{mod.OBIGT}
\title{Functions to Work with the Thermodynamic Database}
\description{
Add or modify species in the thermodynamic database.
}

\usage{
  add.OBIGT(file, species = NULL, force = TRUE)
  mod.OBIGT(..., zap = FALSE)
}

\arguments{
  \item{file}{character, path to a file}
  \item{species}{character, names of species to load from file}
  \item{force}{logical, force replacement of already existing species?}
  \item{...}{character or numeric, properties of species to modify in the thermodynamic database}
  \item{zap}{logical, clear preexisting parameters?}
}

\details{
Note: change made to OBIGT are lost if you reload the database by calling \code{\link{reset}} or \code{\link{OBIGT}} or if you quit the \R session without saving it.

\code{\link{add.OBIGT}} is used to update the thermodynamic database (\code{\link{thermo}$OBIGT}) in the running session.
The format (column names) of the specified file must be the same as the \code{extdata/OBIGT/*.csv} files provided with CHNOSZ.
Note that this includes both the \code{E_units} and \code{model} columns, which were added in versions 1.3.3 and 2.0.0.

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
The name of the species to add or change must be supplied as the first argument of \code{...} or as a named argument (named \samp{name}).
Additional arguments to \code{mod.OBIGT} refer to the name of the property(s) to be updated and are matched to any part of compound column names in \code{\link{thermo}()$OBIGT}.
For instance, either \samp{z} or \samp{T} matches the \samp{z.T} column.
The values provided should also include any order-of-magnitude scaling factors (see \code{\link{thermo}}).

When adding new species, a chemical formula should be included along with the values of any of the thermodynamic properties.
The formula is taken from the \samp{formula} argument, or if that is missing, is taken to be the same as the \samp{name} of the species.
An error occurs if the formula is not valid (i.e. can not be parsed by \code{\link{makeup}}).
For new species, properties that are not specified become NA, except for \samp{state} and \samp{E_units}, which take default values from \code{thermo()$opt}.
These defaults can be overridden by giving a value for \samp{state} or \samp{E_units} in the arguments.

\samp{model}, if missing, is set to \samp{HKF} for \code{state == "aq"} or \samp{CGL} otherwise.
When modifying some existing minerals in OBIGT, \code{model = "CGL"} should be explicitly given in order to override the Berman model.

When modifying species, the parameters indicated by the named arguments of \code{mod.OBIGT} are updated.
Use \code{zap = TRUE} to replace all prexisting parameters (except for \code{state} and \code{model}) with NA values.
}

\value{
The values returned (\code{\link{invisible}}-y) are the indices of the added and/or modified species.
}

\seealso{ \code{\link{thermo}} (description of OBIGT), \code{\link{mod.buffer}} (modify buffer definitions), \code{\link{logB_to_OBIGT}} (fit thermodynamic parameters to formation constants) }

\examples{\dontshow{reset()}
## Modify an existing species (not real properties)
ialanine <- mod.OBIGT("alanine", state = "cr", G = 0, H = 0, S = 0)
# We have made the values of G, H, and S inconsistent
# with the elemental composition of alanine, so the following 
# now produces a message about that
info(ialanine)

## Add an aqueous species (default) with Gibbs energy given in J/mol 
## (the same as the default) and today's date
date <- as.character(Sys.Date())
iCl2O <- mod.OBIGT("Cl2O", date = date, E_units = "J", G = 87738)
info(iCl2O)

## Add a solid species with a name that is different from the formula
mod.OBIGT("lorem-ipsum", formula = "C123", state = "cr", G = -12345678)
# Retrieve the data for this species using either name or formula
info(info("lorem-ipsum"))
info(info("C123"))
# reset database
OBIGT()

## Using add.OBIGT():
# Compare stepwise stability constants of cadmium chloride complexes
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
matplot(P, do.call(cbind, SSH97), type="l")
matplot(P, do.call(cbind, BZA10), type="l", add=TRUE, lwd=2)
legend("topleft", legend=c("", "", "Sverjensky et al., 1997",
  "Bazarkina et al., 2010"), lwd=c(0, 0, 1, 2), bty="n")
# Make reaction labels
y <- c(1.8, 0.2, -0.5, -1)
invisible(lapply(1:4, function(i) {
  text(800, y[i], describe.reaction(subcrt(c(Cdspecies[i], "Cl-",
    Cdspecies[i+1]), c(-1, -1, 1), T=25, P=1)$reaction))
}))
# Restore default database
OBIGT()

# Another use of add.OBIGT()
# Compare Delta G of AABB = UPBB + H2O
# (Figure 9 of Kitadai, 2014)
# Default database has values from Kitadai, 2014
Kit14 <- subcrt(c("[AABB]", "[UPBB]", "H2O"), c(-1, 1, 1), T = seq(0, 300, 10))
# Load superseded parameters for [UPBB] from Dick et al., 2006
mod.OBIGT("[UPBB]", G = -21436, H = -45220, S = 1.62)
DLH06 <- subcrt(c("[AABB]", "[UPBB]", "H2O"), c(-1, 1, 1), T = seq(0, 300, 10))
xlab <- axis.label("T"); ylab <- axis.label("DG", prefix="k")
plot(Kit14$out$T, Kit14$out$G/1000, type = "l", ylim = c(10, 35),
     xlab = xlab, ylab = ylab)
lines(DLH06$out$T, DLH06$out$G/1000, lty = 2)
legend("topleft", c("Dick et al., 2006", "Kitadai, 2014"), lty = c(2, 1))
title(main = "AABB = UPBB + H2O; after Figure 9 of Kitadai, 2014")
# Restore default database
OBIGT()

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
  ylab = axis.label("DG"), ylim = c(-500, 2500))
lines(T, s2$out$G, lty = 2)
# add title and legend
title(main = describe.reaction(s1$reaction))
stxt <- lapply(c("H4SiO4", "SiO2", "SiO2"), expr.species)
legend("top", legend = as.expression(stxt), bty = "n")
legend("topright", c("Stef\u00e1nsson, 2001", "Shock et al., 1989",
  "Apps and Spycher, 2004"), lty = c(0, 1, 2), bty = "n")
OBIGT()
}

\references{
Apps, J. and Spycher, N. (2004) \emph{Data qualification for thermodynamic data used to support THC calculations}. DOC.20041118.0004 ANL-NBS-HS-000043 REV 00. Bechtel SAIC Company, LLC.

Bazarkina, E. F., Zotov, A. V., and Akinfiev, N. N. (2010)  Pressure-dependent stability of cadmium chloride complexes: Potentiometric measurements at 1-1000 bar and 25°C. \emph{Geology of Ore Deposits} \bold{52}, 167--178. \doi{10.1134/S1075701510020054}

Kitadai, N. (2014)  Thermodynamic prediction of glycine polymerization as a function of temperature and pH consistent with experimentally obtained results. \emph{J. Mol. Evol.} \bold{78}, 171--187. \doi{10.1007/s00239-014-9616-1}

Shock, E. L., Helgeson, H. C. and Sverjensky, D. A. (1989) Calculation of the thermodynamic and transport properties of aqueous species at high pressures and temperatures: Standard partial molal properties of inorganic neutral species. \emph{Geochim. Cosmochim. Acta} \bold{53}, 2157--2183. \doi{10.1016/0016-7037(89)90341-4}

Stef{\aacute}nsson, A. (2001) Dissolution of primary minerals of basalt in natural waters. I. Calculation of mineral solubilities from 0\degC to 350\degC. \emph{Chem. Geol.} \bold{172}, 225--250. \doi{10.1016/S0009-2541(00)00263-1}

Sverjensky, D. A., Shock, E. L., and Helgeson, H. C. (1997) Prediction of the thermodynamic properties of aqueous metal complexes to 1000 °C and 5 kbar. \emph{Geochim. Cosmochim. Acta} \bold{61}, 1359--1412. \doi{10.1016/S0016-7037(97)00009-4}
}

\concept{Thermodynamic data}
