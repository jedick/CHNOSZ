\encoding{UTF-8}
\name{nonideal}
\alias{nonideal}
\alias{Bdot}
\title{Activity coefficients of aqueous species}
\description{
Calculate activity coefficients and adjusted (non-ideal) molal properties of aqueous species.
}

\usage{
  nonideal(species, speciesprops, IS, T, P, A_DH, B_DH,
           method=get("thermo")$opt$nonideal)
  Bdot(TC, P, showsplines = "")
}

\arguments{
  \item{species}{names or indices of species for which to calculate nonideal properties}
  \item{speciesprops}{list of dataframes of species properties}
  \item{IS}{numeric, ionic strength(s) used in nonideal calculations, mol kg\eqn{^{-1}}{^-1}}
  \item{T}{numeric, temperature (K)}
  \item{P}{numeric, pressure (bar); required for Helgeson method}
  \item{A_DH}{numeric, A Debye-Huckel coefficient; required for Helgeson method}
  \item{B_DH}{numeric, B Debye-Huckel coefficient; required for Helgeson method}
  \item{method}{character, \samp{Alberty}, \samp{Helgeson}, or \samp{Helgeson0}}
  \item{TC}{numeric, temperature (\degC)}
  \item{showsplines}{character, show isobaric (\samp{T}) or isothermal (\samp{P}) splines}
}

\details{
\code{nonideal} takes a list of dataframes (in \code{speciesprops}) containing the standard molal properties of the identified \code{species}.
The function calculates the *adjusted* properties for given ionic strength (\code{IS}); they are equal to the *standard* values at IS=0.
The function bypasses (leaves unchanged) properties of all species whose charge (determined by the number of Z in their \code{\link{makeup}}) is equal to zero.
The proton (\Hplus) and electron (\eminus) are also bypassed by default; this makes sense if you are setting the pH, i.e. activity of \Hplus, to some value.
To apply the calculations to H+ and/or e-, change \code{\link{thermo}$opt$ideal.H} or \code{ideal.e} to FALSE.
The lengths of \code{IS} and \code{T} supplied in the arguments should be equal to the number of rows of each dataframe in \code{speciesprops}, or length one to use single values throughout.

If \code{method} is \samp{Alberty}, the values of \code{IS} are combined with Alberty's (2003) equation 3.6-1 (extended Debye-Hückel equation) and its derivatives, to calculate adjusted molal properties at the specified ionic strength(s) and temperature(s).
The adjusted molal properties that can be calculated include \samp{G}, \samp{H}, \samp{S} and \samp{Cp}; any columns in the dataframes of \code{speciesprops} with other names are left untouched.

If \code{method} is \samp{Helgeson}, the \dQuote{B-dot} equation is used.
This equation seems to have been originally proposed by Huckel, 1925; parameters were derived for use at high temperature and pressure by Helgeson, 1969; Helgeson et al., 1981; Manning, 2013.
The distance of closest approach (the \dQuote{ion size parameter}) is taken from the UT_SIZES.REF file of the HCh package (Shvarov and Bastrakov, 1992), which is based on Table 2.7 of Garrels and Christ, 1965.
In addition to \code{IS} and \code{T}, this method depends on values of \code{P}, \code{A_DH}, and \code{B_DH} given in the arguments.
The calculation of \dQuote{B-dot}, also used in the equations, is made within \code{nonideal} by calling the \code{Bdot} function.
For some uses, it is desirable to set the \dQuote{B-dot} parameter to zero; this can be done by setting the method to \code{Helgeson0}.
Currently, \samp{G} is the only adjusted molal property that is calculated (but this can be used by \code{\link{subcrt}} to calculate adjusted equilibrium constants).

\code{Bdot} calculates the \dQuote{B-dot} deviation function (Helgeson, 1969) a.k.a. extended term parameter (written as b_gamma; Helgeson et al., 1981) for activity coefficients in NaCl solutions at high temperature and pressure.
Data at \Psat and 0.5 to 5 kb are taken from Helgeson (1969, Table 2 and Figure 3) and Helgeson et al. (1981, Table 27) and extrapolated values at 10 to 30 kb from Manning et al. (2013, Figure 11).
Furthermore, the 10 to 30 kb data were used to generate super-extrapolated values at 40, 50, and 60 kb, which may be encountered using the \code{\link{water.DEW}} calculations.
If all \code{P} correspond to one of the isobaric conditions, the values of \code{Bdot} at \code{T} are calculated by spline fits to the isobaric data.
Otherwise, particular (dependent on the \code{T}) isobaric spline fits are themselves used to construct isothermal splines for the given values of \code{T}; the isothermal splines are then used to generate the values of \code{Bdot} for the given \code{P}.
To see the splines, set \code{showsplines} to \samp{T} to make the first set (isobaric splines) along with the data points, or \samp{P} for examples of isothermal splines at even temperature intervals (here, the symbols are not data, but values generated from the isobaric splines).
This is a crude method of kriging the data, but produces fairly smooth interpolations without adding any external dependencies.
}

\value{
One (\samp{G}) or more (\samp{H}, \samp{S}, \samp{Cp}; currently only with the Alberty method) standard thermodynamic properties (at IS=0) in \code{speciesprops} are replaced by the corresponding adjusted thermodynamic properties (at higher IS).
For all affected species, a column named \code{loggam} (common (base-10) logarithm of gamma, the activity coefficient) is appended to the output dataframe of species properties.
}

\examples{\dontshow{data(thermo)}
### Examples following Alberty, 2003
### (page numbers given below)

## the default method setting is Helgeson;
## change it to Alberty
oldnon <- nonideal("Alberty")

## using nonideal() directly
# p. 273-276: activity coefficient (gamma)
# as a function of ionic strength and temperature
IS <- seq(0, 0.25, 0.005)
T <- c(0, 25, 40)
lty <- 1:3
species <- c("H2PO4-", "HADP-2", "HATP-3", "ATP-4")
col <- rainbow(4)
thermo.plot.new(xlim = range(IS), ylim = c(0, 1),
  xlab = axis.label("IS"), ylab = "gamma")
for(j in 1:3) {
  # use subcrt to generate speciesprops
  speciesprops <- subcrt(species, T = rep(T[j], length(IS)))$out
  # use nonideal to calculate loggamma; this also adjusts G, H, S, Cp,
  # but we don't use them here
  nonidealprops <- nonideal(species, speciesprops, IS = IS, T = convert(T[j], "K"))
  for(i in 1:4) lines(IS, 10^(nonidealprops[[i]]$loggam), lty=lty[j], col=col[i])
}
t1 <- "Activity coefficient (gamma) of -1,-2,-3,-4 charged species"
t2 <- quote("at 0, 25, and 40 "*degree*"C, after Alberty, 2003")
mtitle(as.expression(c(t1, t2)))
legend("topright", lty=c(NA, 1:3), bty="n",
  legend=c(as.expression(axis.label("T")), 0, 25, 40))
legend("top", lty=1, col=col, bty="n",
  legend = as.expression(lapply(species, expr.species)))

## more often, the 'IS' argument of subcrt() is used to compute
## adjusted properties at given ionic strength
# p. 16 Table 1.3: adjusted pKa of acetic acid
# set ideal.H to FALSE to calculate activity coefficients for the proton
# (makes for better replication of the values in Alberty's book)
thermo$opt$ideal.H <<- FALSE
(sres <- subcrt(c("acetate", "H+", "acetic acid"), c(-1, -1, 1),
  IS=c(0, 0.1, 0.25), T=25, property="logK"))
# we're within 0.01 log of Alberty's pK values
Alberty_logK <- c(4.75, 4.54, 4.47)
stopifnot(maxdiff(sres$out$logK, Alberty_logK) < 0.01)
# reset option to default
thermo$opt$ideal.H <<- TRUE

### An example using IS with affinity():
## speciation of phosphate as a function of ionic strength
opar <- par(mfrow=c(2, 1))
basis("CHNOPS+")
Ts <- c(25, 100)
species(c("PO4-3", "HPO4-2", "H2PO4-"))
for(T in Ts) {
  a <- affinity(IS=c(0, 0.14), T=T)
  e <- equilibrate(a)
  if(T==25) diagram(e, ylim=c(-3.0, -2.6), legend.x=NULL)
  else diagram(e, add=TRUE, names=NULL, col="red")
}
title(main="Non-ideality model for phosphate species")
dp <- describe.property(c("pH", "T", "T"), c(7, Ts))
legend("topright", lty=c(NA, 1, 1), col=c(NA, "black", "red"), legend=dp)
text(0.07, -2.76, expr.species("HPO4-2"))
text(0.07, -2.90, expr.species("H2PO4-"))
## phosphate predominance f(IS,pH)
a <- affinity(IS=c(0, 0.14), pH=c(6, 13), T=Ts[1])
d <- diagram(a, fill=NULL)
a <- affinity(IS=c(0, 0.14), pH=c(6, 13), T=Ts[2])
d <- diagram(a, add=TRUE, names=NULL, col="red")
par(opar)


### finished with Alberty equation, let's look at Helgeson
# this is the default setting, but is needed here because
# we set "Alberty" above
nonideal(oldnon) # same as nonideal("Helgeson")

## activity coefficients for monovalent ions at 700 degC, 10 kbar
# after Manning, 2013, Fig. 7
IS <- c(0.001, 0.01, 0.1, 1, 2, 2.79)
# we're above 5000 bar, so need to use IAPWS-95 or DEW
oldwat <- water("DEW")
wprop <- water(c("A_DH", "B_DH"), T=convert(700, "K"), P=10000)
water(oldwat)
# just use an empty table for a single species
speciesprops <- list(data.frame(G=rep(0, length(IS))))
# choose any monovalent species
(nonidealprops <- nonideal("Na+", speciesprops, IS=IS,
  T=convert(700, "K"), P=10000, A_DH=wprop$A_DH, B_DH=wprop$B_DH))
# we get the nonideal Gibbs energy contribution and
# the activity coefficient; check values of the latter
Manning_gamma <- c(0.93, 0.82, 0.65, 0.76, 1.28, 2)
gamma <- 10^nonidealprops[[1]]$loggam
# the error is larger at higher IS
stopifnot(maxdiff(gamma[1], Manning_gamma[1]) < 0.01)
stopifnot(maxdiff(gamma, Manning_gamma) < 0.36)

## data and splines used for calculating B-dot
## (extended term parameter)
Bdot(showsplines = "T")
Bdot(showsplines = "P")
}

\references{
Alberty, R. A. (2003) \emph{Thermodynamics of Biochemical Reactions}, John Wiley & Sons, Hoboken, New Jersey, 397 p. \url{http://www.worldcat.org/oclc/51242181}

Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{http://www.worldcat.org/oclc/517586}

Helgeson, H. C. (1969) Thermodynamics of hydrothermal systems at elevated temperatures and pressures. \emph{Am. J. Sci.} \bold{267}, 729--804. \url{https://doi.org/10.2475/ajs.267.7.729}

Helgeson, H. C., Kirkham, D. H. and Flowers, G. C. (1981) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures. IV. Calculation of activity coefficients, osmotic coefficients, and apparent molal and standard and relative partial molal properties to 600\degC and 5 Kb. \emph{Am. J. Sci.} \bold{281}, 1249--1516. \url{https://doi.org/10.2475/ajs.281.10.1249}

Hückel, E. (1925). The theory of concentrated, aqueous solutions of strong electrolytes. \emph{Physikalische Zeitschrift} \bold{26}, 93--147.

Manning, C. E. (2013) Thermodynamic modeling of fluid-rock interaction at mid-crustal to upper-mantle conditions. \emph{Rev. Mineral. Geochem.} \bold{76}, 135--164. \url{https://doi.org/10.2138/rmg.2013.76.5}

Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \url{https://doi.org/10.2138/rmg.2013.75.5}

Shvarov, Y. and Bastrakov, E. (1999) HCh: A software package for geochemical equilibrium modelling. User's Guide. \emph{Australian Geological Survey Organisation} \bold{1999/25}.
}

\concept{Thermodynamic calculations}
