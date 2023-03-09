\encoding{UTF-8}
\name{nonideal}
\alias{nonideal}
\alias{bgamma}
\title{Activity Coefficients of Aqueous Species}
\description{
Calculate activity coefficients and adjusted molal properties of aqueous species.
}

\usage{
  nonideal(species, speciesprops, IS, T, P, A_DH, B_DH,
           m_star = NULL, method = thermo()$opt$nonideal)
  bgamma(TC, P, showsplines = "")
}

\arguments{
  \item{species}{name of method to use, or names or indices of species for which to calculate nonideal properties}
  \item{speciesprops}{list of dataframes of species properties}
  \item{IS}{numeric, ionic strength(s) used in nonideal calculations, mol kg\eqn{^{-1}}{^-1}}
  \item{T}{numeric, temperature (K)}
  \item{P}{numeric, pressure (bar); required for B-dot or b_gamma equation}
  \item{A_DH}{numeric, A Debye-Huckel coefficient; required for B-dot or b_gamma equation}
  \item{B_DH}{numeric, B Debye-Huckel coefficient; required for B-dot or b_gamma equation}
  \item{m_star}{numeric, total molality of all dissolved species}
  \item{method}{character, \samp{Bdot}, \samp{Bdot0}, \samp{bgamma}, \samp{bgamma0}, or \samp{Alberty}}
  \item{TC}{numeric, temperature (\degC)}
  \item{showsplines}{character, show isobaric (\samp{T}) or isothermal (\samp{P}) splines}
}

\details{
\code{nonideal} calculates activity coefficients and adjusted thermodynamic properties for charged and neutral aqueous species.
At the user level, the main use of this function is to set the method for activity coefficient calculations that gets used by other functions in CHNOSZ.
See the \dQuote{Charged Species} section for a description of the available methods.
Activity coefficient calculations are activated by setting the \code{IS} argument of \code{\link{subcrt}} or \code{\link{affinity}}.
Those functions then call \code{nonideal} with all the arguments needed to perform the calculations.
}

\section{Charged Species}{
The default is to not apply calculations for the proton (\Hplus) and electron (\eminus); this makes sense if you are setting the pH, i.e. activity of \Hplus, to some value.
To apply the calculations to H+ and/or e-, change \code{\link{thermo}()$opt$ideal.H} or \code{ideal.e} to FALSE (see examples).

For the \samp{Alberty} method, the values of \code{IS} are combined with Alberty's (2003) equation 3.6-1 (extended Debye-Hückel equation with an empirical term valid up to 0.25 M ionic strength) and its derivatives (Alberty, 2001), to calculate adjusted molal properties at the specified ionic strength(s) and temperature(s).
The calculations use the equation for the Debye-Hückel constant given by Clarke and Glew, 1980, which is valid between 0 and 150 °C at saturated water vapor pressure (\Psat).

For the \samp{Bdot} method (the default), the \dQuote{B-dot} form of the extended Debye-Hückel equation is used.
This equation is valid at ionic strengths up to approximately 3 mol / kg (Hörbrand et al., 2018).
The distance of closest approach for different ions (the \dQuote{ion size parameter}) is taken from \code{\link{thermo}()$Bdot_acirc}; any species not listed in this file is assigned a value of 4.5 \AA.
The extended term parameter for NaCl-dominated solutions, known as \dQuote{B-dot}, is calculated as a function only of temperature (Helgeson, 1969).
To set the extended term parameter to zero, use the \samp{Bdot0} method.

For the \samp{bgamma} method, the \dQuote{b_gamma} equation is used.
The distance of closest approach is set to a constant (3.72e-8 cm) (e.g., Manning et al., 2013).
The extended term parameter is calculated by calling the \code{bgamma} function.
Alternatively, set the extended term parameter to zero with \samp{bgamma0}.
}

\section{Neutral Species}{
For neutral species, the Setch{\eacute}now equation is used, as described in Shvarov and Bastrakov, 1999.
If \code{\link{thermo}()$opt$Setchenow} is \samp{bgamma0} (the default), the extended term parameter is set to zero and the only non-zero term is the mole fraction to molality conversion factor (using the value of \code{m_star}).
If \code{thermo()$opt$Setchenow} is \samp{bgamma}, the extended term paramter is taken from the setting for the charged species (which can be either \samp{Bdot} or \samp{bgamma}).
Set \code{thermo()$opt$Setchenow} to any other value to disable the calculations for neutral species.
}

\section{Additional Details}{
This information, about the arguments and return values used to perform the calculations, is not normally needed by the user (but the usage is shown in the last example).

For \code{nonideal}, the species can be identified by name or species index in \code{species}.
\code{speciesprops} is a list of dataframes containing the input standard molal properties; normally, at least one column is \samp{G} for Gibbs energy.
The function calculates the *adjusted* properties for given ionic strength (\code{IS}); they are equal to the *standard* values only at IS=0.
The adjusted molal properties that can be calculated include \samp{G}, and (currently only with the Alberty method) \samp{H}, \samp{S} and \samp{Cp}; values of any columns with other names are left untouched.
The lengths of \code{IS} and \code{T} supplied in the arguments should be equal to the number of rows of each dataframe in \code{speciesprops}, or length one to use single values throughout.

In addition to \code{IS} and \code{T}, the \samp{Bdot} and \samp{bgamma} methods depend on values of \code{P}, \code{A_DH}, \code{B_DH}, and \code{m_star} given in the arguments.
\code{m_star}, the total molality of all dissolved species, is used to compute the mole fraction to molality conversion factor.
If \code{m_star} is NULL, it is taken to be equal to \code{IS}, which is an underestimate.
For these methods, \samp{G} is currently the only adjusted molal property that is calculated (but this can be used by \code{\link{subcrt}} to calculate adjusted equilibrium constants).

The return value is the same as the input in \code{speciesprops}, except the input standard thermodynamic properties (at IS=0) are replaced by adjusted properties (at higher IS).
For all affected species, a column named \code{loggam} (common (base-10) logarithm of gamma, the activity coefficient) is appended to the output dataframe of species properties.

\code{bgamma} calculates the extended term parameter (written as b_gamma; Helgeson et al., 1981) for activity coefficients in NaCl-dominated solutions at high temperature and pressure.
Data at \Psat and 0.5 to 5 kb are taken from Helgeson (1969, Table 2 and Figure 3) and Helgeson et al. (1981, Table 27) and extrapolated values at 10 to 30 kb from Manning et al. (2013, Figure 11).
Furthermore, the 10 to 30 kb data were used to generate super-extrapolated values at 40, 50, and 60 kb, which may be encountered using the \code{\link{water.DEW}} calculations.
If all \code{P} correspond to one of the isobaric conditions, the values of \code{Bdot} at \code{T} are calculated by spline fits to the isobaric data.
Otherwise, particular (dependent on the \code{T}) isobaric spline fits are themselves used to construct isothermal splines for the given values of \code{T}; the isothermal splines are then used to generate the values of \code{Bdot} for the given \code{P}.
To see the splines, set \code{showsplines} to \samp{T} to make the first set (isobaric splines) along with the data points, or \samp{P} for examples of isothermal splines at even temperature intervals (here, the symbols are not data, but values generated from the isobaric splines).
This is a basic method of interpolating the data without adding any external dependencies.
}

\examples{\dontshow{reset()}
## Each of the available methods
nonideal("Alberty")
nonideal("bgamma0")
nonideal("bgamma")
nonideal("Bdot0")
nonideal("Bdot") # the default

## What's the activity coefficient of Na+ at
## 25 degC and 1 bar and an ionic strength of 0.7?
sres <- subcrt("Na+", T = 25, IS = 0.7)
# Exponentiate to convert log10(gamma) to gamma
10^sres$out[[1]]$loggam
# Now use a different method
nonideal("bgamma")
sres <- subcrt("Na+", T = 25, IS = 0.7)
10^sres$out[[1]]$loggam

## What are activity coefficients of -3, -2, -1, 0, +1, +2, +3 charged species
## as a function of ionic strength and temperature?
# First choose the method
nonideal("Bdot")
# Define the ionic strength and temperature increments
IS <- c(0.001, 0.01, 0.1, 0.7)
T <- seq(0, 100, 25)
# Use species with charges -3, -2, -1, 0, +1, +2, +3
species <- c("PO4-3", "HPO4-2", "H2PO4-", "H3PO4", "Na+", "Ca+2", "Al+3")
# Initialize empty output table for T (rows) and charge (columns)
gamtab <- matrix(nrow = length(T), ncol = length(species))
rownames(gamtab) <- T
colnames(gamtab) <- -3:3
# Make a list of tables to hold the activity coefficients, one for each IS
gamma <- rep(list(gamtab), length(IS))
names(gamma) <- IS
# Loop over the values of ionic strength
for(i in seq_along(IS)) {
  # Calculate properties of species, including logarithm of activity coefficient
  sres <- subcrt(species, T = T, IS = IS[i])
  # Exponentiate to convert log10(gamma) to gamma, and put the values into the tables
  for(j in seq_along(species)) gamma[[i]][, j] <- 10^sres$out[[j]]$loggam
}
# Print the output and make a plot
print(gamma)
matplot(T, gamma$`0.001`, type = "l")
title(main = "activity coefficients of -3, -2, -1, 0, +1, +2, +3 charged species")

## Alberty, 2003 p. 16 Table 1.3: adjusted pKa of acetic acid
## We use the 'IS' argument in subcrt() to calculate adjusted thermodynamic properties
# Set ideal.H to FALSE to calculate activity coefficients for the proton
# (needed for replication of the values in Alberty's book)
nonideal("Alberty")
thermo("opt$ideal.H" = FALSE)
sres <- subcrt(c("acetate", "H+", "acetic acid"), c(-1, -1, 1),
  IS=c(0, 0.1, 0.25), T=25, property="logK")
# We're within 0.01 of Alberty's pK values
Alberty_logK <- c(4.75, 4.54, 4.47)
# The maximum (absolute) pairwise difference between x and y
max(abs(Alberty_logK - sres$out$logK)) # 0.0072
# Reset option to default
thermo("opt$ideal.H" = TRUE)

## An example using IS with affinity():
## Speciation of phosphate as a function of ionic strength
opar <- par(mfrow = c(2, 1))
basis("CHNOPS+")
Ts <- c(25, 100)
species(c("PO4-3", "HPO4-2", "H2PO4-"))
for(T in Ts) {
  a <- affinity(IS = c(0, 0.14), T = T)
  e <- equilibrate(a)
  if(T==25) diagram(e, ylim = c(-3.0, -2.6), legend.x = NULL)
  else diagram(e, add = TRUE, names = FALSE, col = "red")
}
title(main="Non-ideality model for phosphate species")
dp <- describe.property(c("pH", "T", "T"), c(7, Ts))
legend("topright", lty = c(NA, 1, 1), col = c(NA, "black", "red"), legend = dp)
text(0.07, -2.76, expr.species("HPO4-2"))
text(0.07, -2.90, expr.species("H2PO4-"))
## Phosphate predominance f(IS,pH)
a <- affinity(IS = c(0, 0.14), pH = c(6, 13), T = Ts[1])
d <- diagram(a, fill = NULL)
a <- affinity(IS = c(0, 0.14), pH = c(6, 13), T = Ts[2])
d <- diagram(a, add = TRUE, names = FALSE, col = "red")
par(opar)

## Activity coefficients for monovalent ions at 700 degC, 10 kbar
# After Manning, 2013, Fig. 7
# Here we use the b_gamma equation
nonideal("bgamma")
IS <- c(0.001, 0.01, 0.1, 1, 2, 2.79)
# We're above 5000 bar, so need to use IAPWS-95 or DEW
oldwat <- water("DEW")
sres <- subcrt("Na+", T = 700, P = 10000, IS = IS)
water(oldwat)
# Compare the calculated activity coefficient to values from Manning's figure
gamma <- 10^sres$out[[1]]$loggam
Manning_gamma <- c(0.93, 0.82, 0.65, 0.76, 1.28, 2)
gamma - Manning_gamma

## Plot the data and splines used for calculating b_gamma
## (extended term parameter)
bgamma(showsplines = "T")
bgamma(showsplines = "P")

## A longer example, using nonideal() directly
# Alberty, 2003 p. 273-276: activity coefficient (gamma)
# as a function of ionic strength and temperature
nonideal("Alberty")
IS <- seq(0, 0.25, 0.005)
T <- c(0, 25, 40)
lty <- 1:3
species <- c("H2PO4-", "HADP-2", "HATP-3", "ATP-4")
col <- rainbow(4)
thermo.plot.new(xlim = range(IS), ylim = c(0, 1),
  xlab = axis.label("IS"), ylab = "gamma")
for(j in 1:3) {
  # Use subcrt to generate speciesprops
  speciesprops <- subcrt(species, T = rep(T[j], length(IS)))$out
  # Use nonideal to calculate loggamma; this also adjusts G, H, S, Cp,
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

## Reset method to default
nonideal("Bdot")  # or reset()
}

\references{
Alberty, R. A. (2001) Effect of temperature on standard transformed Gibbs energies of formation of reactants at specified pH and ionic strength and apparent equilibrium constants of biochemical reactions. \emph{J. Phys. Chem. B} \bold{105}, 7865--7870. \doi{10.1021/jp011308v}

Alberty, R. A. (2003) \emph{Thermodynamics of Biochemical Reactions}, John Wiley & Sons, Hoboken, New Jersey, 397 p. \url{https://www.worldcat.org/oclc/51242181}

Clarke, E. C. W. and Glew, D. N. (1980) Evaluation of Debye-Hückel limiting slopes for water between 0 and 150 °C. \emph{J. Chem. Soc. Faraday Trans.} \bold{76}, 1911--1916. \doi{10.1039/f19807601911}

Helgeson, H. C. (1969) Thermodynamics of hydrothermal systems at elevated temperatures and pressures. \emph{Am. J. Sci.} \bold{267}, 729--804. \doi{10.2475/ajs.267.7.729}

Helgeson, H. C., Kirkham, D. H. and Flowers, G. C. (1981) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures. IV. Calculation of activity coefficients, osmotic coefficients, and apparent molal and standard and relative partial molal properties to 600\degC and 5 Kb. \emph{Am. J. Sci.} \bold{281}, 1249--1516. \doi{10.2475/ajs.281.10.1249}

Hörbrand, T., Baumann, T. and Moog, H. C. (2018) Validation of hydrogeochemical databases for problems in deep geothermal energy. \emph{Geotherm. Energy} \bold{6}, 20. \doi{10.1186/s40517-018-0106-3}

Manning, C. E. (2013) Thermodynamic modeling of fluid-rock interaction at mid-crustal to upper-mantle conditions. \emph{Rev. Mineral. Geochem.} \bold{76}, 135--164. \doi{10.2138/rmg.2013.76.5}

Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \doi{10.2138/rmg.2013.75.5}

Shvarov, Y. and Bastrakov, E. (1999) HCh: A software package for geochemical equilibrium modelling. User's Guide. \emph{Australian Geological Survey Organisation} \bold{1999/25}. \url{https://pid.geoscience.gov.au/dataset/ga/25473}
}

\concept{Thermodynamic calculations}
