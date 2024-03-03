\encoding{UTF-8}
\name{EOSregress}
\alias{EOSregress}
\alias{EOSvar}
\alias{EOScalc}
\alias{EOSplot}
\alias{EOSlab}
\alias{EOScoeffs}
\alias{Cp_s_var}
\alias{V_s_var}
\title{Regress equations-of-state parameters for aqueous species}
\description{
Fit experimental volumes and heat capacities using regression equations.
Possible models include the Helgeson-Kirkham-Flowers (HKF) equations of state, or other equations defined using any combination of terms derived from the temperature, pressure and thermodynamic and electrostatic properties of water.
}

\usage{
  EOSregress(exptdata, var = "", T.max = 9999, ...)
  EOSvar(var, T, P, ...)
  EOScalc(coefficients, T, P, ...)
  EOSplot(exptdata, var = NULL, T.max = 9999, T.plot = NULL, 
    fun.legend = "topleft", coefficients = NULL, add = FALSE,
    lty = par("lty"), col=par("col"), ...)
  EOSlab(var, coeff = "")
  EOScoeffs(species, property, P=1)
  Cp_s_var(T = 298.15, P = 1, omega.PrTr = 0, Z = 0)
  V_s_var(T = 298.15, P = 1, omega.PrTr = 0, Z = 0)
}

\arguments{
  \item{exptdata}{dataframe, experimental data}
  \item{var}{character, name(s) of variables in the regression equations}
  \item{T.max}{numeric, maximum temperature for regression, in degrees Kelvin}
  \item{T}{numeric, temperature in Kelvin}
  \item{P}{numeric, pressure in bars}
  \item{...}{arguments specifying additional dependencies of the regression variables}
  \item{T.plot}{numeric, upper limit of temperature range to plot}
  \item{fun.legend}{character, where to place legend on plot}
  \item{coefficients}{dataframe, coefficients to use to make line on plot}
  \item{add}{logical, add lines to an existing plot?}
  \item{lty}{line style}
  \item{col}{color of lines}
  \item{coeff}{numeric, value of equation of state parameter for plot legend}
  \item{species}{character, name of aqueous species}
  \item{property}{character, \samp{Cp} or \samp{V}}
  \item{omega.PrTr}{numeric, value of omega at reference T and P}
  \item{Z}{numeric, charge}
}

\details{
\code{EOSregress} uses a linear model (\code{\link{lm}}) to regress the experimental heat capacity or volume data in \code{exptdata}, which is a data frame with columns \samp{T} (temperature in degrees Kelvin), \samp{P} (pressure in bars), and \samp{Cp} or \samp{V} (heat capacity in cal/mol.K or volume in cm3/mol).
The \samp{Cp} or \samp{V} data must be in the third column.
Only data below the temperature of \code{T.max} are included in the regression.
The regression formula is specified by a vector of names in \code{var}.
The names of the variables can be any combination of the following (listed in the order of search): variables listed in the following table, any available property of \code{\link{water}} (e.g. \samp{V}, \samp{alpha}, \samp{QBorn}), or the name of a function that can be found using \code{\link{get}} in the default environment.
Examples of the latter are \code{Cp_s_var}, \code{V_s_var}, or functions defined by the user in the global environment; the arguments of these functions must include, but are not limited to, \code{T} and \code{P}.

  \tabular{ll}{
    \code{T}           \tab  \eqn{T}{T} (temperature) \cr
    \code{P}           \tab  \eqn{P}{P} (pressure) \cr
    \code{TTheta}      \tab  \eqn{(T-\Theta)}{(T-Theta)} (\eqn{\Theta}{Theta} = 228 K) \cr
    \code{invTTheta}   \tab  \eqn{1/(T-\Theta)}{1/(T-Theta)} \cr
    \code{TTheta2}     \tab  \eqn{(T-\Theta)^2}{(T-Theta)^2} \cr
    \code{invTTheta2}  \tab  \eqn{1/(T-\Theta)^2}{1/(T-Theta)^2} \cr
    \code{invPPsi}     \tab  \eqn{1/(P+\Psi)}{1/(P+Psi)} (\eqn{\Psi}{Psi} = 2600 bar) \cr
    \code{invPPsiTTheta}   \tab  \eqn{1/((P+\Psi)(T-\Theta))}{1/((P+Psi)(T-Theta))} \cr
    \code{TXBorn}          \tab  \eqn{TX}{TX} (temperature times \eqn{X}{X} Born function) \cr
    \code{drho.dT}     \tab  \eqn{d\rho/dT}{drho/dT} (temperature derivative of density of water) \cr
    \code{V.kT}        \tab  \eqn{V\kappa_T}{V.kT} (volume times isothermal compressibility of water) 
  }


\code{EOSvar} calculates the value of the variable named \code{var} (defined as described above) at the specified \code{T} (temperature in degrees Kelvin) and \code{P} (pressure in bars).
This function is used by \code{EOSregress} to get the values of the variables used in the regression.

\code{EOScalc} calculates the predicted heat capacities or volumes using coefficients provided by the result of \code{EOSregress}, at the temperatures and pressures specified by \code{T} and \code{P}.

\code{EOSplot} takes a table of data in \code{exptdata}, runs \code{EOSregress} and \code{EOScalc} and plots the results.
The experimental data are plotted as points, and the calculated values as a smooth line.
The point symbols are filled circles where the calculated value is within 10\% of the experimental value; open circles otherwise.

\code{EOSlab} produces labels for the variables listed above that can be used \code{\link{as.expression}}s in plots.
The value of \code{coeff} is prefixed to the name of the variable (using \code{\link{substitute}}, with a multiplication symbol).
For the properties listed in the table above, and selected properties listed in \code{\link{water}}, the label is formatted using \code{\link{plotmath}} expressions (e.g., with italicized symbols and Greek letters).
If \code{var} is a user-defined function, the function can be given a \samp{label} attribute to provide \code{\link{plotmath}}-style formatting; in this case the appropriate multiplication or division symbol should be specified (see example below).

\code{EOScoeffs} retrieves coefficients in the Helgeson-Kirkham-Flowers equations from the thermodynamic database (\code{\link{thermo}$OBIGT}) for the given aqueous \code{species}.
If the \code{property} is \samp{Cp}, the resulting data frame has column names of \samp{(Intercept)}, \samp{invTTheta2} and \samp{TX}, respectively holding the coefficients \eqn{c_1}{c1}, \eqn{c_2}{c2} and \eqn{\omega}{omega} in the equation \eqn{Cp^\circ = c_1 + c_2/(T-\Theta)^2 + {\omega}TX}{Cp = c1 + c2/(T-Theta)^2 + omega*TX}.
If the \code{property} is \samp{V}, the data frame has column names of \samp{(Intercept)}, \samp{invTTheta} and \samp{Q}, respectively holding the coefficients \eqn{\sigma}{sigma}, \eqn{\xi}{xi} and \eqn{\omega}{omega} in \eqn{V^\circ = \sigma + \xi/(T-\Theta) - {\omega}Q}{V = sigma + xi/(T-Theta) - omega*Q}.
Here, \eqn{\sigma}{sigma} and \eqn{\xi}{xi} are calculated from a1, a2, a3 and a4 in \code{thermo()$OBIGT} at the pressure indicated by \code{P} (default 1 bar).

The original motivation for writing these functions was to explore alternatives or possible modifications to the revised Helgeson-Kirkham-Flowers equations applied to aqueous nonelectrolytes.
As pointed out by Schulte et al., 2001, the functional forms of the equations do not permit retrieving values of the solvation parameter (\eqn{\omega}{omega}) that closely represent the observed trends in both heat capacity and volume at high temperatures (above ca. 200 \degC).

The examples below assume that the \eqn{\omega}{omega} parameter in the HKF functions is a constant (does not depend on T and P), as is appropriate for nonelectrolytes.
For charged species, the variables \code{Cp_s_var} and \code{V_s_var} can be used in the regressions.
They correspond to the solvation contribution to heat capacity or volume, respectively, in the HKF EOS, divided by the value of \eqn{\omega}{omega} at the reference temperature and pressure.
Because these variables are themselves a function of \code{omega.PrTr}, an iterative procedure is needed to perform the regression.

Note that variables \code{QBorn} and \code{V_s_var} are both negated, so that the value of \eqn{\omega}{omega} has its proper sign in the corresponding equations.
}

\value{
For \code{EOSregress}, an object of class \dQuote{lm}.
\code{EOSvar} and \code{EOScalc} both return numeric values.
\code{EOScoeffs} returns a data frame.
}

\seealso{
The vignette \viglink{eos-regress} has more references and examples, including an iterative method to retrieve \code{omega.PrTr}.
}

\references{
  Hnědkovský, L. and Wood, R. H. (1997) Apparent molar heat capacities of aqueous solutions of \CH4, \CO2, \H2S, and \NH3 at temperatures from 304 K to 704 K at a pressure of 28 MPa. \emph{J. Chem. Thermodyn.} \bold{29}, 731--747. \doi{10.1006/jcht.1997.0192}

  Schulte, M. D., Shock, E. L. and Wood, R. H. (1995) The temperature dependence of the standard-state thermodynamic properties of aqueous nonelectrolytes. \emph{Geochim. Cosmochim. Acta} \bold{65}, 3919--3930. \doi{10.1016/S0016-7037(01)00717-7}
}


\examples{
\dontshow{reset()}
## Fit experimental heat capacities of CH4
## using revised Helgeson-Kirkham-Flowers equations
# Read the data from Hnedkovsky and Wood, 1997
f <- system.file("extdata/cpetc/HW97_Cp.csv", package = "CHNOSZ")
d <- read.csv(f)
# Use data for CH4
d <- d[d$species == "CH4", ]
d <- d[, -1]
# Convert J to cal and MPa to bar
d$Cp <- convert(d$Cp, "cal")
d$P <- convert(d$P, "bar")
# Specify the terms in the HKF equations
var <- c("invTTheta2", "TXBorn")
# Perform regression, with a temperature limit
EOSlm <- EOSregress(d, var, T.max = 600)
# Calculate the Cp at some temperature and pressure
EOScalc(EOSlm$coefficients, 298.15, 1)
# Get the database values of c1, c2 and omega for CH4(aq)
CH4coeffs <- EOScoeffs("CH4", "Cp")
## Make plots comparing the regressions
## with the accepted EOS parameters of CH4
opar <- par(mfrow = c(2,2))
EOSplot(d, T.max = 600)
title("Cp of CH4(aq), fit to 600 K")
legend("bottomleft", pch = 1, legend = "Hnedkovsky and Wood, 1997")
EOSplot(d, coefficients = CH4coeffs)
title("Cp from EOS parameters in database")
EOSplot(d, T.max = 600, T.plot = 600)
title("Cp fit to 600 K, plot to 600 K")
EOSplot(d, coefficients = CH4coeffs, T.plot = 600)
title("Cp from EOS parameters in database")
par(opar)

# Continuing from above, with user-defined variables
Theta <- 228  # K 
invTTTheta3 <- function(T, P) (2*T) / (T-T*Theta) ^ 3
invTX <- function(T, P) 1 / T * water("XBorn", T = T, P = P)[,1]
# Print the calculated values of invTTTheta3
EOSvar("invTTTheta3", d$T, d$P)
# Use invTTTheta and invTX in a regression
var <- c("invTTTheta3", "invTX")
EOSregress(d, var)
# Give them a "label" attribute for use in the legend
attr(invTTTheta3, "label") <-
  quote(phantom()\%*\%2 * italic(T) / (italic(T) - italic(T) * Theta) ^ 3)
attr(invTX, "label") <- quote(phantom() / italic(T * X))
# Uncomment the following to make the plot
#EOSplot(d, var)

## Model experimental volumes of CH4
## using HKF equation and an exploratory one
f <- system.file("extdata/cpetc/HWM96_V.csv", package = "CHNOSZ")
d <- read.csv(f)
# Use data for CH4 near 280 bar
d <- d[d$species == "CH4", ]
d <- d[, -1]
d <- d[abs(d$P - 28) < 0.1, ]
d$P <- convert(d$P, "bar")
# The HKF equation
varHKF <- c("invTTheta", "QBorn")
# alpha is the expansivity coefficient of water
varal <- c("invTTheta", "alpha")
opar <- par(mfrow = c(2, 2))
# For both HKF and the expansivity equation,
# we'll fit up to a temperature limit
EOSplot(d, varHKF, T.max = 663, T.plot = 625)
legend("bottomright", pch = 1, legend = "Hnedkovsky et al., 1996")
title("V of CH4(aq), HKF equation")
EOSplot(d, varal, T.max = 663, T.plot = 625)
title("V of CH4(aq), expansivity equation")
EOSplot(d, varHKF, T.max = 663)
title("V of CH4(aq), HKF equation")
EOSplot(d, varal, T.max = 663)
title("V of CH4(aq), expansivity equation")
par(opar)
# Note that the volume regression using the HKF gives
# a result for omega (coefficient on Q) that is
# not consistent with the high-T heat capacities
}

\concept{Extended workflow}
