\encoding{UTF-8}
\name{util.formula}
\alias{util.formula}
\alias{as.chemical.formula}
\alias{mass}
\alias{entropy}
\alias{GHS}
\alias{ZC}
\alias{i2A}
\title{Functions to Work with Chemical Formulas}

\description{
Calculate the standard molal entropy of elements in a compound; calculate the standard molal Gibbs energy or enthalpy of formation, or standard molal entropy, from the other two; list coefficients of selected elements in a chemical formula; calculate the average oxidation state of carbon. Create a stoichiometric matrix for selected species.
}

\usage{
  as.chemical.formula(makeup, drop.zero = TRUE)
  mass(formula)
  entropy(formula)
  GHS(formula, G = NA, H = NA, S = NA, T = 298.15, E_units = "cal")
  ZC(formula)
  i2A(formula)
}

\arguments{
  \item{makeup}{numeric, object returned by \code{makeup}}
  \item{drop.zero}{logical, drop elements with a coefficient of zero?}
  \item{formula}{character, chemical formulas, or numeric, rownumbers in \code{thermo()$OBIGT}}
  \item{G}{numeric, standard molal Gibbs energy of formation from the elements}
  \item{H}{numeric, standard molal enthalpy of formation from the elements}
  \item{S}{numeric, standard molal molal entropy}
  \item{T}{numeric, temperature in Kelvin}
  \item{E_units}{character, energy units (cal or J)}
}

\details{

\code{i2A} returns a stoichiometric matrix representing the elemental composition of the \code{formula}s.
Each column corresponds to an element that is present in at least one of the formulas; some element counts will be zero if not all formula have the same elements.
If a matrix is passed to \code{i2A} it is returned unchanged.

\code{as.chemical.formula} makes a character string representing a chemical formula from a vector of coefficients with names corresponding to the elements (e.g., the output of \code{\link{makeup}}) or from a stoichiometric matrix (output of \code{i2A}).
Each elemental symbol is written followed by its coefficient; negative coefficients are signed.
Any coefficients equal to 1 are not explicitly written, and any charge (indicated by \code{makeup} as \samp{Z}) is shown as a signed number at the end of the formula.
If the formula is uncharged, and the last element has a negative coefficient, +0 is shown at the end of the formula to indicate a charge of zero.

The remaining functions documented here accept vectors of chemical formulas, species indices, or a mixture of both, or stoichiometric matrices with elements on the columns.

\code{mass} and \code{entropy} return the sums of masses or entropies of elements in each of the \code{formula}s.
The masses are calculated using the masses of the elements in their natural isotopic distribution, and the entropies, in cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1}, are calculated using the entropies of the compounds of the pure elements in their stable states at 25 \degC and 1 bar.
The properties of the elements used by this function are taken from \code{\link{thermo}$element}.

\code{GHS} computes one of the standard molal Gibbs energy or enthalpy of formation from the elements, or standard molal entropy, from values of the other two.
The \code{formula}, \code{G}, \code{H} and \code{S} arguments must all have the same length.
The entropies of the elements (\code{Se}) in each \code{formula} are calculated using \code{\link{entropy}}, which gives values in calories.
If \code{E_units} is \samp{J}, the values are converted to Joules.
The equation in effect can be written as \eqn{{\Delta}G^{\circ}={\Delta}H^{\circ}-T{\Delta}S^{\circ}}{delta_G = delta_H - T * delta_S}, where \eqn{{\Delta}S^{\circ}=S-S_e}{delta_S = S - Se} and \eqn{T}{T} is the temperature given in \code{T} (defaults to 298.15 K) (note that \code{G} and \code{H} in the arguments correspond respectively to \eqn{{\Delta}G^{\circ}}{delta_G} and \eqn{{\Delta}H^{\circ}}{delta_H} in the equation).
For each formula, if one of \code{G}, \code{H}, or \code{S} is NA, its value is calculated from the other two.
Otherwise, the values are returned unchanged.

\code{ZC} returns the average oxidation state of carbon (\eqn{Z_{\mathrm{C}}}{Z_C}) calculated from ratios of the elements in the chemical \code{formula}s.
The equation used is
\eqn{Z_{\mathrm{C}}=\frac{Z-n_{\mathrm{H}}+2\left(n_{\mathrm{O}}+n_{\mathrm{S}}\right)+3n_{\mathrm{N}}}{n_{\mathrm{C}}}}{Z_C = (Z - nH + 2*(nO + nS) + 3*nN) / nC}, where the \eqn{n}{n} refer to the number of the indicated element in the formula and \eqn{Z}{Z} is the charge (Dick and Shock, 2011).
The result is \code{NaN} for any formula that does not contain carbon. Elements other than those shown in the equation are not included in the calculation, and produce a warning.
}

\value{
\code{mass}, \code{entropy}, and \code{ZC} return numeric values.
\code{as.chemical.formula} returns a character object.
\code{GHS} returns a matrix with column names \samp{G}, \samp{H} and \samp{S}, and \code{i2A} returns a matrix with column names corresponding to the elements in the formulas.
}

\seealso{
\code{\link{makeup}}, used by \code{mass} and \code{entropy}, and \code{ZC} and \code{i2A} for counting the elements in a formula (the latter two make use of the \code{count.zero} argument).
\code{\link{protein.formula}} has an example of computing ZC for proteins compiled from the RefSeq database.
}

\examples{\dontshow{reset()}
## mass and entropy from chemical formulas
mass("H2O")
entropy("H2O")
mass("-1")  # electron
entropy("-1")

## different ways to get the formula of alanine
iA <- info("alanine")
info(iA)$formula
as.chemical.formula(makeup(iA))

## converting among Gibbs energy, enthalpy, entropy
# calculate the value of G from H and S
GHS("H2O", H=water("H"), S=water("S"))[1, ]
# that not quite equal to the value from water("G");
# probably using different entropies of the elements

## Average oxidation states of carbon
ZC(c("CO2", "CH4", "CHNOSZ"))  # 4, -4, 7
si <- info("LYSC_CHICK")
# Can use species index or formula
ZC(si)
ZC(info(si)$formula)

## calculate the chemical formulas, then
## ZC of all of the proteins in CHNOSZ' database
pf <- protein.formula(thermo()$protein)
range(mass(pf))
# use na.rm=TRUE because we have a "protein" with a formula of H2O
range(ZC(pf), na.rm=TRUE)
}

\references{
  Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLOS One} \bold{6}, e22782. \doi{10.1371/journal.pone.0022782}
}

\concept{Thermodynamic calculations}
