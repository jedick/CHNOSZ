\encoding{UTF-8}
\name{util.expression}
\alias{util.expression}
\alias{expr.species}
\alias{expr.property}
\alias{expr.units}
\alias{axis.label}
\alias{describe.basis}
\alias{describe.property}
\alias{describe.reaction}
\alias{syslab}
\alias{ratlab}
\title{Functions to Express Chemical Formulas and Properties}
\description{Generate expressions suitable for axis labels and plot legends describing chemical species, properties and reactions.}

\usage{
  expr.species(species, state = "aq", value=NULL, log=FALSE, molality=FALSE,
    use.state=FALSE, use.makeup=FALSE)
  expr.property(property, molality = FALSE)
  expr.units(property, prefix = "", per = "mol")
  axis.label(label, units = NULL, basis = get("thermo")$basis, prefix = "",
    molality = FALSE)
  describe.basis(basis = get("thermo")$basis, ibasis = 1:nrow(basis),
    digits = 1, oneline = FALSE, molality = FALSE, use.pH = TRUE)
  describe.property(property, value, digits = 0, oneline = FALSE,
    ret.val = FALSE)
  describe.reaction(reaction, iname = numeric(), states = NULL)
  syslab(system = c("K2O", "Al2O3", "SiO2", "H2O"), dash="-")
  ratlab(ion = "K+", molality = FALSE)
}

\arguments{
  \item{species}{character, formula of a chemical species}
  \item{state}{character, designation of physical state}
  \item{value}{numeric, logarithm of activity or fugacity of species, or value of other property}
  \item{log}{logical, write logarithm of activity/fugacity/molality?}
  \item{molality}{logical, use molality (m) instead of activity (a) for aqueous species?}
  \item{use.state}{logical, include state in expression?}
  \item{use.makeup}{logical, use \code{\link{makeup}} to count the elements?}
  \item{use.pH}{logical, use pH instead of log activity of H+?}
  \item{property}{character, description of chemical property}
  \item{prefix}{character, prefix for units}
  \item{per}{character, denominator in units}
  \item{label}{character, description of species, condition or property}
  \item{units}{character, description of units}
  \item{basis}{data frame, definition of basis species}
  \item{ibasis}{numeric, which basis species to include}
  \item{digits}{numeric, number of digits to show after decimal point}
  \item{oneline}{logical, make descriptions occupy a single line?}
  \item{ret.val}{logical, return only the value with the units?}
  \item{reaction}{data frame, definition of reaction}
  \item{iname}{numeric, show names instead of formulas for these species}
  \item{states}{character, if \samp{all}, show states for all species}
  \item{system}{character, thermodynamic components}
  \item{dash}{character to use for dash between components}
  \item{ion}{character, an ion}
}

\details{

The \code{expr.*} functions create \code{\link{expression}}s using the \code{\link{plotmath}} syntax to describe the names and states and logarithms of activity or fugacity of chemical species, conditions including temperature and pressure and chemical properties such as Gibbs energy and volume.

\code{expr.species} constructs a formatted expression using the formula or name of a single chemical \code{species}.
With no other arguments, the formula is just formatted with the appropriate subscripts and superscripts.
Providing the physical \code{state} adds a variable to the expression (\emph{a} for aqueous species and pure phases, except \emph{f} for gases).
Set \code{molality} to TRUE to write \emph{m} instead of \emph{a} for aqueous species.
The state itself is written in the expression if \code{use.state} is TRUE.
If \code{log} is TRUE, the expression includes a \samp{log} prefix.
Finally, provide a value in \code{value} to write an equation (something like logfO2 = -70), or set it to NA to only write the variable itself (e.g. logfO2).
Set \code{use.makeup} to TRUE to use \code{\link{makeup}} to parse the chemical formula.
This was an older default action that had the undesirable effect of reordering and grouping all the elements, and has been replaced with a different splitting algorithm so that coefficients and charges are sub/superscripted without affecting the intervening text.

\code{expr.property} accepts a description in \code{property} that indicates the chemical property of interest.
Uppercase letters are italicized, and lowercase letters are italicized and subscripted.
Other specific characters are parsed as follows (case-sensitive):

  \tabular{ll}{
     \samp{D} \tab Delta \cr
     \samp{A} \tab bold A (chemical affinity) \cr
     \samp{p} \tab subscript italic P (for isobaric heat capacity) \cr
     \samp{0} \tab degree sign (for a standard-state property) \cr  
     \samp{l} \tab subscript lambda \cr  
     \samp{'} \tab prime symbol \cr  
  }

A \samp{0} gets interpreted as a degree sign only if it does not immediately follow a number (so that e.g. \samp{2.303} can be included in an expression).

Every other character that is one of the \code{\link{letters}} or \code{\link{LETTERS}} in the description of the property is italicized in the expression; other characters such as numerals or mathematical operators are shown without any special formatting.
Special cases for the \code{property} argument (\samp{logK}, \samp{Eh}, \samp{pH}, \samp{pe}, \samp{IS} and \samp{ZC}) are interpreted as simple expressions, and are not parsed according to the above rules.

  \code{expr.units} returns an expression for the units, based on one or more characters appearing in the \code{property}:

  \tabular{ll}{
     \samp{A}, \samp{G}, \samp{H} \tab energy \cr
     \samp{Cp}, \samp{S} \tab energy per Kelvin \cr
     \samp{V} \tab volume \cr
     \samp{E} \tab volume per Kelvin \cr
     \samp{P} \tab pressure \cr
     \samp{T} \tab temperature \cr
     \samp{Eh} \tab electrical potential \cr
     \samp{IS} \tab ionic strength \cr  }

If none of those characters appears in the property, the expression is an empty character (no units). If a \code{prefix} is given, it is added to the expression. The denominator of the units (default \samp{mol}) is taken from the \code{per} argument; it is applied to all units except for \samp{P}, \samp{T}, \samp{Eh}, and \samp{IS}.

\code{axis.label} accepts a generic description of a \code{label}.
If this matches the chemical formula of one of the basis species in the \code{basis} argument, the expression for the label is generated using \code{expr.species} with \code{log} set to the physical state of the basis species.
Otherwise, the expression is built by combining the output of \code{expr.property} with \code{expr.units} (or the value in \code{units}, if it is supplied), placing a comma between the two. 
This function is used extensively in \code{\link{diagram}} and also appears in many of the examples.
Note that \code{\link{diagram}} sets \code{molality} to TRUE if \code{IS} was supplied as an argument to \code{\link{affinity}}.

  \code{describe.basis} makes an expression summarizing the basis species definition (logarithms of activity or fugacity of the basis species) provided in \code{basis}; only the basis species identified by \code{ibasis} are included. 

  \code{describe.property} makes an expression summarizing the properties supplied in \code{property}, along with their \code{value}s. The expressions returned by both functions consist of a property, an equals sign, and a value (with units where appropriate); the expressions have a length equal to the number of property/value pairs. If \code{oneline} is TRUE, the property/value pairs are combined into a single line, separated by commas. The number of digits shown after the decimal point in the values is controlled by \code{digits}. If \code{ret.val} is TRUE, only the values and their units are returned; this is useful for labeling plots with values of temperature.

  \code{describe.reaction} makes an expression summarizing a chemical reaction. The \code{reaction} data frame can be generated using \code{\link{subcrt}}. Based on the sign of their reaction coefficients, species are placed on the reactant (left) or product (right) side of the reaction, where the species with their coefficients are separated by plus signs; the two sides of the reaction are separated by an equals sign. Coefficients equal to 1 are not shown. Chemical formulas of species include a designation of physical state if \code{states} is \samp{all}. Names of species (as provided in \code{reaction}) are shown instead of chemical formulas for the species identified by \code{iname}.

\code{syslab} formats the given thermodynamic components (using \code{expr.species}) and adds intervening en dashes.

\code{ratlab} produces a expression for the activity ratio, viz. (activity of the ion) / [(activity of H+) ^ (charge of the ion)].

}

\seealso{
\code{demo("saturation")} for examples of \code{syslab} and \code{ratlab}.
}

\examples{\dontshow{data(thermo)}
## show descriptions of species and properties on a plot
plot(0, 0, xlim=c(1,5), ylim=c(1,5), xlab="function", ylab="example")
text0 <- function(...) text(..., adj=0)
# species
text0(1, 1, expr.species("CO2"))
text0(1, 2, expr.species("CO2", use.state=TRUE))
text0(1, 3, expr.species("CO2", log=TRUE, use.state=TRUE))
text0(1, 4, expr.species("CO2", log=TRUE))
text0(1, 5, expr.species("CO2", log=TRUE, value=-3))
# properties
text0(2, 1, expr.property("A"))
text0(2, 2, expr.property("DV"))
text0(2, 3, expr.property("DG0f"))
text0(2, 4, expr.property("DCp0,r"))
text0(2, 5, expr.property("T"))
# units
text0(3, 1, expr.units("A", prefix="k"))
text0(3, 2, expr.units("DV"))
text0(3, 3, expr.units("DG0f", prefix="k"))
text0(3, 4, expr.units("DCp0,r"))
text0(3, 5, expr.units("T"))
# axis.label
text0(4, 1, axis.label("DG0f"))
text0(4, 2, axis.label("T"))
text0(4, 3, axis.label("pH"))
text0(4, 4, axis.label("Eh"))
text0(4, 5, axis.label("IS"))
# describe.basis
basis("CHNOS+")
dbasis <- describe.basis(oneline=TRUE, digits=0)
property <- c("P", "T", "Eh", "pH", "IS")
value <- c(1, 42.42, -1, 7, 0.1)
dprop <- describe.property(property, value, oneline=TRUE)
text(3, 1.5, dbasis)
text(3, 2.5, dprop)
dbasis <- describe.basis(ibasis=c(1, 5))
dprop <- describe.property(property[1:2], value[1:2])
legend(2.4, 3.9, legend=c(dbasis, dprop), bty="n")
# describe.reaction
# reaction is automatically balanced since basis species are defined
reaction <- subcrt("glucose", -1)$reaction
text(3, 4.25, describe.reaction(reaction))
text(3, 4.5, describe.reaction(reaction, states="all"))
text(3, 4.75, describe.reaction(reaction, iname=1:4))
title(main="Plot labels for chemical species and thermodynamic properties")
}

\concept{Utility functions}
