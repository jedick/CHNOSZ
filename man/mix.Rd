\encoding{UTF-8}
\name{mix}
\alias{mash}
\alias{rebalance}
\alias{mix}
\title{Combine Diagrams for Multi-Metal Systems}
\description{
  Combine diagrams for different systems by mashing or rebalancing two diagrams or mixing two diagrams with a third.
}

\usage{
  mash(d1, d2)
  rebalance(d1, d2, balance = NULL)
  mix(d1, d2, d3, parts = c(1, 1), .balance = NULL)
}

\arguments{
  \item{d1}{list, output of \code{\link{diagram}} for the first mono-metallic system}
  \item{d2}{list, output of \code{diagram} for the second mono-metallic system}
  \item{balance}{character or numeric, specification of secondary balancing coefficients}
  \item{d3}{list, output of \code{diagram} for the bimetallic system}
  \item{parts}{numeric, amount of each metal (i.e. fixed composition) for the mixed system}
  \item{.balance}{\emph{argument for internal use only}}
}

\details{

These functions make a new \code{\link{affinity}} object from the output of \code{\link{diagram}}.
The result can be used to make a new diagram that shows the combined system.

\code{mash} creates a set of intersecting predominance fields for all possible combinations of species in \code{d1} and \code{d2}.
The new names are formed from the \code{names} used in the source diagrams; for example if "Cp" and "Py" are predominant minerals at the same position in diagrams 1 and 2, the field for the mashed diagram will be labeled "Cp+Py".
The affinities are calculated by summing the formation reactions from the two diagrams to give equal parts of the balancing coefficients in \code{d1} and \code{d2} (that is, equal parts of two different metals).
Note that the actual values of the affinities (and therefore the ratio between the metals) doesn't affect the resulting diagram because the affinities are assigned values of -Inf wherever one of the species is not predominant in the respective single-metal diagram.

\code{mix} is an expanded form of \code{mash} that allows combinations not only between two single-metal diagrams (\code{d1} and \code{d2}) but also between each of those diagrams and third diagram for bimetallic species (\code{d3}).
All combinations of species in all crosses between the diagrams (\code{d1-d2}, \code{d1-d3}, \code{d2-d3}, \code{d3-d3}) are identified.
The mole fractions of species in each combination are computed to satisfy the ratio of metals defined in \code{parts}.
For example, if \code{d1} and \code{d2} are balanced on Fe\S{+2} and VO\s{4}\S{-3}, the species are combined by default to give equal parts of Fe and V.
Note that pairs of bimetallic species in \code{d3} are included as well as single bimetallic species that satisfy the composition in \code{parts} (e.g. FeV for \code{c(1, 1)} or Fe\s{3}V for \code{c(3, 1)}).

From the possible combinations of species, combinations are removed that have a negative mole fraction of any species or that involve any mono-metallic species that has no predominance field in the corresponding single-metal diagram.
The output consists of each unique combination of species, including the combined formation reactions and affinities (in the \code{species} and \code{values} elements of the output list), 
The affinities are assigned values of -Inf wherever one of the species is not predominant in the respective single-metal diagram.
Therefore, either the single-metal diagrams (\code{d1} or \code{d2}) can be recovered by setting \code{parts} to \code{c(1, 0)} or \code{c(0, 1)}, respectively.

NOTE: Unlike the \code{diagram} calls used to make \code{d1} and \code{d2}, which by themselves should produce reasonable diagrams for a single-metal system, the \code{d3} diagram by itself probably has no useful interpretation.
It is only used in \code{mix} as a way to transmit the results of \code{\link{affinity}} for the bimetmal system and the formatted names that are made by \code{diagram}.

\code{rebalance} creates a new set of affinities of reactions \emph{between} species in both systems.
Diagrams for different systems likely use different primary balancing coefficients, such as balancing on different metals.
\code{rebalance} uses \emph{secondary} balancing coefficients, specified acording to \code{balance} (see \code{\link{equilibrate}} for a description of this argument), to determine the reactions between the species in the two systems.
The affinities of these reactions are then used \emph{only} to identify the predominant species at each grid point.
The \emph{returned} value of affinity are carried forward from those used to make the source diagrams (\samp{plotvals} in \code{d1} and \code{d2}), and therefore reflect the primary balancing coefficients.
The returned values are assigned -Inf wherever that species is determined to not predominate according to the secondary balancing.

Because \code{mash} yields finite values of affinity for only a single species at any grid point, the final diagram can be made with any setting of \code{balance}.
\code{mix} gives combinations of species that each have the amount of metals defined in \code{parts}, so it makes no difference whether the final diagram is balanced on either of the metals, or on formula units (\code{balance = 1}).
However, for \code{rebalance}, \code{balance} in the final diagram should be set to \samp{1} to balance on formula units in order to preserve the primary balancing coefficients.

}

\value{
A list object with the same structure as the output from \code{\link{affinity}}, so it can be used as input to \code{diagram}.
}

\seealso{
More examples are in the vignette \viglink{multi-metal}.
}

\examples{
# Don't make all plots because it takes more than
# 10 seconds on CRAN checks 20230313
plot.all <- FALSE
if(plot.all) opar <- par(no.readonly = TRUE)
if(plot.all) par(mfrow = c(2, 2))
# Define basis species with Fe and Cu
basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
xlab <- ratlab("Fe+2", "Cu+")
# Calculate diagram for only Fe-bearing minerals
species(c("pyrite", "pyrrhotite", "magnetite", "hematite"))
aFe <- affinity("Fe+2" = c(0, 12), O2 = c(-40, -16), T = 400, P = 2000)
dFe <- diagram(aFe, xlab = xlab, main = "Fe-S-O-H", plot.it = plot.all)
# Calculate diagram for only Cu-bearing minerals
species(c("covellite", "chalcocite", "tenorite", "cuprite"))
aCu <- affinity(aFe)  # argument recall
dCu <- diagram(aCu, xlab = xlab, main = "Cu-S-O-H", plot.it = plot.all)
### mash() diagram
ac <- mash(dFe, dCu)
diagram(ac, xlab = xlab, main = "Cu-Fe-S-O-H with mash()", plot.it = plot.all)
### rebalance() diagram
ad <- rebalance(dFe, dCu)
diagram(ad, xlab = xlab, balance = 1, main = "Cu-Fe-S-O-H with rebalance()")
db <- describe.basis(3)
leg <- lex(lTP(400, 2000), db)
legend("bottomleft", legend = leg, bty = "n")
if(plot.all) par(opar)}

\concept{Extended workflow}
