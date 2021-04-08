\encoding{UTF-8}
\name{diagram}
\alias{diagram}
\alias{find.tp}
\title{Chemical Activity Diagrams}
\description{
Plot equilibrium chemical activity (1-D speciation) or equal-activity (2-D predominance) diagrams as a function of chemical activities of basis species, temperature and/or pressure.
}

\usage{
  diagram(
    # species affinities or activities
    eout,
    # type of plot
    type = "auto", alpha = FALSE, normalize = FALSE,
    as.residue = FALSE, balance = NULL, groups = as.list(1:length(eout$values)),
    # figure size and sides for axis tick marks
    xrange = NULL, mar = NULL, yline = par("mgp")[1]+0.3, side = 1:4,
    # axis limits and labels
    ylog = TRUE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
    # character sizes
    cex = par("cex"), cex.names = 1, cex.axis = par("cex"),
    # line styles
    lty = NULL, lty.cr = NULL, lty.aq = NULL, lwd = par("lwd"), dotted = NULL,
    spline.method = NULL, contour.method = "edge", levels = NULL,
    # colors
    col = par("col"), col.names = par("col"), fill = NULL,
    fill.NA = "gray80", limit.water = NULL,
    # field and line labels
    names = NULL, format.names = TRUE, bold = FALSE, italic = FALSE, 
    font = par("font"), family = par("family"), adj = 0.5,
    dx = 0, dy = 0, srt = 0, min.area = 0,
    # title and legend
    main = NULL, legend.x = NA,
    # plotting controls
    add = FALSE, plot.it = TRUE, tplot = TRUE, ...)
  find.tp(x)
}

\arguments{
  \item{eout}{list, object returned by \code{\link{equilibrate}} or \code{\link{affinity}}}
  \item{type}{character, type of plot, or name of basis species whose activity to plot}
  \item{alpha}{logical or character (\samp{balance}), for speciation diagrams, plot degree of formation instead of activities?}
  \item{normalize}{logical, divide chemical affinities by balance coefficients (rescale to whole formulas)?}
  \item{as.residue}{logical, divide chemical affinities by balance coefficients (no rescaling)?}
  \item{balance}{character, balancing constraint; see \code{\link{equilibrate}}}
  \item{groups}{list of numeric, groups of species to consider as a single effective species}
  \item{xrange}{numeric, range of x-values between which predominance field boundaries are plotted}
  \item{mar}{numeric, margins of plot frame}
  \item{yline}{numeric, margin line on which to plot the y-axis name}
  \item{side}{numeric, which sides of plot to draw axes}
  \item{xlim}{numeric, limits of x-axis}
  \item{ylim}{numeric, limits of y-axis}
  \item{xlab}{character, label to use for x-axis}
  \item{ylab}{character, label to use for y-axis}
  \item{ylog}{logical, use a logarithmic y-axis (on 1D degree diagrams)?}
  \item{cex}{numeric, character expansion (scaling relative to current)}
  \item{cex.names}{numeric, character expansion factor to be used for names of species on plots}
  \item{cex.axis}{numeric, character expansion factor for names of axes}
  \item{lty}{numeric, line types to be used in plots}
  \item{lty.cr}{numeric, line types for cr-cr boundaries (between two minerals)}
  \item{lty.aq}{numeric, line types for aq-aq boundaries (between two aqueous species)}
  \item{lwd}{numeric, line width}
  \item{dotted}{numeric, how often to skip plotting points on predominance field boundaries (to gain the effect of dotted or dashed boundary lines)}
  \item{spline.method}{character, method used in \code{\link{splinefun}}}
  \item{contour.method}{character, labelling method used in \code{\link{contour}} (use NULL for no labels).}
  \item{levels}{numeric, levels at which to draw contour lines}
  \item{col}{character, color of activity lines (1D diagram) or predominance field boundaries (2D diagram)}
  \item{col.names}{character, colors for labels of species}
  \item{fill}{character, colors used to fill predominance fields}
  \item{fill.NA}{character, color for grid points with NA values}
  \item{limit.water}{NULL or logical, set NA values beyond water stability limits?}
  \item{names}{character, names of species for activity lines or predominance fields}
  \item{format.names}{logical, apply formatting to chemical formulas?}
  \item{bold}{logical, use bold formatting for names?}
  \item{italic}{logical, use italic formatting for names?}
  \item{font}{character, font type for names (has no effect if \code{format.names} is TRUE)}
  \item{family}{character, font family for names}
  \item{adj}{numeric, adjustment for line labels}
  \item{dx}{numeric, x offset for line or field labels}
  \item{dy}{numeric, y offset for line or field labels}
  \item{srt}{numeric, rotation for line labels}
  \item{min.area}{numeric, minimum area of fields that should be labeled, expressed as a fraction of the total plot area}
  \item{main}{character, a main \code{\link{title}} for the plot; \code{NULL} means to plot no title}
  \item{legend.x}{character, description of legend placement passed to \code{\link{legend}}}
  \item{add}{logical, add to current plot?}
  \item{plot.it}{logical, make a plot?}
  \item{tplot}{logical, set up plot with \code{\link{thermo.plot.new}}?}
  \item{x}{matrix, value of the \code{predominant} list element from \code{diagram}}
  \item{...}{additional arguments passed to \code{\link{plot}} or \code{\link{barplot}}}
}

\details{

This function displays diagrams representing either chemical affinities, or equilibrium chemical activities of species.
The first argument is the output from \code{\link{affinity}}, \code{\link{equilibrate}}, or \code{\link{solubility}}.
0-D diagrams, at a single point, are shown as \code{\link{barplot}}s.
1-D diagrams, for a single variable on the x-axis, are plotted as lines.
2-D diagrams, for two variables, are plotted as predominance fields.
The allowed variables are any that \code{\link{affinity}} or the other functions accepts: temperature, pressure, or the chemical activities of the basis species.

A new plot is started unless \code{add} is TRUE.
If \code{plot.it} is FALSE, no plot will be generated but all the intermediate computations will be performed and the results returned.

Line or field labels use the names of the species as provided in \code{eout}; formatting is applied to chemical formulas (see \code{\link{expr.species}}) unless \code{format.names} is FALSE.
Set \code{names} to TRUE or NULL to plot the names, or FALSE, NA, or \code{""} to prevent plotting the names, or a character argument to replace the default species names.
Alternatively, supply a numeric value to \code{names} to indicate a subset of default names that should or shouldn't be plotted (positive and negative indices, respectively).
Use \code{col.names} and \code{cex.names} to change the colors and size of the labels.
Use \code{cex} and \code{cex.axis} to adjust the overall character expansion factors (see \code{\link{par}}) and those of the axis labels.
The x- and y-axis labels are automatically generated unless they are supplied in \code{xlab} and \code{ylab}. 

If \code{groups} is supplied, the activities of the species identified in each numeric element of this list are multiplied by the balance coefficients of the species, then summed together.
The names of the list are used to label the lines or fields for the summed activities of the resulting groups.

}

\section{\code{type} argument}{
This paragraph describes the effect of the \code{type} argument when the output from \code{affinity} is being used.
For \code{type} set to \samp{auto}, and with 0 or 1 variables defined in \code{\link{affinity}}, the property computed by \code{affinity} for each species is plotted.
This is usually the affinity of the formation reactions, but can be set to some other property (using the \code{property} argument of \code{affinity}), such as the equilibrium constant (\samp{logK}).
For two variables, equilibrium predominance (maximum affinity) fields are plotted.
This \dQuote{maximum affinity method} (Dick, 2019) uses balancing coefficients that are specified by the \code{balance} argument.
If \code{type} is \samp{saturation}, the function plots the line for each species where the affinity of formation equals zero (see \code{demo("saturation")} for an example).
If for a given species no saturation line is possible or the range of the diagram does not include the saturation line, the function prints a message instead.
If \code{type} is the name of a basis species, then the equilibrium activity of the selected basis species in each of the formation reactions is plotted (see the \CO2-acetic acid example in \code{\link{buffer}}).
In the case of 2-D diagrams, both of these options use \code{\link{contour}} to draw the lines, with the method specified in \code{contour.method}.

This paragraph describes the effect of the \code{type} argument when the output from \code{solubility} is being used.
For one mineral or gas, if \code{type} set to \samp{auto}, the equilibrium activities of each aqueous species are plotted.
If \code{type} is \samp{loga.balance}, the activity of the balancing basis species (i.e. total solubility) is plotted; this is represented by contours on a 2-D diagram.
For two or more minerals or gases, if \code{type} set to \samp{auto}, the values of \samp{loga.balance} (overall minimum solubility) are plotted.
If \code{type} is \samp{loga.equil}, the solubilities of the individual minerals and gases are plotted.
For examples that use these features, see \code{\link{solubility}} and various \code{\link{demos}}: \samp{DEW}, \samp{contour}, \samp{gold}, \samp{solubility}, \samp{sphalerite}.
}

\section{1-D diagrams}{
For 1-D diagrams, the default setting for the y-axis is a logarithmic scale (unless \code{alpha} is TRUE) with limits corresponding to the range of logarithms of activities (or 0,1 if \code{alpha} is TRUE); these actions can be overridden by \code{ylog} and \code{ylim}.
If \code{legend.x} is NA (the default), the lines are labeled with the names of the species near the maximum value.
Otherwise, a \code{\link{legend}} is placed at the location identified by \code{legend.x}, or omitted if \code{legend.x} is NULL.

If \code{alpha} is TRUE, the fractional degrees of formation (ratios of activities to total activity) are plotted.
Or, setting \code{alpha} to \samp{balance} allows the activities to be multiplied by the number of the balancing component; this is useful for making \dQuote{percent carbon} diagrams where the species differ in carbon number.
The line type and line width can be controlled with \code{lty} and \code{lwd}, respectively.
To connect the points with splines instead of lines, set \code{spline.method} to one of the methods in \code{\link{splinefun}}.
}

\section{2-D diagrams}{
On 2-D diagrams, the fields represent the species with the highest equilibrium activity.
\code{fill} determines the color of the predominance fields, \code{col} that of the boundary lines.
The default of NULL for \code{fill} uses a light blue, light tan, and darker tan color for fields with aqueous species, one solid, or two solids.
These correspond to the web colors "aliceblue", "antiquewhite", and "burlywood" with some transparency added; see \viglink{multi-metal} for an example with two solids produced using \code{\link{mix}}.
If all the species in the diagram have the same state, or if the \code{fill} argument is NA or a 0-length value, the predominance fields are transparent, i.e. no fill color is used.
Otherwise, \code{fill} can be any \code{\link{colors}}, or the word \samp{rainbow}, \samp{heat}, \samp{terrain}, \samp{topo}, or \samp{cm}, indicating a palette from \pkg{grDevices}.
Starting with \R version 3.6.0, \code{fill} can be the name of any available HCL color palette, matched in the same way as the \code{palette} argument of \code{hcl.colors}.

\code{fill.NA} gives the color for empty fields, i.e. points at which NA values are present for any species.
This may occur when there are missing thermodynamic data or the temperature or pressure are not in the range of the equations of state.
To make overlay diagrams easier to construct, the default for \code{fill.NA} is automatically changed to \samp{transparent} when \code{add} is TRUE.

If \code{limit.water} is TRUE, the diagram is clipped to the the water stability region on Eh-pH (and some other) diagrams.
That is, predominance fields are shown only where water is stable, and \code{fill.NA} is used for areas where \H2O is not stable.
The default of NULL for \code{limit.water} does not clip the main diagram but instead overlays it on the water stability fields.
Change \code{limit.water} to FALSE to not show the water stability regions at all; this is automatically done if \code{limit.water} is missing and \code{add} is TRUE.

The default line-drawing algorithm uses \code{\link{contourLines}} to obtain smooth-looking diagonal and curved lines, at the expense of not coinciding exactly with the rectangular grid that is used for drawing colors.
\code{lty}, \code{col}, and \code{lwd} can be specified, but limiting the lines via \code{xrange} is not currently supported.
Set \code{lty.cr} or \code{lty.aq} to 0 to suppress boundary lines between minerals or aqueous species.

To go back to the old behavior for drawing lines, set \code{dotted} to \samp{0}.
The old behavior does not respect \code{lty}; instead, the style of the boundary lines on 2-D diagrams can be altered by supplying one or more non-zero integers in \code{dotted}, which indicates the fraction of line segments to omit; a value of \samp{1} or NULL for \code{dotted} has the effect of not drawing the boundary lines.

\code{normalize} and \code{as.residue} apply only to the 2-D diagrams, and only when \code{eout} is the output from \code{affinity}.
With \code{normalize}, the activity boundaries are calculated between the residues of the species (the species divided by the balance coefficients), then the activities are rescaled to the whole species formulas.
With \code{as.residue}, the activity boundaries are calculated between the residues of the species, and no rescaling is performed.
}

\section{Activity Coefficients}{
The wording in this page and names of variables in functions refer exclusively to \samp{activities} of aqueous species.
However, if activity coefficients are calculated (using the \code{IS} argument in \code{\link{affinity}}), then these variables are effectively transformed to molalities (see \code{tests/testthat/} \code{test-logmolality.R}).
So that the labels on diagrams are adjusted accordingly, \code{\link{diagram}} sets the \code{molality} argument of \code{\link{axis.label}} to TRUE if \code{IS} was supplied as an argument to \code{\link{affinity}}.
The labeling as molality takes effect even if \code{IS} is set to 0; this way, by including (or not) the \code{IS = 0} argument to \code{affinity}, the user decides whether to label aqueous species variables as molality (or activity) for calculations at zero ionic strength (where molality = activity).
}

\section{Other Functions}{
\code{find.tp} finds the locations in a matrix of integers that are surrounded by the greatest number of different values.
The function counts the unique values in a 3x3 grid around each point and returns a matrix of indices (similar to \code{\link{which}(..., arr.ind = TRUE)}) for the maximum count (ties result in more than one pair of indices).
It can be used with the output from \code{diagram} for calculations in 2 dimensions to approximately locate the triple points on the diagram.
}

\section{Value}{
  \code{diagram} returns an \code{\link{invisible}} list containing, first, the contents of \code{eout}, i.e. the output of \code{\link{affinity}} or \code{\link{equilibrate}} supplied in the function call.
  To this are added the names of the plotted variable in \code{plotvar}, the labels used for species (which may be \code{\link{plotmath}} expressions if \code{format.names} is TRUE) in \code{names}, and the values used for plotting in a list named \code{plotvals}.
  For 1-D diagrams, \code{plotvals} usually corresponds to the chemical activities of the species (i.e. \code{eout$loga.equil}), or, if \code{alpha} is \code{TRUE}, their mole fractions (degrees of formation).
  For 2-D diagrams, \code{plotvals} corresponds to the values of affinity (from \code{eout$values}) divided by the respective balancing coefficients for each species.
  For 2-D diagrams, the output also contains the matrices \code{predominant}, which identifies the predominant species in \code{eout$species} at each grid point, and \code{predominant.values}, which has the affinities of the predominant species divided by the balancing coefficients (if \code{eout} is the output of \code{affinity}) or the activities of the predominant species (if \code{eout} is the output of \code{equilibrate}).
  The rows and columns of these matrices correspond to the \emph{x} and \emph{y} variables, respectively.
}

\seealso{ 
Other examples are present in the help for \code{\link{buffer}}, and even more can be found in \code{\link{demos}}.
}

\examples{
\dontshow{reset()}
## calculate the equilibrium logarithm of activity of a 
## basis species in different reactions
basis("CHNOS")
species(c("ethanol", "lactic acid", "deoxyribose", "ribose"))
a <- affinity(T=c(0, 150))
diagram(a, type="O2", legend.x="topleft", col=rev(rainbow(4)), lwd=2)
title(main="Equilibrium logfO2 for 1e-3 mol/kg of CO2 and ... ")

### 1-D diagrams: logarithms of activities

## Degrees of formation of ionized forms of glycine
## After Fig. 1 of Aksu and Doyle, 2001 
basis("CHNOS+")
species(ispecies <- info(c("glycinium", "glycine", "glycinate")))
a <- affinity(pH=c(0, 14))
e <- equilibrate(a)
diagram(e, alpha=TRUE, lwd=1)
title(main=paste("Degrees of formation of aqueous glycine species\n",
  "after Aksu and Doyle, 2001"))

## Degrees of formation of ATP species as a function of 
## temperature, after LaRowe and Helgeson, 2007, Fig. 10b
# to make a similar diagram, activity of Mg+2 here is set to
# 10^-4, which is different from LH07, who used 10^-3 total molality
basis(c("CO2", "NH3", "H2O", "H3PO4", "O2", "H+", "Mg+2"),
  c(999, 999, 999, 999, 999, -5, -4))
species(c("HATP-3", "H2ATP-2", "MgATP-2", "MgHATP-"))
a <- affinity(T=c(0, 120, 25))
e <- equilibrate(a)
diagram(e, alpha=TRUE)  
title(main=paste("Degrees of formation of ATP species,\n",
  "pH=5, log(aMg+2)=-3. After LaRowe and Helgeson, 2007"),
  cex.main=0.9)

### 2-D diagrams: predominance diagrams
### these use the maximum affinity method

## Fe-S-O at 200 deg C, after Helgeson, 1970
basis(c("Fe", "oxygen", "S2"))
species(c("iron", "ferrous-oxide", "magnetite",
  "hematite", "pyrite", "pyrrhotite"))
# the calculations include the phase transitions of
# pyrrhotite; no additional step is needed
a <- affinity(S2=c(-50, 0), O2=c(-90, -10), T=200)
diagram(a, fill="heat")
title(main=paste("Fe-S-O, 200 degrees C, 1 bar",
  "After Helgeson, 1970", sep="\n"))

## pe-pH diagram for hydrated iron sulfides,
## goethite and pyrite, after Majzlan et al., 2006
basis(c("Fe+2", "SO4-2", "H2O", "H+", "e-"), 
  c(0, log10(3), log10(0.75), 999, 999))
species(c("rhomboclase", "ferricopiapite", "hydronium jarosite",
  "goethite", "melanterite", "pyrite"))
a <- affinity(pH=c(-1, 4, 256), pe=c(-5, 23, 256))
d <- diagram(a, main="Fe-S-O-H, after Majzlan et al., 2006")
water.lines(d, lwd=2)
text(3, 22, describe.basis(thermo()$basis[2:3,], digits=2, oneline=TRUE))
text(3, 21, describe.property(c("T", "P"), c(25, 1), oneline=TRUE))

## aqueous Al species, after Tagirov and Schott, 2001
basis(c("Al+3", "F-", "H+", "O2", "H2O"))
AlOH <- c("Al(OH)4-", "Al(OH)3", "Al(OH)2+", "AlOH+2")
Al <- "Al+3"
AlF <- c("AlF+2", "AlF2+", "AlF3", "AlF4-")
AlOHF <- c("Al(OH)2F2-", "Al(OH)2F", "AlOHF2")
species(c(AlOH, Al, AlF, AlOHF), "aq")
res <- 300
a <- affinity(pH = c(0.5, 6.5, res), `F-` = c(-2, -9, res), T = 200)
diagram(a, fill = "terrain")
dprop <- describe.property(c("T", "P"), c(200, "Psat"))
legend("topright", legend = dprop, bty = "n")
mtitle(c("Aqueous aluminum species",
         "After Tagirov and Schott, 2001 Fig. 4d"), cex = 0.95)

## Temperature-Pressure: kayanite-sillimanite-andalusite
# cf. Fig. 49 of Helgeson et al., 1978
# this is a system of one component (Al2SiO5), however:
# - number of basis species must be the same as of elements
# - avoid using H2O or other aqueous species because of
#     T/P limits of the water() calculations;
basis(c("corundum", "quartz", "oxygen"))
species(c("kyanite", "sillimanite", "andalusite"))
# database has transition temperatures of kyanite and andalusite
# at 1 bar only, so we permit calculation at higher temperatures
a <- affinity(T=c(200, 900, 99), P=c(0, 9000, 101), exceed.Ttr=TRUE)
d <- diagram(a, fill=NULL)
slab <- syslab(c("Al2O3", "SiO2", "H2O"))
mtitle(c(as.expression(slab), "after Helgeson et al., 1978"))
# find the approximate position of the triple point
tp <- find.tp(d$predominant)
Ttp <- a$vals[[1]][tp[1, 2]]
Ptp <- rev(a$vals[[2]])[tp[1, 1]]
points(Ttp, Ptp, pch=10, cex=5)
}

\references{

Aksu, S. and Doyle, F. M. (2001) Electrochemistry of copper in aqueous glycine solutions. \emph{J. Electrochem. Soc.} \bold{148}, B51--B57.

Dick, J. M. (2019) CHNOSZ: Thermodynamic calculations and diagrams for geochemistry. \emph{Front. Earth Sci.} \bold{7}:180. \doi{10.3389/feart.2019.00180}

Helgeson, H. C. (1970) A chemical and thermodynamic model of ore deposition in hydrothermal systems. \emph{Mineral. Soc. Amer. Spec. Pap.} \bold{3}, 155--186. \url{https://www.worldcat.org/oclc/583263}

Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K. (1978) Summary and critique of the thermodynamic properties of rock-forming minerals. \emph{Am. J. Sci.} \bold{278-A}, 1--229. \url{https://www.worldcat.org/oclc/13594862}

LaRowe, D. E. and Helgeson, H. C. (2007) Quantifying the energetics of metabolic reactions in diverse biogeochemical systems: electron flow and ATP synthesis. \emph{Geobiology} \bold{5}, 153--168. \doi{10.1111/j.1472-4669.2007.00099.x}

Majzlan, J., Navrotsky, A., McClesky, R. B. and Alpers, C. N. (2006) Thermodynamic properties and crystal structure refinement of ferricopiapite, coquimbite, rhomboclase, and Fe\s{2}(SO\s{4})\s{3}(H\s{2}O)\s{5}. \emph{Eur. J. Mineral.} \bold{18}, 175--186. \doi{10.1127/0935-1221/2006/0018-0175}

Tagirov, B. and Schott, J. (2001) Aluminum speciation in crustal fluids revisited. \emph{Geochim. Cosmochim. Acta} \bold{65}, 3965--3992. \doi{10.1016/S0016-7037(01)00705-0}

}

\concept{Main workflow}
