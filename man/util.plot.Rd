\encoding{UTF-8}
\name{util.plot}
\alias{util.plot}
\alias{thermo.plot.new}
\alias{thermo.axis}
\alias{label.plot}
\alias{usrfig}
\alias{label.figure}
\alias{water.lines}
\alias{mtitle}
\alias{ZC.col}
\alias{add.alpha}
\title{Functions to Create and Modify Plots}

\description{
  Initialize a new plot window using preset parameters, add an axis or title to a plot, generate labels for axes and subplots, add stability lines for water, get colors for a set of numeric values.
}

\usage{
  thermo.plot.new(xlim, ylim, xlab, ylab, cex = par("cex"),
    mar = NULL, lwd = par("lwd"), side = c(1,2,3,4), 
    mgp = c(1.7, 0.3, 0), cex.axis = par("cex"), col = par("col"),
    yline = NULL, axs = "i", plot.box = TRUE, las = 1,
    xline = NULL, grid = "", col.grid = "gray", ...)
  thermo.axis(lab = NULL, side = 1:4, line = 1.5, cex = par("cex"),
    lwd = par("lwd"), col = par("col"), grid = "", col.grid = "gray",
    plot.line = FALSE)
  label.plot(x, xfrac = 0.07, yfrac = 0.93, paren = FALSE,
    italic = FALSE, ...)
  usrfig()
  label.figure(x, xfrac = 0.05, yfrac = 0.95, paren = FALSE,
    italic = FALSE, ...)
  water.lines(eout, which = c("oxidation","reduction"),
    lty = 2, lwd = 1, col = par("fg"), plot.it = TRUE)
  mtitle(main, line = 0, spacing = 1, ...)
  ZC.col(z)
  add.alpha(col, alpha)
}

\arguments{
  \item{xlim}{numeric, limits of the \eqn{x}{x}-axis}
  \item{ylim}{numeric, limits of the \eqn{y}{y}-axis}
  \item{xlab}{character, \eqn{x}{x}-axis label}
  \item{ylab}{character, \eqn{y}{y}-axis label}
  \item{cex}{numeric, character expansion factor for labels}
  \item{mar}{numeric, width (number of lines) of margins on each side of plot}
  \item{lwd}{numeric, line width}
  \item{side}{numeric, which sides of plot to draw axes}
  \item{mgp}{numeric, sizes of margins of plot}
  \item{cex.axis}{numeric, character expansion factor for names of axes}
  \item{col}{character, color}
  \item{yline}{numeric, margin line on which to plot \eqn{y}{y}-axis name}
  \item{axs}{character, setting for axis limit calculation}
  \item{plot.box}{logical, draw a box around the plot?}
  \item{las}{numeric, style for axis labels}
  \item{xline}{numeric, margin line on which to plot \eqn{x}{x}-axis name}
  \item{grid}{character, type of grid (\samp{major}, \samp{minor}, or \samp{both})}
  \item{col.grid}{character, color of the grid lines}
  \item{plot.line}{logical, draw axis lines?}
  \item{...}{further arguments passed to \code{par} or \code{mtext}}
  \item{lab}{character, axis label}
  \item{line}{numeric, margin line on which to place axis label or plot title}
  \item{x}{character, label to place on plot}
  \item{xfrac}{numeric, fractional location on \eqn{x}{x}-axis for placement of label}
  \item{yfrac}{numeric, fractional location on \eqn{y}{y}-axis for placement of label}
  \item{paren}{logical, add parentheses around label text?}
  \item{italic}{logical, italicize label text?}
  \item{eout}{data frame, output of \code{\link{affinity}}, \code{\link{equilibrate}}, or \code{\link{diagram}}}
  \item{which}{character, which of oxidation/reduction lines to plot}
  \item{lty}{numeric, line type}
  \item{plot.it}{logical, plot the lines?}
  \item{main}{character, text for plot title}
  \item{spacing}{numeric, spacing between multiple lines}
  \item{z}{numeric, set of values}
  \item{alpha}{character, hexadecimal value of color transparency (alpha)}
}

\details{

\code{thermo.plot.new} sets parameters for a new plot, creates a new plot using \code{\link{plot.new}}, and adds the axes tick marks to the plot.
Plot parameters (see \code{\link{par}}) including \code{cex}, \code{mar}, \code{lwd}, \code{mgp} and \code{axs} can be given, as well as a numeric vector in \code{side} identifying which sides of the plot receive tick marks.
\code{yline}, if present, denotes the margin line (default \code{\link{par}('mgp')[1]}) where the y-axis name is plotted.
\code{thermo.axis} is the function that actually adds the axes, including inward-pointing major and minor tick marks (often used for thermodynamic property diagrams).

Use \code{grid} to add a grid to the plot, corresponding to either the major ticks (solid lines), minor ticks (dashed lines), or both.
The grid can be made by adding \code{grid} argument to \code{\link{diagram}}, or by calling \code{thermo.axis} after \code{diagram} (see example).

\code{water.lines} plots lines representing the oxidation and reduction stability limits of water on Eh/pe/\logfO2/\logfH2 vs pH/\T/\P diagrams.
The x- and y-variables and their ranges are taken from \code{eout}.
Values of \T, \P, pH, and \logaH2O, not corresponding to either axis, are also taken from \code{eout}.
\code{which} controls which lines are drawn (\samp{oxidation}, \samp{reduction}, or both (the default)).
The value of \code{swapped} in the output reflects whether pH, \T, or \P is on the x-axis (TRUE) or y-axis (FALSE).
\code{NA} is returned for any diagram for variables that can not be processed (including diagrams with more than 2 variables).

\code{label.plot} and \code{label.figure} add identifying text within the plot region and figure region.
The value given for \code{x} is made into a label, optionally italicized and with parentheses (like \ifelse{latex}{\eqn{(a)}}{\ifelse{html}{\out{(<i>a</i>)}}{(a)}}).
The location of the label is controlled by \code{xfrac} and \code{yfrac} (the fractional coordinates of either the plot or figure region), and \code{...} can include other parameters such as \code{cex} and \code{adj} that are passed to \code{\link{text}}.

\code{usrfig} returns the limits of the figure region in \dQuote{user} coordinates (i.e. the limits of the plot region, from \code{\link{par}("usr")}).
It is a supporting function for \code{label.figure} but is also useful for other circumstances where information must be added at a particular location in a figure.

\code{mtitle} can be used to add a multi-line title to a plot.
It loops over each element of \code{main} and places it on a separate margin line using \code{\link{mtext}}.
The spacing of the last (bottom) line from the edge of the plot is specified by \code{line}.
This function exists to facilitate using \code{\link{expression}}s in multiline titles.

\code{ZC.col} uses \CRANpkg{colorspace}) to generate colors from a diverging palette (red - light grey - blue) corresponding to the values in \code{z}.
Red is associated with lower values of \code{z}.
This function is intended to generate colors for distinguishing average oxidation state of carbon \code{\link{ZC}}, but any numeric values can be supplied.

\code{add.alpha} adds transparency to a color by appending the value of \code{alpha} to the hexadecimal representation of the color given in \code{col}.
}

\seealso{
\code{\link{diagram}} uses \code{thermo.plot.new} to set up a new plot, unless the argument \code{tplot} is set to FALSE in \code{diagram}.
}

\examples{
basis(c("H2S", "H2O", "H+", "e-"))
species(c("HS-", "H2S", "HSO4-", "SO4-2"))
a <- affinity(pH = c(0, 12), Eh = c(-1, 1), T = 200)
opar <- par(mfrow=c(2, 2))
diagram(a, grid = "both")
title(main = 'diagram(a, grid = "both")')
diagram(a, grid = "major")
title(main = 'diagram(a, grid = "major")')
diagram(a, grid = "minor")
title(main = 'diagram(a, grid = "minor")')
diagram(a, fill = "terrain")
thermo.axis(grid = "major", col.grid = "slategray")
title(main = 'thermo.axis(grid = "major")')
par(thermo()$opar)
par(opar)
}

\concept{Utility functions}
