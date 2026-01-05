\encoding{UTF-8}
\name{util.plot}
\alias{util.plot}
\alias{label.plot}
\alias{usrfig}
\alias{label.figure}
\alias{water.lines}
\alias{mtitle}
\alias{ZC.col}
\alias{hyphen.in.pdf}
\title{Utility functions for plots}

\description{
  Label plots and figures, generate labels for axes and subplots, add stability lines for water, or get colors for a set of numeric values.
}

\usage{
  label.plot(label, xfrac = 0.07, yfrac = 0.93, paren = FALSE,
    italic = FALSE, ...)
  usrfig()
  label.figure(label, xfrac = 0.05, yfrac = 0.95, paren = FALSE,
    italic = FALSE, ...)
  water.lines(eout, which = c("oxidation","reduction"),
    lty = 2, lwd = 1, col = par("fg"), plot.it = TRUE)
  mtitle(main, line = 0, spacing = 1, ...)
  ZC.col(z)
  hyphen.in.pdf(x)
}

\arguments{
  \item{label}{character, label to place on plot}
  \item{xfrac}{numeric, fractional location on \eqn{x}{x}-axis for placement of label}
  \item{yfrac}{numeric, fractional location on \eqn{y}{y}-axis for placement of label}
  \item{paren}{logical, add parentheses around label text?}
  \item{...}{further arguments passed to \code{text} or \code{mtext}}
  \item{italic}{logical, italicize label text?}
  \item{eout}{data frame, output of \code{\link{affinity}}, \code{\link{equilibrate}}, or \code{\link{diagram}}}
  \item{which}{character, which of oxidation/reduction lines to plot}
  \item{lty}{numeric, line type}
  \item{lwd}{numeric, line width}
  \item{col}{character, color}
  \item{plot.it}{logical, plot the lines?}
  \item{main}{character, text for plot title}
  \item{line}{numeric, margin line on which to place plot title}
  \item{spacing}{numeric, spacing between multiple lines}
  \item{z}{numeric, set of values}
  \item{x}{character, any text object}
}

\details{

\code{water.lines} plots lines representing the oxidation and reduction stability limits of water on Eh/pe/\logfO2/\logfH2 vs pH/\T/\P diagrams.
The x- and y-variables and their ranges are taken from \code{eout}.
Values of \T, \P, pH, and \logaH2O, not corresponding to either axis, are also taken from \code{eout}.
\code{which} controls which lines are drawn (\samp{oxidation}, \samp{reduction}, or both (the default)).
The value of \code{swapped} in the output reflects whether pH, \T, or \P is on the x-axis (TRUE) or y-axis (FALSE).
\code{NA} is returned for any diagram for variables that can not be processed (including diagrams with more than 2 variables).

\code{label.plot} and \code{label.figure} add identifying text within the plot region and figure region.
The value given for \code{label} is used, optionally italicized and with parentheses (like \ifelse{latex}{\eqn{(a)}}{\ifelse{html}{\out{(<i>a</i>)}}{(a)}}).
The location of the label is controlled by \code{xfrac} and \code{yfrac} (the fractional coordinates of either the plot or figure region),
and \code{...} can include other parameters such as \code{cex} and \code{adj} that are passed to \code{\link{text}}.

\code{usrfig} returns the limits of the figure region in \dQuote{user} coordinates (i.e. the limits of the plot region, from \code{\link{par}("usr")}).
It is a supporting function for \code{label.figure} but is also useful for other circumstances where information must be added at a particular location in a figure.

\code{mtitle} can be used to add a multi-line title to a plot.
It loops over each element of \code{main} and places it on a separate margin line using \code{\link{mtext}}.
The spacing of the last (bottom) line from the edge of the plot is specified by \code{line}.
This function exists to facilitate using \code{\link{expression}}s in multiline titles.

\code{ZC.col} uses \CRANpkg{colorspace} to generate colors from a diverging palette (red - light grey - blue) corresponding to the values in \code{z}.
Red is associated with lower values of \code{z}.
This function is intended to generate colors for distinguishing average oxidation state of carbon \code{\link{ZC}}, but any numeric values can be supplied.

If a PDF device is active, \code{hyphen.in.pdf} replaces all instances of ASCII hyphen (\verb{-}) in \code{x} with the Unicode hyphen (\verb{\uad}).
This is done so that hyphens don't get displayed as minus signs in R graphics output, but it only works for the PDF device.
}

\examples{
basis(c("H2S", "H2O", "H+", "e-"))
species(c("H2S", "HS-", "S3-", "SO2", "HSO4-", "SO4-2"))
a <- affinity(pH = c(0, 12), Eh = c(-1, 1), T = 200)
opar <- par(mfrow = c(2, 2))
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
