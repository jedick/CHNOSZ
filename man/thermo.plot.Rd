\encoding{UTF-8}
\name{thermo.plot}
\alias{thermo.plot}
\alias{thermo.plot.new}
\alias{thermo.axis}
\title{Functions to create and modify plots}

\description{
  Initialize a new plot window using preset parameters or add an axis to a plot.
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
  \item{...}{further arguments passed to \code{par}}
  \item{lab}{character, axis label}
  \item{line}{numeric, margin line on which to place axis label}
  \item{plot.line}{logical, draw axis lines?}
}

\details{

\code{thermo.plot.new} sets parameters for a new plot, creates a new plot using \code{\link{plot.new}}, and adds the axes tick marks to the plot.
Plot parameters (see \code{\link{par}}) including \code{cex}, \code{mar}, \code{lwd}, \code{mgp} and \code{axs} can be given, as well as a numeric vector in \code{side} identifying which sides of the plot receive tick marks.
\code{yline}, if present, denotes the margin line (default \code{\link{par}('mgp')[1]}) where the y-axis name is plotted.
\code{thermo.axis} is the function that actually adds the axes, including inward-pointing major and minor tick marks (often used for thermodynamic property diagrams).

Use \code{grid} to add a grid to the plot, corresponding to either the major ticks (solid lines), minor ticks (dashed lines), or both.
The grid can be made by adding \code{grid} argument to \code{\link{diagram}}, or by calling \code{thermo.axis} after \code{diagram} (see example).

}

\seealso{
\code{\link{diagram}} uses \code{thermo.plot.new} to set up a new plot, unless the argument \code{tplot} is set to FALSE in \code{diagram}.
}

\examples{
# Start a new plot with specific settings
thermo.plot.new(c(0, 10), c(0, 100), xlab = "X variable", ylab = "Y variable")

# Add major and minor tick marks to an existing plot
plot(1:10)
thermo.axis()
}

\concept{Utility functions}
