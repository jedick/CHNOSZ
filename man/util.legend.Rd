\encoding{UTF-8}
\name{util.legend}
\alias{util.legend}
\alias{lNaCl}
\alias{lS}
\alias{lT}
\alias{lP}
\alias{lTP}
\title{Functions to make legend text}
\description{Generate expressions describing system conditions that can be added to legends of plots.}

\usage{
  lNaCl(x, digits = 2)
  lS(x, digits = 3)
  lT(x, digits = 0)
  lP(x, digits = 0)
  lTP(x, y, digits = 0)
}

\arguments{
  \item{x}{numeric, value of the property}
  \item{digits}{numeric, digits for rounding}
  \item{y}{numeric, value of pressure}
}

\details{
These functions are used to make expressions for common chemical system variables that can be used in plot legends.
\code{lNaCl} describe the molality of NaCl, and \code{lS} the total molality of sulfur.
\code{lT} and \code{lP} describe the temperature and pressure.
\code{lTP} describe the temperature and pressure together, separated by a comma.
}

\seealso{
  \code{\link{util.expression}}
}

\examples{
plot.new()
l <- c(lTP(100, "Psat"), lNaCl(1), lS(1e-3))
legend("center", legend = l)
}

\concept{Utility functions}
