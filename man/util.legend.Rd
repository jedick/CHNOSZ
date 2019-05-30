\encoding{UTF-8}
\name{util.legend}
\alias{util.legend}
\alias{lNaCl}
\alias{lS}
\alias{lT}
\alias{lP}
\alias{lTP}
\alias{lex}
\title{Functions to Make Legend Text}
\description{Generate expressions suitable plot legends describing system conditions.}

\usage{
  lNaCl(x, digits = 2)
  lS(x, digits = 3)
  lT(x, digits = 0)
  lP(x, digits = 0)
  lTP(x, y, digits = 0)
  lex(...)
}

\arguments{
  \item{x}{numeric, value of the property}
  \item{digits}{numeric, digits for rounding}
  \item{y}{numeric, value of pressure}
  \item{...}{language, objects to combine in an expression}
}

\details{
These functions are used to make expressions for common chemical system variables that can be used in plot legends.
\code{lNaCl} describe the molality of NaCl, and \code{lS} the total molality of sulfur.
\code{lT} and \code{lP} describe the temperature and pressure.
\code{lTP} describe the temperature and pressure together, separated by a comma.

The above functions return language objects, which can be combined with \code{lex} to make an expression that when used in \code{\link{legend}} appears on multiple lines.

}

\seealso{
  \code{\link{util.expression}}
}

\examples{
plot.new()
l <- lex(lTP(100, "Psat"), lNaCl(1), lS(1e-3))
legend("center", l)
}

\concept{Utility functions}
