\encoding{UTF-8}
\name{util.list}
\alias{util.list}
\alias{which.pmax}
\title{Functions to Work with Lists}

\description{
  Combine lists or perform arithmetic operations on elements of lists.
}

\usage{
  which.pmax(elts, na.rm = FALSE, pmin = FALSE)
}

\arguments{
  \item{elts}{list, numeric vectors for which to find maximum values (in parallel) (\code{which.pmax}).}
  \item{na.rm}{logical, remove missing values?}
  \item{pmin}{logical, find minimum values instead of maximum ones?}
}

\details{
\code{which.pmax} takes a list of equal-length numeric vectors (or objects that can be coerced to numeric) in \code{elts} and returns the index of the vector holding the maximum value at each position.
If \code{na.rm} is \code{TRUE}, values of \code{NA} are removed; if \code{pmin} is \code{TRUE} the function finds locations of the minimum values instead.
}

\concept{Utility functions}
