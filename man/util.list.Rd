\encoding{UTF-8}
\name{util.list}
\alias{util.list}
\alias{which.pmax}
\title{Functions to Work with Lists}

\description{
  Identify list elements that have the maximum (or minimum) values.
}

\usage{
  which.pmax(x, maximum = TRUE)
}

\arguments{
  \item{x}{list of numeric vectors}
  \item{maximum}{logical, find maximum values or minimum values?}
}

\details{
\code{which.pmax} takes a list of equal-length numeric vectors or equal-dimension arrays in \code{x} and returns the index of the list element that has the maximum value at each point.
Change \code{maximum} to \code{FALSE} to find the minimum values instead.
}

\concept{Utility functions}
