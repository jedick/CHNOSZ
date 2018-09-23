\encoding{UTF-8}
\name{util.test}
\alias{util.test}
\alias{maxdiff}
\alias{expect_maxdiff}
\title{Functions for Writing Tests}
\description{
Functions modelled after the \code{expect_} functions in \CRANpkg{testthat}.
}

\usage{
  maxdiff(x, y)
  expect_maxdiff(object, expected, maxdiff = 0)
}

\arguments{
  \item{x}{numeric object}
  \item{y}{numeric object}
  \item{object}{numeric, object to test}
  \item{expected}{numeric, expected value}
  \item{maxdiff}{numeric, maximum pairwise difference between object and expected value}
}

\details{
\code{maxdiff} computes the maximum (absolute) pairwise difference between x and y, i.e. \code{max(abs(y - x))}.

\code{expect_maxdiff} tests that the maximum of the pairwise differences between two objects is less than the value of the argument \code{maxdiff}.
The function uses \code{\link[testthat]{expect}} to generate an expectation in the \CRANpkg{testthat} framework.
}

\concept{Utility functions}
