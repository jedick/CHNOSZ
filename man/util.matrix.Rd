\encoding{UTF-8}
\name{util.matrix}
\alias{util.matrix}
\alias{invertible.combs}
\title{Functions for Various Matrix Operations}
\description{
  Find rows of a matrix that form invertible (linearly independent) combinations.
}

\usage{
  invertible.combs(A, nmax=20)
}

\arguments{
  \item{A}{A matrix, with at least as many rows as columns.}
  \item{nmax}{The maximum number of rows to consider.}
}

\details{
  Given a matrix \code{A}, with number of rows equal to or greater than the number of columns, return the combinations of row numbers that constitute invertible square matrices. Consider only the first \code{nmax} rows of the original matrix (to save time for large systems).
}

\examples{\dontshow{reset()}
## what combinations of the 20 common amino acids have
## a linearly independent stoichiometry with five elements?
# the names of the amino acids
aanames <- aminoacids("")
# their species indices
iaa <- suppressMessages(info(aanames))
# the full stoichiometric matrix
A <- i2A(iaa)
# the invertible combinations
icA <- invertible.combs(A)
stopifnot(nrow(icA)==6067)
# that's a bit less than 40% of all possible combinations
nrow(icA) / ncol(combn(20, 5))
# count the occurrences of each amino acid
counts <- table(icA)
names(counts) <- aminoacids(1)
(sc <- sort(counts))
# the two sulfur-containing ones show up most frequently
stopifnot(tail(names(sc), 2)==c("C", "M"))
}

\concept{Utility functions}
