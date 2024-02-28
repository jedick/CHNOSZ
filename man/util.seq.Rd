\encoding{UTF-8}
\name{util.seq}
\alias{util.seq}
\alias{aminoacids}

\title{Functions to Work with Sequence Data}

\description{
  Return names or one- or three-letter abbreviations of amino acids.
}

\usage{
  aminoacids(nchar = 1, which = NULL)
}

\arguments{
  \item{nchar}{numeric, \eqn{1} to return one-letter, \eqn{3} to return three-letter abbreviations for amino acids}
  \item{which}{character, which amino acids to name}
}

\details{

\code{aminoacids} returns the one-letter abbreviations (\code{nchar}=\samp{1}) or the three-letter abbreviations (\code{nchar}=\samp{3}) or the names of the neutral amino acids (\code{nchar}=\samp{""}) or the names of the amino acids with ionized side chains (\code{nchar}=\samp{"Z"}).
The output includes 20 amino acids in alphabetic order by 1-letter abbreviation (the order used in \code{thermo()$protein}), unless \code{which} is provided, indicating the desired amino acids (either as 1- or 3-letter abbreviations or names of the neutral amino acids).

}

\seealso{
  \code{\link{protein.formula}} for calculating the chemical formulas of proteins.
}

\examples{
# Three-letter abbreviations of amino acids
aminoacids(3)
}

\concept{Protein properties}
