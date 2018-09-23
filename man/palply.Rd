\encoding{UTF-8}
\name{palply}
\alias{palply}

\title{Conditional Parallel Processing}
\description{
  Use multiple processors for large calculations.
}

\usage{
  palply(varlist, X, FUN, ...)
}

\arguments{
  \item{...}{equivalent to the same argument in \code{\link[parallel]{parLapply}}}
  \item{varlist}{character, names of variables to export using \code{\link[parallel]{clusterExport}}}
  \item{X}{vector, argument for \code{\link{lapply}} or \code{parLapply}}
  \item{FUN}{function, argument for \code{lapply} or \code{parLapply}}
}

\details{
\code{palply} is a wrapper function to run \code{\link[parallel]{parLapply}} if length of \code{X} > \code{\link{thermo}$opt$paramin} and multiple cores are available, otherwise it runs \code{\link{lapply}}.
Note that \code{parLapply} is called with \code{methods} set to FALSE.
If lots of package startup messages are created when running \code{\link[parallel]{makeCluster}} (which is called by \code{palply}), it can probably be stopped by adding a test for \code{\link{interactive}} sessions around any \code{\link{library}} commands in the \code{\link{Rprofile}}.
}

\seealso{
\code{\link{read.fasta}}, \code{\link{count.aa}}, \code{\link{affinity}}, \code{\link{equil.boltzmann}} and \code{\link{equil.reaction}} for functions that use \code{palply}.
Tests are in \file{tests/test-util.program.R}, and a \dQuote{real world} example is in \file{demos/density.R}.
}

\concept{Utility functions}
