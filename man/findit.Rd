\encoding{UTF-8}
\name{findit}
\alias{findit}
\alias{plot_findit}
\title{Gridded Search to Optimize Objective Functions}
\description{
  Use a gridded search to find a combination of one or more of chemical activities of basis species, temperature and/or pressure that maximize or minimize a objective function of the metastable equilibrium chemical activities of the species of interest.
}

\usage{
  findit (lims = list(), objective = "CV", niter = NULL, iprotein = NULL, 
    plot.it = TRUE, T = 25, P = "Psat", res = NULL, labcex = 0.6, 
    loga2 = NULL, loga.balance = 0, rat = NULL,
    balance = NULL, normalize = FALSE)
  plot_findit(x, which=NULL, mar=c(3.5,5,2,2), xlab="iteration", \dots)
}

\arguments{
  \item{lims}{list, specification of search limits}
  \item{objective}{character, name of objective function to optimize}
  \item{niter}{numeric, number of iterations}
  \item{res}{numeric, grid resolution (number of points on one edge)}
  \item{iprotein}{numeric, indices of proteins}
  \item{plot.it}{logical, make a plot?}
  \item{T}{numeric, temperature}
  \item{P}{numeric, pressure; or character, "Psat"}
  \item{labcex}{numeric, character expansion for plot labels}
  \item{loga2}{numeric, reference logarithms of activity of species}
  \item{loga.balance}{numeric, logarithm of total activity of balanced quantity (passed to \code{\link{diagram}})}
  \item{rat}{numeric, ratio of edge length in successive iterations}
  \item{balance}{character or numeric, balanced quantity (passed to \code{\link{diagram}})}
  \item{normalize}{logical, normalize chemical formulas by the balance vector? (passed to \code{\link{diagram}})}
  \item{x}{list, object of class \code{findit}}
  \item{which}{numeric, which of the parameters to plot}
  \item{mar}{numeric, plot margin specification}
  \item{xlab}{character, x-axis label}
  \item{\dots}{additional arguments passed to \code{plot}}
}

\details{
\code{findit} implements a gridded optimization to find the minimum or maximum value of an \code{\link{objective}} function.
The variables are one or more of the chemical activities, temperature and/or pressure whose ranges are listed in \code{lims}.
Generally, the system (\code{\link{basis}} species and \code{\link{species}} of interest) must be set up before calling this function.
If \code{iprotein} is supplied, indicating a set of proteins to use in the calculation, the definition of the \code{species} is not required. 
\code{lims} is a list, each element of which is vector having a name that is the formula of one of the basis species, \samp{T} or \samp{P} and a pair of values indicating the range of the named parameter.
The values are the logarithms of activities of the basis species, or temperature or pressure (in the user's units; see \code{\link{util.units}}).
If either \samp{T} or \samp{P} is missing from the list in \code{lims}, the calculations are performed at isothermal and/or isobaric conditions indicated by \code{T} and \code{P} arguments.

Taking \eqn{nd} as the number of dimensions (number of variables in \code{lims}), default values of \code{niter} and \code{res} come from the following table.
These settings have been selected to be able to run the function quickly in the higher dimensions.
Detailed studies of a system might have to use more iterations and/or higher resolutions.

\tabular{rrrrr}{
  \code{nd} \tab \code{niter} \tab \code{res} \tab grid points (\code{res^nd}) \tab \code{rat} \cr
  1 \tab 4 \tab 128 \tab 128 \tab 0.7 \cr
  2 \tab 6 \tab 64 \tab 4096 \tab 0.7 \cr
  3 \tab 6 \tab 16 \tab 4096 \tab 0.8 \cr
  4 \tab 8 \tab 8 \tab  4096 \tab 0.9 \cr
  5 \tab 12 \tab 6 \tab 7776 \tab 0.9 \cr
  6 \tab 12 \tab 4 \tab 4096 \tab 0.95 \cr
  7 \tab 12 \tab 4 \tab 16384 \tab 0.95 \cr
}

The function performs \code{niter} iterations.
At first, the limits of the parameters given in \code{lims} define the extent of a \eqn{nd}-dimensional box around the space of interest.
The value of \code{objective} is calculated at each of the \eqn{res^{nd}}{res^nd} grid points and and optimum value located (see \code{\link{revisit}}).
In the next iteration the new search box is centered on the location of the optimum value, and the edges are shrunk so their length is \code{rat} * the length in the previous step.
If the limits of any of the parameters extend beyond those in \code{lims}, they are pushed in to fit (preserving the difference between them).

\code{plot_findit} plots the values of the parameters and the objective function as a function of the number of iterations.
}

\value{
\code{findit} returns a list having class \code{findit} with elements \code{value} (values of the parameters, and value of the objective function, at each iteration), \code{lolim} (lower limits of the parameters) and \code{hilim} (upper limits of the parameters).
}

\seealso{ \code{demo("findit")} and \code{test-findit.R} for examples. }

\concept{Extended workflow}
