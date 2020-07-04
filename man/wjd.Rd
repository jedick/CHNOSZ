\encoding{UTF-8}
\name{wjd}
\alias{wjd}
\alias{element.potentials}
\alias{is.near.equil}
\alias{guess}
\alias{run.wjd}
\alias{run.guess}
\alias{equil.potentials}
\title{Gibbs Energy Minimization by Steepest Descent}
\description{
  Find the quantities of chemical species, subject to constant elemental bulk composition of the system, that minimize the Gibbs energy of the system.
}

\usage{
  wjd(  
    A = matrix(c(
      1,2,2,0,0,1,0,0,0,1,
      0,0,0,1,2,1,1,0,0,0,
      0,0,1,0,0,0,1,1,2,1),ncol=3,
      dimnames=list(NULL,c("H","N","O"))),
    G0.RT = c(
      -10.021,-21.096,-37.986,-9.846,-28.653,
      -18.918,-28.032,-14.640,-30.594,-26.111),
    Y = c(0.1,0.35,0.5,0.1,0.35,0.1,0.1,0.1,0.1,0.1),
    P = 51,
    nlambda = 101,
    imax = 10,
    Gfrac = 1e-7
  )
  element.potentials(w, plot.it=FALSE, iplot=1:ncol(w$A))
  is.near.equil(w, tol=0.01, quiet=FALSE)
  guess(  
    A = matrix(c(
      1,2,2,0,0,1,0,0,0,1,
      0,0,0,1,2,1,1,0,0,0,
      0,0,1,0,0,0,1,1,2,1),ncol=3,
      dimnames=list(NULL,c("H","N","O"))),
    B = c(2,1,1), method="stoich", minX = 0.001, iguess = 1, ic = NULL
  )
  run.wjd(ispecies, B = NULL, method = "stoich",
    Y = run.guess(ispecies, B, method), P=1, T=25, nlambda=101, imax = 10,
    Gfrac = 1e-7, tol = 0.01)
  run.guess(ispecies, B = NULL, method = "stoich", iguess = NULL)
  equil.potentials(w, tol=0.01, T=25)
}

\arguments{
  \item{A}{matrix, chemical formulas of the species (elements on columns)}
  \item{G0.RT}{numeric, the Gibbs energies / RT, at a single temperature (length equal to number of species)}
  \item{Y}{numeric, initial solution, a positive set of values (numbers of moles, length equal to number of species)}
  \item{P}{numeric, pressure in atmospheres}
  \item{nlambda}{numeric, number of values of fractional distance change (\eqn{\lambda}{lambda}) tested at each step.}
  \item{imax}{numeric, maximum number of iterations}
  \item{Gfrac}{numeric, Gibbs energy change of system, as a fraction of total system energy in previous step, below which iterations will stop}
  \item{w}{list, output from \code{wjd}}
  \item{plot.it}{logical, make a plot?}
  \item{iplot}{numeric, which elements for which to make plots}
  \item{tol}{numeric, maximum difference in chemical potentials that counts as equilibrium}
  \item{quiet}{logical, don't output messages?}
  \item{B}{numeric, numbers of moles of the elements}
  \item{method}{character, method used for generating an initial solution}
  \item{minX}{numeric, minimum mole number for 'central' method}
  \item{iguess}{numeric, which guess to return}
  \item{ic}{numeric, which combination(s) of variable species to use (NULL for all)}
  \item{ispecies}{numeric, species indices (rownumbers of \code{thermo}$OBIGT)}
  \item{T}{numeric, temperature in \degC}
}

\details{

\code{wjd} implements the steepest descent algorithm for Gibbs energy minimization in a closed system described by White et al., 1958. 
\dQuote{Gibbs energy} (G) referred to here is the same as the \dQuote{free energy} (F) used by White et al., 1958. 
\code{wjd} itself is independent of other functions or datasets in CHNOSZ, but \code{run.wjd} and \code{run.guess} are provided to make access to the thermodynamic database in CHNOSZ easier.

The default values of \code{A}, \code{G0.RT}, \code{Y} and \code{P} correspond to the example problem described by White et al., 1958 for gases in the H, N, O system at 3500 K. 
Note that for this, and for any other equilibrium problem that can be simulated using this function, the mole quantities in \code{Y} must all be positive numbers. 
Operationally, this vector defines not only the \dQuote{initial solution} but also the bulk composition of the system; it is not possible to define the bulk composition using mole numbers of elements alone. 
The \code{dimnames} attribute in the default value for \code{A} gives the names of the elements; this attribute is not necessary for the function to operate, but is used in the examples below to help label the plots.

  White et al., 1958 describe in detail the computation of the direction of steepest descent by means of Lagrange multipliers. They propose an iterative solution to the energy minimization problem, where at each step the mole numbers of species are recomputed and a new steepest descent direction calculated from there. However, the authors only give general guidelines for computing the value of \eqn{\lambda}{lambda}, that is, the fraction of the total distance the system actually moves in the direction of steepest descent from the current point at each iteration. The two constraints given for determining the value of \eqn{\lambda}{lambda} are that all mole numbers of species are positive, and the Gibbs energy of the system actually decreases (the minimum point is not passed). In the code described here, at each iteration the minimum value of lambda, not exceeding unity, that violates the first condition is computed directly (a value of one is assigned if the mole numbers remain positive through the entire range). With the default setting of \code{nlambda}, 101 values of \eqn{\lambda}{lambda} at even intervals from 0 to this maximum permissible value are tested for the second condition at each iteration, and the highest conforming value is selected. If a value of 0 occurs, it means that the algorithm has reached an endpoint independently of the iteration and convergence tests (\code{rho} and \code{Gfrac}; see below). If this occurs, the value of \code{nlambda} might have to be increased depending on the user's needs.

  The number of iterations is controlled by \code{imax} and \code{Gfrac}. The maximum number of iterations is set by \code{imax}; it can even be zero, though such a setting might only be useful in combination with \code{element.potentials} to characterize the initial state of a minimization problem. Within the limit of \code{imax}, the iterations continue until the magnitutde of the change in total Gibbs energy of the system, as a fraction of the system's energy in the previous iteration, is lower than the value specified in \code{Gfrac}. For the first example below, the default setting of \code{Gfrac} causes the algorithm to stop after six iterations.

  Using the output of \code{wjd}, provided in the argument \code{w}, \code{element.potentials} calculates the chemical potentials of the elements in the system. It does so by combining the values of \code{G0.RT} of species with the inverses of stoichiometric matrices of combinations of species whose elemental compositions are linearly independent from each other. These possible combinations are constructed using the function \code{\link{invertible.combs}}. The value returned by \code{element.potentials} is a matrix, with each column corresponding to a different element and each row to a different combination of species. The entries in the matrix are the chemical potentials of the elements divided by \eqn{RT}{RT}. If \code{plot.it} is set to TRUE, the chemical potentials of the elements are plotted as a function of species combination number, with as many plots as elements, unless \code{iplot} is set to another value (e.g. \samp{c(1,3)} for only elements 1 and 3). In the first example below, the number of unique combinations of species is 120, but only 76 of these combinations provide stoichiometric independence.

There is no guarantee that the algorithm will converge on a global (or even be close to a local) minimum.
However, some tests are available to help assess the likelihood that a solution is close to equilibrium.
A necessary condition of equilibrium is that the chemical potentials of the elements be independent of the particular combination of species used to compute them.
\code{is.near.equil} compares the chemical potentials for each element, computed using \code{element.potentials}, with the value of \code{tol}.
If, for each element, the range of potentials/RT (difference between minimum and maximum) is less than \code{tol}, the result is TRUE, otherwise the function returns FALSE, and prints a message unless \code{quiet} is TRUE.
The default value of \code{tol} corresponds to an energy of 0.01 * 1.9872 * 298.15 = ca. 6 cal/mol at 25 \degC.

One of the constraints of the algorithm coded in \code{wjd} is that the initial solution, and all iterations, require positive (non-zero) numbers of moles of every species.
Often, when investigating an equilibrium problem, the stoichiometric constraints are expressed most readily in terms of bulk composition -- numbers of moles of each element.
\code{guess} is a function to make initial guesses about the numbers of moles of all species in the system subject to the positivity constraints.
Its system-specific arguments are the stoichiometric matrix \code{A} (as defined above for \code{wjd}) and the bulk composition vector \code{B}, giving the number of moles of each element, in the same order that the elements appear in \code{A}.
The first \code{method} available in \code{guess} generates the \sQuote{central} solution of the system of linear equations using the \code{\link[limSolve]{xranges}} function from \CRANpkg{limSolve}. The central solution is the mean of ranges of unknowns. The inequality constraint, or minimum number of moles of any species, is given by \code{minX}.

The second method for \code{guess} \sQuote{stoich} (and the default for \code{run.guess} and \code{run.wjd}) is to test successive combinations of species whose elemental compositions are linearly independent.
The linearly independent combinations tested are all those from \code{\link{invertible.combs}} if \code{ic} is NULL, or only those identified by \code{ic}.
Each combination is referred to as the \sQuote{variable} species; the moles of all \sQuote{other} species are set to a single value.
This value is determined by the constraint that the greatest proportion, relative to the bulk composition in \code{B}, of any element contributed by all the \sQuote{other} species is equal to a value in \code{max.frac} (see code).
The function tests nine hard-coded values of \code{max.frac} from 0.01 to 0.99, at each one solving for the moles of the \sQuote{variable} species that make up the difference in numbers of moles of elements.
If the numbers of moles of all the \sQuote{variable} species is positive, the guess is accepted.
The first accepted guess is returned if \code{iguess} is 1 (the default); other values of \code{iguess} indicate which guess to return.
If \code{iguess} is NULL, all results are returned in a list, with non-successful guesses indicated by NA.
In the first example below, of the 76 combinations of species whose elemental compositions are linearly independent, 32 yield guesses following this method.

\code{run.wjd} is a wrapper function to run \code{wjd}, provided the \code{ispecies} in the thermodynamic database (\code{\link{thermo}$OBIGT}), the chemical composition of the system in \code{B}, and the pressure \code{P} and temperature \code{T}; the latter are passed to \code{\link{subcrt}} (with \code{exceed.Ttr = TRUE}) to generate the values of \code{G0.RT} for \code{wjd}.
Alternatively to \code{B}, the initial guess of numbers of moles of species can be provided in \code{Y}; otherwise as many combinations of \code{Y} as returned from \code{run.guess} are tested until one is found that \code{is.near.equil}.
The function gives an error if none of the guesses in \code{Y} is near equilibrium, within the tolerance set by \code{tol}.

\code{run.guess} is a wrapper function to call \code{guess} using the stoichiometric matrix \code{A} built from the \code{ispecies} indices in the thermodynamic database.

\code{equil.potentials} returns the average (\code{\link{colMeans}}) of \code{element.potentials(w)}, or NULL if \code{is.near.equil(w, tol=tol)} is FALSE.
The output of this function can be used as the \code{emu} argument for \code{\link{basis.logact}} to calculate the corresponding activities of the basis species.
}

\value{
\code{wjd} returns a list with the problem definition and results: elements \code{A}, \code{G0.RT}, \code{Y}, and \code{P} are as supplied in the arguments; the results are in \code{X} (final mole numbers of species), \code{F.Y} (Gibbs energy of the system at initial conditions and after each iteration), \code{lambda} (value used for \eqn{\lambda}{lambda} at each iteration), and \code{elements} (matrix of moles of elements at initial conditions and after each iteration; iterations on the columns and elements on the rows).
}

\seealso{
\code{demo("wjd")} for a longer example, and \code{\link{invertible.combs}}, used by \code{element.potentials} to find combinations of species that are compositionally independent.
}

\examples{\dontshow{reset()}
## run the function with default settings to reproduce
## the example problem in White et al., 1958
w <- wjd()
# the mole fractions are very close to those shown in the
# last column of Table III in the paper
print(w$X)
# the Gibbs energy of the system decreased, 
# and by a smaller amount, at each iteration
print(diff(w$F.Y))
# there are 76 unique combinations of species that can be
# used to calculate the chemical potentials of the elements
stopifnot(nrow(invertible.combs(w$A))==76)
# what the scatter in those chemical potentials looks like
\dontshow{opar <- par(no.readonly=TRUE)}
ep <- element.potentials(w, plot.it=TRUE)
# the differences in chemical potentials / RT are all less than 0.01
is.near.equil(w)  # TRUE

## run the algorithm for only one iteration
w <- wjd(imax=1)
# the scatter in chemical potentials is much greater
ep <- element.potentials(w, plot.it=TRUE)
# and we're pretty far from equilibrium
is.near.equil(w)  # FALSE
\dontshow{par(opar)}

## test all of the guesses of inititial mole quantities
## provided by guess() using default bulk composition (H2NO)
# 9 of them are not is.near.equil with the tolerance lowered to 0.0001
sapply( 1:32, function(i) 
  is.near.equil(wjd(Y=guess(method = "stoich", iguess=i)), tol=0.0001) )

## using run.wjd(): 20 crystalline amino acids
iaa <- info(aminoacids(""), "cr")
# starting with one mole of each amino acid
w <- run.wjd(iaa, Y=rep(1, 20), imax=20)
# the following is TRUE (FALSE if tol is left at default)
is.near.equil(w, tol=0.012)
# in this assemblage, what are the amino acids 
# in order of increasing abundance?
aminoacids()[order(w$X)]
# because the elements are redistributed among the species,
# the total number of moles of species does not remain constant
sum(w$X)  # <20
}

\references{
White, W. B., Johnson, S. M. and Dantzig, G. B. (1958) Chemical equilibrium in complex mixtures. \emph{J. Chem. Phys.} \bold{28}, 751--755. \url{https://doi.org/10.1063/1.1744264}
}

\concept{Extended workflow}
