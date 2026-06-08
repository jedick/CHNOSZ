\encoding{UTF-8}
\name{JANAF.to.OBIGT}
\alias{JANAF.to.OBIGT}
\title{Import thermodynamic data from NIST-JANAF tables}
\description{
Parse JANAF table for standard thermodynamic properties, fit a Cp equation, and add the parameters to OBIGT.
Only gases are currently supported.
}

\usage{
  JANAF.to.OBIGT(file, abbrv = NULL, T_max = 1500, MAE_max = 1, plot_Cp = FALSE)
}

\arguments{
  \item{file}{character, path to JANAF txt file (can be a URL)}
  \item{abbrv}{character, abbreviation (or alternative name) of species for \code{abbrv} column in OBIGT}
  \item{T_max}{numeric, maximum T (K) for fitting Cp equation}
  \item{MAE_max}{numeric, maximum mean absolute error allowed for Cp equation}
  \item{plot_Cp}{logical, whether to make a plot of actual and predicted Cp values}
}

\details{

This function parses the \code{file} to obtain values of the species name, formula, state,
standard thermodynamic properties at 25 \degC and 1 bar (G, H, S, Cp), and Cp coefficients.
A five-parameter Cp equation is used (see \code{\link{thermo}}) to obtain the coefficients (a, b, c, d, e).
The parameters are then added to OBIGT using \code{\link{mod.OBIGT}}.
If a species with the same name and state already exists, its parameters are replaced.

An error is produced for species that are not gases, or if the MAE of the Cp fit is greater than \code{MAE_max}.

}

\section{Note}{

Values of standard G and H calculated by \code{\link{subcrt}} at temperatures other than 25 \degC differ from those listed in the JANAF tables.
The reason is that \code{subcrt} returns \emph{apparent} Gibbs energies and enthalpies of formation,
which neglects the temperature effects for elements in their reference states
(see discussion of the Benson-Helgeson convention on p. 60 of Anderson, 2005).

}

\value{
The species index of the new species in OBIGT.
}

\examples{
# Test with a local file
file <- system.file("extdata/misc/C-127.txt", package = "CHNOSZ")
# Uncomment to use a URL instead
#file <- "https://janaf.nist.gov/tables/C-127.txt"
# Add parameters to OBIGT and make a Cp plot
inew <- JANAF.to.OBIGT(file, abbrv = "acetylene", plot_Cp = TRUE)
# Show calculated properties
subcrt(inew)
}

\references{
Anderson GM. 2005. \emph{Thermodynamics of Natural Systems}. 2nd ed. Cambridge, UK: Cambridge University Press. 648 p. \url{https://worldcat.org/title/474880901}

NIST. 1998. \emph{NIST-JANAF Thermochemical Tables}. \doi{10.18434/T42S31}
}
