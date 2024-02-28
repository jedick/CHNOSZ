\encoding{UTF-8}
\name{CHNOSZ-package}
\alias{CHNOSZ-package}
\docType{package}
\title{Thermodynamic Calculations and Diagrams for Geochemistry}
\description{
CHNOSZ is a package for thermodynamic calculations, primarily with applications in geochemistry and compositional biology.
It can be used to calculate the standard molal thermodynamic properties and chemical affinities of reactions relevant to geobiochemical processes, and to visualize the equilibrium activities of species on chemical speciation and predominance diagrams.
}

\section{Warm Tips}{
  \itemize{
    \item To view the manual, run \code{help.start()} then select \sQuote{Packages} and \sQuote{CHNOSZ}. Examples in the function help pages can be run by pasting the code block into the \R console.
    \item Also check out the vignette \viglink{anintro} (\emph{An Introduction to CHNOSZ}).
    \item Run the command \code{\link{examples}()} to run all of the examples provided in CHNOSZ. This should take about a minute.
  }
}

\section{Getting Help}{
Each help page (other than this one) has been given one of the following \dQuote{concept index entries}:
\itemize{
  \item Main workflow: \code{\link{info}}, \code{\link{subcrt}}, \code{\link{basis}}, \code{\link{species}}, \code{\link{affinity}}, \code{\link{equilibrate}}, \code{\link{diagram}}
  \item Extended workflow: \code{\link{swap.basis}}, \code{\link{buffer}}, \code{\link{mosaic}}, \code{\link{EOSregress}}
  \item Thermodynamic data: \code{\link{data}}, \code{\link{extdata}}, \code{\link{add.OBIGT}}, \code{\link{util.data}}
  \item Thermodynamic calculations: \code{\link{util.formula}}, \code{\link{makeup}}, \code{\link{util.units}}, \code{\link{Berman}}, \code{\link{nonideal}}, \code{\link{util.misc}}
  \item Water properties: \code{\link{water}}, \code{\link{util.water}}, \code{\link{DEW}}, \code{\link{IAPWS95}}
  \item Protein properties: \code{\link{protein.info}}, \code{\link{add.protein}}, \code{\link{util.protein}}, \code{\link{util.seq}}, \code{\link{ionize.aa}}
  \item Other tools: \code{\link{examples}}, \code{\link{taxonomy}}
  \item Utility functions: \code{\link{util.expression}}, \code{\link{util.plot}}, \code{\link{util.array}}, \code{\link{util.list}}, \code{\link{palply}}
}
These concept entries are visible to \code{\link{help.search}} (aka \code{??}).
For example, help pages related to thermodynamic data can be listed using \code{??"thermodynamic data"}.
}

\section{Warning}{
All thermodynamic data and examples are provided on an as-is basis.
It is up to you to check not only the accuracy of the data, but also the \emph{suitability of the data AND computational techniques} for your problem.
By combining data taken from different sources, it is possible to build an inconsistent and/or nonsensical calculation.
An attempt has been made to provide a default database (OBIGT) that is internally consistent, but no guarantee can be made.
If there is any doubt about the accuracy or suitability of data for a particular problem, please consult the primary sources (see \code{\link{thermo.refs}}).
}

\section{Acknowledgements}{
  This package would not exist without the scientific influence and friendship of the late Professor Harold C. Helgeson.
  The \file{src/H2O92D.f} file with Fortran code for calculating the thermodynamic and electrostatic properties of \H2O is modified from the \acronym{SUPCRT92} package (Johnson et al., 1992).

  Work on CHNOSZ at U.C. Berkeley from ca. 2003 to 2008 was supported in part by research grants to HCH from the U.S. National Science Foundation and Department of Energy.
  In 2009--2011, development of this package was partially supported by NSF grant EAR-0847616 to JMD.
}

\references{
Johnson, J. W., Oelkers, E. H. and Helgeson, H. C. (1992) SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\degC. \emph{Comp. Geosci.} \bold{18}, 899--947. \doi{10.1016/0098-3004(92)90029-Q}
}

\keyword{package}
