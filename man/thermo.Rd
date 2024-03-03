\encoding{UTF-8}
\name{thermo}
\alias{reset}
\alias{OBIGT}
\alias{thermo}
\alias{CHNOSZ}
\title{Thermodynamic database and system settings}
\description{
Run \code{reset()} to reset all of the data used in CHNOSZ to default values.
This includes the computational settings, thermodynamic database, and system settings (chemical species).

The system settings are changed using \code{\link{basis}} and \code{\link{species}}.
To clear the system settings (the default, i.e. no species loaded), run \code{basis("")}; to clear only the formed species, run \code{species(delete = TRUE)}

The thermodynamic database is changed using \code{\link{add.OBIGT}} and \code{\link{mod.OBIGT}}.
To restore the default database without altering the species settings, run \code{OBIGT()}.

The computational settings are changed using \code{\link{water}}, \code{\link{P.units}}, \code{\link{T.units}}, \code{\link{E.units}}, and some other commands (e.g. \code{\link{mod.buffer}}).

All the data are stored in the \code{thermo} data object in an environment named \code{CHNOSZ}.
\code{thermo()} is a convenience function to access or modify parts of this object, in particular some computational settings, for example, \code{thermo("opt$ideal.H" = FALSE)} (see \code{\link{nonideal}}).

The source data are provided with CHNOSZ as \code{*.csv} files in the \code{extdata/thermo} and \code{extdata/OBIGT} directories of the package.
These files are used to build the \code{thermo} object, which is described below.
}

\usage{
  reset()
  OBIGT(no.organics = FALSE)
  thermo(...)
}

\arguments{
  \item{no.organics}{logical, load the database without data files for organic species (NOTE: \CH4 is listed as an \dQuote{inorganic} species)?}
  \item{...}{list, one or more arguments whose names correspond to the setting to modify}
}

\format{

  \itemize{
     
    \item \code{thermo()$opt} 
    List of computational settings.
    Square brackets indicate default values.
    Note that the units of \code{G.tol} and \code{Cp.tol} depend on the \code{E_units} for each species in \code{thermo()$OBIGT}.
    Therefore, species with \code{E_units} of \samp{J} have a lower absolute tolerance for producing messages (because 4.184 J = 1 cal).
    \tabular{lll}{
      \code{E.units} \tab character \tab The user's units of energy ([\samp{J}] or \samp{cal}) (see \code{\link{subcrt}})\cr
      \code{T.units} \tab character \tab The user's units of temperature ([\samp{C}] or \samp{K})\cr
      \code{P.units} \tab character \tab The user's units of pressure ([\samp{bar}] or \samp{MPa})\cr
      \code{state} \tab character \tab The default physical state for searching species [\samp{aq}] (see \code{\link{info}})\cr
      \code{water} \tab character \tab Computational option for properties of water ([\samp{SUPCRT}] or \samp{IAPWS}; see \code{\link{water}})\cr
      \code{G.tol} \tab numeric \tab Difference in G above which \code{\link{check.GHS}} produces a message (cal mol\ifelse{latex}{\eqn{^{-1}}}{\ifelse{html}{\out{<sup>-1</sup>}}{^-1}}) [100]\cr
      \code{Cp.tol} \tab numeric \tab Difference in Cp above which \code{\link{check.EOS}} produces a message (cal K\ifelse{latex}{\eqn{^{-1}}}{\ifelse{html}{\out{<sup>-1</sup>}}{^-1}} mol\ifelse{latex}{\eqn{^{-1}}}{\ifelse{html}{\out{<sup>-1</sup>}}{^-1}}) [1]\cr
      \code{V.tol} \tab numeric \tab Difference in V above which \code{\link{check.EOS}} produces a message (cm\ifelse{latex}{\eqn{^{3}}}{\ifelse{html}{\out{<sup>3</sup>}}{^3}} mol\ifelse{latex}{\eqn{^{-1}}}{\ifelse{html}{\out{<sup>-1</sup>}}{^-1}}) [1]\cr
      \code{varP} \tab logical \tab Use variable-pressure standard state for gases? [\code{FALSE}] (see \code{\link{subcrt}})\cr
      \code{IAPWS.sat} \tab character \tab State of water for saturation properties [\samp{liquid}] (see \code{\link{util.water}})\cr
      \code{paramin} \tab integer \tab Minimum number of calculations to launch parallel processes [1000] (see \code{\link{palply}}) \cr
      \code{ideal.H} \tab logical \tab Should \code{\link{nonideal}} ignore the proton? [\code{TRUE}] \cr
      \code{ideal.e} \tab logical \tab Should \code{\link{nonideal}} ignore the electron? [\code{TRUE}] \cr
      \code{nonideal} \tab character \tab Option for charged species in \code{\link{nonideal}} [\code{Bdot}] \cr
      \code{Setchenow} \tab character \tab Option for neutral species in \code{\link{nonideal}} [\code{bgamma0}] \cr
      \code{Berman} \tab character \tab User data file for mineral parameters in the Berman equations [\code{NA}] \cr
      \code{maxcores} \tab numeric \tab Maximum number of cores for parallel calculations with \code{\link{palply}} [\code{2}] \cr
      \code{ionize.aa} \tab numeric \tab Calculate properties of ionized proteins when \Hplus is in basis species (see \code{\link{affinity}}) [\code{TRUE}]
}

    \item \code{thermo()$element}
  Dataframe containing the thermodynamic properties of elements taken from Cox et al., 1989, Wagman et al., 1982, and (for Am, Pu, Np, Cm) Thoenen et al., 2014.
  The standard molal entropy (\eqn{S}(\code{Z})) at 25 \degC and 1 bar for the \dQuote{element} of charge (\code{Z}) was calculated from \eqn{S}(H2,g) + 2\eqn{S}(\code{Z}) =  2\eqn{S}(H+), where the standard molal entropies of H2,g and H+ were taken from Cox et al., 1989.
  The mass of \code{Z} is taken to be zero.
  Accessing this data frame using \code{\link{mass}} or \code{\link{entropy}} will select the first entry found for a given element; i.e., values from Wagman et al., 1982 will only be retrieved if the properties of the element are not found from Cox et al., 1989.
      \tabular{lll}{
      \code{element}  \tab character  \tab Symbol of element\cr
      \code{state}  \tab character \tab Stable state of element at 25 \degC and 1 bar\cr
      \code{source} \tab character \tab Source of data\cr
      \code{mass}  \tab numeric \tab Mass of element (in natural isotopic distribution;\cr
      \tab \tab referenced to a mass of 12 for \S{12}C)\cr
      \code{s}   \tab numeric \tab Entropy of the compound of the element in its stable\cr
      \tab \tab state at 25 \degC and 1 bar (cal K\S{-1} mol\S{-1})\cr
      \code{n}  \tab numeric \tab Number of atoms of the element in its stable\cr
      \tab \tab compound at 25 \degC and 1 bar
    }

    \item \code{thermo()$OBIGT}

  This dataframe is a thermodynamic database of standard molal thermodynamic properties and equations of state parameters of species.
  \dQuote{OrganoBioGeoTherm} is the name of a Windows GUI interface to SUPCRT92 that was produced in Harold C. Helgeson's Laboratory of Theoretical Geochemistry and Biogeochemistry at the University of California, Berkeley.
  The \acronym{OBIGT} database was originally distributed with that program, and was the starting point for the database in CHNOSZ.

  Note the following database conventions:
  \itemize{
     \item The combination of \code{name} and \code{state} defines a species in \code{thermo()$OBIGT}. A species cannot be duplicated (this is checked when running \code{reset()}).
     \item English names of inorganic gases are used only for the gas state. The dissolved species is named with the chemical formula. Therefore, \code{info("oxygen")} refers to the gas, and \code{info("O2")} refers to the aqueous species.
  } 

  Each entry is referenced to one or two literature sources listed in \code{thermo()$refs}.
  Use \code{\link{thermo.refs}} to look up the citation information for the references.
  See the vignette \viglink{OBIGT} for a complete description of the sources of data.

The identifying characteristics of species and their standard molal thermodynamic properties at 25 \degC and 1 bar are located in the first 13 columns of \code{thermo()$OBIGT}:

    \tabular{lll}{
      \code{name}     \tab character \tab Species name\cr
      \code{abbrv}    \tab character \tab Species abbreviation\cr
      \code{formula}  \tab character \tab Species formula\cr
      \code{state}    \tab character \tab Physical state\cr
      \code{ref1}  \tab character \tab Primary source\cr
      \code{ref2}  \tab character \tab Secondary source\cr
      \code{date}     \tab character \tab Date of data entry (\href{https://en.wikipedia.org/wiki/ISO_8601}{ISO 8601} extended format)\cr
      \code{model}    \tab character \tab Model for thermodynamic properties of the species\cr
      \code{E_units}     \tab character \tab Units of energy: \samp{J} for Joules or \samp{cal} for calories\cr
      \code{G}        \tab numeric   \tab Standard molal Gibbs energy of formation\cr
      \tab \tab from the elements (J|cal mol\S{-1})\cr
      \code{H}        \tab numeric   \tab Standard molal enthalpy of formation\cr
      \tab \tab from the elements (J|cal mol\S{-1})\cr
      \code{S}        \tab numeric   \tab Standard molal entropy (J|cal mol\S{-1} K\S{-1})\cr
      \code{Cp}       \tab numeric   \tab Standard molal isobaric heat capacity (J|cal mol\S{-1} K\S{-1})\cr
      \code{V}	      \tab numeric   \tab Standard molal volume (cm\S{3} mol\S{-1})
    }

\code{model} must be one of \samp{H2O}, \samp{HKF}, \samp{CGL}, \samp{Berman}, \samp{AD}, or \samp{DEW}.
\samp{H2O} is reserved for liquid water, the properties of which are calculated using one of several available models (see \code{\link{water}}).
Most aqueous species use \samp{HKF} (the revised Helgeson-Kirkham-Flowers model).
Properties of aqueous species with \code{model} set to \samp{AD} are calculated using the Akinfiev-Diamond model, and those with \samp{DEW} are calculated using the \code{\link{DEW}} model.
Many minerals in the default database use the \samp{Berman} model (see \code{\link{Berman}}).
All other species use \samp{CGL} (general crystalline, gas, liquid model).
Properties of these species are calculated using a heat capacity function with up to six terms; the exponent on the final term can be defined in the database (see below).

The meanings of the remaining columns depend on the model for each species.
The names of these columns are compounded from those of the parameters in the HKF equations of state and general heat capacity equation; for example, column 13 is named \code{a1.a}.
Scaling of the values by integral powers of ten (i.e., orders of magnitude; OOM) for the HKF parameters (this also includes the DEW model) is based on the usual (but by no means universal) convention in the literature.

Columns 15-22 for aqueous species (parameters in the revised HKF equations of state).
NOTE: Most older papers use units of calories for these parameters, so \samp{cal} is listed here; the actual units for each species are set in the \code{E_units} column.

    \tabular{lll}{
      \code{a1} \tab numeric \tab \a1 * 10 (cal mol\S{-1} bar\S{-1})\cr
      \code{a2} \tab numeric \tab \a2 * 10\S{-2} (cal mol\S{-1})\cr
      \code{a3} \tab numeric \tab \a3 (cal K mol\S{-1} bar\S{-1})\cr
      \code{a4} \tab numeric \tab \a4 * 10\S{-4} (cal mol\S{-1} K)\cr
      \code{c1} \tab numeric \tab \c1 (cal mol\S{-1} K\S{-1})\cr
      \code{c2} \tab numeric \tab \c2 * 10\S{-4} (cal mol\S{-1} K)\cr
      \code{omega} \tab numeric \tab \omega * 10\S{-5} (cal mol\S{-1})\cr
      \code{Z}  \tab numeric \tab Charge
    }

Columns 15-22 for crystalline, gas and liquid species (\Cp = \eqn{a} + \eqn{bT} + \eqn{cT}\S{-2} + \eqn{dT}\S{-0.5} + \eqn{eT}\S{2} + \eqn{fT}\S{lambda}).
NOTE: As of CHNOSZ 2.0.0, OOM scaling for heat capacity coefficients has been removed, and new entries use units of Joules unless indicated by setting \code{E_units} to \samp{cal}.

    \tabular{lll}{
      \code{a} \tab numeric \tab \eqn{a} (J K\S{-1} mol\S{-1})\cr
      \code{b} \tab numeric \tab \eqn{b} (J K\S{-2} mol\S{-1})\cr
      \code{c} \tab numeric \tab \eqn{c} (J K mol\S{-1})\cr
      \code{d} \tab numeric \tab \eqn{d} (J K\S{-0.5} mol\S{-1})\cr
      \code{e} \tab numeric \tab \eqn{e} (J K\S{-3} mol\S{-1})\cr
      \code{f} \tab numeric \tab \eqn{f} (J K\S{-lambda-1} mol\S{-1})\cr
      \code{lambda} \tab numeric \tab \lambda (exponent on the \eqn{f} term)\cr
      \code{T} \tab numeric \tab Positive value: Temperature (K) of polymorphic transition or phase stability limit\cr
      \code{T} \tab numeric \tab Negative value: Opposite of temperature (K) limit of \Cp equation (see \viglink{FAQ} for details)\cr
    }

Columns 15-17 for aqueous species using the Akinfiev-Diamond model. Note that the \code{c} column is used to store the \eqn{\xi}{xi} parameter. Columns 18-22 are not used.

    \tabular{lll}{
      \code{a} \tab numeric \tab \eqn{a} (cm\S{3} g\S{-1}) \cr
      \code{b} \tab numeric \tab \eqn{b} (cm\S{3} K\S{0.5} g\S{-1}) \cr
      \code{c} \tab numeric \tab \eqn{\xi}{xi} \cr
      \code{d} \tab numeric \tab \eqn{XX1} NA \cr
      \code{e} \tab numeric \tab \eqn{XX2} NA \cr
      \code{f} \tab numeric \tab \eqn{XX3} NA \cr
      \code{lambda} \tab numeric \tab \eqn{XX4} NA \cr
      \code{Z} \tab numeric \tab \eqn{Z} NA \cr
    }

    \item \code{thermo()$refs}
    References for thermodynamic data.
    \tabular{lll}{
      \code{key} \tab character \tab Source key\cr
      \code{author} \tab character \tab Author(s)\cr
      \code{year} \tab character \tab Year\cr
      \code{citation} \tab character \tab Citation (journal title, volume, and article number or pages; or book or report title)\cr
      \code{note} \tab character \tab Short description of the compounds or species in this data source\cr
      \code{URL} \tab character \tab URL\cr
    }

    \item \code{thermo()$buffers}

    Dataframe which contains definitions of buffers of chemical activity. Each named buffer can be composed of one or more species, which may include any species in the thermodynamic database and/or any protein. The calculations provided by \code{\link{buffer}} do not take into account polymorphic transitions of minerals, so individual polymorphs of such minerals must be specified in the buffers.
    \tabular{lll}{
      \code{name} \tab character \tab Name of buffer\cr
      \code{species} \tab character \tab Name of species\cr
      \code{state} \tab character \tab Physical state of species\cr
      \code{logact} \tab numeric \tab Logarithm of activity (fugacity for gases)
    }

    \item \code{thermo()$protein}
    Data frame of amino acid compositions of selected proteins. Most of the compositions were taken from the SWISS-PROT/UniProt online database (Boeckmann et al., 2003) and the protein and organism names usually follow the conventions adopted there. In some cases different isoforms of proteins are identified using modifications of the protein names; for example, \samp{MOD5.M} and \code{MOD5.N} proteins of \samp{YEAST} denote the mitochondrial and nuclear isoforms of this protein. See \code{\link{pinfo}} to search this data frame by protein name, and other functions to work with the amino acid compositions.
    \tabular{lll}{
      \code{protein} \tab character \tab Identification of protein\cr
      \code{organism} \tab character \tab Identification of organism\cr
      \code{ref} \tab character \tab Reference key for source of compositional data\cr
      \code{abbrv} \tab character \tab Abbreviation or other ID for protein\cr
      \code{chains} \tab numeric \tab Number of polypeptide chains in the protein\cr
      \code{Ala}\dots\code{Tyr} \tab numeric \tab Number of each amino acid in the protein
    }

    \item \code{thermo()$groups}
    This is a dataframe with 22 columns for the amino acid sidechain, backbone and protein backbone groups ([Ala]..[Tyr],[AABB],[UPBB]) whose rows correspond to the elements C, H, N, O, S. It is used to quickly calculate the chemical formulas of proteins that are selected using the \code{iprotein} argument in \code{\link{affinity}}.

    \item \code{thermo()$basis}
    Initially \code{NULL}, reserved for a dataframe written by \code{\link{basis}} upon definition of the basis species. The number of rows of this dataframe is equal to the number of columns in \dQuote{...} (one for each element).
     \tabular{lll}{
        \code{...} \tab numeric \tab One or more columns of stoichiometric\cr
        \tab \tab coefficients of elements in the basis species\cr
        \code{ispecies} \tab numeric \tab Rownumber of basis species in \code{thermo()$OBIGT}\cr
        \code{logact} \tab numeric \tab Logarithm of activity or fugacity of basis species\cr
        \code{state} \tab character \tab Physical state of basis species\cr
     }

    \item \code{thermo()$species}
    Initially \code{NULL}, reserved for a dataframe generated by \code{\link{species}} to define the species of interest. The number of columns in \dQuote{...} is equal to the number of basis species (i.e., rows of \code{thermo()$basis}).
    \tabular{lll}{
       \code{...} \tab numeric \tab One or more columns of stoichiometric\cr
       \tab \tab coefficients of basis species in the species of interest\cr
       \code{ispecies} \tab numeric \tab Rownumber of species in \code{thermo()$OBIGT}\cr
       \code{logact} \tab numeric \tab Logarithm of activity or fugacity of species\cr
       \code{state} \tab character \tab Physical state of species\cr
       \code{name} \tab character \tab Name of species\cr
    }

    \item \code{thermo()$stoich}
    A precalculated stoichiometric matrix for the default database. This is a matrix, not a data frame, and as such can accept duplicated row names, corresponding to chemical formulas of the species. See \code{\link{retrieve}}, and the first test in \code{inst/tinytest/test-retrieve.R} for how to update this.
    \tabular{lll}{
       \code{rownames} \tab character \tab Chemical formulas from \code{thermo()$OBIGT}\cr
       \code{...} \tab numeric \tab Stoichiometry, one column for each element present in any species
    }

    \item \code{thermo()$Bdot_acirc}
    Values of ion size parameter (\AA) for species, taken from the UT_SIZES.REF file of the HCh package (Shvarov and Bastrakov, 1999), which is based on Table 2.7 of Garrels and Christ, 1965.
    This is used in \code{\link{nonideal}} with the default \samp{Bdot} method.
    Custom ion size parameters can be added to this vector; to override a default value for a species, either replace the numeric value for that species or prepend a named numeric value (for duplicated species, the first value is used).
    See \code{demo("yttrium")} for an example of adding and overriding species.

    \item \code{thermo()$Berman}
    A data frame with thermodynamic parameters for minerals in the Berman equations, assembled from files in \samp{extdata/Berman} and used in \code{\link{Berman}}.

  }  % end of itemize with long descriptions

} % end of format

\section{Note on polymorphic transitions}{
To enable the calculation of thermodynamic properties of polymorphic transitions, higher-temperature polymorphs of minerals are listed in OBIGT with states \samp{cr2}, \samp{cr3}, etc.
The standard thermodynamic properties of high-temperature polymorphs at 25 \degC and 1 bar are apparent values that are consistent with given values of enthalpy of transition (where available) at the transition temperature (\eqn{T_{tr}}{Ttr}).
See the \viglink{FAQ} question \dQuote{How can minerals with polymorphic transitions be added to the database?} for details of the retrieval of standard thermodynamic properties of polymorphs used in OBIGT.
}

\seealso{
Other data files, including those supporting the examples and vignettes, are documented separately at \code{\link{extdata}}. 
}

\examples{
## Where are the data files in CHNOSZ?
system.file("extdata", package = "CHNOSZ")
# What files make up OBIGT?
# Note: file names with _aq, _cr, _gas, or _liq
# are used in the default database
dir(system.file("extdata/OBIGT", package = "CHNOSZ"))

## Exploring thermo()$OBIGT
# What physical states are present
unique(thermo()$OBIGT$state)
# Formulas of ten random species
n <- nrow(thermo()$OBIGT)
thermo()$OBIGT$formula[runif(10)*n]

## Adding an element
old <- thermo()$element
# Element symbol, state, source (can be anything),
# mass, entropy, and number in compound
Xprops <- data.frame(element = "X", state = "cr",
  source = "user", mass = 100, s = 100, n = 1)
new <- rbind(old, Xprops)
thermo(element = new)
# Now "X" is recognized as an element in other functions
mass("X10")
# Restore default settings to remove X
reset()
}

\references{
Cox, J. D., Wagman, D. D. and Medvedev, V. A., eds. (1989) \emph{CODATA Key Values for Thermodynamics}. Hemisphere Publishing Corporation, New York, 271 p. \url{https://www.worldcat.org/oclc/18559968}

Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{https://www.worldcat.org/oclc/517586}

Thoenen, T., Hummel, W., Berner, U. and Curti, E. (2014) \emph{The PSI/Nagra Chemical Thermodynamic Database 12/07}. Paul Scherrer Institut. \url{https://www.psi.ch/en/les/database}

Wagman, D. D., Evans, W. H., Parker, V. B., Schumm, R. H., Halow, I., Bailey, S. M., Churney, K. L. and Nuttall, R. L. (1982) The NBS tables of chemical thermodynamic properties. Selected values for inorganic and C\s{1} and C\s{2} organic substances in SI units. \emph{J. Phys. Chem. Ref. Data} \bold{11} (supp. 2), 1--392. \url{https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf}
}

\concept{Thermodynamic data}
