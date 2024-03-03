\encoding{UTF-8}
\name{util.units}
\alias{util.units}
\alias{P.units}
\alias{T.units}
\alias{E.units}
\alias{convert}
\title{Functions to convert units}
\description{
These functions convert values between units and set the user's preferred units.
}

\usage{
  P.units(units = NULL)
  T.units(units = NULL)
  E.units(units = NULL)
  convert(value, units, T = 298.15, P = 1, pH = 7, logaH2O = 0)
}

\arguments{
  \item{units}{character, name of units to set or convert to/from}
  \item{value}{numeric, value(s) to be converted}
  \item{T}{numeric, temperature (Kelvin), used in \samp{G}-\samp{logK}, \samp{pe}-\samp{Eh} and \samp{logfO2}-\samp{E0} conversions}
  \item{P}{numeric, pressure (bar), used in \samp{logfO2}-\samp{E0} conversions}
  \item{pH}{numeric, pH, used in \samp{logfO2}-\samp{E0} conversions}
  \item{logaH2O}{numeric, logarithm of activity of water, used in \samp{logfO2}-\samp{E0} conversions}
}

\details{
The units settings are used by \code{\link{subcrt}}, \code{\link{affinity}}, and \code{\link{diagram}} to accept input in or convert output to the units desired by the user.
The settings, which can be queried or changed with \code{T.units}, \code{E.units} and \code{P.units}, refer to the units of temperature (\code{C} or \code{K}), energy (\code{J} or \code{cal}), and pressure (\code{bar} or \code{MPa}).
(The first value in each of those pairs refers to the default units).

The actual units conversions are handled by \code{convert}, through which \code{values} are transformed into destination \code{units} (names not case sensitive).
The possible conversions and settings for the \code{units} argument are shown in the following table.
Note that \samp{Eh} and \samp{E0} both stand for the value of Eh (oxidation-reduction potential in volts); they have different names so that one can choose to convert between Eh and either \samp{pe} or \samp{logfO2}. 

  \tabular{lll}{
     property \tab units \tab setting of \code{units} argument \cr
     temperature \tab \degC, K \tab \code{C}, \code{K} \cr
     pressure \tab bar, MPa \tab \code{bar}, \code{MPa} \cr
     energy \tab cal, J \tab \code{cal}, \code{J} \cr
     energy \tab J, cm\eqn{^3}{^3} bar \tab \code{joules}, \code{cm3bar} \cr
     energy \tab J, [none] \tab \code{G}, \code{logK} \cr
     oxidation potential \tab volt, [none] \tab \code{Eh}, \code{pe} \cr
     oxidation potential \tab volt, [none] \tab \code{E0}, \code{logfO2} \cr
  }

Another use of the function is to convert the results from \code{\link{solubility}} into parts per billion, million, or thousand.
These destination units are specified by \samp{ppb}, \samp{ppm}, or \samp{ppt}.
Additionally, the logarithms can be chosen with \samp{logppb}, \samp{logppm}, and \samp{logppt}.
See \code{\link{demo}("contour")} and \code{\link{demo}("sphalerite")} for examples.
}

\examples{\dontshow{reset()}
## Direct usage of convert
# Temperature (Kelvin) to degrees C
convert(273.15, "C")
# Temperature (degrees C) to Kelvin
convert(100, "K")
# Gibbs energy (J mol-1) to/from logK
convert(1000, "logK")		
convert(1000, "logK", T = 373.15)
convert(1, "G") 
# Eh (volt) to pe 
convert(-1, "pe")		
convert(-1, "pe", T = 373.15)
# logfO2 to E0 (volt)
convert(-80, "E0")             
convert(-80, "E0", pH = 5)
convert(-80, "E0", pH = 5, logaH2O = -5)
# Convert from calories to Joules
convert(1, "J")    # 1 cal = 4.184 J
# Convert from Joules to calories
convert(1, "cal")  # 1 J = 0.239 cal
# Convert cm3bar to Joules
convert(10, "joules") # 10 cm3.bar = 1 J

## Setting the units
# Make K the units for temperature arguments to subcrt() and affinity()
T.units("K") 
# Return to default - degrees C
T.units("C")
}

\concept{Thermodynamic calculations}
