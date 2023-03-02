# CHNOSZ/DEW.R
# R functions for the Deep Earth Water model
# 20170924 jmd

# Most code here was translated from VBA macros in the DEW spreadsheet (DEW_May19_2017_11.0.2 .xlsm)
# Comments starting with ' were transferred from the DEW spreadsheet
# In the original, functions return zero for invalid input; here, NA is used instead

# Equations were selected using default options in the DEW spreadsheet:
# Density of water equation    (1) - Zhang & Duan (2005)
# Dielectric Constant Equation (4) - Sverjensky et al. (2014)
# Water Free Energy Equation   (2) - Integral of Volume

# 'Returns the density of water at the input pressure and temperature, in units of g/cm^3.
# 'pressure       - The pressure to calculate the density of water at, in bars
# 'temperature    - The temperature to calculate the density of water at, in Celsius
# 'error          - The density returned will calculate a pressure which differs from the input pressure by the value of "error" or less.
# The default value of error is taken from equations in DEW spreadsheet (table "Calculations")
calculateDensity <- function(pressure, temperature, error = 0.01) {
  myfunction <- function(pressure, temperature) {
    minGuess <- 1E-05
    guess <- 1E-05
    equation <- 1 # 'The maxGuess is dependent on the value of "equation"
    maxGuess <- 7.5 * equation - 5  
    calcP <- 0
    # 'Loop through and find the density
    calculateDensity <- NA
    for(i in 1:50) {
      # 'Calculates the pressure using the specified equation
      calcP <- calculatePressure(guess, temperature)
      # 'If the calculated pressure is not equal to input pressure, this determines a new
      # 'guess for the density based on current guess and how the calculated pressure
      # 'relates to the input pressure. In effect, this a form of a bisection method.
      if(abs(calcP - pressure) > error) {
        if(calcP > pressure) {
          maxGuess <- guess
          guess <- (guess + minGuess) / 2
        } else if(calcP < pressure) {
          minGuess <- guess
          guess <- (guess + maxGuess) / 2
        }
      } else {
        calculateDensity <- guess
        break
      }
    }
    calculateDensity
  }
  # Make input pressure and temperature the same length
  if(length(pressure) < length(temperature)) pressure <- rep(pressure, length.out=length(temperature))
  if(length(temperature) < length(pressure)) temperature <- rep(temperature, length.out=length(pressure))
  # Use a loop to process vectorized input
  sapply(1:length(pressure), function(i) myfunction(pressure[i], temperature[i]))
}

# 'Returns the Gibbs Free Energy of water in units of cal/mol.
# 'pressure           - The pressure to calculate the Gibbs Free Energy at, in bars
# 'temperature        - The temperature to calculate the Gibbs Free Energy at, in Celsius
calculateGibbsOfWater <- function(pressure, temperature) {

  # 'Equation created by Brandon Harrison. This models data for the Gibbs energy at 1 kb as a function of temperature,
  # 'then defines the gibbs free energy as the integral over the volume as a function of temperature.
  myfunction <- function(pressure, temperature) {
    # 'Gibbs Free Energy of water at 1 kb. This equation is a polynomial fit to data as a function of temperature.
    # 'It is valid in the range of 100 to 1000 C.
    GAtOneKb <- 2.6880734E-09 * temperature^4 + 6.3163061E-07 * temperature^3 - 
      0.019372355 * temperature^2 - 16.945093 * temperature - 55769.287

    if(pressure < 1000) {                 # 'Simply return zero, this method only works at P >= 1000 bars
      integral <- NA
    } else if(pressure == 1000) {         # 'Return the value calculated above from the polynomial fit
      integral <- 0
    } else if(pressure > 1000) {          # 'Integrate from 1 kb to P over the volume
      integral <- 0
      # 'Integral is sum of rectangles with this width. This function in effect limits the spacing
      # 'to 20 bars so that very small pressures do not have unreasonably small widths. Otherwise the width
      # 'is chosen such that there are always 500 steps in the numerical integration. This ensures that for very
      # 'high pressures, there are not a huge number of steps calculated which is very computationally taxing.
      spacing <- ifelse((pressure - 1000) / 500 < 20, 20, (pressure - 1000) / 500)
      for(i in seq(1000, pressure, by = spacing)) {
        # 'This integral determines the density only down to an error of 100 bars
        # 'rather than the standard of 0.01. This is done to save computational
        # 'time. Tests indicate this reduces the computation by about a half while
        # 'introducing little error from the standard of 0.01.
        integral <- integral + (18.01528 / calculateDensity(i, temperature, 100) / 41.84) * spacing
      }
    }
    GAtOneKb + integral
  }
  # Make input pressure and temperature the same length
  if(length(pressure) < length(temperature)) pressure <- rep(pressure, length.out=length(temperature))
  if(length(temperature) < length(pressure)) temperature <- rep(temperature, length.out=length(pressure))
  # Use a loop to process vectorized input
  sapply(1:length(pressure), function(i) myfunction(pressure[i], temperature[i]))
}

# 'Returns the Dielectric constant of water at the given density and temperature.
# 'density        - The density of water to use in calculating epsilon, in g/cm^3
# 'temperature    - The temperature to calculate epsilon with, in Celsius
calculateEpsilon <- function(density, temperature) {
  # 'Power Function - Created by Dimitri Sverjensky and Brandon Harrison
  # 'Relevant parameters
  a1 <- -0.00157637700752506
  a2 <- 0.0681028783422197
  a3 <- 0.754875480393944
  b1 <- -8.01665106535394E-05
  b2 <- -0.0687161761831994
  b3 <- 4.74797272182151

  A <- a1 * temperature + a2 * sqrt(temperature) + a3
  B <- b1 * temperature + b2 * sqrt(temperature) + b3

  exp(B) * density ^ A
}

# 'Outputs the value of Q in units of bar^-1
# 'pressure           - The pressure to calculate Q at, in bars
# 'temperature        - The temperature to calculate Q at, in Celsius
# 'density            - The density at the input pressure and temperature, input simply to save time, in g/cm^3
calculateQ <- function(density, temperature) {
  eps <- calculateEpsilon(density, temperature)
  depsdrho <- calculate_depsdrho(density, temperature)
  drhodP <- calculate_drhodP(density, temperature)

  depsdrho * drhodP / eps ^2
}

### Unexported functions ###

# 'Returns the pressure of water corresponding to the input density and temperature, in units of bars.
# 'density        - The density to use in finding a pressure, in g/cm^3
# 'temperature    - The temperature to use in finding a pressure, in Celsius
calculatePressure <- function(density, temperature) {
  m <- 18.01528         # 'Molar mass of water molecule in units of g/mol
  ZD05_R <- 83.144      # 'Gas Constant in units of cm^3 bar/mol/K
  ZD05_Vc <- 55.9480373 # 'Critical volume in units of cm^3/mol
  ZD05_Tc <- 647.25     # 'Critical temperature in units of Kelvin

  TK <- temperature + 273.15   # 'Temperature must be converted to Kelvin
  Vr <- m / density / ZD05_Vc
  Tr <- TK / ZD05_Tc

  B <- 0.349824207 - 2.91046273 / (Tr * Tr) + 2.00914688 / (Tr * Tr * Tr)
  C <- 0.112819964 + 0.748997714 / (Tr * Tr) - 0.87320704 / (Tr * Tr * Tr)
  D <- 0.0170609505 - 0.0146355822 / (Tr * Tr) + 0.0579768283 / (Tr * Tr * Tr)
  E <- -0.000841246372 + 0.00495186474 / (Tr * Tr) - 0.00916248538 / (Tr * Tr * Tr)
  f <- -0.100358152 / Tr
  g <- -0.00182674744 * Tr

  delta <- 1 + B / Vr + C / (Vr * Vr) + D / Vr^4 + E / Vr^5 + (f / (Vr * Vr) + g / Vr^4) * exp(-0.0105999998 / (Vr * Vr))

  ZD05_R * TK * density * delta / m
}

# 'Calculates the partial derivative of density with respect to pressure, i.e. (d(rho)/dP)_T, in units of g^3/cm^3/bar
# 'density        - The density of water, in g/cm^3
# 'temperature    - The temperature of water, in Celsius
calculate_drhodP <- function(density, temperature) {
  m <- 18.01528          # 'Molar mass of water molecule in units of g/mol
  ZD05_R <- 83.144       # 'Gas Constant in units of cm^3 bar/mol/K
  ZD05_Vc <- 55.9480373  # 'Critical volume in units of cm^3/mol
  ZD05_Tc <- 647.25      # 'Critical temperature in units of Kelvin

  TK <- temperature + 273.15       # 'temperature must be converted to Kelvin
  Tr <- TK / ZD05_Tc
  cc <- ZD05_Vc / m                # 'This term appears frequently in the equation and is defined here for convenience
  Vr <- m / (density * ZD05_Vc)

  B <- 0.349824207 - 2.91046273 / (Tr * Tr) + 2.00914688 / (Tr * Tr * Tr)
  C <- 0.112819964 + 0.748997714 / (Tr * Tr) - 0.87320704 / (Tr * Tr * Tr)
  D <- 0.0170609505 - 0.0146355822 / (Tr * Tr) + 0.0579768283 / (Tr * Tr * Tr)
  E <- -0.000841246372 + 0.00495186474 / (Tr * Tr) - 0.00916248538 / (Tr * Tr * Tr)
  f <- -0.100358152 / Tr
  g <- 0.0105999998 * Tr

  delta <- 1 + B / Vr + C / (Vr^2) + D / Vr^4 + E / Vr^5 + (f / (Vr^2) + g / Vr^4) * exp(-0.0105999998 / Vr^2)

  kappa <- B * cc + 2 * C * (cc^2) * density + 4 * D * cc^4 * density^3 + 5 * E * cc^5 * density^4 +
    (2 * f * (cc^2) * density + 4 * g * cc^4 * density^3 - (f / (Vr^2) + g / Vr^4) * (2 * 0.0105999998 * (cc^2) * density)) * exp(-0.0105999998 / (Vr^2))

  m / (ZD05_R * TK * (delta + density * kappa))
}

# 'Returns the partial derivative of the dielectric constant with respect to density in units of cm^3/g.
# 'density        - The density of water to calculate with, in g/cm^3
# 'temperature    - The temperature to calculate with, in Celsius
calculate_depsdrho <- function(density, temperature) {
  # 'Power Function - Created by Dimitri Sverjensky and Brandon Harrison
  # 'Relevant parameters
  a1 <- -0.00157637700752506
  a2 <- 0.0681028783422197
  a3 <- 0.754875480393944
  b1 <- -8.01665106535394E-05
  b2 <- -0.0687161761831994
  b3 <- 4.74797272182151

  A <- a1 * temperature + a2 * sqrt(temperature) + a3
  B <- b1 * temperature + b2 * sqrt(temperature) + b3

  A * exp(B) * density ^ (A - 1)
}

### Testing functions ###
# These unexported functions are included for testing purposes only.
# In CHNOSZ, the g function and omega(P,T) are calculated via hkf().

# 'Returns the value of omega at the input P and T.
# The value returned is 'in units of cal/mol and NOT multiplied by 10^-5.
#'pressure   - Pressure to calculate at, in bars
#'temperature- Temperature to calculate at, in Celsius
#'density    - Density of water to calculate omega at, in g/cm^3.
#'wref       - The value of omega at standard pressure and temperature, in units of cal/mol.
#'Z          - The charge of the species
calculateOmega <- function(pressure, temperature, density, wref, Z) {
  # 'These equations are given by Shock et al. (1992)
  eta <- 166027       # 'Value in units of Angstroms cal mol^-1
  # 'Defines the electrostatic radius at reference pressure and temperature
  reref <- Z * Z / (wref / eta + Z / 3.082)
  # 'This represents the pressure and temperature dependent solvent function
  g <- calculateG(pressure, temperature, density)
  # 'Defines the electrostatic radius at the input P and T
  re <- reref + abs(Z) * g
  omega <- eta * (Z * Z / re - Z / (3.082 + g))
  # 'If species is hydrogen, the species is neutral, or the pressure is above 6 kb,
  # 'this equation is not necessary because omega is very close to wref.
  if(Z==0) omega[] <- wref
  omega[pressure > 6000] <- wref
}

# 'Returns the value of the g function. If the density is greater than 1 g/cm^3, then zero is returned.
# 'pressure   - The pressure to calculate at, in bars
# 'temperature- The temperature to calculate at, in celsius
# 'density    - The density of water at which to calculate g at, in g/cm^3
calculateG <- function(pressure, temperature, density) {
  T <- temperature
  P <- pressure
  a_g <- -2.037662 + 0.005747 * T - 6.557892E-06 * T * T
  b_g <- 6.107361 - 0.01074377 * T + 1.268348E-05 * T * T
  # 'Calculates the difference function in the case where we need to calculate at Psat conditions
  f <- (((T - 155) / 300)^4.8 + 36.66666 * ((T - 155) / 300)^16) *
      (-1.504956E-10 * (1000 - P)^3 + 5.017997E-14 * (1000 - P)^4)
  f[P > 1000 | T < 155 | T > 355] <- 0
  g <- a_g * (1 - density)^b_g - f
  # Use g = 0 for density >= 1
  g[density >= 1] <- 0
  g
}
