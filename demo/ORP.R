# CHNOSZ/demo/ORP.R 
# First version 20100715 jmd

library(CHNOSZ)

# Calculate the temperature dependence of 
# potentials vs. standard hydrogen electrode (SHE) of various electrodes (Ag/AgCl)
# and ORP standards (ZoBell, Light's, (tri)iodide) 

# Use the extended Debye-Huckel equation with Alberty's parameters
oldnon <- nonideal("Alberty")

# Bard et al.'s fit to the potential
# (Bard, Parson, Jordan, Standard Potentials In Aqueous Solution, 1985)
AgAgCl.Bard <- function(T, high.T = TRUE) {
  # we use the corrected high-T formula from wikipedia
  if(high.T) return(0.23737 - 5.3783e-4 * T - 2.3728e-6 * T^2 - 2.2671e-9 * (T+273))
  else return(0.23695 - 4.8564e-4 * T - 3.4205e-6 * T^2 - 5.869e-9 * (T+273))
}

# Function to calculate the potential of Ag/AgCl vs. SHE
# Ag(s) + Cl- = AgCl(s) + e-
# logK = -pe - logaCl
# pe = -logK - logmCl - loggamCl
# ORP = RT/F * (logK - logmCl - loggamCl)
AgAgCl <- function(T, mKCl = 4) {
  # mKCl is the molality of KCl in the electrolyte
  # We take it as a first approximation to be equal to
  # the molality of Cl- (and to the ionic strength)
  logmCl <- log10(mKCl)
  # Get the logK for the reaction
  logK <- subcrt(c("Ag", "Cl-", "AgCl", "e-"), c(-1, -1, 1, 1), c("cr", "aq", "cr", "aq"), T = T)$out$logK
  # Get the activity coefficient for Cl-
  loggamCl <- log(10) * subcrt("Cl-", T = T, IS = mKCl)$out[[1]]$loggam
  # Get the pe for the solution
  pe <- -logK - logmCl - loggamCl
  # Convert that to Eh
  Eh <- convert(pe, "Eh", T = convert(T, "K"))
  return(Eh)
}

ZoBell <- function(T) {
  # Doesn't work very well because we ignore the
  # ferricyanide and ferrocyanide complexes
  # Fe+2 = Fe+3 + e-
  # logK = logaFe3 - logaFe2 - pe
  # Get the logK for the reaction
  logK <- subcrt(c("Fe+2", "Fe+3", "e-"), c(-1, 1, 1), T = T)$out$logK
  # We use the recipe from standard methods (table 2580:II)
  # 1.4080 g K4Fe(CN)6.3H2O -> 0.0033333 mol Fe+2
  # 1.0975 g K3Fe(CN)6      -> 0.0033333 mol Fe+3
  # 7.4555 g KCl            -> 0.1 mol Cl-
  logmFe2 <- logmFe3 <- log10(0.0033333)
  # Get the loggam for the iron species
  loggamFe2 <- log(10) * subcrt("Fe+2", T = T, IS = 1)$out[[1]]$loggam
  loggamFe3 <- log(10) * subcrt("Fe+3", T = T, IS = 1)$out[[1]]$loggam
  # Get the pe for the solution
  pe <- -logK + logmFe3 + loggamFe3 - logmFe2 - loggamFe2
  # Convert to Eh
  Eh <- convert(pe, "Eh", T = convert(T, "K"))
  return(Eh)
}

ZoBell.table <- function(T = NULL, which = NULL) {
  # Oxidation-reduction potential of ZoBell's solution 
  # from Standard Methods for Water and Wastewater or YSI 
  # (interpolated and/or extrapolated as necessary)
  # Standard methods (1997) table 2580:I
  Eh.T.SMW <- 1:30
  Eh.SMW <- c(0.481, 0.479, 0.476, 0.474, 0.472, 0.47, 0.468, 0.465, 0.463, 0.461, 
  0.459, 0.457, 0.454, 0.452, 0.45, 0.448, 0.446, 0.443, 0.441, 0.439, 0.437, 
  0.435, 0.432, 0.43, 0.428, 0.426, 0.424, 0.421, 0.419, 0.417)
  # From YSI (2005):
  # Measuring ORP on YSI 6-Series Sondes: Tips, Cautions and Limitations
  # NOTE: these values are vs. Ag/AgCl (4 M KCl)
  Eh.T.YSI <- seq(-5, 50, by = 5)
  Eh.YSI <- c(267.0, 260.5, 254.0, 247.5, 241.0, 234.5, 228.0, 221.5, 215.0, 208.5, 202.0, 195.5)/1000
  # Spline function for each of the tables
  SMW <- splinefun(Eh.T.SMW, Eh.SMW)
  YSI <- splinefun(Eh.T.YSI, Eh.YSI)
  # Just one of the tables
  Eh.fun <- get(which)
  Eh.T <- get(paste("Eh.T", which, sep = "."))
  if(is.null(T)) T <- Eh.T
  return(data.frame(T = T, Eh = Eh.fun(T)))
}

Light <- function(T) {
  # This is going to look something like
  # Fe+2 = Fe+3 + e-
  # logK = logaFe3 - logaFe2 - pe
  # Get the logK for the reaction
  logK <- subcrt(c("Fe+2", "Fe+3", "e-"), c(-1, 1, 1), T = T)$out$logK
  # We use the recipe from standard methods (table 2580:II)
  # 39.21 g Fe(NH4)2(SO4)2(H2O)6 -> 0.1 mol Fe+2
  # 48.22 g Fe(NH4)(SO4)2(H2O)12 -> 0.1 mol Fe+3
  logmFe2 <- logmFe3 <- log10(0.1)
  # Get the loggam for the iron species
  loggamFe2 <- log(10) * subcrt("Fe+2", T = T, IS = 0.2)$out[[1]]$loggam
  loggamFe3 <- log(10) * subcrt("Fe+3", T = T, IS = 0.2)$out[[1]]$loggam
  # Get the pe for the solution
  pe <- -logK + logmFe3 + loggamFe3 - logmFe2 - loggamFe2
  # Convert to Eh
  Eh <- convert(pe, "Eh", T = convert(T, "K"))
  return(Eh)
}

Iodide.table <- function(T=NULL) {
  # Oxidation-reduction potential of Thermo's iodide solution
  # from thermo instruction sheet 255218-001 (articlesFile_18739)
  T.Iodide <- seq(0, 50, 5)
  Eh.Iodide <- c(438, 435, 431, 428, 424, 420, 415, 411, 406, 401, 396)/1000
  Iodide <- splinefun(T.Iodide, Eh.Iodide)
  if(is.null(T)) T <- T.Iodide
  return(data.frame(T = T, Eh = Iodide(T)))
}

Iodide <- function(T) {
  # This is going to look something like
  # 3I- = I3- + 2e-
  # logK = -2pe + logaI3 - 3logaI
  # Get the logK for the reaction
  logK <- subcrt(c("I-", "I3-", "e-"), c(-3, 1, 2), T = T)$out$logK
  # Could the activities be 0.1 M ... or something else?
  logmI <- log10(2)
  logmI3 <- log10(0.01)
  # Get the loggam for the iodine species
  loggamI <- log(10) * subcrt("I-", T = T, IS = 0.2)$out[[1]]$loggam
  loggamI3 <- log(10) * subcrt("I3-", T = T, IS = 0.2)$out[[1]]$loggam
  # Get the pe for the solution
  pe <- ( -logK + logmI3 + loggamI3 - 3 * (logmI - loggamI) ) / 2
  # Convert to Eh
  Eh <- convert(pe, "Eh", T = convert(T, "K"))
  return(Eh)
}

figure <- function() {
  # Make some figures
  # Temperatures in degrees C
  T <- seq(0, 100, 5)
  # Temperature-Eh diagram for various electrodes
  thermo.plot.new(ylim = c(0, 0.8), xlim = c(0, 100), 
    ylab = axis.label("Eh"), xlab = axis.label("T"))
  # Ag/AgCl electrode (Bard et al. fit)
  points(T, AgAgCl.Bard(T), pch = 0)
  # Ag/AgCl electrode (equilibrium calculations) 
  lines(T, AgAgCl(T))
  # ZoBell's solution (SMW table 2580)
  SMW <- ZoBell.table(which = "SMW")
  points(SMW$T, SMW$Eh, pch = 1)
  # ZoBell's solution (YSI tech report table)
  YSI <- ZoBell.table(which = "YSI")
  # Make these values referenced to SHE instead of Ag/AgCl
  Eh.YSI <- YSI$Eh + AgAgCl(YSI$T)
  points(YSI$T, Eh.YSI, pch = 2)
  # Light's solution (equilibrium values)
  lines(T, Light(T))
  # Light's solution (at 25 degrees only)
  points(25, 0.475 + 0.200, pch = 3)
  # Thermo's I-/I3- solution
  Thermo <- Iodide.table()
  points(Thermo$T, Thermo$Eh, pch = 4)
  # Calculated I-/I3- values
  lines(T, Iodide(T))
  # Add some labels
  text(c(30, 30, 30, 50), c(0.72, 0.5, 0.35, 0.25), 
    c("Light", "ZoBell", "(Tri)Iodide", "Ag/AgCl"))
  title(main = "Potentials vs standard hydrogen electrode (SHE)")
}

# Finally, make the plot
figure()

# Reset the nonideality method to default
nonideal(oldnon)
