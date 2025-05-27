# CHNOSZ/demo/MgATP.R
# Speciation of ATP with H+ and Mg+2
# 20250423 Moved from anintro.Rmd

library(CHNOSZ)

# Speciation plots, similar to Figures 1.2-1.5 in Alberty (2003).
# These figures show the distribution of differently charged species of
# adenosine triphosphate (ATP) as a function of pH,
# and the average number of H+ and Mg+2 bound to ATP in solution as a function of pH or pMg.

par(mfrow = c(2, 2), mar = c(3.1, 3.6, 2.1, 1.6), mgp = c(1.8, 0.5, 0))

# For the following calculations, we change the nonideality method to Alberty.
# This is a simpler formulation than the default,
# with parameters that are suitable for biochemical species at relatively low temperatures.

oldnon <- nonideal("Alberty")

# Use the following commands to set the basis species,
# add the variously protonated ATP species,
# calculate the affinities of the formation reactions,
# equilibrate the system, and make a degree of formation (alpha) or mole fraction diagram.
# This is similar to Figure 1.3 of Alberty (2003), but is calculated for I = 0 M and T = 100 degrees C.

basis("MgCHNOPS+")
species(c("ATP-4", "HATP-3", "H2ATP-2", "H3ATP-", "H4ATP"))
T <- 100
a <- affinity(pH = c(3, 9), T = T)
e <- equilibrate(a)
d <- diagram(e, alpha = TRUE, tplot = FALSE)
title(main = describe.property("T", T))

# Note that we have saved the numeric results of diagram(), i.e. the degrees of formation of the species (alpha).
# With that, we can calculate and plot the average number of protons bound per ATP molecule.
# To do so, we use R's rbind() and do.call() to turn alpha into a matrix,
# then multiply by the number of protons bound to each species,
# and sum the columns to get the total (i.e. average proton number, N_H+).

alphas <- do.call(rbind, d$plotvals)
nH <- alphas * 0:4
Hlab <- substitute(italic(N)[H^`+`])
plot(a$vals[[1]], colSums(nH), type = "l", xlab = "pH", ylab=Hlab, lty=2, col=2)

# Adding the IS argument to affinity(), we can now plot N_H+ at the given ionic strength.
# Here we set plot.it = FALSE diagram() because we use the computed alpha to make our own plot.
# This is similar to Figure 1.3 of Alberty (2003), but at higher temperature.

a <- affinity(pH = c(3, 9), IS = 0.25, T = T)
e <- equilibrate(a)
d <- diagram(e, alpha = TRUE, plot.it = FALSE)
alphas <- do.call(rbind, d$plotvals)
nH <- alphas * 0:4
lines(a$vals[[1]], colSums(nH))
legend("topright", legend = c("I = 0 M", "I = 0.25 M"), lty = 2:1, col = 2:1, cex = 0.8)
ATP.H <- substitute("ATP and H"^`+`)
title(main = ATP.H)

# Next, we add the Mg+2-complexed ATP species:

species(c("MgATP-2", "MgHATP-", "MgH2ATP", "Mg2ATP"), add = TRUE)

# Here is a function to calculate and plot N_H+ for a given pMg:

Hplot <- function(pMg, IS = 0.25) {
  basis("Mg+2", -pMg)
  a <- affinity(pH = c(3, 9), IS = IS, T = T)
  e <- equilibrate(a)
  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
  alphas <- do.call(rbind, d$plotvals)
  NH <- alphas * c(0:4, 0, 1, 2, 0)
  lines(a$vals[[1]], colSums(NH), lty = 7 - pMg, col = 7 - pMg)
}

# With that function in hand, we plot the lines corresponding to pMg = 2 to 6.
# This is similar to Figure 1.4 of Alberty (2003):

plot(c(3, 9), c(0, 2), type = "n", xlab = "pH", ylab = Hlab)
lapply(2:6, Hplot)
legend("topright", legend = paste("pMg = ", 2:6), lty = 5:1, col = 5:1, cex = 0.8)
ATP.H.Mg <- substitute("ATP and H"^`+`~"and Mg"^`+2`)
title(main = ATP.H.Mg)

# The next function calculates and plots the average number of Mg+2 bound to ATP (N_Mg+2) for a given pH.
# Here we multiply alpha by the number of Mg+2 in each species, and negate loga_Mg+2 (the variable used in affinity()) to get pMg.

Mgplot <- function(pH, IS = 0.25) {
  basis("pH", pH)
  a <- affinity(`Mg+2` = c(-2, -7), IS = IS, T = T)
  e <- equilibrate(a)
  d <- diagram(e, alpha = TRUE, plot.it = FALSE)
  alphas <- do.call(rbind, d$plotvals)
  NMg <- alphas * species()$`Mg+`
  lines(-a$vals[[1]], colSums(NMg), lty = 10 - pH, col = 10 - pH)
}

# Using that function, we plot the lines corresponding to pH = 3 to 9.
# This is similar to Figure 1.5 of Alberty (2003):

Mglab <- substitute(italic(N)[Mg^`+2`])
plot(c(2, 7), c(0, 1.2), type = "n", xlab = "pMg", ylab = Mglab)
lapply(3:9, Mgplot)
legend("topright", legend = paste("pH = ", 3:9), lty = 7:1, col = 7:1, cex = 0.8)
title(main = ATP.H.Mg)

# Now that we're finished, we can reset the nonideality method to the default.

nonideal(oldnon)
