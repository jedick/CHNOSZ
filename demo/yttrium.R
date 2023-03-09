# CHNOSZ/demo/yttrium.R
# Make diagrams similar to Figure 11 of Guan et al. (2020)
# https://doi.org/10.1016/j.gca.2020.04.015
# 20221122 jmd version 1

library(CHNOSZ)

# Function to add Y(III)-Cl species using logB values from Table 9 of Guan et al. (2020)
add.Y.species <- function(P, plot.it = FALSE) {
  # Temperature
  T <- c(25, seq(50, 500, 50))
  if(P == 800) {
    logB <- list(
      c(1.04, 0.48, 0.13, 0.43, 1.12, 2.06, 3.21, 4.58, 6.26, 8.39, 11.14),
      c(-9.14, -7.34, -4.14, -1.37, 1.12, 3.46, 5.78, 8.24, 11.05, 14.48, 18.77),
      c(-14, -11.48, -7.06, -3.25, 0.14, 3.32, 6.45, 9.74, 13.47, 18.01, 23.65),
      c(-15.94, -13.2, -8.39, -4.27, -0.61, 2.79, 6.15, 9.67, 13.63, 18.45, 24.41)
    )
  } else if(P == 1000) {
    logB <- list(
      c(1.13, 0.54, 0.16, 0.43, 1.09, 2, 3.1, 4.39, 5.9, 7.69, 9.85),
      c(-9.33, -7.51, -4.3, -1.55, 0.9, 3.18, 5.4, 7.68, 10.15, 12.94, 16.16),
      c(-14.24, -11.71, -7.27, -3.49, -0.14, 2.95, 5.95, 9.01, 12.3, 16, 20.25),
      c(-16.19, -13.43, -8.62, -4.52, -0.91, 2.41, 5.62, 8.89, 12.4, 16.33, 20.84)
    )
  } else stop("logB values for P =", P, "are not available here")
  # Define species and coefficients in formation reactions
  species <- list(
    c("Y+3", "Cl-", "YCl+2"),
    c("Y+3", "Cl-", "YCl2+"),
    c("Y+3", "Cl-", "YCl3"),
    c("Y+3", "Cl-", "YCl4-")
  )
  coeffs <- list(
    c(-1, -1, 1),
    c(-1, -2, 1),
    c(-1, -3, 1),
    c(-1, -4, 1)
  )
  # Fit the formation constants to thermodynamic parameters and add them to OBIGT
  for(i in 1:4) logB.to.OBIGT(logB[[i]], species[[i]], coeffs[[i]], T = T, P = P, tolerance = 0.6, npar = 5)
  # Plot the given and fitted values
  if(plot.it) {
    par(mfrow = c(2, 2))
    for(i in 1:4) {
      sres <- subcrt(species[[i]], coeffs[[i]], T = T, P = P)
      plot(T, sres$out$logK, type = "l", xlab = axis.label("T"), ylab = axis.label("logK"))
      points(T, logB[[i]], pch = 19)
      title(describe.reaction(sres$reaction))
      if(i == 1) legend("topleft", c("Guan et al. (2020)", "logB.to.OBIGT()"), pch = c(19, NA), lty = c(NA, 1))
      legend("bottomright", paste(P, "bar"), bty = "n")
    }
  }
}

# Function to plot distribution of Y(III) chloride species at T and P
Y_Cl <- function() {

  # Define total molality of NaCl
  # Start at 0.1 because we can't use 0 in the logarithmic value needed for affinity()
  mNaCl <- seq(0.1, 4.9, 0.2)

  # Define T, P, and pH values
  Ts <- c(200, 350, 500)
  Ps <- c(800, 800, 1000)
  pHs <- c(3, 0.3)

  # Setup plot
  par(mfrow = c(3, 2))
  par(cex = 0.9)

  # Loop over T and P
  for(i in 1:3) {
    T <- Ts[i]
    P <- Ps[i]

    # Add new species
    add.Y.species(P)
    # Setup chemical system
    basis(c("Y+3", "Cl-", "e-"))
    species(c("Y+3", "YCl+2", "YCl2+", "YCl3", "YCl4-"))

    # Loop over pH
    for(j in 1:2) {
      pH <- pHs[j]

      # Calculate molality of Cl- and ionic strength
      NaCl_props <- suppressMessages(lapply(mNaCl, NaCl, T = T, P = P, pH = pH))
      # Turn the list into a data frame
      NaCl_props <- do.call(rbind, lapply(NaCl_props, as.data.frame))
      # Calculate affinity of formation reactions
      a <- affinity("Cl-" = log10(NaCl_props$m_Cl), IS = NaCl_props$IS, T = T, P = P)
      # Calculate species distribution for total Y(III) equal to 0.01 m
      m_Y <- 0.01
      e <- equilibrate(a, loga.balance = log10(m_Y))
      # Make x-axis show total m(NaCl) instead of logm(Cl-) 20221208
      e$vals[[1]] <- mNaCl
      # Set colors similar to Guan et al. (2020)
      col <- 2:6
      # Only label lines above 1/20 = 0.05 mol fraction
      mol <- 10^do.call(cbind, lapply(e$loga.equil, as.data.frame))
      molfrac <- mol / rowSums(mol)
      ilab <- apply(molfrac > 0.05, 2, any)
      names <- e$species$name
      names[!ilab] <- ""
      d <- diagram(e, alpha = TRUE, xlim = c(0, 5), ylim = c(0, 1),
        xlab = expr.species("NaCl", molality = TRUE), ylab = "Fraction of Y-Cl species",
        names = names, col = col, lty = 1, lwd = 2, mar = c(3, 3.5, 2.5, 3.5))
      # Calculate and plot coordination number
      CN <- 1 * d$plotvals[[2]] + 2 * d$plotvals[[3]] + 3 * d$plotvals[[4]] + 4 * d$plotvals[[5]]
      # Rescale to y-axis limits [0, 1]
      CN_scaled <- CN / 4
      lines(mNaCl, CN_scaled, lty = 3, lwd = 3, col = "darkorange")
      # Add ticks for CN
      axis(4, seq(0, 1, 0.25), labels = 0:4, lwd.ticks = 2, tcl = -0.5, mgp = c(1.7, 0.8, 0))
      mtext("Cl coordination number", 4, 1.7, las = 0, cex = par("cex"))
      # Make title
      lab <- lTP(T, P)
      lab[[4]] <- bquote(pH == .(pH))
      title(lab, cex.main = 1)

    }

    if(i==1) mtext("After Guan et al. (2020)", line = 0.8, adj = -2.7, cex = 0.9, font = 2)
  }
}

# Use non-default ion size parameters 20230309
Bdot_acirc <- thermo()$Bdot_acirc
# Cl- and Y+3 override the defaults, and YCl+2 is a new species
Bdot_acirc <- c("Cl-" = 4, "Y+3" = 5, "YCl+2" = 4, Bdot_acirc)
thermo("Bdot_acirc" = Bdot_acirc)

# Run the functions to make plots for the demo
opar <- par(no.readonly = TRUE)
add.Y.species(800, plot.it = TRUE)
add.Y.species(1000, plot.it = TRUE)
Y_Cl()

# Restore plot settings and CHNOSZ settings
par(opar)
reset()
