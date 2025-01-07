# CHNOSZ/demo/sphalerite.R
# Sphalerite solubility after Akinfiev and Tagirov, 2014, Fig. 13
# 20190526 jmd initial version
library(CHNOSZ)
opar <- par(no.readonly = TRUE)

# Set up chemical system
basis(c("Zn+2", "Cl-", "H2S", "H2O", "O2", "H+"))
basis("H2S", log10(0.05))
species("ZnS")
iaq <- retrieve("Zn", c("O", "H", "Cl", "S"), "aq")

# A function to make a single plot
plotfun <- function(T = 400, P = 500, m_tot = 0.1, pHmin = 4, logppmmax = 3) {
  # Get pH values
  res <- 100
  pH <- seq(pHmin, 10, length.out = res)
  # Calculate speciation in NaCl-H2O system at given pH
  NaCl <- NaCl(m_tot = m_tot, T = T, P = P, pH = pH)

  # Calculate solubility with mosaic (triggered by bases argument) to account for HS- and H2S speciation
  s <- solubility(iaq, bases = c("H2S", "HS-"), pH = pH, "Cl-" = log10(NaCl$m_Cl), T = T, P = P, IS = NaCl$IS)

  # Convert log activity to log ppm
  sp <- convert(s, "logppm")
  diagram(sp, type = "loga.equil", ylim = c(-5, logppmmax))
  diagram(sp, add = TRUE, lwd = 2, col = "green3")

  # Add water neutrality line
  pKw <- - subcrt(c("H2O", "OH-", "H+"), c(-1, 1, 1), T = T, P = P)$out$logK
  abline(v = pKw / 2, lty = 2, lwd = 2, col = "blue1")

  # Add legend
  l <- lex(lNaCl(m_tot), lTP(T, P))
  legend("topright", legend = l, bty = "n")
}

plotfun()
title(main = ("Solubility of sphalerite, after Akinfiev and Tagirov, 2014, Fig. 13"), font.main = 1)

par(opar)

### The following code for making multiple plots is not used in the demo ###

# A function to make a page of plots
pagefun <- function() {
  # Set the values of temperature, pressure, and total NaCl
  T <- c(400, 400, 250, 250, 100, 100)
  # Use a list to be able to mix numeric and character values for P
  P <- list(500, 500, "Psat", "Psat", "Psat", "Psat")
  m_tot <- c(0.1, 1, 0.1, 1, 0.1, 1)
  # The plots have differing limits
  pHmin <- c(4, 4, 2, 2, 2, 2)
  logppmmax <- c(3, 3, 2, 2, 0, 0)
  # Make the plots
  par(mfrow = c(3, 2))
  for(i in 1:6) plotfun(T = T[i], P = P[[i]], m_tot = m_tot[i], pHmin = pHmin[i], logppmmax = logppmmax[i])
}

# A function to make a png file with all the plots
pngfun <- function() {
  png("sphalerite.png", width = 1000, height = 1200, pointsize = 24)
  pagefun()
  # Add an overall title
  par(xpd = NA)
  text(1, 14, "Solubility of sphalerite, after Akinfiev and Tagirov, 2014, Fig. 13", cex = 1.5)
  par(xpd = FALSE)
  dev.off()
}

# We don't run these functions in the demo
#pagefun()
#pngfun()
