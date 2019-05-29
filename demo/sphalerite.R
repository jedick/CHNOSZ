# CHNOSZ/demo/sphalerite.R
# sphalerite solubility after Akinfiev and Tagirov, 2014, Fig. 13
# 20190526 jmd initial version

# set up chemical system
basis(c("ZnS", "Cl-", "H2S", "H2O", "O2", "H+"))
iZn <- retrieve("Zn", c("O", "H", "Cl", "S"), "aq")
species(iZn)

# a function to make a single plot
plotfun <- function(T = 400, P = 500, m_tot = 0.1, pHmin = 4, logppmmax = 3) {
  # calculate NaCl speciation from simplified model
  NaCl <- NaCl(T = T, P = P, m_tot = m_tot)
  basis("Cl-", log10(NaCl$m_Cl))
  basis("H2S", log10(0.05))

  # use mosaic to account for HS- and H2S speciation
  m <- mosaic(c("H2S", "HS-"), pH = c(pHmin, 10), T = T, P = P, IS = NaCl$IS)
  s <- solubility(m$A.species)

  # convert log activity to log ppm
  sp <- convert(s, "logppm")
  diagram(sp, ylim = c(-5, logppmmax))
  diagram(sp, type = "loga.balance", add = TRUE, lwd = 2, col = "green3")

  # add water neutrality line
  pKw <- - subcrt(c("H2O", "OH-", "H+"), c(-1, 1, 1), T = T, P = P)$out$logK
  abline(v = pKw / 2, lty = 2, lwd = 2, col = "blue1")

  # add legend
  lNaCl <- substitute(italic(m)[NaCl] == x~mol~kg^-1, list(x = m_tot))
  lT <- substitute(x*degree*C, list(x = T))
  if(identical(P, "Psat")) lP <- quote(italic(P)[sat]) else lP <- substitute(x~bar, list(x = P))
  lTP <- substitute(list(x, y), list(x = lT, y = lP))
  ltxt <- as.expression(c(lNaCl, lTP))
  legend("topright", legend = ltxt, bty = "n")
}

plotfun()
title(main = ("Solubility of sphalerite, after Akinfiev and Tagirov, 2014, Fig. 13"), font.main = 1)

### the following code for making multiple plots is not used in the demo ###

# a function to make a page of plots
pagefun <- function() {
  # set the values of temperature, pressure, and total NaCl
  T <- c(400, 400, 250, 250, 100, 100)
  # use a list to be able to mix numeric and character values for P
  P <- list(500, 500, "Psat", "Psat", "Psat", "Psat")
  m_tot <- c(0.1, 1, 0.1, 1, 0.1, 1)
  # the plots have differing limits
  pHmin <- c(4, 4, 2, 2, 2, 2)
  logppmmax <- c(3, 3, 2, 2, 0, 0)
  # make the plots
  par(mfrow = c(3, 2))
  for(i in 1:6) plotfun(T = T[i], P = P[[i]], m_tot = m_tot[i], pHmin = pHmin[i], logppmmax = logppmmax[i])
}

# a function to make a png file with all the plots
pngfun <- function() {
  png("sphalerite.png", width = 1000, height = 1200, pointsize = 24)
  pagefun()
  # add an overall title
  par(xpd = NA)
  text(1, 14, "Solubility of sphalerite, after Akinfiev and Tagirov, 2014, Fig. 13", cex = 1.5)
  par(xpd = FALSE)
  dev.off()
}

# we don't run these functions in the demo
#pagefun()
#pngfun()
