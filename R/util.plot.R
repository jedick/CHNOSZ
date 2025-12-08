# CHNOSZ/util.plot.R
# Utility functions for plots

label.plot <- function(label, xfrac = 0.07, yfrac = 0.93, paren = FALSE, italic = FALSE, ...) {
  # Make a text label e.g., "(a)" in the corner of a plot
  # xfrac, yfrac: fraction of axis where to put label (default top right)
  # paren: put a parenthesis around the text, and italicize it?
  if(italic) label <- bquote(italic(.(label)))
  if(paren) label <- bquote(group("(", .(label), ")"))
  pu <- par('usr')
  x <- pu[1]+xfrac*(pu[2]-pu[1])
  y <- pu[3]+yfrac*(pu[4]-pu[3])
  # Conversion for logarithmic axes
  if(par("xlog")) x <- 10^x
  if(par("ylog")) y <- 10^y
  text(x, y, labels = label, ...)
}

usrfig <- function() {
  # Function to get the figure limits in user coordinates
  # Get plot limits in user coordinates (usr) and as fraction [0,1] of figure region (plt)
  xusr <- par('usr')[1:2]; yusr <- par('usr')[3:4]
  xplt <- par('plt')[1:2]; yplt <- par('plt')[3:4]
  # Linear model to calculate figure limits in user coordinates
  xlm <- lm(xusr ~ xplt); ylm <- lm(yusr ~ yplt)
  xfig <- predict.lm(xlm, data.frame(xplt = c(0, 1)))
  yfig <- predict.lm(ylm, data.frame(yplt = c(0, 1)))
  return(list(x = xfig, y = yfig))
}

label.figure <- function(label, xfrac = 0.05, yfrac = 0.95, paren = FALSE, italic = FALSE, ...) {
  # Function to add labels outside of the plot region  20151020
  f <- usrfig()
  # Similar to label.plot(), except we have to set xpd here
  opar <- par(xpd = NA)
  if(italic) label <- bquote(italic(.(label)))
  if(paren) label <- bquote(group("(", .(label), ")"))
  # Calculate location for label
  x <- f$x[1] + xfrac * (f$x[2] - f$x[1])
  y <- f$y[1] + yfrac * (f$y[2] - f$y[1])
  # Conversion for logarithmic axes
  if(par("xlog")) x <- 10^x
  if(par("ylog")) y <- 10^y
  text(x, y, labels = label, ...)
  par(opar)
}

water.lines <- function(eout, which = c('oxidation','reduction'),
  lty = 2, lwd = 1, col = par('fg'), plot.it = TRUE) {

  # Draw water stability limits for Eh-pH, logfO2-pH, logfO2-T or Eh-T diagrams
  # (i.e. redox variable is on the y axis)

  # Get axes, T, P, and xpoints from output of affinity() or equilibrate()
  if(missing(eout)) stop("'eout' (the output of affinity(), equilibrate(), or diagram()) is missing")
  # Number of variables used in affinity()
  nvar1 <- length(eout$vars)
  # If these were on a transect, the actual number of variables is less
  dim <- dim(eout$loga.equil[[1]]) # for output from equilibrate()
  if(is.null(dim)) dim <- dim(eout$values[[1]]) # for output from affinity()
  nvar2 <- length(dim)
  # We only work on diagrams with 1 or 2 variables
  if(!nvar1 %in% c(1, 2) | !nvar2 %in% c(1, 2)) return(NA)

  # If needed, swap axes so redox variable is on y-axis
  # Also do this for 1-D diagrams 20200710
  if(is.na(eout$vars[2])) eout$vars[2] <- "nothing"
  swapped <- FALSE
  if(eout$vars[2] %in% c("T", "P", "nothing")) {
    eout$vars <- rev(eout$vars)
    eout$vals <- rev(eout$vals)
    swapped <- TRUE
  }
  xaxis <- eout$vars[1]
  yaxis <- eout$vars[2]
  xpoints <- eout$vals[[1]]
  # Make xaxis "nothing" if it is not pH, T, or P 20201110
  # (so that horizontal water lines can be drawn for any non-redox variable on the x-axis)
  if(!identical(xaxis, "pH") & !identical(xaxis, "T") & !identical(xaxis, "P")) xaxis <- "nothing"

  # T and P are constants unless they are plotted on one of the axes
  T <- eout$T
  if(eout$vars[1] == "T") T <- envert(xpoints, "K")
  P <- eout$P
  if(eout$vars[1] == "P") P <- envert(xpoints, "bar")
  # logaH2O is 0 unless given in eout$basis
  iH2O <- match("H2O", rownames(eout$basis))
  if(is.na(iH2O)) logaH2O <- 0 else logaH2O <- as.numeric(eout$basis$logact[iH2O])
  # pH is 7 unless given in eout$basis or plotted on one of the axes
  iHplus <- match("H+", rownames(eout$basis))
  if(eout$vars[1] == "pH") pH <- xpoints
  else if(!is.na(iHplus)) {
    minuspH <- eout$basis$logact[iHplus]
    # Special treatment for non-numeric value (happens when a buffer is used, even for another basis species)
    if(can.be.numeric(minuspH)) pH <- -as.numeric(minuspH) else pH <- NA
  }
  else pH <- 7

  # O2state is gas unless given in eout$basis
  iO2 <- match("O2", rownames(eout$basis))
  if(is.na(iO2)) O2state <- "gas" else O2state <- eout$basis$state[iO2]
  # H2state is gas unles given in eout$basis
  iH2 <- match("H2", rownames(eout$basis))
  if(is.na(iH2)) H2state <- "gas" else H2state <- eout$basis$state[iH2]

  # Where the calculated values will go
  y.oxidation <- y.reduction <- NULL
  if(xaxis %in% c("pH", "T", "P", "nothing") & yaxis %in% c("Eh", "pe", "O2", "H2")) {
    # Eh/pe/logfO2/logaO2/logfH2/logaH2 vs pH/T/P
    if('reduction' %in% which) {
      logfH2 <- logaH2O # usually 0
      if(yaxis == "H2") {
        logK <- suppressMessages(subcrt(c("H2", "H2"), c(-1, 1), c("gas", H2state), T = T, P = P, convert = FALSE))$out$logK
        # This is logfH2 if H2state == "gas", or logaH2 if H2state == "aq"
        logfH2 <- logfH2 + logK
        y.reduction <- rep(logfH2, length.out = length(xpoints))
      } else {
        logK <- suppressMessages(subcrt(c("H2O", "O2", "H2"), c(-1, 0.5, 1), c("liq", O2state, "gas"), T = T, P = P, convert = FALSE))$out$logK 
        # This is logfO2 if O2state == "gas", or logaO2 if O2state == "aq"
        logfO2 <- 2 * (logK - logfH2 + logaH2O)
        if(yaxis == "O2") y.reduction <- rep(logfO2, length.out = length(xpoints))
        else if(yaxis == "Eh") y.reduction <- convert(logfO2, 'E0', T = T, P = P, pH = pH, logaH2O = logaH2O)
        else if(yaxis == "pe") y.reduction <- convert(convert(logfO2, 'E0', T = T, P = P, pH = pH, logaH2O = logaH2O), "pe", T = T)
      }
    }
    if('oxidation' %in% which) {
      logfO2 <- logaH2O # usually 0
      if(yaxis == "H2") {
        logK <- suppressMessages(subcrt(c("H2O", "O2", "H2"), c(-1, 0.5, 1), c("liq", "gas", H2state), T = T, P = P, convert = FALSE))$out$logK 
        # This is logfH2 if H2state == "gas", or logaH2 if H2state == "aq"
        logfH2 <- logK - 0.5*logfO2 + logaH2O
        y.oxidation <- rep(logfH2, length.out = length(xpoints))
      } else {
        logK <- suppressMessages(subcrt(c("O2", "O2"), c(-1, 1), c("gas", O2state), T = T, P = P, convert = FALSE))$out$logK 
        # This is logfO2 if O2state == "gas", or logaO2 if O2state == "aq"
        logfO2 <- logfO2 + logK
        if(yaxis == "O2") y.oxidation <- rep(logfO2, length.out = length(xpoints))
        else if(yaxis == "Eh") y.oxidation <- convert(logfO2, 'E0', T = T, P = P, pH = pH, logaH2O = logaH2O)
        else if(yaxis == "pe") y.oxidation <- convert(convert(logfO2, 'E0', T = T, P = P, pH = pH, logaH2O = logaH2O), "pe", T = T)
      }
    }
  } else return(NA)

  # Now plot the lines
  if(plot.it) {
    if(swapped) {
      if(nvar1 == 1 | nvar2 == 2) {
        # Add vertical lines on 1-D diagram 20200710
        abline(v = y.oxidation[1], lty = lty, lwd = lwd, col = col)
        abline(v = y.reduction[1], lty = lty, lwd = lwd, col = col)
      } else {
        # xpoints above is really the ypoints
        lines(y.oxidation, xpoints, lty = lty, lwd = lwd, col = col)
        lines(y.reduction, xpoints, lty = lty, lwd = lwd, col = col)
      }
    } else {
      lines(xpoints, y.oxidation, lty = lty, lwd = lwd, col = col)
      lines(xpoints, y.reduction, lty = lty, lwd = lwd, col = col)
    }
  }
  # Return the values
  return(invisible(list(xpoints = xpoints, y.oxidation = y.oxidation, y.reduction = y.reduction, swapped = swapped)))
}

mtitle <- function(main, line = 0, spacing = 1, ...) {
  # Make a possibly multi-line plot title 
  # Useful for including expressions on multiple lines 
  # 'line' is the margin line of the last (bottom) line of the title
  len <- length(main)
  for(i in 1:len) mtext(main[i], line = line + (len - i)*spacing, ...)
}

# Get colors for range of ZC values 20170206
ZC.col <- function(z) {
  # Scale values to [1, 1000]
  z <- z * 999/diff(range(z))
  z <- round(z - min(z)) + 1
  # Diverging (blue - light grey - red) palette
  # dcol <- colorspace::diverge_hcl(1000, c = 100, l = c(50, 90), power = 1)
  # Use precomputed values
  file <- system.file("extdata/misc/bluered.txt", package = "CHNOSZ")
  dcol <- read.table(file, as.is = TRUE)[[1]]
  # Reverse the palette so red is at lower ZC (more reduced)
  rev(dcol)[z]
}

# To get hyphen instead of minus sign 20220630
# https://stackoverflow.com/questions/10438398/any-way-to-disable-the-minus-hack-in-pdf-poscript-output
hyphen.in.pdf <- function(x) {
  # We only want to make the substitution in a pdf device (won't work in png, e.g. for knitr vignettes)
  if(identical(names(dev.cur()), "pdf")) gsub("-", "\uad", x, fixed = TRUE) else x
}
