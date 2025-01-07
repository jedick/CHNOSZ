# CHNOSZ/diagram.R

# Plot equilibrium chemical activity and predominance diagrams 
# 20061023 jmd v1
# 20120927 work with output from either equil() or affinity(), 
#   gather plotvals independently of plot parameters (including nd),
#   single return statement

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("equilibrate.R")
#source("util.plot.R")
#source("util.character.R")
#source("util.misc.R")

diagram <- function(
  # Species affinities or activities
  eout, 
  # Type of plot
  type = "auto", alpha = FALSE, normalize = FALSE, as.residue = FALSE, balance = NULL,
  groups = as.list(1:length(eout$values)), xrange = NULL,
  # Figure size and sides for axis tick marks
  mar = NULL, yline = par("mgp")[1]+0.3, side = 1:4,
  # Axis limits and labels
  ylog = TRUE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, 
  # Sizes
  cex = par("cex"), cex.names = 1, cex.axis = par("cex"),
  # Line styles
  lty = NULL, lty.cr = NULL, lty.aq = NULL, lwd = par("lwd"), dotted = NULL,
  spline.method = NULL, contour.method = "edge", levels = NULL,
  # Colors
  col = par("col"), col.names = par("col"), fill = NULL,
  fill.NA = "gray80", limit.water = NULL,
  # Field and line labels
  names = NULL, format.names = TRUE, bold = FALSE, italic = FALSE,
  font = par("font"), family = par("family"), adj = 0.5, dx = 0, dy = 0, srt = 0,
  min.area = 0,
  # Title and legend
  main = NULL, legend.x = NA,
  # Plotting controls
  add = FALSE, plot.it = TRUE, tplot = TRUE, ...
) {

  ### Argument handling ###

  ## Check that eout is a valid object
  efun <- eout$fun
  if(length(efun) == 0) efun <- ""
  if(!(efun %in% c("affinity", "rank.affinity", "equilibrate") | grepl("solubilit", efun)))
    stop("'eout' is not the output from one of these functions: affinity, rank.affinity, equilibrate, or solubility")
  # For solubility(), default type is loga.balance 20210303
  if(grepl("solubilit", efun) & missing(type)) type <- "loga.balance"
  # Check balance argument for rank.affinity() 20220416
  if(efun == "rank.affinity") {
    if(!identical(balance, 1)) {
      if(!is.null(balance)) stop("balance = 1 or NULL is required for plotting output of rank.affinity()")
      if(is.null(balance)) message("diagram: setting balance = 1 for plotting output of rank.affinity()")
      balance <- 1
    }
  }

  ## 'type' can be:
  #    'auto'                - values returned by affinity() (aout)
  #    'loga.equil'          - activities of formed species (eout from equilibrate() or solubility())
  #    name of basis species - activity of a basis species (aout)
  #    'saturation'          - affinity=0 line for each species (aout)
  #    'loga.balance'        - total activities of formed species (eout from solubility())
  eout.is.aout <- FALSE
  plot.loga.basis <- FALSE
  if(type %in% c("auto", "saturation")) {
    if(!"loga.equil" %in% names(eout)) {
      eout.is.aout <- TRUE
      # Get the balancing coefficients
      if(type == "auto") {
        bal <- balance(eout, balance)
        n.balance <- bal$n.balance
        balance <- bal$balance
      } else n.balance <- rep(1, length(eout$values))
      # In case all coefficients are negative (for bimetal() examples) 20200713
      #   e.g. H+ for minerals with basis(c("Fe+2", "Cu+", "hydrogen sulfide", "oxygen", "H2O", "H+"))
      if(all(n.balance < 0)) n.balance <- -n.balance
    }
  } else if(type %in% rownames(eout$basis)) {
    # To calculate the loga of basis species at equilibrium
    if(!missing(groups)) stop("can't plot equilibrium activities of basis species for grouped species")
    if(isTRUE(alpha) | is.character(alpha)) stop("equilibrium activities of basis species not available with alpha = TRUE")
    plot.loga.basis <- TRUE
  } else if(type == "loga.equil" & !"loga.equil" %in% names(eout)) stop("'eout' is not the output from equil()") 
  else if(!type %in% c("loga.equil", "loga.balance")) stop(type, " is not a valid diagram type")

  ## Consider a different number of species if we're grouping them together
  ngroups <- length(groups)

  ## Keep the values we plot in one place so they can be modified, plotted and eventually returned
  # Unless something happens below, we'll plot the loga.equil from equilibrate()
  plotvals <- eout$loga.equil
  plotvar <- "loga.equil"

  ## Handle loga.balance here (i.e. solubility calculations)
  if(type == "loga.balance") {
    plotvals <- list(eout$loga.balance)
    plotvar <- "loga.balance"
  }

  ## Modify default arguments with add = TRUE
  if(add) {
    # Change missing fill.NA to transparent
    if(missing(fill.NA)) fill.NA <- "transparent"
    # Change missing limit.water to FALSE 20200819
    if(missing(limit.water)) limit.water <- FALSE
  }

  ## Number of dimensions (T, P or chemical potentials that are varied)
  # length(eout$vars) - the number of variables = the maximum number of dimensions
  # length(dim(eout$values[[1]])) - nd = 1 if it was a transect along multiple variables
  nd <- min(length(eout$vars), length(dim(eout$values[[1]])))

  ## Deal with output from affinity()
  if(eout.is.aout) {
    # Plot property from affinity(), divided by balancing coefficients
    plotvals <- lapply(1:length(eout$values), function(i) {
      # We divide by the balancing coefficients if we're working with affinities
      # This is not normalizing the formulas! it's balancing the reactions...
      #   normalizing the formulas is done below
      eout$values[[i]] / n.balance[i]
    })
    names(plotvals) <- names(eout$values)
    plotvar <- eout$property
    if(efun == "rank.affinity") {
      plotvar <- "rank.affinity"
      message(paste("diagram: plotting average affinity ranking for", length(plotvals), "groups"))
    } else if(plotvar == "A") {
      # We change 'A' to 'A/(2.303RT)' so the axis label is made correctly
      # 20171027 use parentheses to avoid ambiguity about order of operations
      plotvar <- "A/(2.303RT)"
      if(nd == 2 & type == "auto") message("diagram: using maximum affinity method for 2-D diagram")
      else if(nd == 2 & type == "saturation") message("diagram: plotting saturation lines for 2-D diagram")
      else message("diagram: plotting A/(2.303RT) / n.balance")
    } else message(paste("diagram: plotting", plotvar, " / n.balance"))
  }

  ## Use molality instead of activity if the affinity calculation includes ionic strength 20171101
  molality <- "IS" %in% names(eout)

  ## When can normalize and as.residue be used
  if(normalize | as.residue) {
    if(normalize & as.residue) stop("'normalize' and 'as.residue' can not both be TRUE")
    if(!eout.is.aout) stop("'normalize' or 'as.residue' can be TRUE only if 'eout' is the output from affinity()")
    if(nd != 2) stop("'normalize' or 'as.residue' can be TRUE only for a 2-D (predominance) diagram")
    if(normalize) message("diagram: using 'normalize' in calculation of predominant species")
    else message("diagram: using 'as.residue' in calculation of predominant species")
  }

  ## Sum affinities or activities of species together in groups 20090524
  # Use lapply/Reduce 20120927
  if(!missing(groups)) {
    # Loop over the groups
    plotvals <- lapply(groups, function(ispecies) {
      # Remove the logarithms ...
      if(eout.is.aout) act <- lapply(plotvals[ispecies], function(x) 10^x)
      # and, for activity, multiply by n.balance 20170207
      else act <- lapply(seq_along(ispecies), function(i) eout$n.balance[ispecies[i]] * 10^plotvals[[ispecies[i]]])
      # then, sum the activities
      return(Reduce("+", act))
    })
    # Restore the logarithms
    plotvals <- lapply(plotvals, function(x) log10(x))
    # Combine the balancing coefficients for calculations using affinity
    if(eout.is.aout) n.balance <- sapply(groups, function(ispecies) sum(n.balance[ispecies]))
  }

  ## Calculate the equilibrium logarithm of activity of a basis species
  ## (such that affinities of formation reactions are zero)
  if(plot.loga.basis) {
    ibasis <- match(type, rownames(eout$basis))
    # The logarithm of activity used in the affinity calculation
    is.loga.basis <- can.be.numeric(eout$basis$logact[ibasis])
    if(!is.loga.basis) stop(paste("the logarithm of activity for basis species", type, "is not numeric - was a buffer selected?"))
    loga.basis <- as.numeric(eout$basis$logact[ibasis])
    # The reaction coefficients for this basis species
    nu.basis <- eout$species[, ibasis]
    # The logarithm of activity where affinity = 0
    plotvals <- lapply(1:length(eout$values), function(x) {
      # NB. eout$values is a strange name for affinity ... should be named something like eout$affinity ...
      loga.basis - eout$values[[x]]/nu.basis[x]
    })
    plotvar <- type
  }

  ## alpha: plot fractional degree of formation
  # Scale the activities to sum = 1  ... 20091017
  # Allow scaling by balancing component 20171008
  if(isTRUE(alpha) | is.character(alpha)) {
    # Remove the logarithms
    act <- lapply(plotvals, function(x) 10^x)
    if(identical(alpha, "balance")) for(i in 1:length(act)) act[[i]] <- act[[i]] * eout$n.balance[i]
    # Sum the activities
    sumact <- Reduce("+", act)
    # Divide activities by the total
    alpha <- lapply(act, function(x) x/sumact)
    plotvals <- alpha
    plotvar <- "alpha"
  }

  ## Identify predominant species
  predominant <- NA
  H2O.predominant <- NULL
  # Whether we're considering multiple species, based on the plotting variable
  pv_multi <- plotvar %in% c("loga.equil", "alpha", "A/(2.303RT)", "rank.affinity") & type != "saturation"
  if(pv_multi) {
    pv <- plotvals
    # Some additional steps for affinity values, but not for equilibrated activities
    if(eout.is.aout) {
      for(i in 1:length(pv)) {
        # TODO: The equilibrium.Rmd vignette shows predominance diagrams using
        #   'normalize' and 'as.residue' settings that are consistent with equilibrate(),
        #   but what is the derivation of these equations?
        if(normalize) pv[[i]] <- (pv[[i]] + eout$species$logact[i] / n.balance[i]) - log10(n.balance[i])
        else if(as.residue) pv[[i]] <- pv[[i]] + eout$species$logact[i] / n.balance[i]
      }
    }
    if(grepl("solubilities", efun)) {
      # For solubilites of multiple minerals, find the minimum value (most stable mineral) 20210321
      mypv <- Map("-", pv)
      predominant <- which.pmax(mypv)
    } else {
      # For all other diagrams, we want the maximum value (i.e. maximum affinity method)
      predominant <- which.pmax(pv)
    }

    # Show water stability region
    if((is.null(limit.water) | isTRUE(limit.water)) & nd == 2) {
      wl <- water.lines(eout, plot.it = FALSE)
      # Proceed if water.lines produced calculations for this plot
      if(!identical(wl, NA)) {
        H2O.predominant <- predominant
        # For each x-point, find the y-values that are outside the water stability limits
        for(i in seq_along(wl$xpoints)) {
          ymin <- min(c(wl$y.oxidation[i], wl$y.reduction[i]))
          ymax <- max(c(wl$y.oxidation[i], wl$y.reduction[i]))
          if(!wl$swapped) {
            # The actual calculation
            iNA <- eout$vals[[2]] < ymin | eout$vals[[2]] > ymax
            # Assign NA to the predominance matrix
            H2O.predominant[i, iNA] <- NA
          } else {
            # As above, but x- and y-axes are swapped
            iNA <- eout$vals[[1]] < ymin | eout$vals[[1]] > ymax
            H2O.predominant[iNA, i] <- NA
          }
        }
      }
    }
  }

  ## Create some names for lines/fields if they are missing
  is.pname <- FALSE
  onames <- names
  if(identical(names, FALSE) | identical(names, NA)) names <- ""
  else if(!is.character(names)) {
    # Properties of basis species or reactions?
    if(eout$property %in% c("G.basis", "logact.basis")) names <- rownames(eout$basis)
    else {
      if(!missing(groups)) {
        if(is.null(names(groups))) names <- paste("group", 1:length(groups), sep = "")
        else names <- names(groups)
      }
      else names <- as.character(eout$species$name)
      # Remove non-unique organism or protein names
      if(all(grepl("_", names))) {
        is.pname <- TRUE
        # Everything before the underscore (the protein)
        pname <- gsub("_.*$", "", names)
        # Everything after the underscore (the organism)
        oname <- gsub("^.*_", "", names)
        # If the pname or oname are all the same, use the other one as identifying name
        if(length(unique(pname)) == 1) names <- oname
        if(length(unique(oname)) == 1) names <- pname
      }
      # Append state to distinguish ambiguous species names
      isdup <- names %in% names[duplicated(names)]
      if(any(isdup)) names[isdup] <- paste(names[isdup],
        " (", eout$species$state[isdup], ")", sep = "")
    }
  }
  # Numeric values indicate a subset 20181007
  if(all(is.numeric(onames))) {
    if(isTRUE(all(onames > 0))) names[-onames] <- ""
    else if(isTRUE(all(onames < 0))) names[-onames] <- ""
    else stop("numeric 'names' should be all positive or all negative")
  }

  ## Apply formatting to chemical formulas 20170204
  if(all(grepl("_", names))) is.pname <- TRUE
  if(format.names & !is.pname) {
    # Check if names are a deparsed expression (used in mix()) 20200718
    parsed <- FALSE
    if(any(grepl("paste\\(", names))) {
      exprnames <- parse(text = names)
      if(length(exprnames) != length(names)) stop("parse()-ing names gives length not equal to number of names")
      parsed <- TRUE
    } else {
      exprnames <- as.expression(names)
      # Get formatted chemical formulas
      for(i in seq_along(exprnames)) {
        # Don't try to format the names if they have "+" followed by a character
        #   (created by mix()ing diagrams with format.names = FALSE);
        #   expr.species() can't handle it 20200722
        if(!grepl("\\+[a-zA-Z]", names[i])) exprnames[[i]] <- expr.species(exprnames[[i]])
      }
    }
    # Apply bold or italic
    bold <- rep(bold, length.out = length(exprnames))
    italic <- rep(italic, length.out = length(exprnames))
    for(i in seq_along(exprnames)) {
      if(bold[i]) exprnames[[i]] <- substitute(bold(a), list(a = exprnames[[i]]))
      if(italic[i]) exprnames[[i]] <- substitute(italic(a), list(a = exprnames[[i]]))
    }
    # Only use the expression if it's different from the unformatted names
    if(parsed | !identical(as.character(exprnames), names)) names <- exprnames
  }

  ## Where we'll put extra output for predominance diagrams (namesx, namesy)
  out2D <- list()

  ### Now we're getting to the plotting ###

  if(plot.it) {

    ### General plot parameters ###

    ## Handle line type/width/color arguments
    if(is.null(lty)) {
      if(type == "loga.balance" | nd == 2) lty <- 1
      else lty <- 1:ngroups
    }
    lty <- rep(lty, length.out = length(plotvals))
    lwd <- rep(lwd, length.out = length(plotvals))
    col <- rep(col, length.out = length(plotvals))
    
    # Function to get label for i'th variable 20230809
    # (uses custom labels from 'labels' list element added by mosaic)
    getlabel <- function(ivar) {
      label <- eout$vars[ivar]
      if(!is.null(eout$labels)) {
        if(label %in% names(eout$labels)) {
          label <- eout$labels[[label]]
        }
      }
      label
    }

    if(nd == 0) {

      ### 0-D diagram - bar graph of properties of species or reactions
      # Plot setup
      if(missing(ylab)) ylab <- axis.label(plotvar, units = "", molality = molality)
      barplot(unlist(plotvals), names.arg = names, ylab = ylab, cex.names = cex.names, col = col, ...)
      if(!is.null(main)) title(main = main)

    } else if(nd == 1) {

      ### 1-D diagram - lines for properties or chemical activities
      xvalues <- eout$vals[[1]]
      if(missing(xlim)) xlim <- c(xvalues[1], rev(xvalues)[1])
      # Initialize the plot
      if(!add) {
        if(missing(xlab)) xlab <- axis.label(getlabel(1), basis = eout$basis, molality = molality)
        if(missing(ylab)) {
          ylab <- axis.label(plotvar, units = "", molality = molality)
          if(plotvar == "rank.affinity") ylab <- "Average affinity ranking"
          # Use ppb, ppm, ppt (or log ppb etc.) for converted values of solubility 20190526
          if(grepl("solubility.", eout$fun, fixed = TRUE)) {
            ylab <- strsplit(eout$fun, ".", fixed = TRUE)[[1]][2]
            ylab <- gsub("log", "log ", ylab)
          }
        }
        # To get range for y-axis, use only those points that are in the xrange
        if(is.null(ylim)) {
          isx <- xvalues >= min(xlim) & xvalues <= max(xlim)
          xfun <- function(x) x[isx]
          myval <- sapply(plotvals, xfun)
          ylim <- extendrange(myval)
        }
        if(tplot) thermo.plot.new(xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, cex = cex, mar = mar, yline = yline, side = side, ...)
        else plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
      }
      # Draw the lines
      spline.n <- 256 # the number of values at which to calculate splines
      if(is.null(spline.method) | length(xvalues) > spline.n) {
        for(i in 1:length(plotvals)) lines(xvalues, plotvals[[i]], col = col[i], lty = lty[i], lwd = lwd[i])
      } else {
        # Plot splines instead of lines connecting the points 20171116
        spline.x <- seq(xlim[1], xlim[2], length.out = spline.n)
        for(i in 1:length(plotvals)) {
          spline.y <- splinefun(xvalues, plotvals[[i]], method = spline.method)(spline.x)
          lines(spline.x, spline.y, col = col[i], lty = lty[i], lwd = lwd[i])
        }
      }
      if(type %in% c("auto", "loga.equil") & !is.null(legend.x)) {
        # 20120521: use legend.x = NA to label lines rather than make legend
        if(is.na(legend.x)) {
          maxvals <- do.call(pmax, pv)
          # Label placement and rotation
          dx <- rep(dx, length.out = length(plotvals))
          dy <- rep(dy, length.out = length(plotvals))
          srt <- rep(srt, length.out = length(plotvals))
          cex.names <- rep(cex.names, length.out = length(plotvals))
          # Don't assign to adj becuase that messes up the missing test below
          alladj <- rep(adj, length.out = length(plotvals))
          for(i in 1:length(plotvals)) {
            # y-values for this line
            myvals <- as.numeric(plotvals[[i]])
            # Don't take values that lie close to or above the top of plot
            myvals[myvals > ylim[1] + 0.95*diff(ylim)] <- ylim[1]
            # If we're adding to a plot, don't take values that are above the top of this plot
            if(add) {
              this.ylim <- par("usr")[3:4]
              myvals[myvals > this.ylim[1] + 0.95*diff(this.ylim)] <- this.ylim[1]
            }
            # The starting x-adjustment
            thisadj <- alladj[i]
            # If this line has any of the overall maximum values, use only those values
            # (useful for labeling straight-line affinity comparisons 20170221)
            is.max <- myvals == maxvals
            if(any(is.max) & plotvar != "alpha") {
              # Put labels on the median x-position
              imax <- median(which(is.max))
            } else {
              # Put labels on the maximum of the line
              # (useful for labeling alpha plots)
              imax <- which.max(myvals)
              # Avoid the sides of the plot; take care of reversed x-axis
              if(missing(adj)) {
                if(sign(diff(xlim)) > 0) {
                  if(xvalues[imax] > xlim[1] + 0.8*diff(xlim)) thisadj <- 1
                  if(xvalues[imax] < xlim[1] + 0.2*diff(xlim)) thisadj <- 0
                } else {
                  if(xvalues[imax] > xlim[1] + 0.2*diff(xlim)) thisadj <- 0
                  if(xvalues[imax] < xlim[1] + 0.8*diff(xlim)) thisadj <- 1
                }
              }
            }
            # Also include y-offset (dy) and y-adjustment (labels bottom-aligned with the line)
            # ... and srt (string rotation) 20171127
            text(xvalues[imax] + dx[i], plotvals[[i]][imax] + dy[i], labels = names[i], adj = c(thisadj, 0),
              cex = cex.names[i], srt = srt[i], font = font, family = family)
          }
        } else legend(x = legend.x, lty = lty, legend = names, col = col, cex = cex.names, lwd = lwd, ...)
      }
      # Add a title
      if(!is.null(main)) title(main = main)

    } else if(nd == 2) {

      ### 2-D diagram - fields indicating species predominance, or contours for other properties

      ### Functions for constructing predominance area diagrams
      ## Color fill function
      fill.color <- function(xs, ys, out, fill, nspecies) {
        # Handle min/max reversal
        if(xs[1] > xs[length(xs)]) {
          tc <- out
          t <- numeric()
          for(i in 1:length(xs)) {
            t[i] <- xs[length(xs)+1-i]
            tc[, i] <- out[, length(xs)+1-i]
          }
          out <- tc
          xs <- t
        }
        if(ys[1] > ys[length(ys)]) {
          tc <- out
          t <- numeric()
          for(i in 1:length(ys)) {
            t[i] <- ys[length(ys)+1-i]
            tc[i, ] <- out[length(ys)+1-i, ]
          }
          out <- tc
          ys <- t
        }
        # The z values
        zs <- out
        for(i in 1:nrow(zs)) zs[i,] <- out[nrow(zs)+1-i,]
        zs <- t(zs)
        breaks <- c(-1, 0, 1:nspecies) + 0.5
        # Use fill.NA for NA values
        zs[is.na(zs)] <- 0
        image(x = xs, y = ys, z = zs, col = c(fill.NA, fill), add = TRUE, breaks = breaks, useRaster = TRUE)
      }

      ## Curve plot function
      # 20091116 replaced plot.curve with plot.line; different name, same functionality, *much* faster
      plot.line <- function(out, xlim, ylim, dotted, col, lwd, xrange) {
        # Plot boundary lines between predominance fields
        vline <- function(out, ix) {
          ny <- nrow(out)
          xs <- rep(ix, ny*2+1)
          ys <- c(rep(ny:1, each = 2), 0)
          y1 <- out[, ix]
          y2 <- out[, ix+1]
          # No line segment inside a stability field
          iy <- which(y1 == y2)
          ys[iy*2] <- NA
          # No line segment at a dotted position
          iyd <- rowSums(sapply(dotted, function(y) ys%%y == 0)) > 0
          ys[iyd] <- NA
          return(list(xs = xs, ys = ys))
        }
        hline <- function(out, iy) {
          nx <- ncol(out)
          ys <- rep(iy, nx*2+1)
          xs <- c(0, rep(1:nx, each = 2))
          x1 <- out[iy, ]
          x2 <- out[iy+1, ]
          # No line segment inside a stability field
          ix <- which(x1 == x2)
          xs[ix*2] <- NA
          # No line segment at a dotted position
          ixd <- rowSums(sapply(dotted, function(x) xs%%x == 0)) > 0
          xs[ixd] <- NA
          return(list(xs = xs, ys = ys))
        }
        clipfun <- function(z, zlim) {
          if(zlim[2] > zlim[1]) {
            z[z>zlim[2]] <- NA
            z[z<zlim[1]] <- NA
          } else {
            z[z>zlim[1]] <- NA
            z[z<zlim[2]] <- NA
          }
          return(z)
        }
        rx <- (xlim[2] - xlim[1]) / (ncol(out) - 1)
        ry <- (ylim[2] - ylim[1]) / (nrow(out) - 1)
        # Vertical lines
        xs <- ys <- NA
        for(ix in 1:(ncol(out)-1)) {
          vl <- vline(out,ix)
          xs <- c(xs,vl$xs,NA)
          ys <- c(ys,vl$ys,NA)
        }
        xs <- xlim[1] + (xs - 0.5) * rx
        ys <- ylim[1] + (ys - 0.5) * ry
        ys <- clipfun(ys, ylim)
        if(!is.null(xrange)) xs <- clipfun(xs, xrange)
        lines(xs, ys, col = col, lwd = lwd)
        # Horizontal lines
        xs <- ys <-NA
        for(iy in 1:(nrow(out)-1)) {
          hl <- hline(out, iy)
          xs <- c(xs, hl$xs, NA)
          ys <- c(ys, hl$ys, NA)
        }
        xs <- xlim[1] + (xs - 0.5) * rx
        ys <- ylim[2] - (ys - 0.5) * ry
        xs <- clipfun(xs, xlim)
        if(!is.null(xrange)) xs <- clipfun(xs, xrange)
        lines(xs, ys, col = col, lwd = lwd)
      }

      ## New line plotting function 20170122
      contour.lines <- function(predominant, xlim, ylim, lty, col, lwd) {
        # The x and y values
        xs <- seq(xlim[1], xlim[2], length.out = dim(predominant)[1])
        ys <- seq(ylim[1], ylim[2], length.out = dim(predominant)[2])
        # Reverse any axis that has decreasing values
        if(diff(xlim) < 0) {
          predominant <- predominant[nrow(predominant):1, ]
          xs <- rev(xs)
        }
        if(diff(ylim) < 0) {
          predominant <- predominant[, ncol(predominant):1]
          ys <- rev(ys)
        }
	# The categories (species/groups/etc) on the plot
	zvals <- na.omit(unique(as.vector(predominant)))
        
        # Initialize list and counter for line x,y values 20240615
        linesout <- list()
        iout <- 1

        if(is.null(lty.aq) & is.null(lty.cr)) {

          # DEFAULT method: loop over species
          for(i in 1:(length(zvals)-1)) {
            # Get the "z" values
            # Use + 0 trick to convert T/F to 1/0 20220524
            z <- (predominant == zvals[i]) + 0
            z[is.na(z)] <- 0
            # Use contourLines() instead of contour() in order to get line coordinates 20181029
            cLines <- contourLines(xs, ys, z, levels = 0.5)
            if(length(cLines) > 0) {
              # Loop in case contourLines returns multiple lines
              for(k in 1:length(cLines)) {
                # Draw the lines
                lines(cLines[[k]][2:3], lty = lty[zvals[i]], col = col[zvals[i]], lwd = lwd[zvals[i]])

                # Store the x and y values (list components 2 and 3)
                linesout[[iout]] <- cLines[[k]][[2]]
                names(linesout)[iout] <- paste0("x", k, "_", zvals[i])
                linesout[[iout+1]] <- cLines[[k]][[3]]
                names(linesout)[iout+1] <- paste0("y", k, "_", zvals[i])
                iout <- iout + 2
              }
            }
            # Mask species to prevent double-plotting contour lines
            predominant[z == zvals[i]] <- NA
          }

        } else {

          # ALTERNATE (slower) method to handle lty.aq and lty.cr: take each possible pair of species
          # Reinstated on 20210301
          for(i in 1:(length(zvals)-1)) {
            for(j in (i+1):length(zvals)) {
              z <- predominant
              # Draw contours only for this pair
              z[!z %in% c(zvals[i], zvals[j])] <- NA
              # Give them neighboring values (so we get one contour line)
              z[z == zvals[i]] <- 0
              z[z == zvals[j]] <- 1
              # Use contourLines() instead of contour() in order to get line coordinates 20181029
              cLines <- contourLines(xs, ys, z, levels = 0.5)
              if(length(cLines) > 0) {
                # Loop in case contourLines returns multiple lines
                for(k in 1:length(cLines)) {
                  # Draw the lines
                  mylty <- lty[zvals[i]]
                  if(!is.null(lty.cr)) {
                    # Use lty.cr for cr-cr boundaries 20190530
                    if(all(grepl("cr", eout$species$state[c(zvals[i], zvals[j])]))) mylty <- lty.cr
                  }
                  if(!is.null(lty.aq)) {
                    # Use lty.aq for aq-aq boundaries 20190531
                    if(all(grepl("aq", eout$species$state[c(zvals[i], zvals[j])]))) mylty <- lty.aq
                  }
                  lines(cLines[[k]][2:3], lty = mylty, col = col[zvals[i]], lwd = lwd[zvals[i]])

                  # Store the x and y values (list components 2 and 3)
                  linesout[[iout]] <- cLines[[k]][[2]]
                  names(linesout)[iout] <- paste0("x", k, "_", zvals[i], ".", zvals[j])
                  linesout[[iout+1]] <- cLines[[k]][[3]]
                  names(linesout)[iout+1] <- paste0("y", k, "_", zvals[i], ".", zvals[j])
                  iout <- iout + 2

                }
              }
            }
          }

        }

        # Return x,y coordinates of lines padded to equal length
        # https://stackoverflow.com/questions/34570860/adding-na-to-make-all-list-elements-equal-length 20181029
        # For compatibility with R 3.1.0, don't use lengths() here 20190302
        lapply(linesout, `length<-`, max(sapply(linesout, length)))

      }

      ## To add labels
      plot.names <- function(out, xs, ys, xlim, ylim, names, srt, min.area) {
        # Calculate coordinates for field labels
        # Revisions: 20091116 for speed, 20190223 work with user-specified xlim and ylim
        namesx <- namesy <- rep(NA, length(names))
        # Even if 'names' is NULL, we run the loop in order to generate namesx and namesy for the output 20190225
        area.plot <- length(xs) * length(ys)
        for(i in seq_along(groups)) {
          this <- which(out == i, arr.ind = TRUE)
          if(length(this) == 0) next
          xsth <- xs[this[, 2]]
          ysth <- rev(ys)[this[, 1]]
          # Use only values within the plot range
          rx <- range(xlim)
          ry <- range(ylim)
          xsth <- xsth[xsth >= rx[1] & xsth <= rx[2]]
          ysth <- ysth[ysth >= ry[1] & ysth <= ry[2]]
          if(length(xsth) == 0 | length(ysth) == 0) next
          # Skip plotting names if the fields are too small 20200720
          area <- max(length(xsth), length(ysth))
          frac.area <- area / area.plot
          if(!frac.area >= min.area) next
          namesx[i] <- mean(xsth)
          namesy[i] <- mean(ysth)
        }
        # Fields that really exist on the plot
        if(!is.null(names)) {
          cex <- rep(cex.names, length.out = length(names))
          col <- rep(col.names, length.out = length(names))
          font <- rep(font, length.out = length(names))
          family <- rep(family, length.out = length(names))
          srt <- rep(srt, length.out = length(names))
          dx <- rep(dx, length.out = length(names))
          dy <- rep(dy, length.out = length(names))
          for(i in seq_along(names)) {
            if(!(identical(col[i], 0)) & !is.na(col[i]))
              text(namesx[i] + dx[i], namesy[i] + dy[i], labels = names[i], cex = cex[i], col = col[i], font = font[i], family = family[i], srt = srt[i])
          }
        }
        return(list(namesx = namesx, namesy = namesy))
      }

      ### Done with supporting functions
      ### Now we really get to make the diagram itself

      # Colors to fill predominance fields
      if(is.null(fill)) {
        if(length(unique(eout$species$state)) == 1) fill <- "transparent"
        else fill <- ifelse(grepl("cr,cr", eout$species$state), "#DEB88788", ifelse(grepl("cr", eout$species$state), "#FAEBD788", "#F0F8FF88"))
      } else if(identical(fill, NA) | length(fill) == 0) fill <- "transparent"
      else if(isTRUE(fill[1] == "rainbow")) fill <- rainbow(ngroups)
      else if(isTRUE(fill[1] %in% c("heat", "terrain", "topo", "cm"))) fill <- get(paste0(fill[1], ".colors"))(ngroups)
      else if(getRversion() >= "3.6.0" & length(fill) == 1) {
        # Choose an HCL palette 20190411
        # Matching adapted from hcl.colors()
        fx <- function(x) tolower(gsub("[-, _, \\,, (, ), \\ , \\.]", "", x))
        p <- charmatch(fx(fill), fx(hcl.pals()))
        if(!is.na(p)) {
          if(!p < 1L) {
            fill <- hcl.colors(ngroups, fill)
          }
        }
      }
      fill <- rep(fill, length.out = ngroups)
      # The x and y values 
      xs <- eout$vals[[1]]
      ys <- eout$vals[[2]]
      # The limits of the calculation; they aren't necessarily increasing, so don't use range()
      xlim.calc <- c(xs[1], tail(xs, 1))
      ylim.calc <- c(ys[1], tail(ys, 1))
      # Add if(is.null) to allow user-specified limits 20190223
      if(is.null(xlim)) {
        if(add) xlim <- par("usr")[1:2]
        else xlim <- xlim.calc
      }
      if(is.null(ylim)) {
        if(add) ylim <- par("usr")[3:4]
        else ylim <- ylim.calc
      }
      # Initialize the plot
      if(!add) {
        if(is.null(xlab)) xlab <- axis.label(getlabel(1), basis = eout$basis, molality = molality)
        if(is.null(ylab)) ylab <- axis.label(getlabel(2), basis = eout$basis, molality = molality)
        if(tplot) thermo.plot.new(xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
          cex = cex, cex.axis = cex.axis, mar = mar, yline = yline, side = side, ...)
        else plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
        # Add a title
        if(!is.null(main)) title(main = main)
      }

      # Start with NA value for x,y locations of lines
      linesout <- NA

      if(identical(predominant, NA)) {

        # No predominance matrix, so we're contouring properties
        contour.method <- rep(contour.method, length.out = length(plotvals))
        drawlabels <- TRUE
        if(identical(contour.method, NULL) | identical(contour.method[1], NA) | identical(contour.method[1], "")) drawlabels <- FALSE
        if(type == "saturation") {
          # For saturation plot, contour affinity = 0 for all species
          for(i in 1:length(plotvals)) {
            zs <- plotvals[[i]]
            # Skip plotting if this species has no possible saturation line, or a line outside the plot range
            if(length(unique(as.numeric(zs))) == 1) {
              message("diagram: no saturation line possible for ", names[i])
              next
            }
            if(all(zs < 0) | all(zs > 0)) {
              message("diagram: beyond range for saturation line of ", names[i])
              next
            }
            if(drawlabels) contour(xs, ys, zs, add = TRUE, col = col, lty = lty, lwd = lwd, labcex = cex, levels = 0, labels = names[i], method = contour.method[i])
            else contour(xs, ys, zs, add = TRUE, col = col, lty = lty, lwd = lwd, labcex = cex, levels = 0, labels = names[i], drawlabels = FALSE)
          }
        } else {
          # Contour solubilities (loga.balance), or properties using first species only
          if(length(plotvals) > 1) warning("showing only first species in 2-D property diagram")
          zs <- plotvals[[1]]
          if(is.null(levels)) {
            if(drawlabels) contour(xs, ys, zs, add = TRUE, col = col, lty = lty, lwd = lwd, labcex = cex, method = contour.method[1])
            else contour(xs, ys, zs, add = TRUE, col = col, lty = lty, lwd = lwd, labcex = cex, drawlabels = FALSE)
          } else {
            if(drawlabels) contour(xs, ys, zs, add = TRUE, col = col, lty = lty, lwd = lwd, labcex = cex, method = contour.method[1], levels = levels)
            else contour(xs, ys, zs, add = TRUE, col = col, lty = lty, lwd = lwd, labcex = cex, levels = levels, drawlabels = FALSE)
          }
        }
        # Keep the x,y coordinates of the names around to add to the output
        pn <- list(namesx = NULL, namesy = NULL)

      } else {

        # Put predominance matrix in the right order for image()
        zs <- t(predominant[, ncol(predominant):1])
        if(!is.null(H2O.predominant)) {
          zsH2O <- t(H2O.predominant[, ncol(H2O.predominant):1])
          # This colors in fields and greyed-out H2O non-stability regions
          if(!is.null(fill)) fill.color(xs, ys, zsH2O, fill, ngroups)
          # Clip entire diagram to H2O stability region?
          if(isTRUE(limit.water)) zs <- zsH2O
        }
        # This colors in fields (possibly a second time, to overlay on H2O regions)
        if(!is.null(fill)) fill.color(xs, ys, zs, fill, ngroups)

        # Plot field labels
        pn <- plot.names(zs, xs, ys, xlim, ylim, names, srt, min.area)
        # Only draw the lines if there is more than one field  20180923
        #   (to avoid warnings from contour, which seem to be associated with weird
        #   font metric state and subsequent errors adding e.g. subscripted text to plot)
        if(length(na.omit(unique(as.vector(zs)))) > 1) {
          if(!is.null(dotted)) plot.line(zs, xlim.calc, ylim.calc, dotted, col, lwd, xrange = xrange)
          else linesout <- contour.lines(predominant, xlim.calc, ylim.calc, lty = lty, col = col, lwd = lwd)
        }
        # Re-draw the tick marks and axis lines in case the fill obscured them
        has.color <- FALSE
        if(!identical(unique(fill), "transparent")) has.color <- TRUE
        if(!is.null(H2O.predominant) & !identical(fill.NA, "transparent")) has.color <- TRUE
        if(any(is.na(zs)) & !identical(fill.NA, "transparent")) has.color <- TRUE
        if(tplot & !add & has.color) {
          thermo.axis()
          box()
        }
      }

      # Done with the 2D plot!
      out2D <- list(namesx = pn$namesx, namesy = pn$namesy, linesout = linesout)

    } # end if(nd == 2)
  } # end if(plot.it)

  # Even if plot = FALSE, return the diagram clipped to the water stability region (for testing) 20200719
  if(isTRUE(limit.water) & !is.null(H2O.predominant)) predominant <- H2O.predominant

  # Make a matrix with the affinities or activities of predominant species 20200724
  # (for calculating affinities of metastable species - multi-metal.Rmd example)
  predominant.values <- NA
  if(!identical(predominant, NA)) {
    predominant.values <- plotvals[[1]]
    predominant.values[] <- NA
    for(ip in na.omit(unique(as.vector(predominant)))) {
      ipp <- predominant == ip
      ipp[is.na(ipp)] <- FALSE
      # 20201219 Change eout$values to plotvals (return normalized affinities or equilibrium activities)
      predominant.values[ipp] <- plotvals[[ip]][ipp]
    }
  }

  outstuff <- list(plotvar = plotvar, plotvals = plotvals, names = names, predominant = predominant, predominant.values = predominant.values)
  # Include the balance name and coefficients if we diagrammed affinities 20200714
  if(eout.is.aout) outstuff <- c(list(balance = balance, n.balance = n.balance), outstuff)
  out <- c(eout, outstuff, out2D)
  invisible(out)
}

find.tp <- function(x) {
  # Find triple points in an matrix of integers  20120525 jmd
  #   these are the locations closest to the greatest number of different values
  # Rearrange the matrix in the same way that diagram() does for 2-D predominance diagrams
  x <- t(x[, ncol(x):1])
  # All of the indexes for the matrix
  id <- which(x > 0, arr.ind = TRUE)
  # We'll do a brute-force count at each position
  n <- sapply(1:nrow(id), function(i) {
    # Row and column range to look at (3x3 except at edges)
    r1 <- max(id[i, 1] - 1, 0)
    r2 <- min(id[i, 1] + 1, nrow(x))
    c1 <- max(id[i, 2] - 1, 0)
    c2 <- min(id[i, 2] + 1, ncol(x))
    # The number of unique values
    return(length(unique(as.numeric(x[r1:r2, c1:c2]))))
  })
  # Which positions have the most counts?
  imax <- which(n == max(n))
  # Return the indices
  return(id[imax, ])
}
