# CHNOSZ/diagram.R
# plot equilibrium chemical activity and predominance diagrams 
# 20061023 jmd v1
# 20120927 work with output from either equil() or affinity(), 
#   gather plotvals independently of plot parameters (including nd),
#   single return statement

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("equilibrate.R")
#source("util.plot.R")
#source("util.character.R")
#source("util.misc.R")

diagram <- function(
  # species affinities or activities
  eout, 
  # type of plot
  type="auto", alpha=FALSE, normalize=FALSE, as.residue=FALSE, balance=NULL,
  groups=as.list(1:length(eout$values)), xrange=NULL,
  # figure size and sides for axis tick marks
  mar=NULL, yline=par("mgp")[1]+0.3, side=1:4,
  # axis limits and labels
  ylog=TRUE, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
  # sizes
  cex=par("cex"), cex.names=1, cex.axis=par("cex"),
  # line styles
  lty=NULL, lwd=par("lwd"), dotted=NULL, spline.method=NULL, contour.method="edge",
  # colors
  col=par("col"), col.names=par("col"), fill=NULL,
  fill.NA="slategray1", limit.water=TRUE,
  # field and line labels
  names=NULL, format.names=TRUE, bold=FALSE, italic=FALSE,
  font=par("font"), family=par("family"), adj=0.5, dy=0, srt=0,
  # title and legend
  main=NULL, legend.x=NA,
  # plotting controls
  add=FALSE, plot.it=TRUE, tplot=TRUE, ...
) {

  ### argument handling ###

  ## check that eout is valid input
  if(!"sout" %in% names(eout)) stop("'eout' does not look like output from equil() or affinity()")

  ## 'type' can be:
  #    'auto'                - property from affinity() (1D) or maximum affinity (affinity 2D) (aout) or loga.equil (eout)
  #    'loga.equil'          - equilibrium activities of species of interest (eout)
  #    name of basis species - equilibrium activity of a basis species (aout)
  #    'saturation'          - affinity=0 line for each species (2D)
  eout.is.aout <- FALSE
  plot.loga.basis <- FALSE
  if(type %in% c("auto", "saturation")) {
    if(!"loga.equil" %in% names(eout)) {
      eout.is.aout <- TRUE
      # get the balancing coefficients
      if(type=="auto") n.balance <- balance(eout, balance)$n.balance
      else n.balance <- rep(1, length(eout$values))
    }
  } else if(type %in% rownames(eout$basis)) {
    # to calculate the loga of basis species at equilibrium
    if(!missing(groups)) stop("can't plot equilibrium activities of basis species for grouped species")
    if(isTRUE(alpha) | is.character(alpha)) stop("equilibrium activities of basis species not available with alpha=TRUE")
    plot.loga.basis <- TRUE
  } else if(type=="loga.equil" & !"loga.equil" %in% names(eout)) stop("'eout' is not the output from equil()") 
  else if(type!="loga.equil") stop(type, " is not a valid diagram type")

  ## consider a different number of species if we're grouping them together
  ngroups <- length(groups)

  ## keep the values we plot in one place so they can be modified, plotted and eventually returned
  # unless something happens below, we'll plot the loga.equil from equilibrate()
  plotvals <- eout$loga.equil
  plotvar <- "loga.equil"

  ## number of dimensions (T, P or chemical potentials that are varied)
  # length(eout$vars) - the number of variables = the maximum number of dimensions
  # length(dim(eout$values[[1]])) - nd=1 if it was a transect along multiple variables
  nd <- min(length(eout$vars), length(dim(eout$values[[1]])))

  ## deal with output from affinity()
  if(eout.is.aout) {
    # plot property from affinity(), divided by balancing coefficients
    plotvals <- lapply(1:length(eout$values), function(i) {
      # we divide by the balancing coefficients if we're working with affinities
      # this is not normalizing the formulas! it's balancing the reactions...
      # normalizing the formulas is done below
      eout$values[[i]] / n.balance[i]
    })
    plotvar <- eout$property
    # we change 'A' to 'A/(2.303RT)' so the axis label is made correctly
    # 20171027 use parentheses to avoid ambiguity about order of operations
    if(plotvar=="A") {
      plotvar <- "A/(2.303RT)"
      if(nd==2 & type=="auto") message("diagram: using maximum affinity method for 2-D diagram")
      else if(nd==2 & type=="saturation") message("diagram: plotting saturation lines for 2-D diagram")
      else message("diagram: plotting A/(2.303RT) / n.balance")
    } else message(paste("diagram: plotting", plotvar, " / n.balance"))
  }

  ## use molality instead of activity if the affinity calculation include ionic strength 20171101
  use.molality <- "IS" %in% names(eout)

  ## when can normalize and as.residue be used
  if(normalize | as.residue) {
    if(normalize & as.residue) stop("'normalize' and 'as.residue' can not both be TRUE")
    if(!eout.is.aout) stop("'normalize' or 'as.residue' can be TRUE only if 'eout' is the output from affinity()")
    if(nd!=2) stop("'normalize' or 'as.residue' can be TRUE only for a 2-D (predominance) diagram")
    if(normalize) message("diagram: using 'normalize' in calculation of predominant species")
    else message("diagram: using 'as.residue' in calculation of predominant species")
  }

  ## sum affinities or activities of species together in groups 20090524
  # using lapply/Reduce 20120927
  if(!missing(groups)) {
    # loop over the groups
    plotvals <- lapply(groups, function(ispecies) {
      # remove the logarithms
      if(eout.is.aout) act <- lapply(plotvals[ispecies], function(x) 10^x)
      # and, for activity, multiply by n.balance 20170207
      else act <- lapply(seq_along(ispecies), function(i) eout$n.balance[ispecies[i]] * 10^plotvals[[ispecies[i]]])
      # sum the activities
      return(Reduce("+", act))
    })
    # restore the logarithms
    plotvals <- lapply(plotvals, function(x) log10(x))
    # we also combine the balancing coefficients for calculations using affinity
    if(eout.is.aout) n.balance <- sapply(groups, function(ispecies) sum(n.balance[ispecies]))
  }

  ## calculate the equilibrium logarithm of activity of a basis species
  ## (such that affinities of formation reactions are zero)
  if(plot.loga.basis) {
    ibasis <- match(type, rownames(eout$basis))
    # the logarithm of activity used in the affinity calculation
    is.loga.basis <- can.be.numeric(eout$basis$logact[ibasis])
    if(!is.loga.basis) stop(paste("the logarithm of activity for basis species", type, "is not numeric - was a buffer selected?"))
    loga.basis <- as.numeric(eout$basis$logact[ibasis])
    # the reaction coefficients for this basis species
    nu.basis <- eout$species[, ibasis]
    # the logarithm of activity where affinity = 0
    plotvals <- lapply(1:length(eout$values), function(x) {
      # eout$values is a strange name for affinity ... should be named something like eout$affinity ...
      loga.basis - eout$values[[x]]/nu.basis[x]
    })
    plotvar <- type
  }

  ## alpha: plot fractional degree of formation
  # scale the activities to sum=1  ... 20091017
  # allow scaling by balancing component 20171008
  if(isTRUE(alpha) | is.character(alpha)) {
    # remove the logarithms
    act <- lapply(plotvals, function(x) 10^x)
    if(identical(alpha, "balance")) for(i in 1:length(act)) act[[i]] <- act[[i]] * eout$n.balance[i]
    # sum the activities
    sumact <- Reduce("+", act)
    # divide activities by the total
    alpha <- lapply(act, function(x) x/sumact)
    plotvals <- alpha
    plotvar <- "alpha"
  }

  ## identify predominant species
  predominant <- NA
  linesout <- NA
  if(plotvar %in% c("loga.equil", "alpha", "A/(2.303RT)") & type!="saturation") {
    pv <- plotvals
    # some additional steps for affinity values, but not for equilibrated activities
    if(eout.is.aout) {
      for(i in 1:length(pv)) {
        # TODO: see vignette for an explanation for how this is normalizing
        # the formulas in a predominance calculation
        if(normalize) pv[[i]] <- (pv[[i]] + eout$species$logact[i] / n.balance[i]) - log10(n.balance[i])
        else if(as.residue) pv[[i]] <- pv[[i]] + eout$species$logact[i] / n.balance[i]
      }
    }
    predominant <- which.pmax(pv)
    # clip plot to water stability region
    if(limit.water) {
      wl <- water.lines(eout, plot.it=FALSE)
      # proceed if water.lines produced calculations for this plot
      if(!identical(wl, NA)) {
        # for each x-point, find the y-values that are outside the water stability limits
        for(i in seq_along(wl$xpoints)) {
          ymin <- min(c(wl$y.oxidation[i], wl$y.reduction[i]))
          ymax <- max(c(wl$y.oxidation[i], wl$y.reduction[i]))
          if(!wl$swapped) {
            # the actual calculation
            iNA <- eout$vals[[2]] < ymin | eout$vals[[2]] > ymax
            # assign NA to the predominance matrix
            predominant[i, iNA] <- NA
          } else {
            # as above, but x- and y-axes are swapped
            iNA <- eout$vals[[1]] < ymin | eout$vals[[1]] > ymax
            predominant[iNA, i] <- NA
          }
        }
      }
    }
  }

  ## where we'll put extra output for predominance diagrams (namesx, namesy, inames)
  out2D <- list()

  ### now on to the plotting ###

  if(plot.it) {

    ### general plot parameters ###

    ## handle line type/width/color arguments
    if(is.null(lty)) lty <- 1:ngroups
    lty <- rep(lty, length.out=ngroups)
    lwd <- rep(lwd, length.out=ngroups)
    col <- rep(col, length.out=ngroups)
    col.names <- rep(col.names, length.out=ngroups)

    ## make up some names for lines/fields if they are missing
    is.pname <- FALSE
    onames <- names
    if(missing(names) | all(is.numeric(names))) {
      # properties of basis species or reactions?
      if(eout$property %in% c("G.basis", "logact.basis")) names <- rownames(eout$basis)
      else {
        if(!missing(groups)) {
          if(is.null(names(groups))) names <- paste("group", 1:length(groups), sep="")
          else names <- names(groups)
        }
        else names <- as.character(eout$species$name)
        # remove non-unique organism or protein names
        if(all(grepl("_", names))) {
          is.pname <- TRUE
          # everything before the underscore (the protein)
          pname <- gsub("_.*$", "", names)
          # everything after the underscore (the organism)
          oname <- gsub("^.*_", "", names)
          # if the pname or oname are all the same, use the other one as identifying name
          if(length(unique(pname))==1) names <- oname
          if(length(unique(oname))==1) names <- pname
        }
        # append state to distinguish ambiguous species names
        isdup <- names %in% names[duplicated(names)]
        if(any(isdup)) names[isdup] <- paste(names[isdup],
          " (", eout$species$state[isdup], ")", sep="")
      }
    }
    # numeric values indicate a subset 20181007
    if(all(is.numeric(onames))) {
      if(all(onames > 0)) names[-onames] <- ""
      else if(all(onames < 0)) names[-onames] <- ""
      else stop("numeric 'names' should be all positive or all negative")
    }

    ## apply formatting to chemical formulas 20170204
    if(all(grepl("_", names))) is.pname <- TRUE
    if(format.names & !is.pname) {
      exprnames <- as.expression(names)
      for(i in seq_along(exprnames)) {
        exprnames[[i]] <- expr.species(exprnames[[i]])
        if(bold) exprnames[[i]] <- substitute(bold(a), list(a=exprnames[[i]]))
        if(italic) exprnames[[i]] <- substitute(italic(a), list(a=exprnames[[i]]))
      }
      names <- exprnames
    }

    if(nd==0) {

      ### 0-D diagram - bar graph of properties of species or reactions
      # plot setup
      if(missing(ylab)) ylab <- axis.label(plotvar, units="", use.molality=use.molality)
      barplot(unlist(plotvals), names.arg=names, ylab=ylab, cex.names=cex.names, col=col, ...)
      if(!is.null(main)) title(main=main)

    } else if(nd==1) {

      ### 1-D diagram - lines for properties or chemical activities
      xvalues <- eout$vals[[1]]
      if(missing(xlim)) xlim <- range(xvalues)  # TODO: this is backward if the vals are not increasing
      # initialize the plot
      if(!add) {
        if(missing(xlab)) xlab <- axis.label(eout$vars[1], basis=eout$basis, use.molality=use.molality)
        if(missing(ylab)) ylab <- axis.label(plotvar, units="", use.molality=use.molality)
        # to get range for y-axis, use only those points that are in the xrange
        if(is.null(ylim)) {
          isx <- xvalues >= min(xlim) & xvalues <= max(xlim)
          xfun <- function(x) x[isx]
          myval <- sapply(plotvals, xfun)
          ylim <- extendrange(myval)
        }
        if(tplot) thermo.plot.new(xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, cex=cex, mar=mar, yline=yline, side=side, ...)
        else plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
      }
      # draw the lines
      spline.n <- 256 # the number of values at which to calculate splines
      if(is.null(spline.method) | length(xvalues) > spline.n) {
        for(i in 1:length(plotvals)) lines(xvalues, plotvals[[i]], col=col[i], lty=lty[i], lwd=lwd[i])
      } else {
        # plot splines instead of lines connecting the points 20171116
        spline.x <- seq(xlim[1], xlim[2], length.out=spline.n)
        for(i in 1:length(plotvals)) {
          spline.y <- splinefun(xvalues, plotvals[[i]], method=spline.method)(spline.x)
          lines(spline.x, spline.y, col=col[i], lty=lty[i], lwd=lwd[i])
        }
      }
      if(!add & !is.null(legend.x)) {
        # 20120521: use legend.x=NA to label lines rather than make legend
        if(is.na(legend.x)) {
          maxvals <- do.call(pmax, pv)
          dy <- rep(dy, length.out=length(plotvals))
          srt <- rep(srt, length.out=length(plotvals))
          # don't assign to adj becuase that messes up the missing test below
          alladj <- rep(adj, length.out=length(plotvals))
          for(i in 1:length(plotvals)) {
            # y-values for this line
            myvals <- as.numeric(plotvals[[i]])
            # don't take values that lie close to or above the top of plot
            myvals[myvals > ylim[1] + 0.95*diff(ylim)] <- ylim[1]
            # the starting x-adjustment
            thisadj <- alladj[i]
            # if this line has any of the overall maximum values, use only those values
            # (useful for labeling straight-line affinity comparisons 20170221)
            is.max <- myvals==maxvals
            if(any(is.max) & plotvar != "alpha") {
              # put labels on the median x-position
              imax <- median(which(is.max))
            } else {
              # put labels on the maximum of the line
              # (useful for labeling alpha plots)
              imax <- which.max(myvals)
              # try to avoid the sides of the plot; take care of reversed x-axis
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
            # also include y-offset (dy) and y-adjustment (labels bottom-aligned with the line)
            # .. and srt (string rotation) 20171127
            text(xvalues[imax], plotvals[[i]][imax] + dy[i], labels=names[i], adj=c(thisadj, 0), cex=cex.names, srt=srt[i], font=font, family=family)
          }
        } else legend(x=legend.x, lty=lty, legend=names, col=col, cex=cex.names, lwd=lwd, ...)
      }
      # add a title
      if(!is.null(main)) title(main=main)

    } else if(nd==2) {

      ### 2-D diagram - fields indicating species predominance, or contours for other properties

      ### functions for constructing predominance area diagrams
      ## color fill function
      fill.color <- function(xs, ys, out, fill, nspecies) {
        # handle min/max reversal
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
        # the z values
        zs <- out
        for(i in 1:nrow(zs)) zs[i,] <- out[nrow(zs)+1-i,]
        zs <- t(zs)
        breaks <- c(-1, 0, 1:nspecies) + 0.5
        # use fill.NA for NA values
        zs[is.na(zs)] <- 0
        image(x=xs, y=ys, z=zs, col=c(fill.NA, fill), add=TRUE, breaks=breaks, useRaster=TRUE)
      }
      ## curve plot function
      # 20091116 replaced plot.curve with plot.line; different
      # name, same functionality, *much* faster
      plot.line <- function(out, xlim, ylim, dotted, col, lwd, xrange) {
        # plot boundary lines between predominance fields
        vline <- function(out, ix) {
          ny <- nrow(out)
          xs <- rep(ix, ny*2+1)
          ys <- c(rep(ny:1, each=2), 0)
          y1 <- out[, ix]
          y2 <- out[, ix+1]
          # no line segment inside a stability field
          iy <- which(y1==y2)
          ys[iy*2] <- NA
          # no line segment at a dotted position
          iyd <- rowSums(sapply(dotted, function(y) ys%%y==0)) > 0
          ys[iyd] <- NA
          return(list(xs=xs, ys=ys))
        }
        hline <- function(out, iy) {
          nx <- ncol(out)
          ys <- rep(iy, nx*2+1)
          xs <- c(0, rep(1:nx, each=2))
          x1 <- out[iy, ]
          x2 <- out[iy+1, ]
          # no line segment inside a stability field
          ix <- which(x1==x2)
          xs[ix*2] <- NA
          # no line segment at a dotted position
          ixd <- rowSums(sapply(dotted, function(x) xs%%x==0)) > 0
          xs[ixd] <- NA
          return(list(xs=xs, ys=ys))
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
        # vertical lines
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
        lines(xs, ys, col=col, lwd=lwd)
        # horizontal lines
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
        lines(xs, ys, col=col, lwd=lwd)
      }
      ## new line plotting function 20170122
      contour.lines <- function(predominant, xlim, ylim, lty, col, lwd) {
        # the x and y values
        xs <- seq(xlim[1], xlim[2], length.out=dim(predominant)[1])
        ys <- seq(ylim[1], ylim[2], length.out=dim(predominant)[2])
        # reverse any axis that has decreasing values
        if(diff(xlim) < 0) {
          predominant <- predominant[nrow(predominant):1, ]
          xs <- rev(xs)
        }
        if(diff(ylim) < 0) {
          predominant <- predominant[, ncol(predominant):1]
          ys <- rev(ys)
        }
	# the categories (species/groups/etc) on the plot
	zvals <- na.omit(unique(as.vector(predominant)))
        # initialize list to output line coordinates
        linesout <- list()
        iout <- 1
	# take each possible pair of species
	for(i in 1:(length(zvals)-1)) {
	  for(j in (i+1):length(zvals)) {
	    z <- predominant
	    # draw contours only for this pair
	    z[!z %in% c(zvals[i], zvals[j])] <- NA
	    # give them neighboring values (so we get one contour line)
	    z[z==zvals[i]] <- 0
	    z[z==zvals[j]] <- 1
            # use contourLines() instead of contour() in order to get line coordinates 20181029
	    cLines <- contourLines(xs, ys, z, levels=0.5)
            if(length(cLines) > 0) {
              # loop in case contourLines returns multiple lines
              for(k in 1:length(cLines)) {
                # draw the lines
                lines(cLines[[k]][2:3], lty=lty, col=col, lwd=lwd)
                # keep the x and y values (list components 2 and 3)
                linesout[[iout]] <- cLines[[k]][[2]]
                names(linesout)[iout] <- paste0("x", k, "_", i, ".", j)
                linesout[[iout+1]] <- cLines[[k]][[3]]
                names(linesout)[iout+1] <- paste0("y", k, "_", i, ".", j)
                iout <- iout + 2
              }
            }
	  }
	}
        # https://stackoverflow.com/questions/34570860/adding-na-to-make-all-list-elements-equal-length 20181029
        lapply(linesout, `length<-`, max(lengths(linesout)))
      }
      ## label plot function
      # calculate coordinates for field labels
      plot.names <- function(out, xs, ys, names) {
        ll <- ngroups
        namesx <- numeric(ll); namesy <- numeric(ll); n <- numeric(ll)
        for(j in nrow(out):1) {
          # 20091116 for speed, loop over ngroups instead of k (columns)
          for(i in 1:ll) {
            k <- which(out[j,]==i)
            if(length(k)==0) next
            namesx[i] <- namesx[i] + sum(xs[k])
            namesy[i] <- namesy[i] + length(k)*ys[nrow(out)+1-j]
            n[i] <- n[i] + length(k)
          }
        }
        namesx <- namesx[n!=0]
        namesy <- namesy[n!=0]
        inames <- n!=0
        n <- n[n!=0]
        namesx <- namesx/n
        namesy <- namesy/n
        # plot field labels
        # the cex argument in this function specifies the character 
        # expansion of the labels relative to the current
        if(!is.null(names) & any(inames)) text(namesx, namesy, labels=names[inames], cex=cex.names, col=col.names[inames], font=font, family=family)
        return(list(namesx=namesx, namesy=namesy, inames=which(inames)))
      }

      ### done with predominance diagram functions
      ### now on to the diagram itself

      # colors to fill predominance fields
      # default to heat colors if we're on screen, or to transparent if we're adding to a plot
      if(missing(fill)) {
        if(add) fill <- "transparent"
        else if(any(grepl(names(dev.cur()), c("X11cairo", "quartz", "windows")))) fill <- "heat"
      }
      if(is.null(fill)) fill <- "transparent"
      else if(isTRUE(fill[1]=="rainbow")) fill <- rainbow(ngroups)
      else if(isTRUE(fill[1] %in% c("heat", "terrain", "topo", "cm"))) fill <- get(paste0(fill[1], ".colors"))(ngroups)
      fill <- rep(fill, length.out=ngroups)
      # modify the default for fill.NA
      if(add & missing(fill.NA)) fill.NA <- "transparent"
      # the x and y values 
      xs <- eout$vals[[1]]
      ys <- eout$vals[[2]]
      # the limits; they aren't necessarily increasing, so don't use range()
      xlim <- c(xs[1], tail(xs, 1))
      ylim <- c(ys[1], tail(ys, 1))
      # initialize the plot
      if(!add) {
        if(is.null(xlab)) xlab <- axis.label(eout$vars[1], basis=eout$basis, use.molality=use.molality)
        if(is.null(ylab)) ylab <- axis.label(eout$vars[2], basis=eout$basis, use.molality=use.molality)
        if(tplot) thermo.plot.new(xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
          cex=cex, cex.axis=cex.axis, mar=mar, yline=yline, side=side, ...)
        else plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
        # add a title
        if(!is.null(main)) title(main=main)
      }
      # colors and curves (predominance), or contours (properties)
      if(identical(predominant, NA)) {
        if(type=="saturation") {
          # for saturation plot, contour affinity=0 for all species
          for(i in 1:length(plotvals)) {
            zs <- plotvals[[i]]
            # skip plotting if this species has no possible saturation line, or a line outside the plot range
            if(length(unique(as.numeric(zs)))==1) {
              message("diagram: no saturation line possible for ", names[i])
              next
            }
            if(all(zs < 0) | all(zs > 0)) {
              message("diagram: beyond range for saturation line of ", names[i])
              next
            }
            if(identical(contour.method, NULL) | identical(contour.method, NA) | identical(contour.method, ""))
              contour(xs, ys, zs, add=TRUE, col=col, lty=lty, lwd=lwd, labcex=cex, levels=0, labels=names[i], drawlabels=FALSE)
            else contour(xs, ys, zs, add=TRUE, col=col, lty=lty, lwd=lwd, labcex=cex, levels=0, labels=names[i], method=contour.method)
          }
        } else {
          # otherwise, make contours of properties using first species only
          if(length(plotvals) > 1) warning("showing only first species in 2-D property diagram")
          zs <- plotvals[[1]]
          contour(xs, ys, zs, add=TRUE, col=col, lty=lty, lwd=lwd, labcex=cex, method=contour.method)
        }
        pn <- list(namesx=NULL, namesy=NULL, inames=NULL)
      } else {
        # put predominance matrix in the right order for image() etc
        zs <- t(predominant[, ncol(predominant):1])
        if(!is.null(fill)) fill.color(xs, ys, zs, fill, ngroups)
        pn <- plot.names(zs, xs, ys, names)
        # only draw the lines if there is more than one field  20180923
        # (to avoid warnings from contour, which seem to be associated with weird
        # font metric state and subsequent errors adding e.g. subscripted text to plot)
        if(length(na.omit(unique(as.vector(zs)))) > 1) {
          if(!is.null(dotted)) plot.line(zs, xlim, ylim, dotted, col, lwd, xrange=xrange)
          else linesout <- contour.lines(predominant, xlim, ylim, lty=lty, col=col, lwd=lwd)
        }
        # re-draw the tick marks and axis lines in case the fill obscured them
        if(tplot & !identical(fill, "transparent")) thermo.axis()
      } # done with the 2D plot!
      out2D <- list(namesx=pn$namesx, namesy=pn$namesy, inames=pn$inames)
    } # end if(nd==2)
  } # end if(plot.it)
  out <- c(eout, list(plotvar=plotvar, plotvals=plotvals, names=names, predominant=predominant))
  if(!identical(predominant, NA) & is.null(dotted)) out <- c(out, list(lines=linesout))
  out <- c(out, out2D)
  return(invisible(out))
}


strip <- function(affinity, ispecies=NULL, col=NULL, ns=NULL,
  xticks=NULL, ymin=-0.2, xpad=1, cex.names=0.7) {
  # make strip chart(s) showing the degrees of formation
  # of species as color bars of varying widths
  # extracted from bison/plot.R 20091102 jmd
  # figure out defaults
  a <- affinity
  xlim <- range(a$vals[[1]])
  xlab <- axis.label(a$vars[1])
  if(is.null(ispecies)) ispecies <- list(1:nrow(a$species))
  if(!is.list(ispecies)) ispecies <- list(ispecies)
  if(!is.null(ns) & !is.list(ns)) ns <- list(ns)
  if(is.null(col)) col <- rainbow(length(ispecies[[1]]))
  # how many strip charts on this plot
  # determined by the length of the ispecies list
  nstrip <- length(ispecies)
  # start up the plot
  plot.xlim <- c(xlim[1]-xpad,xlim[2]+xpad)
  ymax <- nstrip+0.3
  thermo.plot.new(xlim=plot.xlim,ylim=c(ymin,ymax),xlab=xlab,ylab="",
      side=c(1,3),mar=par('mar'),do.box=FALSE)
  if(!is.null(xticks)) {
    # mark the positions of the sites on the x-axis
    for(i in 1:5) lines(rep(xticks[i],2),c(ymin,ymin+0.1),lwd=6,col=col[i])
    for(i in 1:5) lines(rep(xticks[i],2),c(ymax,ymax-0.1),lwd=6,col=col[i])
  }
  for(j in 1:nstrip) {
    # get the degrees of formation
    e <- equilibrate(a, normalize=TRUE, ispecies=ispecies[[j]])
    # depict the relative stabilities of the proteins as color bars
    # make vertical color bars sizes proportional to activities
    xs <- seq(xlim[1],xlim[2],length.out=length(e$loga.equil[[1]]))
    # total height of the stack
    ly <- 0.5  
    # where to start plotting on the y-axis
    y0 <- j-ly
    for(i in 1:length(e$loga.equil[[1]])) {
      # create a vector of activities
      loga <- numeric()
      for(k in 1:length(ispecies[[j]])) loga <- c(loga, e$loga.equil[[k]][i])
      # turn the log activities into alphas
      act <- 10^loga
      alpha <- act/sum(act)
      alpha.order <- order(alpha, decreasing=FALSE)
      # plot the bars, least abundant first 
      dy <- y0    # keep track of the height of the stack
      for(k in 1:length(ispecies[[j]])) {
        y1 <- dy
        y2 <- y1 + alpha[alpha.order[k]] * ly
        dy <- y2
        # note: lwd should be lowered here for very high-resolution plots
        lines(c(xs[i],xs[i]),c(y1,y2),col=col[alpha.order[k]],lwd=3,lend="butt")
      }
    }
    # label the color bar
    text(xlim[1],j-ly*1.4,names(ispecies[j]),adj=0,cex=cex.names)
    # add inset plot showing the relative numbers of species
    if(!is.null(ns)) {
      ys1 <- y0 - 0.85*ly
      ys2 <- y0 - 0.15*ly
      xmax <- xlim[2]
      xmin <- xlim[1] + 17/22*diff(xlim)
      xss <- seq(xmin,xmax,length.out=length(ispecies[[j]]))
      yss <- numeric()
      for(i in 1:length(ispecies[[j]])) yss <- c(yss,ys1+(ys2-ys1)*(ns[[j]][i]/max(ns[[j]])))
      points(xss,yss,pch=20,cex=0.5)
      lines(xss,yss)
    }
  }
}

find.tp <- function(x) {
  # find triple points in an matrix of integers  20120525 jmd
  # these are the locations closest to the greatest number of different values
  # rearrange the matrix in the same way that diagram() does for 2-D predominance diagrams
  x <- t(x[, ncol(x):1])
  # all of the indexes for the matrix
  id <- which(x > 0, arr.ind=TRUE)
  # we'll do a brute-force count at each position
  n <- sapply(1:nrow(id), function(i) {
    # row and column range to look at (3x3 except at edges)
    r1 <- max(id[i, 1]-1, 0)
    r2 <- min(id[i, 1]+1, nrow(x))
    c1 <- max(id[i, 2]-1, 0)
    c2 <- min(id[i, 2]+1, ncol(x))
    # the number of unique values
    return(length(unique(as.numeric(x[r1:r2, c1:c2]))))
  })
  # now which positions have the most counts?
  imax <- which(n==max(n))
  # return the indices
  return(id[imax, ])
}

