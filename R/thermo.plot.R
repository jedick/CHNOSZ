# CHNOSZ/thermo.plot.R
# Functions to create and modify plots

# Start a new plot with some customized settings
thermo.plot.new <- function(xlim,ylim,xlab,ylab,cex = par('cex'),mar = NULL,lwd = par('lwd'),side = c(1,2,3,4),
  mgp = c(1.7,0.3,0),cex.axis = par('cex'),col = par('col'),yline = NULL,axs = 'i',plot.box = TRUE,
  las = 1,xline = NULL, grid = "", col.grid = "gray", ...) {
  thermo <- get("thermo", CHNOSZ)
  # 20120523 store the old par in thermo()$opar
  if(is.null(thermo$opar)) {
    thermo$opar <- par(no.readonly = TRUE)
    assign("thermo", thermo, CHNOSZ)
  }
  # 20090324 mar handling: NULL - a default setting; NA - par's setting
  # 20090413 changed mar of top side from 2 to 2.5
  marval <- c(3, 3.5, 2.5, 1)
  if(identical(mar[1], NA)) marval <- par("mar")
  # 20181007 get mar from the current device (if it exists) and par("mar") is not the default
  if(!is.null(dev.list())) {
    if(!identical(par("mar"), c(5.1, 4.1, 4.1, 2.1))) marval <- par("mar")
  }
  # Assign marval to mar if the latter is NULL or NA
  if(!is.numeric(mar)) mar <- marval
  par(mar = mar,mgp = mgp,tcl = 0.3,las = las,xaxs = axs,yaxs = axs,cex = cex,lwd = lwd,col = col,fg = col, ...)
  plot.new()
  plot.window(xlim = xlim,ylim = ylim)
  if(plot.box) box()
  # Labels
  if(is.null(xline)) xline <- mgp[1]
  thermo.axis(xlab,side = 1,line = xline,cex = cex.axis,lwd = NULL)
  if(is.null(yline)) yline <- mgp[1]
  thermo.axis(ylab,side = 2,line = yline,cex = cex.axis,lwd = NULL)
  # (optional) tick marks
  if(1 %in% side) thermo.axis(NULL,side = 1,lwd = lwd, grid = grid, col.grid = col.grid, plot.line = !plot.box)
  if(2 %in% side) thermo.axis(NULL,side = 2,lwd = lwd, grid = grid, col.grid = col.grid, plot.line = !plot.box)
  if(3 %in% side) thermo.axis(NULL,side = 3,lwd = lwd, plot.line = !plot.box)
  if(4 %in% side) thermo.axis(NULL,side = 4,lwd = lwd, plot.line = !plot.box)
}

# Function to add axes and axis labels to plots,
#   with some default style settings (rotation of numeric labels)
# With the default arguments (no labels specified), it plots only the axis lines and tick marks
#   (used by diagram() for overplotting the axis on diagrams filled with colors).
thermo.axis <- function(lab = NULL, side = 1:4,line = 1.5, cex = par('cex'), lwd = par('lwd'),
  col = par('col'), grid  =  "", col.grid = "gray", plot.line = FALSE) {
  if(!is.null(lwd)) {
    for(thisside in side) {

      ## Get the positions of major tick marks
      at <- axis(thisside,labels = FALSE,tick = FALSE) 
      # Get nicer divisions for axes that span exactly 15 units 20200719
      if(thisside %in% c(1,3)) lim <- par("usr")[1:2]
      if(thisside %in% c(2,4)) lim <- par("usr")[3:4]
      if(abs(diff(lim)) == 15) at <- seq(lim[1], lim[2], length.out = 6)
      if(abs(diff(lim)) == 1.5) at <- seq(lim[1], lim[2], length.out = 4)
      # Make grid lines
      if(grid %in% c("major", "both") & thisside == 1) abline(v = at, col = col.grid)
      if(grid %in% c("major", "both") & thisside == 2) abline(h = at, col = col.grid)
      ## Plot major tick marks and numeric labels
      do.label <- TRUE
      if(missing(side) | (missing(cex) & thisside %in% c(3,4))) do.label <- FALSE
      # col and col.ticks: plot the tick marks but no line (we make it with box() in thermo.plot.new()) 20190416
      # mat: don't plot ticks at the plot limits 20190416
      if(thisside %in% c(1, 3)) pat <- par("usr")[1:2]
      if(thisside %in% c(2, 4)) pat <- par("usr")[3:4]
      mat <- setdiff(at, pat)
      if(plot.line) axis(thisside, at = mat, labels = FALSE, tick = TRUE, lwd = lwd, col.axis = col, col = col)
      else axis(thisside, at = mat, labels = FALSE, tick = TRUE, lwd = lwd, col.axis = col, col = NA, col.ticks = col)
      # Plot only the labels at all major tick points (including plot limits) 20190417
      if(do.label) axis(thisside, at = at, tick = FALSE, col = col)

      ## Plot minor tick marks
      # The distance between major tick marks
      da <- abs(diff(at[1:2]))
      # Distance between minor tick marks
      di <- da / 4
      if(!da %% 3) di <- da / 3
      else if(da %% 2 | !(da %% 10)) di <- da / 5
      # Number of minor tick marks
      if(thisside %in% c(1,3)) {
        ii <- c(1,2) 
        myasp <- par('xaxp')
      } else {
        ii <- c(3,4)
        myasp <- par('yaxp')
      }
      myusr <- par('usr')[ii]
      daxis <- abs(diff(myusr))
      nt <- daxis / di + 1
      ## If nt isn't an integer, it probably
      ## means the axis limits don't correspond
      ## to major tick marks (expect problems)
      ##at <- seq(myusr[1],myusr[2],length.out = nt)
      # Start from (bottom/left) of axis?
      bl <- 1
      #if(myasp[2] == myusr[2]) bl <- 2
      # Is forward direction (top/right)?
      tr <- 1
      if(xor(myusr[2] < myusr[1] , bl == 2)) tr <- -1
      #at <- myusr[bl] + tr * di * seq(0:(nt-1))
      # Well all of that doesn't work in a lot of cases,
      # where none of the axis limits correspond to
      # major tick marks. perhaps the following will work
      at <- myusr[1] + tr * di * (0:(nt-1))
      # Apply an offset
      axt <- axTicks(thisside)[1]
      daxt <- (axt - myusr[1])/di
      daxt <- (daxt-round(daxt))*di
      at <- at + daxt
      ## Get the positions of major tick marks and make grid lines
      if(grid %in% c("minor", "both") & thisside == 1) abline(v = at, col=col.grid, lty = 3)
      if(grid %in% c("minor", "both") & thisside == 2) abline(h = at, col=col.grid, lty = 3)
      tcl <- par('tcl') * 0.5
      at <- setdiff(at, pat)
      if(plot.line) axis(thisside, labels = FALSE, tick = TRUE, lwd = lwd, col.axis = col, at = at, tcl = tcl, col = col)
      else axis(thisside, labels = FALSE, tick = TRUE, lwd = lwd, col.axis = col, at = at, tcl = tcl, col = NA, col.ticks = col)
    }
  }
  # Rotate labels on side axes
  for(thisside in side) {
    if(thisside %in% c(2,4)) las <- 0 else las <- 1
    if(!is.null(lab)) mtext(lab, side = thisside, line = line, cex = cex, las = las)
  }
}

