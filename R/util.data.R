# CHNOSZ/util.data.R
# Check entries in the thermodynamic database

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.formula.R")
#source("util.data.R")
#source("util.character.R")

thermo.refs <- function(key=NULL, keep.duplicates=FALSE) {
  ## Return references for thermodynamic data.
  ## 20110615 browse.refs() first version
  ## 20170212 thermo.refs() remove browsing (except for table of all sources)
  # 'key' can be
  # NULL: show a table of all sources in a browser
  # character: return data for each listed source key
  # numeric: open one or two web pages for each listed species
  # list: the output of subcrt()
  ## First retrieve the sources table
  thermo <- get("thermo", CHNOSZ)
  x <- thermo$refs[order(thermo$refs$note), ]
  ## Show a table in the browser if 'key' is NULL 
  if(is.null(key)) {
    # Create the html links
    cite <- x$citation
    x$citation <- sprintf("<a href='%s' target='_blank'>%s</a>", x$URL, cite)
    notlinked <- x$URL=="" | is.na(x$URL)
    x$citation[notlinked] <- cite[notlinked]
    # Remove the last (URL) component
    #x$URL <- NULL
    x <- x[1:5]
    # Count the number of times each source is cited in thermo()$OBIGT
    # e.g. if key is "Kel60" we match "Kel60 [S92]" but not "Kel60.1 [S92]"
    # http://stackoverflow.com/questions/6713310/how-to-specify-space-or-end-of-string-and-space-or-start-of-string
    # We also have to escape keys with "+" signs
    ns1 <- sapply(x$key, function(x) sum(grepl(gsub("+", "\\+", paste0(x, "($|\\s)"), fixed=TRUE), thermo$OBIGT$ref1)) )
    ns2 <- sapply(x$key, function(x) sum(grepl(gsub("+", "\\+", paste0(x, "($|\\s)"), fixed=TRUE), thermo$OBIGT$ref2)) )
    number <- ns1 + ns2
    number[number==0] <- ""
    # Now that we're using the sortTable() from w3schools.com, numbers are sorted like text
    # Add leading zeros to make the numbers sortable 20170317
    # (the zeros disappear somewhere in the rendering of the page)
    number <- formatC(number, width = 3, format = "d", flag = "0")
    # append the counts to the table to be shown
    x <- c(list(number=number), x)
    # Title to display for web page
    title <- "References for thermodynamic data in CHNOSZ"
    ### The following is adapted from print.findFn in package 'sos'
    f0 <- tempfile()
    File <- paste(f0, ".html", sep="")
    #Dir <- dirname(File)
    #js <- system.file("extdata/js", "sorttable.js", package = "CHNOSZ")
    #file.copy(js, Dir)
    ## Sundar's original construction:
    con <- file(File, "wt")
    on.exit(close(con))
    .cat <- function(...)
      cat(..., "\n", sep = "", file = con, append = TRUE)
    ## Start
    cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
                "http://www.w3.org/TR/html4/strict.dtd">\n', file = con)
    .cat("<html>")
    .cat("<head>")
    .cat("<title>", title, "</title>")
    # sorttable.js is "Blocked for security reasons" in Gmail 20170317
    #.cat("<script src=sorttable.js type='text/javascript'></script>")
    # https://www.w3schools.com/howto/howto_js_sort_table.asp
    .cat('<script type="text/javascript">
	  function sortTable(n) {
	    var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
	    table = document.getElementById("thermorefs");
	    switching = true;
	    dir = "asc";
	    while (switching) {
	      switching = false;
	      rows = table.getElementsByTagName("TR");
	      for (i = 1; i < (rows.length - 1); i++) {
		shouldSwitch = false;
		x = rows[i].getElementsByTagName("TD")[n];
		y = rows[i + 1].getElementsByTagName("TD")[n];
		if (dir == "asc") {
		  if (x.innerHTML > y.innerHTML) {
		    shouldSwitch= true;
		    break;
		  }
		} else if (dir == "desc") {
		  if (x.innerHTML < y.innerHTML) {
		    shouldSwitch= true;
		    break;
		  }
		}
	      }
	      if (shouldSwitch) {
		rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
		switching = true;
		switchcount ++;
	      } else {
		if (switchcount == 0 && dir == "asc") {
		  dir = "desc";
		  switching = true;
		}
	      }
	    }
	  }
</script>')
    .cat("</head>")
    ### Boilerplate text
    .cat("<body>")
    .cat(paste0('<h1>References for thermodynamic data in <a href="https://chnosz.net"><font color="red">CHNOSZ</font></a> ',
               packageDescription("CHNOSZ")$Version), " (", packageDescription("CHNOSZ")$Date, ')</h1>')
    .cat("<h3>Click on a column header to sort, or on a citation to open the URL in new window.</h3>")
    .cat("<h4>Column 'number' gives the number of times each reference appears in thermo()$OBIGT.</h4>")
    .cat('<p>See also the vignette <a href="https://chnosz.net/vignettes/OBIGT.html">OBIGT thermodynamic database</a>.</p>')
    ### Start table and headers
    .cat("<table id='thermorefs' border='1'>")
    .cat("<tr>")
    #.cat(sprintf("  <th>%s</th>\n</tr>",
    #             paste(names(x), collapse = "</th>\n  <th>")))
    for(i in 1:length(x)) .cat(sprintf('  <th onclick="sortTable(%s)">%s</th>', i-1, names(x)[i]))
    .cat("</tr>")
    ### Now preparing the body of the table
    paste.list <- c(lapply(x, as.character), sep = "</td>\n  <td>")
    tbody.list <- do.call("paste", paste.list)
    tbody <- sprintf("<tr>\n  <td>%s</td>\n</tr>", tbody.list)
    tbody <- sub("<td><a", "<td class=link><a", tbody, useBytes = TRUE)
    .cat(tbody)
    ### Finish it!
    .cat("</table></body></html>")
    ### End adaptation from print.findFn
    # Show table in browser
    browseURL(File)
    cat("thermo.refs: table of references is shown in browser\n")
  } else if(is.character(key)) {
    # Return citation information for the given source(s)
    # We omit the [S92] in "HDNB78 [S92]" etc.
    key <- gsub("\ .*", "", key)
    ix <- match(key, x$key)
    ina <- is.na(ix)
    if(any(is.na(ix))) message(paste("thermo.refs: reference key(s)",
      paste(key[ina], collapse = ","), "not found"))
    return(x[ix, ])
  } else if(is.numeric(key)) {
    # Get the source keys for the indicated species
    sinfo <- suppressMessages(info(key, check.it = FALSE))
    if(keep.duplicates) {
      # Output a single reference for each species 20180927
      # (including duplicated references, and not including ref2)
      mysources <- sinfo$ref1
    } else {
      mysources <- unique(c(sinfo$ref1, sinfo$ref2))
      mysources <- mysources[!is.na(mysources)]
    }
    return(thermo.refs(mysources))
  } else if(is.list(key)) {
    if("species" %in% names(key)) ispecies <- key$species$ispecies
    else if("reaction" %in% names(key)) ispecies <- key$reaction$ispecies
    else stop("list does not appear to be a result from subcrt()")
    if(is.null(ispecies)) stop("list does not appear to be a result from subcrt()")
    return(thermo.refs(ispecies))
  }
}

check.EOS <- function(eos, model, prop, ret.diff = TRUE) {
  # Compare calculated properties from thermodynamic parameters with given (database) values.
  # Print message and return the calculated value if tolerance is exceeded
  # or NA if the difference is within the tolerance.
  # 20110808 jmd
  thermo <- get("thermo", CHNOSZ)
  # Get calculated value based on EOS
  Theta <- 228  # K
  if(model %in% c("HKF", "DEW")) {
    # Run checks for aqueous species
    if(prop=="Cp") {
      ## Value of X consistent with IAPWS95
      #X <- -2.773788E-7
      # We use the value of X consistent with SUPCRT
      X <- -3.055586E-7
      refval <- eos$Cp
      calcval <- eos$c1 + eos$c2/(298.15-Theta)^2 + eos$omega*298.15*X
      tol <- thermo$opt$Cp.tol
      units <- paste(eos$E_units, "K-1 mol-1")
    } else if(prop=="V") {
      ## Value of Q consistent with IAPWS95
      #Q <- 0.00002483137
      # Value of Q consistent with SUPCRT92
      Q <- 0.00002775729
      refval <- eos$V
      calcval <- 41.84*eos$a1 + 41.84*eos$a2/2601 + 
        (41.84*eos$a3 + 41.84*eos$a4/2601) / (298.15-Theta) - Q * eos$omega
      isJoules <- eos$E_units == "J"
      if(any(isJoules)) calcval[isJoules] <- convert(calcval[isJoules], "cal")
      tol <- thermo$opt$V.tol
      units <- "cm3 mol-1"
    }
  } else {
    # Run checks for non-aqueous species (i.e., CGL)
    if(prop=="Cp") {
      refval <- eos$Cp
      Tr <- 298.15
      calcval <- eos$a + eos$b*Tr + eos$c*Tr^-2 + eos$d*Tr^-0.5 + eos$e*Tr^2 + eos$f*Tr^eos$lambda
      tol <- thermo$opt$Cp.tol
      units <- paste(eos$E_units, "K-1 mol-1")
    }
  }
  # Calculate the difference
  diff <- calcval - refval
  if(ret.diff) return(diff)
  else {
    # Return the calculated value if the difference is greater than tol
    if(!is.na(calcval)) {
      if(!is.na(refval)) {
        if(abs(diff) > tol) {
          message(paste("check.EOS: calculated ", prop, " of ", eos$name, "(", eos$state,
            ") differs by ", round(diff,2), " ", units, " from database value", sep=""))
          return(calcval)
        }
      } else return(calcval)
    }
  }
  # Return NA in most cases
  return(NA)
}

check.GHS <- function(ghs, ret.diff = TRUE) {
  # Compare G calculated from H and S with given (database) values
  # Print message and return the calculated value if tolerance is exceeded
  # or NA if the difference is within the tolerance
  # 20110808 jmd
  thermo <- get("thermo", CHNOSZ)
  # Get calculated value based on H and S
  Se <- entropy(as.character(ghs$formula))
  iscalories <- ghs$E_units == "cal"
  if(any(iscalories)) Se[iscalories] <- convert(Se[iscalories], "cal")
  refval <- ghs$G
  DH <- ghs$H
  S <- ghs$S
  Tr <- 298.15
  calcval <- DH - Tr * (S - Se)
  # Now on to the comparison
  # Calculate the difference
  diff <- calcval - refval
  if(ret.diff) return(diff)
  else if(!is.na(calcval)) {
    if(!is.na(refval)) {
      diff <- calcval - refval
      if(abs(diff) > thermo$opt$G.tol) {
        message(paste("check.GHS: calculated G of ", ghs$name, "(", ghs$state,
          ") differs by ", round(diff), " ", ghs$E_units, " mol-1 from database value", sep=""))
        return(calcval)
      }
    } else return(calcval)
  } else {
    # Calculating a value of G failed, perhaps because of missing elements
    return(NULL)
  }
  # Return NA in most cases
  return(NA)
}

check.OBIGT <- function() {
  # Function to check self-consistency between
  # values of Cp and V vs. EOS parameters
  # and among G, H, S values
  # 20110808 jmd replaces 'check=TRUE' argument of info()
  checkfun <- function(what) {
    message(paste("check.OBIGT: checking", what))
    # Looking at thermo$OBIGT
    if(what=="OBIGT") tdata <- get("thermo", CHNOSZ)$OBIGT
    else if(what=="DEW") tdata <- read.csv(system.file("extdata/OBIGT/DEW.csv", package = "CHNOSZ"), as.is = TRUE)
    else if(what=="SLOP98") tdata <- read.csv(system.file("extdata/OBIGT/SLOP98.csv", package = "CHNOSZ"), as.is = TRUE)
    else if(what=="SUPCRT92") tdata <- read.csv(system.file("extdata/OBIGT/SUPCRT92.csv", package = "CHNOSZ"), as.is = TRUE)
    else if(what=="AS04") tdata <- read.csv(system.file("extdata/OBIGT/AS04.csv", package = "CHNOSZ"), as.is = TRUE)
    else if(what=="AD") tdata <- read.csv(system.file("extdata/OBIGT/AD.csv", package = "CHNOSZ"), as.is = TRUE)
    else if(what=="GEMSFIT") tdata <- read.csv(system.file("extdata/OBIGT/GEMSFIT.csv", package = "CHNOSZ"), as.is = TRUE)
    ntot <- nrow(tdata)
    # Where to keep the results
    DCp <- DV <- DG <- rep(NA,ntot)
    # First get the species that use HKF equations
    isHKF <- tdata$model %in% c("HKF", "DEW")
    if(any(isHKF)) {
      eos.HKF <- OBIGT2eos(tdata[isHKF,], "HKF")
      DCp.HKF <- check.EOS(eos.HKF, "HKF", "Cp")
      DV.HKF <- check.EOS(eos.HKF, "HKF", "V")
      cat(paste("check.OBIGT: GHS for", sum(isHKF), "species with HKF model in", what, "\n"))
      DG.HKF <- check.GHS(eos.HKF)
      # Store the results
      DCp[isHKF] <- DCp.HKF
      DV[isHKF] <- DV.HKF
      DG[isHKF] <- DG.HKF
    }
    # Then other species, if they are present
    if(sum(!isHKF) > 0) {
      eos.cgl <- OBIGT2eos(tdata[!isHKF,], "cgl")
      DCp.cgl <- check.EOS(eos.cgl, "CGL", "Cp")
      cat(paste("check.OBIGT: GHS for", sum(!isHKF), "cr,gas,liq species in", what, "\n"))
      DG.cgl <- check.GHS(eos.cgl)
      DCp[!isHKF] <- DCp.cgl
      DG[!isHKF] <- DG.cgl
    }
    # Put it all together
    out <- data.frame(table = what, ispecies = 1:ntot, name = tdata$name, state = tdata$state, E_units = tdata$E_units, DCp = DCp, DV = DV, DG = DG)
    return(out)
  }
  # Check default database (OBIGT)
  out <- checkfun("OBIGT")
  # Check optional data
  out <- rbind(out, checkfun("DEW"))
  out <- rbind(out, checkfun("SLOP98"))
  out <- rbind(out, checkfun("SUPCRT92"))
  out <- rbind(out, checkfun("AS04"))
  out <- rbind(out, checkfun("GEMSFIT"))
  # Set differences within a tolerance to NA
  out$DCp[abs(out$DCp) < 1] <- NA
  out$DV[abs(out$DV) < 1] <- NA
  out$DG[abs(out$DG) < 500] <- NA
  # Take out species where all reported differences are NA
  ina <- is.na(out$DCp) & is.na(out$DV) & is.na(out$DG)
  out <- out[!ina,]
  # Round the values
  out$DCp <- round(out$DCp,2)
  out$DV <- round(out$DV,2)
  out$DG <- round(out$DG)
  # Return the results
  return(out)
}

RH2OBIGT <- function(compound=NULL, state="cr", file=system.file("extdata/adds/RH98_Table15.csv", package="CHNOSZ")) {
  # Get thermodynamic properties and equations of state parameters using 
  #   group contributions from Richard and Helgeson, 1998   20120609 jmd
  # Read the compound names, physical states, chemical formulas and group stoichiometry from the file
  # We use check.names=FALSE because the column names are the names of the groups,
  #   and are not syntactically valid R names, and stringsAsFactors=FALSE
  #   so that formulas are read as characters (for checking with as.chemical.formula)
  dat <- read.csv(file, check.names=FALSE, stringsAsFactors=FALSE)
  # "compound" the compound names and states from the file
  comate.arg <- comate.dat <- paste(dat$compound, "(", dat$state, ")", sep="")
  # "compound" the compound names and states from the arguments
  if(!is.null(compound)) comate.arg <- paste(compound, "(", state, ")", sep="")
  # Identify the compounds
  icomp <- match(comate.arg, comate.dat)
  # Check if all compounds were found
  ina <- is.na(icomp)
  if(any(ina)) stop(paste("compound(s)", paste(comate.arg[ina], collapse=" "), "not found in", file))
  # Initialize output data frame
  out <- get("thermo", CHNOSZ)$OBIGT[0, ]
  # Loop over the compounds
  for(i in icomp) {
    # The group stoichiometry for this compound
    thisdat <- dat[i, ]
    # Take out groups that are NA or 0
    thisdat <- thisdat[, !is.na(thisdat)]
    thisdat <- thisdat[, thisdat!=0]
    # Identify the groups in this compound
    igroup <- 4:ncol(thisdat)
    ispecies <- info(colnames(thisdat)[igroup], state=thisdat$state)
    # Check if all groups were found
    ina <- is.na(ispecies)
    if(any(ina)) stop(paste("group(s)", paste(colnames(thisdat)[igroup][ina], collapse=" "), "not found in", thisdat$state, "state"))
    # Group additivity of properties and parameters: add contributions from all groups
    thiseos <- t(colSums(get("thermo", CHNOSZ)$OBIGT[ispecies, 10:22] * as.numeric(thisdat[, igroup])))
    # Group additivity of chemical formula
    formula <- as.chemical.formula(colSums(i2A(ispecies) * as.numeric(thisdat[, igroup])))
    # Check if the formula is the same as in the file
    if(!identical(formula, thisdat$formula)) 
      stop(paste("formula", formula, "of", comate.dat[i], "(from groups) is not identical to", thisdat$formula, "(listed in file)" ))
    # Build the front part of OBIGT data frame
    thishead <- data.frame(name=thisdat$compound, abbrv=NA, formula=formula, state=thisdat$state, 
      ref1=NA, ref2=NA, date=as.character(Sys.Date()), E_units = "cal", stringsAsFactors=FALSE)
    # Insert the result into the output
    out <- rbind(out, cbind(thishead, thiseos))
  }
  return(out)
}

# Dump all thermodynamic data in default and optional OBIGT files 20171121
dumpdata <- function(file = NULL) {
  # The default database (OBIGT)
  dat <- get("thermo", CHNOSZ)$OBIGT
  OBIGT <- cbind(source = "OBIGT", dat)
  # Optional data
  dat <- read.csv(system.file("extdata/OBIGT/DEW.csv", package = "CHNOSZ"), as.is = TRUE)
  DEW <- cbind(source = "DEW", dat)
  dat <- read.csv(system.file("extdata/OBIGT/SLOP98.csv", package = "CHNOSZ"), as.is = TRUE)
  SLOP98 <- cbind(source = "SLOP98", dat)
  dat <- read.csv(system.file("extdata/OBIGT/SUPCRT92.csv", package = "CHNOSZ"), as.is = TRUE)
  SUPCRT92 <- cbind(source = "SUPCRT92", dat)
  # More optional data 20220929
  dat <- read.csv(system.file("extdata/OBIGT/AS04.csv", package = "CHNOSZ"), as.is = TRUE)
  AS04 <- cbind(source = "AS04", dat)
  dat <- read.csv(system.file("extdata/OBIGT/AD.csv", package = "CHNOSZ"), as.is = TRUE)
  AD <- cbind(source = "AD", dat)
  dat <- read.csv(system.file("extdata/OBIGT/GEMSFIT.csv", package = "CHNOSZ"), as.is = TRUE)
  GEMSFIT <- cbind(source = "GEMSFIT", dat)
  # Put it all together
  out <- rbind(OBIGT, SUPCRT92, SLOP98, AS04, AD, DEW, GEMSFIT)
  # Quote columns 2 (name) and 3 (abbrv) because they have commas for some entries
  if(!is.null(file)) write.csv(out, file, row.names = FALSE, quote = c(2, 3))
  else(return(out))
}

### Unexported functions ###

# Take a data frame in the format of thermo()$OBIGT of one or more rows,
#   remove scaling factors from equations-of-state parameters,
#   and apply new column names depending on the state.
# If fixGHS is TRUE a missing one of G, H or S for any species is calculated
#   from the other two and the chemical formula of the species.
# If toJoules is TRUE, convert parameters to Joules 20220325
# This function is used by both info() and subcrt() when retrieving entries from the thermodynamic database.
OBIGT2eos <- function(OBIGT, state, fixGHS = FALSE, toJoules = FALSE) {

  # Figure out the model for each species 20220929
  model <- OBIGT$model
  model[is.na(model)] <- ""
  isCGL <- model == "CGL"
  isHKF <- model == "HKF"
  isDEW <- model == "DEW"
  isAD <- model == "AD"
  # Remove scaling factors for the HKF and DEW species
  #   protect this by an if statement to workaround error in subassignment to empty subset of data frame in R < 3.6.0
  #   (https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17483) 20190302
  is_HKF <- isHKF | isDEW
  if(any(is_HKF)) OBIGT[is_HKF, 15:22] <- t(t(OBIGT[is_HKF, 15:22]) * 10 ^ c(-1, 2, 0, 4, 0, 4, 5, 0))
  # For AD species, set NA values in unused columns
  if(any(isAD)) OBIGT[isAD, 18:21] <- NA
  # Change column names depending on the model
  if(all(isAD)) colnames(OBIGT)[15:22] <- c("a", "b", "xi", "XX1", "XX2", "XX3", "XX4", "Z") 
  else if(all(isHKF | isAD | isDEW)) colnames(OBIGT)[15:22] <- c("a1", "a2", "a3", "a4", "c1", "c2", "omega", "Z") 
  else colnames(OBIGT)[15:22] <- c("a", "b", "c", "d", "e", "f", "lambda", "T")

  if(toJoules) {
    # Convert parameters from calories to Joules 20220325
    # [Was: convert parameters from Joules to calories 20190530]
    iscal <- OBIGT$E_units == "cal"
    if(any(iscal)) {
      OBIGT[iscal, c(10:13, 15:20)] <- convert(OBIGT[iscal, c(10:13, 15:20)], "J")
      # We only convert the last column for aqueous species (HKF parameter: omega), not for CGL species (arbitrary exponent: lambda)  20190903
      isaq <- OBIGT$state == "aq"
      if(any(isaq)) OBIGT[iscal & isaq, 21] <- convert(OBIGT[iscal & isaq, 21], "J")
      # Also update the E_units column 20220325
      OBIGT$E_units[iscal] <- "J"
    }
  }

  if(fixGHS) {
    # Fill in one of missing G, H, S;
    #   for use esp. by subcrt() because NA for one of G, H or S precludes calculations at high T
    # Which entries are missing just one
    imiss <- which(rowSums(is.na(OBIGT[, 10:12])) == 1)
    if(length(imiss) > 0) {
      for(i in 1:length(imiss)) {
        # Calculate the missing value from the others
        ii <- imiss[i]
        GHS <- as.numeric(GHS(as.character(OBIGT$formula[ii]), G = OBIGT[ii, 10], H = OBIGT[ii, 11], S = OBIGT[ii, 12],
                              E_units = ifelse(toJoules, "J", OBIGT$E_units[ii])))
        icol <- which(is.na(OBIGT[ii, 10:12]))
        OBIGT[ii, icol + 9] <- GHS[icol]
      }
    }
  }

  OBIGT
}
