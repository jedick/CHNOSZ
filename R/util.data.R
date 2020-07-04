# CHNOSZ/util.data.R
# check entries in the thermodynamic database

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.formula.R")
#source("util.data.R")
#source("util.character.R")

thermo.refs <- function(key=NULL, keep.duplicates=FALSE) {
  ## return references for thermodynamic data.
  ## 20110615 browse.refs() first version
  ## 20170212 thermo.refs() remove browsing (except for table of all sources)
  # 'key' can be
  # NULL: show a table of all sources in a browser
  # character: return data for each listed source key
  # numeric: open one or two web pages for each listed species
  # list: the output of subcrt()
  ## first retrieve the sources table
  thermo <- get("thermo", CHNOSZ)
  x <- thermo$refs[order(thermo$refs$note), ]
  ## show a table in the browser if 'key' is NULL 
  if(is.null(key)) {
    # create the html links
    cite <- x$citation
    x$citation <- sprintf("<a href='%s' target='_blank'>%s</a>", x$URL, cite)
    notlinked <- x$URL=="" | is.na(x$URL)
    x$citation[notlinked] <- cite[notlinked]
    # remove the last (URL) component
    #x$URL <- NULL
    x <- x[1:5]
    # count the number of times each source is cited in thermo$OBIGT
    # e.g. if key is "Kel60" we match "Kel60 [S92]" but not "Kel60.1 [S92]"
    # http://stackoverflow.com/questions/6713310/how-to-specify-space-or-end-of-string-and-space-or-start-of-string
    # we also have to escape keys with "+" signs
    ns1 <- sapply(x$key, function(x) sum(grepl(gsub("+", "\\+", paste0(x, "($|\\s)"), fixed=TRUE), thermo$OBIGT$ref1)) )
    ns2 <- sapply(x$key, function(x) sum(grepl(gsub("+", "\\+", paste0(x, "($|\\s)"), fixed=TRUE), thermo$OBIGT$ref2)) )
    number <- ns1 + ns2
    number[number==0] <- ""
    # now that we're using the sortTable() from w3schools.com, numbers are sorted like text
    # add leading zeros to make the numbers sortable 20170317
    # (the zeros disappear somewhere in the rendering of the page)
    number <- formatC(number, width = 3, format = "d", flag = "0")
    # append the counts to the table to be shown
    x <- c(list(number=number), x)
    # title to display for web page
    title <- "References for thermodynamic data in CHNOSZ"
    ### the following is adapted from print.findFn in package 'sos'
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
    ## start
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
    ### boilerplate text
    .cat("<body>")
    .cat('<h1>References for thermodynamic data in <a href="http://chnosz.net"><font color="red">CHNOSZ</font></a></h1>')
    .cat("<h3>Click on a column header to sort, or on a citation to open the URL in new window.</h3>")
    .cat("<h4>Column 'number' gives the number of times each reference appears in thermo$OBIGT.</h4>")
    .cat('<p>See also the vignette <a href="http://chnosz.net/vignettes/OBIGT.html">Thermodynamic data in CHNOSZ</a>.</p>')
    ### start table and headers
    .cat("<table id='thermorefs' border='1'>")
    .cat("<tr>")
    #.cat(sprintf("  <th>%s</th>\n</tr>",
    #             paste(names(x), collapse = "</th>\n  <th>")))
    for(i in 1:length(x)) .cat(sprintf('  <th onclick="sortTable(%s)">%s</th>', i-1, names(x)[i]))
    .cat("</tr>")
    ### now preparing the body of the table
    paste.list <- c(lapply(x, as.character), sep = "</td>\n  <td>")
    tbody.list <- do.call("paste", paste.list)
    tbody <- sprintf("<tr>\n  <td>%s</td>\n</tr>", tbody.list)
    tbody <- sub("<td><a", "<td class=link><a", tbody, useBytes = TRUE)
    .cat(tbody)
    ### finish it!
    .cat("</table></body></html>")
    ### end adaptation from print.findFn
    # show table in browser
    browseURL(File)
    cat("thermo.refs: table of references is shown in browser\n")
  } else if(is.character(key)) {
    # return citation information for the given source(s)
    # we omit the [S92] in "HDNB78 [S92]" etc.
    key <- gsub("\ .*", "", key)
    ix <- match(key, x$key)
    ina <- is.na(ix)
    if(any(is.na(ix))) message(paste("thermo.refs: reference key(s)",
      paste(key[ina], collapse = ","), "not found"))
    return(x[ix, ])
  } else if(is.numeric(key)) {
    # get the source keys for the indicated species
    sinfo <- suppressMessages(info(key))
    if(keep.duplicates) {
      # output a single reference for each species 20180927
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

checkEOS <- function(eos, state, prop, ret.diff=FALSE) {
  # compare calculated properties from equation-of-state
  # parameters with reference (tabulated) values
  # print message and return the calculated value
  # if tolerance is exceeded
  # or NA if the difference is within the tolerance
  # 20110808 jmd
  thermo <- get("thermo", CHNOSZ)
  # get calculated value based on EOS
  Theta <- 228  # K
  if(identical(state, "aq")) {
    if(prop=="Cp") {
      # value of X consistent with IAPWS95
      X <- -2.773788E-7
      # we use the value of X consistent with SUPCRT
      X <- -3.055586E-7
      refval <- eos$Cp
      calcval <- eos$c1 + eos$c2/(298.15-Theta)^2 + eos$omega*298.15*X
      tol <- thermo$opt$Cp.tol
      units <- paste(eos$E_units, "K-1 mol-1")
    } else if(prop=="V") {
      # value of Q consistent with IAPWS95
      Q <- 0.00002483137
      # value of Q consistent with SUPCRT92
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
    # all other states
    if(prop=="Cp") {
      refval <- eos$Cp
      Tr <- 298.15
      calcval <- eos$a + eos$b*Tr + eos$c*Tr^-2 + eos$d*Tr^-0.5 + eos$e*Tr^2 + eos$f*Tr^eos$lambda
      tol <- thermo$opt$Cp.tol
      units <- paste(eos$E_units, "K-1 mol-1")
    }
  }
  # calculate the difference
  diff <- calcval - refval
  if(ret.diff) return(diff)
  else {
    # return the calculated value
    # if the difference is higher than tol
    if(!is.na(calcval)) {
      if(!is.na(refval)) {
        if(abs(diff) > tol) {
          message(paste("checkEOS: ", prop, " of ", eos$name, " ", eos$state, " (", rownames(eos),
            ") differs by ", round(diff,2), " ", units, " from tabulated value", sep=""))
          return(calcval)
        }
      } else return(calcval)
    }
  }
  # return NA in most cases
  return(NA)
}

checkGHS <- function(ghs, ret.diff=FALSE) {
  # compare calculated G from H and S with reference (tabulated) values
  # print message and return the calculated value if tolerance is exceeded
  # or NA if the difference is within the tolerance
  # 20110808 jmd
  thermo <- get("thermo", CHNOSZ)
  # get calculated value based on H and S
  Se <- entropy(as.character(ghs$formula))
  isJoules <- ghs$E_units == "J"
  if(any(isJoules)) Se[isJoules] <- convert(Se[isJoules], "J")
  refval <- ghs$G
  DH <- ghs$H
  S <- ghs$S
  Tr <- 298.15
  calcval <- DH - Tr * (S - Se)
  # now on to the comparison
  # calculate the difference
  diff <- calcval - refval
  if(ret.diff) return(diff)
  else if(!is.na(calcval)) {
    if(!is.na(refval)) {
      diff <- calcval - refval
      if(abs(diff) > thermo$opt$G.tol) {
        message(paste("checkGHS: G of ", ghs$name, " ", ghs$state, " (", rownames(ghs),
          ") differs by ", round(diff), " ", ghs$E_units, " mol-1 from tabulated value", sep=""))
        return(calcval)
      }
    } else return(calcval)
  } else {
    # calculating a value of G failed, perhaps because of missing elements
    return(NULL)
  }
  # return NA in most cases
  return(NA)
}

check.OBIGT <- function() {
  # function to check self-consistency between
  # values of Cp and V vs. EOS parameters
  # and among G, H, S values
  # 20110808 jmd replaces 'check=TRUE' argument of info()
  checkfun <- function(what) {
    # looking at thermo$OBIGT
    if(what=="OBIGT") tdata <- get("thermo", CHNOSZ)$OBIGT
    else if(what=="DEW") tdata <- read.csv(system.file("extdata/OBIGT/DEW.csv", package="CHNOSZ"), as.is=TRUE)
    else if(what=="SLOP98") tdata <- read.csv(system.file("extdata/OBIGT/SLOP98.csv", package="CHNOSZ"), as.is=TRUE)
    else if(what=="SUPCRT92") tdata <- read.csv(system.file("extdata/OBIGT/SUPCRT92.csv", package="CHNOSZ"), as.is=TRUE)
    else if(what=="OldAA") tdata <- read.csv(system.file("extdata/OBIGT/OldAA.csv", package="CHNOSZ"), as.is=TRUE)
    else if(what=="AS04") tdata <- read.csv(system.file("extdata/OBIGT/AS04.csv", package="CHNOSZ"), as.is=TRUE)
    else if(what=="AkDi") tdata <- read.csv(system.file("extdata/OBIGT/AkDi.csv", package="CHNOSZ"), as.is=TRUE)
    ntot <- nrow(tdata)
    # where to keep the results
    DCp <- DV <- DG <- rep(NA,ntot)
    # first get the aqueous species
    isaq <- tdata$state=="aq"
    if(any(isaq)) {
      eos.aq <- OBIGT2eos(tdata[isaq,], "aq")
      DCp.aq <- checkEOS(eos.aq, "aq", "Cp", ret.diff = TRUE)
      DV.aq <- checkEOS(eos.aq, "aq", "V", ret.diff = TRUE)
      cat(paste("check.OBIGT: GHS for", sum(isaq), "aq species in", what, "\n"))
      DG.aq <- checkGHS(eos.aq, ret.diff = TRUE)
      # store the results
      DCp[isaq] <- DCp.aq
      DV[isaq] <- DV.aq
      DG[isaq] <- DG.aq
    }
    # then other species, if they are present
    if(sum(!isaq) > 0) {
      eos.cgl <- OBIGT2eos(tdata[!isaq,], "cgl")
      DCp.cgl <- checkEOS(eos.cgl, "cgl", "Cp", ret.diff = TRUE)
      cat(paste("check.OBIGT: GHS for", sum(!isaq), "cr,gas,liq species in", what, "\n"))
      DG.cgl <- checkGHS(eos.cgl, ret.diff = TRUE)
      DCp[!isaq] <- DCp.cgl
      DG[!isaq] <- DG.cgl
    }
    # put it all together
    out <- data.frame(table = what, ispecies = 1:ntot, name = tdata$name, state = tdata$state, E_units = tdata$E_units, DCp = DCp, DV = DV, DG = DG)
    return(out)
  }
  # check default database (OBIGT)
  out <- checkfun("OBIGT")
  # check optional data
  out <- rbind(out, checkfun("DEW"))
  out <- rbind(out, checkfun("SLOP98"))
  out <- rbind(out, checkfun("SUPCRT92"))
  out <- rbind(out, checkfun("OldAA"))
  # set differences within a tolerance to NA
  out$DCp[abs(out$DCp) < 1] <- NA
  out$DV[abs(out$DV) < 1] <- NA
  out$DG[abs(out$DG) < 500] <- NA
  # take out species where all reported differences are NA
  ina <- is.na(out$DCp) & is.na(out$DV) & is.na(out$DG)
  out <- out[!ina,]
  # round the values
  out$DCp <- round(out$DCp,2)
  out$DV <- round(out$DV,2)
  out$DG <- round(out$DG)
  # how to make the file at extdata/thermo/OBIGT_check.csv
  # write.csv(out,"OBIGT_check.csv",na="",row.names=FALSE)
  # return the results
  return(out)
}

RH2OBIGT <- function(compound=NULL, state="cr", file=system.file("extdata/adds/RH98_Table15.csv", package="CHNOSZ")) {
  # get thermodynamic properties and equations of state parameters using 
  # group contributions from Richard and Helgeson, 1998   20120609 jmd
  # read the compound names, physical states, chemical formulas and group stoichiometry from the file
  # we use check.names=FALSE because the column names are the names of the groups,
  # and are not syntactically valid R names, and stringsAsFactors=FALSE
  # so that formulas are read as characters (for checking with as.chemical.formula)
  dat <- read.csv(file, check.names=FALSE, stringsAsFactors=FALSE)
  # "compound" the compound names and states from the file
  comate.arg <- comate.dat <- paste(dat$compound, "(", dat$state, ")", sep="")
  # "compound" the compound names and states from the arguments
  if(!is.null(compound)) comate.arg <- paste(compound, "(", state, ")", sep="")
  # identify the compounds
  icomp <- match(comate.arg, comate.dat)
  # check if all compounds were found
  ina <- is.na(icomp)
  if(any(ina)) stop(paste("compound(s)", paste(comate.arg[ina], collapse=" "), "not found in", file))
  # initialize output data frame
  out <- get("thermo", CHNOSZ)$OBIGT[0, ]
  # loop over the compounds
  for(i in icomp) {
    # the group stoichiometry for this compound
    thisdat <- dat[i, ]
    # take out groups that are NA or 0
    thisdat <- thisdat[, !is.na(thisdat)]
    thisdat <- thisdat[, thisdat!=0]
    # identify the groups in this compound
    igroup <- 4:ncol(thisdat)
    ispecies <- info(colnames(thisdat)[igroup], state=thisdat$state)
    # check if all groups were found
    ina <- is.na(ispecies)
    if(any(ina)) stop(paste("group(s)", paste(colnames(thisdat)[igroup][ina], collapse=" "), "not found in", thisdat$state, "state"))
    # group additivity of properties and parameters: add contributions from all groups
    thiseos <- t(colSums(get("thermo", CHNOSZ)$OBIGT[ispecies, 9:21] * as.numeric(thisdat[, igroup])))
    # group additivity of chemical formula
    formula <- as.chemical.formula(colSums(i2A(ispecies) * as.numeric(thisdat[, igroup])))
    # check if the formula is the same as in the file
    if(!identical(formula, thisdat$formula)) 
      stop(paste("formula", formula, "of", comate.dat[i], "(from groups) is not identical to", thisdat$formula, "(listed in file)" ))
    # build the front part of OBIGT data frame
    thishead <- data.frame(name=thisdat$compound, abbrv=NA, formula=formula, state=thisdat$state, 
      ref1=NA, ref2=NA, date=today(), E_units = "cal", stringsAsFactors=FALSE)
    # insert the result into the output
    out <- rbind(out, cbind(thishead, thiseos))
  }
  return(out)
}

# dump all thermodynamic data in default and optional OBIGT files 20171121
dumpdata <- function(file=NULL) {
  # default database (OBIGT)
  dat <- get("thermo", CHNOSZ)$OBIGT
  OBIGT <- cbind(source="OBIGT", dat)
  # optional data
  dat <- read.csv(system.file("extdata/OBIGT/DEW.csv", package="CHNOSZ"), as.is=TRUE)
  DEW <- cbind(source="DEW", dat)
  dat <- read.csv(system.file("extdata/OBIGT/SLOP98.csv", package="CHNOSZ"), as.is=TRUE)
  SLOP98 <- cbind(source="SLOP98", dat)
  dat <- read.csv(system.file("extdata/OBIGT/SUPCRT92.csv", package="CHNOSZ"), as.is=TRUE)
  SUPCRT92 <- cbind(source="SUPCRT92", dat)
  # put it all together
  out <- rbind(OBIGT, DEW, SLOP98, SUPCRT92)
  # quote columns 2 (name) and 3 (abbrv) because they have commas for some entries
  if(!is.null(file)) write.csv(out, file, row.names=FALSE, quote=c(2, 3))
  else(return(out))
}

### unexported functions ###

# Take a data frame in the format of thermo$OBIGT of one or more rows,
#   remove scaling factors from equations-of-state parameters,
#   and apply new column names depending on the state.
# And convert energy units from J to cal (used by subcrt()) 20190530
# If fixGHS is TRUE a missing one of G, H or S for any species is calculated
#   from the other two and the chemical formula of the species.
# This function is used by both info and subcrt when retrieving entries from the thermodynamic database.
OBIGT2eos <- function(OBIGT, state, fixGHS = FALSE, tocal = FALSE) {
  # remove scaling factors from EOS parameters
  # and apply column names depending on the EOS
  if(identical(state, "aq")) {
    # species in the Akinfiev-Diamond model (AkDi) have NA for Z 20190219
    isAkDi <- is.na(OBIGT[, 21])
    # remove scaling factors for the HKF species, but not for the AkDi species
    # protect this by an if statement to workaround error in subassignment to empty subset of data frame in R < 3.6.0
    # (https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17483) 20190302
    if(any(!isAkDi)) OBIGT[!isAkDi, 14:21] <- t(t(OBIGT[!isAkDi, 14:21]) * 10^c(-1,2,0,4,0,4,5,0))
    # for AkDi species, set NA values in remaining columns (for display only)
    if(any(isAkDi)) OBIGT[isAkDi, 17:20] <- NA
    # if all of the species are AkDi, change the variable names
    if(all(isAkDi)) colnames(OBIGT)[14:21] <- c('a','b','xi','XX1','XX2','XX3','XX4','Z') 
    else colnames(OBIGT)[14:21] <- c('a1','a2','a3','a4','c1','c2','omega','Z') 
  } else {
    OBIGT[,14:21] <- t(t(OBIGT[,14:21]) * 10^c(0,-3,5,0,-5,0,0,0))
    colnames(OBIGT)[14:21] <- c('a','b','c','d','e','f','lambda','T')
  }
  if(tocal) {
    # convert values from Joules to calories 20190530
    iJ <- OBIGT$E_units=="J"
    if(any(iJ)) {
      # we only convert column 20 for aqueous species (omega), not for cgl species (lambda)  20190903
      if(identical(state, "aq")) OBIGT[iJ, c(9:12, 14:20)] <- convert(OBIGT[iJ, c(9:12, 14:20)], "cal")
      else OBIGT[iJ, c(9:12, 14:19)] <- convert(OBIGT[iJ, c(9:12, 14:19)], "cal")
    }
  }
  if(fixGHS) {
    # fill in one of missing G, H, S
    # for use esp. by subcrt because NA for one of G, H or S 
    # will preclude calculations at high T
    # which entries are missing just one
    imiss <- which(rowSums(is.na(OBIGT[,9:11]))==1)
    if(length(imiss) > 0) {
      for(i in 1:length(imiss)) {
        # calculate the missing value from the others
        ii <- imiss[i]
        GHS <- as.numeric(GHS(as.character(OBIGT$formula[ii]), G=OBIGT[ii,9], H=OBIGT[ii,10], S=OBIGT[ii,11],
                              E_units = ifelse(tocal, "cal", OBIGT$E_units[ii])))
        icol <- which(is.na(OBIGT[ii,9:11]))
        OBIGT[ii,icol+8] <- GHS[icol]
      }
    }
  }
  return(OBIGT)
}
