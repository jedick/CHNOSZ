# CHNOSZ/makeup.R

## If this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.formula.R")
#source("util.character.R")

makeup <- function(formula, multiplier=1, sum=FALSE, count.zero=FALSE) {
  # Return the elemental makeup (counts) of a chemical formula
  # that may contain suffixed and/or parenthetical subformulas and/or charge
  # and negative or positive, fractional coefficients
  # A matrix is processed by rows
  if(is.matrix(formula)) 
    return(lapply(seq_len(nrow(formula)), 
      function(i) makeup(formula[i, ])))
  # A named numeric object is returned untouched
  # (needed for recursive operation of the function?
  #  - if this is not done it messes up the HOX example in protein.info.Rd)
  if(!is.null(names(formula)) & is.numeric(formula)) return(formula)
  # A list of named objects is also returned untouched
  if(is.list(formula) & !is.null(names(formula[[1]]))) return(formula)
  # Prepare to multiply the formula by the multiplier, if given
  if(length(multiplier) > 1 & length(multiplier) != length(formula))
    stop("multiplier does not have length = 1 or length = number of formulas")
  multiplier <- rep(multiplier, length(formula))
  # If the formula argument has length > 1, apply the function over each formula
  if(length(formula) > 1) {
    # Get formulas for any species indices in the argument
    formula <- get.formula(formula)
    out <- lapply(seq_along(formula), function(i) {
      makeup(formula[i], multiplier[i])
    })
    # If sum is TRUE, take the sum of all formulas
    if(sum) {
      out <- unlist(out)
      out <- tapply(out, names(out), sum)
    } else if(count.zero) {
      # If count.zero is TRUE, all elements appearing in any 
      # of the formulas are counted for each species
      # First construct elemental makeup showing zero of each element
      em0 <- unlist(out)
      # Exclude NA from the element names 20190802
      em0 <- em0[!is.na(em0)]
      em0 <- tapply(em0, names(em0), sum)
      em0[] <- 0
      # Create NA matrix for NA formulas 20190802
      emNA <- em0
      emNA[] <- NA
      # Then sum each formula and the zero vector,
      # using tapply to group the elements
      out <- lapply(out, function(x) {
        if(anyNA(x)) emNA else {
          xem <- c(x, em0)
          tapply(xem, names(xem), sum)
        }
      })
    }
    return(out)
  }
  # If the formula argument is numeric,
  # and if the thermo object is available,
  # get the formula of that numbered species from thermo()$OBIGT
  if("CHNOSZ" %in% .packages()) {
    thermo <- get("thermo", CHNOSZ)
    if(is.numeric(formula)) formula <- thermo$OBIGT$formula[formula]
  }
  # First deal with charge
  cc <- count.charge(formula)
  # count.elements doesn't know about charge so we need
  # to explicate the elemental symbol for it
  formula <- cc$uncharged
  if(cc$Z != 0 ) formula <- paste(formula, "Z", cc$Z, sep="")
  # Now "Z" will be counted
  # If there are no subformulas, just use count.elements
  if(length(grep("(\\(|\\)|\\*|\\:)", formula))==0) {
    out <- count.elements(formula)
  } else {
    # Count the subformulas
    cf <- count.formulas(formula)
    # Count the elements in each subformula
    ce <- lapply(names(cf), count.elements)
    # Multiply elemental counts by respective subformula counts
    mcc <- lapply(seq_along(cf), function(i) ce[[i]]*cf[i])
    # Unlist the subformula counts and sum them together by element
    um <- unlist(mcc)
    out <- unlist(tapply(um, names(um), sum, simplify=FALSE))
  }
  # All done with the counting, now apply the multiplier
  out <- out * multiplier
  # Complain if there are any elements that look strange
  if("CHNOSZ" %in% .packages()) {
    are.elements <- names(out) %in% thermo$element$element
    if(!all(are.elements)) warning(paste("element(s) not in thermo()$element:", 
      paste(names(out)[!are.elements], collapse=" ") ))
  }
  # Done!
  return(out)
}

count.elements <- function(formula) {
  # Count the elements in a chemical formula   20120111 jmd
  # This function expects a simple formula,
  # no charge or parenthetical or suffixed subformulas
  # Regular expressions inspired by an answer on
  # http://stackoverflow.com/questions/4116786/parsing-a-chemical-formula-from-a-string-in-c
  if(is.na(formula)) return(NA)
  #elementRegex <- "([A-Z][a-z]*)([0-9]*)"
  elementSymbol <- "([A-Z][a-z]*)"
  # Here, element coefficients can be signed (+ or -) and have a decimal point
  elementCoeff <- "((\\+|-|\\.|[0-9])*)"
  elementRegex <- paste(elementSymbol, elementCoeff, sep="")
  # Stop if it doesn't look like a chemical formula 
  validateRegex <- paste("^(", elementRegex, ")+$", sep="")
  if(length(grep(validateRegex, formula)) == 0)
    stop(paste("'",formula,"' is not a simple chemical formula", sep="", collapse="\n"))
  # Where to put the output
  element <- character()
  count <- numeric()
  # From now use "f" for formula to make writing the code easier
  f <- formula
  # We want to find the starting positions of all the elemental symbols
  # Make substrings, starting at every position in the formula
  fsub <- sapply(1:nchar(f), function(i) substr(f, i, nchar(f)))
  # Get the numbers (positions) that start with an elemental symbol
  # i.e. an uppercase letter
  ielem <- grep("^[A-Z]", fsub)
  # For each elemental symbol, j is the position before the start of the next 
  # symbol (or the position of the last character of the formula)
  jelem <- c(tail(ielem - 1, -1), nchar(f))
  # Assemble the stuff: each symbol-coefficient combination
  ec <- sapply(seq_along(ielem), function(i) substr(f, ielem[i], jelem[i]))
  # Get the individual element symbols and coefficients
  myelement <- gsub(elementCoeff, "", ec)
  mycount <- as.numeric(gsub(elementSymbol, "", ec))
  # Any missing coefficients are unity
  mycount[is.na(mycount)] <- 1
  # Append to the output
  element <- c(element, myelement)
  count <- c(count, mycount)
  # In case there are repeated elements, sum all of their counts
  # (tapply hint from https://stat.ethz.ch/pipermail/r-help/2011-January/265341.html)
  # use simplify=FALSE followed by unlist to get a vector, not array 20171005
  out <- unlist(tapply(count, element, sum, simplify=FALSE))
  # tapply returns alphabetical sorted list. keep the order appearing in the formula
  out <- out[match(unique(element), names(out))]
  return(out)
}

### Unexported functions ###

# Returns a list with named elements
#   `Z` - the numeric value of the charge
#   `uncharged` - the original formula string excluding the charge
count.charge <- function(formula) {
  # Count the charge in a chemical formula   20120113 jmd
  # Everything else is counted as the uncharged part of the formula
  # Examples:
  # C3H6NO2-    # charge is -1, uncharged part is C3H6NO2
  # XYZ-1.2     # charge is -1.2, uncharged part is XYZ
                # ( makeup() treats this equivalently to XY-0.2 )
  # C-1.567     # charge is -1.567, uncharged part is C1 (one carbon)
  # C-1.567+0   # charge is zero, uncharged part C-1.567 (-1.567 carbons)
  # Most of the time the formula isn't charged
  Z <- 0
  uncharged <- formula
  # We have charge if there is a plus or minus sign that is 
  # followed by only digits (possibly also a decimal point)
  # at the end of the formula
  charged <- grep("(\\+|-)(\\.|[0-9])*$", formula)
  if(length(charged) > 0) {
    fsplit <- unlist(strsplit(formula, ""))
    # The position of the sign of charge
    isign <- tail(grep("(\\+|-)", fsplit), 1)
    # The sign as +1 or -1
    if(fsplit[isign]=="+") sign <- 1 else sign <- -1
    # If the +/- symbol is at the end by itself thats a +/- 1
    # otherwise we multiply the magnitude by the sign
    nf <- nchar(formula)
    if(isign==nf) Z <- sign 
    else Z <- sign * as.numeric(substr(formula, isign+1 ,nf))
    # All set ... Zzzz this is wearing me out!
    # got our charge, so remove it from the formula
    uncharged <- substr(formula, 1, isign-1)
  } 
  # Assemble the output
  out <- list(Z=Z, uncharged=uncharged)
  return(out)
}

# Returns a numeric vector with names refering to each of the subformulas or the whole formula if there are no subformulas
count.formulas <- function(formula) {
  # Count the subformulas in a chemical formula   20120112 jmd
  # The formula may include multiple unnested parenthetical subformulas
  # and/or one suffixed subformula (separated by * or :, with no parentheses in suffix)
  # e.g. formula for sepiolite, Mg4Si6O15(OH)2(H2O)2*4H2O gives
  #           count
  # OH            2
  # H2O           6
  # Mg4Si6O15     1
  # Other formulas that are able to be parsed
  # NaCa2(Al5Si13)O36*14H2O   # multiple-digit coefficient on suffix
                              # and no explicit coefficient on parenthetical term
  # C6H14N2O2:HCl             # a suffix with no numbers at all
  # Where to keep the output
  subform <- character()
  count <- numeric()
  # Deal with parentheses if we have them
  # First split the characters
  fsplit <- strsplit(formula, "")[[1]]
  # Then identify the positions of the parentheses
  iopen <- grep("\\(", fsplit)
  iclose <- grep("\\)", fsplit)
  # Do we have parentheses?
  if(length(iopen) > 0 | length(iclose) > 0) {
    # Are all parentheses paired?
    if(length(iopen) != length(iclose)) stop("formula has unpaired parentheses")
    # Are the parentheses unnested?
    iparen <- as.numeric(matrix(c(iopen, iclose), nrow=2, byrow=TRUE))
    if(any(diff(iparen) < 0)) stop(paste("formula has nested parentheses: ", formula))
    # iend will be the last position including coefficients
    # (e.g. the position of 2 in (OH)2)
    iend <- iclose
    # inum are all the positions with any part of a number
    # (including sign and decimal point)
    inum <- grep("\\+|-|\\.|[0-9]", fsplit)
    # ichar are all other positions
    nf <- nchar(formula)
    ichar <- (1:nf)[-inum]
    # Loop over the parentheses
    for(i in seq_along(iopen)) {
      # If the position after close paren is a number, parse it
      if((iclose[i]+1) %in% inum) {
        # iend is the position of the next character, minus one
        ic <- head(ichar[ichar-iclose[i] > 0], 1)
        # If it's not present, then we're at the end of the formula
        if(length(ic)==0) iend[i] <- nf else iend[i] <- ic - 1
        # All the stuff after the iclose up to iend is the coefficient
        count <- c(count, as.numeric(substr(formula,iclose[i]+1, iend[i])))
      } else count <- c(count, 1)
      # Everything between iopen and iclose is the subformula
      subform <- c(subform, substr(formula,iopen[i]+1, iclose[i]-1))
    }
    # We can remove all of the paranthetical terms and their coefficients
    for(i in seq_along(iopen)) fsplit[iopen[i]:iend[i]] <- ""
    formula <- paste(fsplit, collapse="")
  }
  # Deal with a suffixed subformula if we have one
  # We assume that there is at most one suffix, so
  # at most two strings after both of these splitting tests
  fsplit <- unlist(strsplit(formula, "*", fixed=TRUE))
  fsplit <- unlist(strsplit(fsplit, ":", fixed=TRUE))
  formula <- fsplit[1]
  if(length(fsplit) > 1) {
    # Parse the coefficient; this time it's in front
    f2 <- fsplit[2]
    f2split <- unlist(strsplit(f2, ""))
    # The positions of the numbers
    inum <- grep("\\+|-|\\.|[0-9]", f2split)
    if(length(inum)==0) {
      # No numbers, we have one of the subformula
      mycount <- 1
      mysub <- f2
    } else {
      # The position of the first character
      ic <- head((seq_along(f2split))[-inum], 1)
      # If the first character is before the first number,
      # the coefficient is one
      if(ic < inum[1]) mycount <- 1
      else mycount <- as.numeric(substr(f2, inum[1], ic-1))
      # We also have the subformula
      mysub <- substr(f2, ic, nchar(f2))
    }
    # Add to existing subformula (as with H2O in sepiolite)
    # or append the coefficient and subformula
    subform <- c(subform, mysub)
    count <- c(count, mycount)
  }
  # Add in the remaining formula, if there is any left
  if(!is.null(formula)) {
    subform <- c(subform, formula)
    count <- c(count, 1)
  }
  # Assemble the counts 
  out <- tapply(count, subform, sum)
  return(out)
}

