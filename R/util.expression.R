# CHNOSZ/util.expression.R
# write descriptions of chemical species, properties, reactions, conditions
# modified from describe(), axis.label()  20120121 jmd

## if this file is interactively sourced, the following are also needed to provide unexported functions:
#source("util.character.R")

expr.species <- function(species, state = "aq", value = NULL, log = FALSE, molality = FALSE, use.state = FALSE, use.makeup = FALSE) {
  # make plotting expressions for chemical formulas
  # that include subscripts, superscripts (if charged)
  # and optionally designations of states +/- loga or logf prefix
  if(length(species) > 1) (stop("more than one species"))
  # convert to character so that "1", "2", etc. don't get converted to chemical formulas via makeup()
  species <- as.character(species)
  if(use.makeup) {
    # the counts of elements in the species:
    # here we don't care too much if an "element" is a real element
    # (listed in thermo()$element), so we suppress warnings
    elements <- suppressWarnings(try(makeup(species), TRUE))
  } else elements <- split.formula(species)
  # if species can't be parsed as a chemical formula, we don't do the formula formatting
  if(inherits(elements, "try-error") | !is.numeric(elements)) expr <- species
  else {
    # where we'll put the expression
    expr <- ""
    # loop over elements
    for(i in 1:length(elements)) {
      if(names(elements)[i] != 'Z') {
        # append the elemental symbol
        expr <- substitute(paste(a, b), list(a=expr, b=names(elements)[i]))
        # recover the coefficient
        coeff <- elements[i]
        if(coeff!=1) {
          # append the coefficient
          expr <- substitute(a[b], list(a=expr, b=as.character(coeff)))
        }
      } else {
        # for charged species, don't show "Z" but do show e.g. "+2"
        coeff <- elements[i]
        if(coeff==-1) coeff <- "-"
        else if(coeff==1) coeff <- "+"
        else if(coeff > 0) coeff <- paste("+", as.character(coeff), sep="")
        else coeff <- as.character(coeff)
        # append the coefficient as a superscript
        expr <- substitute(a^b, list(a=expr, b=coeff))
      }
    }
  }
  # write the physical state
  if(use.state) {
    # subscript it if we're not giving the value
    if(is.null(value)) expr <- substitute(a[group('(',italic(b),')')],list(a=expr, b=state))
    else expr <- substitute(a*group('(',italic(b),')'),list(a=expr, b=state))
  }
  # write a variable and value if given
  if(!is.null(value) | log | molality) {
    # write [logarithm of] activity or fugacity (or molality 20171101)
    var <- "a"
    if(molality) var <- "m"
    if(state %in% c("g", "gas")) var <- "f"
    expr <- substitute(italic(a)[b], list(a = var, b = expr))
    # use the logarithm?
    if(log) expr <- substitute(log ~ a, list(a = expr))
    # write the value if not NULL or NA
    if(!is.null(value)) {
      if(!is.na(value)) expr <- substitute(a == b, list(a = expr, b = value))
    }
  }
  return(expr)
}

expr.property <- function(property, molality=FALSE) {
  # a way to make expressions for various properties
  # e.g. expr.property('DG0r') for standard molal Gibbs 
  # energy change of reaction
  propchar <- s2c(property)
  expr <- ""
  # some special cases
  if(is.na(property)) return("")
  if(property=="logK") return(quote(log~italic(K)))
  # grepl here b/c diagram() uses "loga.equil" and "loga.basis"
  if(grepl("loga", property)) {
    if(molality) return(quote(log~italic(m)))
    else return(quote(log~italic(a)))
  }
  if(property=="alpha") return(quote(alpha))
  if(property=="Eh") return("Eh")
  if(property=="pH") return("pH")
  if(property=="pe") return("pe")
  if(property=="IS") return(quote(italic(I)))
  if(property=="ZC") return(quote(italic(Z)[C]))
  # process each character in the property abbreviation
  prevchar <- character()
  for(i in 1:length(propchar)) {
    if(i > 1) prevchar <- thischar
    thischar <- propchar[i]
    # unless indicated below, uppercase letters are italicized
    # and lowercase letters are italicized and subscripted
    # (includes f for property of formation and r for property of reaction)
    if(thischar %in% LETTERS) thisexpr <- substitute(italic(a), list(a=thischar))
    else if(thischar %in% letters) thisexpr <- substitute(""[italic(a)], list(a=thischar))
    else thisexpr <- substitute(a, list(a=thischar))
    # D for greek Delta
    # p for subscript italic P (in Cp)
    # 0 for degree sign (but not immediately following a number e.g. 2.303)
    # l for subscript small lambda
    # ' for prime symbol (like "minute")
    if(thischar=='D') thisexpr <- substitute(Delta)
    if(thischar=='p') thisexpr <- substitute(a[italic(P)], list(a=""))
    if(thischar=='0' & !can.be.numeric(prevchar)) thisexpr <- substitute(degree)
    if(thischar=='l') thisexpr <- substitute(a[lambda], list(a=""))
    if(thischar=="'") thisexpr <- substitute(minute)
    # put it together
    expr <- substitute(a*b, list(a=expr, b=thisexpr))
  }
  return(expr)
}

expr.units <- function(property, prefix="", per="mol") {
  # make an expression describing units
  # unless we have match below, there will be no units
  # (e.g., logK, pH, pe)
  expr <- ""
  # A, G, H - energy
  if(grepl("A", property)) expr <- substitute(a, list(a=E.units()))
  if(grepl("G", property)) expr <- substitute(a, list(a=E.units()))
  if(grepl("H", property) & !grepl("pH", property)) expr <- substitute(a, list(a=E.units()))
  # Cp, S - energy (per K)
  if(grepl("Cp", property)) expr <- substitute(a~K^-1, list(a=E.units()))
  if(grepl("S", property)) expr <- substitute(a~K^-1, list(a=E.units()))
  # V - volume
  if(grepl("V", property)) expr <- substitute(a^3, list(a="cm"))
  # E - volume (per K)
  if(grepl("E", property)) expr <- substitute(a^3~K^-1, list(a="cm"))
  # P - pressure
  if(grepl("P", property)) expr <- substitute(a, list(a=P.units()))
  # T - temperature
  if(grepl("T", property)) {
    expr <- substitute(a, list(a=T.units()))
    # add a degree sign for Celsius
    if(T.units()=="C") expr <- substitute(degree*a, list(a=expr))
  }
  # Eh - electrical potential
  if(grepl("Eh", property)) expr <- substitute(a, list(a="volt"))
  # IS - ionic strength
  if(grepl("IS", property)) expr <- substitute(a, list(a=mol~kg^-1))
  if(!expr=="") {
    # add prefix if appropriate
    if(!prefix=="") expr <- substitute(a*b, list(a=prefix, b=expr))
    # add mol^-1 if appropriate
    if(!any(sapply(c("P", "T", "Eh", "IS"), function(x) grepl(x, property))))
      expr <- substitute(a~b^-1, list(a=expr, b=per))
  }
  return(expr)
}

axis.label <- function(label, units=NULL, basis=thermo()$basis, prefix="", molality=FALSE) {
  # make a formatted axis label from a generic description
  # it can be a chemical property, condition, or chemical activity in the system
  # if the label matches one of the basis species
  # or if the state is specified, it's a chemical activity
  # 20090826: just return the argument if a comma is already present
  # (it's good for custom labels that shouldn't be italicized)
  if(grepl(",", label)) return(label)
  if(label %in% rownames(basis)) {
    # 20090215: the state this basis species is in
    state <- basis$state[match(label, rownames(basis))]
    # get the formatted label
    desc <- expr.species(label, state = state, log = TRUE, molality = molality)
  } else {
    # the label is for a chemical property or condition
    # make the label by putting a comma between the property and the units
    property <- expr.property(label, molality=molality)
    if(is.null(units)) units <- expr.units(label, prefix=prefix)
    # no comma needed if there are no units
    if(units=="") desc <- substitute(a, list(a=property))
    else desc <- substitute(list(a, b), list(a=property, b=units))
  }
  # done!
  return(desc)
}

describe.basis <- function(basis = thermo()$basis, ibasis = 1:nrow(basis),
  digits = 1, oneline = FALSE, molality = FALSE, use.pH = TRUE) {
  # make expressions for the chemical activities/fugacities of the basis species
  propexpr <- valexpr <- character()
  for(i in ibasis) {
    if(can.be.numeric(basis$logact[i])) {
      # we have an as.numeric here in case the basis$logact is character
      # (by inclusion of a buffer for one of the other basis species)
      if(rownames(basis)[i]=="H+" & use.pH) {
        propexpr <- c(propexpr, "pH")
        valexpr <- c(valexpr, format(round(-as.numeric(basis$logact[i]), digits), nsmall=digits))
      } else {
        # propexpr is logarithm of activity or fugacity
        propexpr <- c(propexpr, expr.species(rownames(basis)[i], state=basis$state[i], log = TRUE, molality = molality))
        valexpr <- c(valexpr, format(round(as.numeric(basis$logact[i]), digits), nsmall=digits))
      }
    } else {
      # a non-numeric value is the name of a buffer
      valexpr <- c(valexpr, basis$logact[i])
      # propexpr is pH, activity or fugacity
      if(rownames(basis)[i]=="H+" & use.pH) propexpr <- c(propexpr, "pH")
      else propexpr <- c(propexpr, expr.species(rownames(basis)[i], state=basis$state[i], value = NA, log = FALSE, molality = molality))
    }
  }
  # write an equals sign between the property and value
  desc <- character()
  for(i in seq_along(propexpr)) {
    thisdesc <- substitute(a==b, list(a=propexpr[[i]], b=valexpr[[i]]))
    if(oneline) {
      # put all the property/value equations on one line, separated by commas
      if(i==1) desc <- substitute(a, list(a=thisdesc))
      else desc <- substitute(list(a, b), list(a=desc, b=thisdesc))
    } else desc <- c(desc, thisdesc)
  }
  return(as.expression(desc))
}

describe.property <- function(property=NULL, value=NULL, digits=0, oneline=FALSE, ret.val=FALSE) {
  # make expressions for pressure, temperature, other conditions
  if(is.null(property) | is.null(value)) stop("property or value is NULL")
  propexpr <- valexpr <- character()
  for(i in 1:length(property)) {
    propexpr <- c(propexpr, expr.property(property[i]))
    if(is.na(value[i])) thisvalexpr <- ""
    else if(value[i]=="Psat") thisvalexpr <- quote(italic(P)[sat])
    else {
      thisvalue <- format(round(as.numeric(value[i]), digits), nsmall=digits)
      thisunits <- expr.units(property[i])
      thisvalexpr <- substitute(a~b, list(a=thisvalue, b=thisunits))
    }
    valexpr <- c(valexpr, as.expression(thisvalexpr))
  } 
  # with ret.val=TRUE, return just the value with the units (e.g. 55 degC)
  if(ret.val) return(valexpr)
  # write an equals sign between the property and value
  desc <- character()
  for(i in seq_along(propexpr)) {
    if(is.na(value[i])) thisdesc <- propexpr[[i]]
    else thisdesc <- substitute(a==b, list(a=propexpr[[i]], b=valexpr[[i]]))
    if(oneline) {
      # put all the property/value equations on one line, separated by commas
      if(i==1) desc <- substitute(a, list(a=thisdesc))
      else desc <- substitute(list(a, b), list(a=desc, b=thisdesc))
    } else desc <- c(desc, thisdesc)
  }
  return(as.expression(desc))
}

describe.reaction <- function(reaction, iname=numeric(), states=NULL) {
  # make an expression describing the reaction that is
  # the 'reaction' part of subcrt() output
  reactexpr <- prodexpr <- character()
  # loop over the species in the reaction
  for(i in 1:nrow(reaction)) {
    # get the name or the chemical formula of the species
    if(i %in% iname) species <- reaction$name[i]
    else {
      # should the chemical formula have a state?
      if(identical(states,"all") | i %in% states) species <- expr.species(reaction$formula[i], state=reaction$state[i], use.state=TRUE)
      else species <- expr.species(reaction$formula[i], state=reaction$state[i])
    }
    # get the absolute value of the reaction coefficient
    abscoeff <- abs(reaction$coeff[i])
    # put the coefficient in if it's not 1
    if(abscoeff==1) coeffspec <- species
    else {
      # we put in some space if the coefficient comes before a name
      if(i %in% iname) coeffspec <- substitute(a~b, list(a=abscoeff, b=species))
      else coeffspec <- substitute(a*b, list(a=abscoeff, b=species))
    }
    # is it a reactant or product?
    if(reaction$coeff[i] < 0) {
      if(length(reactexpr)==0) reactexpr <- substitute(a, list(a=coeffspec))
      else reactexpr <- substitute(a+b, list(a=reactexpr, b=coeffspec))
    } else {
      # it's a product
      if(length(prodexpr)==0) prodexpr <- substitute(a, list(a=coeffspec))
      else prodexpr <- substitute(a+b, list(a=prodexpr, b=coeffspec))
    }
  }
  # put an equals sign between reactants and products
  # change this to unicode for the reaction double-arrow 20190218 \u21cc
  # go back to equals - the arrow only works out-of-the-box on linux 20190817
  desc <- substitute(a ~ "=" ~ b, list(a=reactexpr, b=prodexpr))
  return(desc)
}

# make formatted text for activity ratio 20170217
# allow changing the bottom ion 20200716
ratlab <- function(top = "K+", bottom = "H+", molality = FALSE) {
  # the charges
  Ztop <- makeup(top)["Z"]
  Zbottom <- makeup(bottom)["Z"]
  # the text for the exponents
  exp.bottom <- as.character(Ztop)
  exp.top <- as.character(Zbottom)
  if(exp.top=="1") exp.top <- ""
  if(exp.bottom=="1") exp.bottom <- ""
  # the expression for the top and bottom
  expr.top <- expr.species(top)
  expr.bottom <- expr.species(bottom)
  # with molality, change a to m
  a <- ifelse(molality, "m", "a")
  # the final expression
  substitute(log~(italic(a)[expr.top]^exp.top / italic(a)[expr.bottom]^exp.bottom),
             list(a = a, expr.top = expr.top, exp.top = exp.top, expr.bottom = expr.bottom, exp.bottom = exp.bottom))
}

# make formatted text for thermodynamic system 20170217
syslab <- function(system = c("K2O", "Al2O3", "SiO2", "H2O"), dash="-") {
  for(i in seq_along(system)) {
    expr <- expr.species(system[i])
    # use en dash here
    if(i==1) lab <- expr else lab <- substitute(a*dash*b, list(a=lab, dash=dash, b=expr))
  }
  lab
}

### unexported function ###

split.formula <- function(formula) {
  ## like makeup(), but split apart the formula based on
  ## numbers (subscripts); don't scan for elemental symbols 20171018
  # if there are no numbers or charge, return the formula as-is
  # Change [0-9]? to [\\.0-9]* (recognize decimal point in numbers, recognize charge longer than one digit) 20220608
  if(! (grepl("[\\.0-9]", formula) | grepl("\\+[\\.0-9]*$", formula) | grepl("-[\\.0-9]*$", formula))) return(formula)
  # first split off charge
  # (assume that no subscripts are signed)
  Z <- 0
  hascharge <- grepl("\\+[\\.0-9]*$", formula) | grepl("-[\\.0-9]*$", formula)
  if(hascharge) {
    # for charge, we match + or - followed by zero or more numbers at the end of the string
    if(grepl("\\+[\\.0-9]*$", formula)) {
      fsplit <- strsplit(formula, "+", fixed=TRUE)[[1]]
      if(is.na(fsplit[2])) Z <- 1 else Z <- as.numeric(fsplit[2])
    }
    if(grepl("-[\\.0-9]*$", formula)) {
      fsplit <- strsplit(formula, "-")[[1]]
      # for formula=="H-citrate-2", unsplit H-citrate
      if(length(fsplit) > 2) {
        f2 <- tail(fsplit, 1)
        f1 <- paste(head(fsplit, -1), collapse="-")
        fsplit <- c(f1, f2)
      }
      if(is.na(fsplit[2])) Z <- -1 else Z <- -as.numeric(fsplit[2])
    }
    formula <- fsplit[1]
  }
  # to get strings, replace all numbers with placeholder (#), then split on that symbol
  # the outer gsub is to replace multiple #'s with one
  numhash <- gsub("#+", "#", gsub("[\\.0-9]", "#", formula))
  strings <- strsplit(numhash, "#")[[1]]
  # to get coefficients, replace all characters (non-numbers) with placeholder, then split
  charhash <- gsub("#+", "#", gsub("[^\\.0-9]", "#", formula))
  coeffs <- strsplit(charhash, "#")[[1]]
  # if the first coefficient is empty, remove it
  if(coeffs[1]=="") coeffs <- tail(coeffs, -1) else {
    # if the first string is empty, treat the first coefficient as a leading string (e.g. in 2-octanone)
    if(strings[1]=="") {
      strings[2] <- paste0(coeffs[1], strings[2])
      coeffs <- tail(coeffs, -1)
      strings <- tail(strings, -1)
    }
  }
  # if we're left with no coefficients, just return the string
  if(length(coeffs)==0 & Z==0) return(strings)
  # if we're missing a coefficient, append one
  if(length(coeffs) < length(strings)) coeffs <- c(coeffs, 1)
  # use strings as names for the numeric coefficients
  coeffs <- as.numeric(coeffs)
  names(coeffs) <- strings
  # include charge if it is not 0
  if(Z!=0) coeffs <- c(coeffs, Z=Z)
  return(coeffs)
}
