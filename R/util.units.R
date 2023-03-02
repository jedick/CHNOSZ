# CHNOSZ/util.units.R
# Set units and convert values between units

P.units <- function(units=NULL) {
  ## Change units of pressure or list the current one
  # Show the current units, if none is specified
  if(is.null(units)) return(get("thermo", CHNOSZ)$opt$P.units)
  # Argument handling
  units <- tolower(units)
  if(!units %in% c("bar","mpa")) stop("units of pressure must be either bar or MPa")
  # Set the units and return them
  if(units=="bar") with(CHNOSZ, thermo$opt$P.units <- "bar")
  if(units=="mpa") with(CHNOSZ, thermo$opt$P.units <- "MPa")
  message("changed pressure units to ", get("thermo", CHNOSZ)$opt$P.units)
}

T.units <- function(units=NULL) {
  ## Change units of temperature or list the current one
  # Show the current units, if none is specified
  if(is.null(units)) return(get("thermo", CHNOSZ)$opt$T.units)
  # Argument handling
  units <- tolower(units)
  if(!units %in% c("c","k")) stop("units of temperature must be either C or K")
  # Set the units and return them
  if(units=="c") with(CHNOSZ, thermo$opt$T.units <- "C")
  if(units=="k") with(CHNOSZ, thermo$opt$T.units <- "K")
  message("changed temperature units to ", get("thermo", CHNOSZ)$opt$T.units)
}

E.units <- function(units=NULL) {
  ## Change units of energy or list the current one
  # Show the current units, if none is specified
  if(is.null(units)) return(get("thermo", CHNOSZ)$opt$E.units)
  # Argument handling
  units <- tolower(units)
  if(!units %in% c("cal","j")) stop("units of energy must be either cal or J")
  # Set the units and return them
  if(units=="cal") with(CHNOSZ, thermo$opt$E.units <- "cal")
  if(units=="j") with(CHNOSZ, thermo$opt$E.units <- "J")
  message("changed energy units to ", get("thermo", CHNOSZ)$opt$E.units)
}

convert <- function(value, units, T=298.15,
  P=1, pH=7, logaH2O=0) {
  # Converts value(s) to the specified units

  # Process a list value if it's the output from solubility 20190525
  if(is.list(value) & !is.data.frame(value)) {
    if(!isTRUE(value$fun %in% c("solubility", "solubilities"))) stop("'value' is a list but is not the output from solubility()")
    if(!is.character(units)) stop("please specify a character argument for the destination units (e.g. ppm or logppm)")
    # Determine the element from 'balance' or 'in.terms.of', if it's available
    element <- value$in.terms.of
    if(is.null(element)) element <- value$balance
    grams.per.mole <- mass(element)
    message(paste("solubility: converting to", units, "by weight using the mass of", element))
    ppfun <- function(loga, units, grams.per.mole) {
      # Exponentiate loga to get molality
      moles <- 10^loga
      # Convert moles to mass (g)
      grams <- moles * grams.per.mole
      # Convert grams to ppb, ppm, or ppt
      ppx <- NULL
      # 1 ppt = 1 g / kg H2O
      # 1 ppm = 1 mg / kg H2O
      if(grepl("ppt", units)) ppx <- grams * 1e0
      if(grepl("ppm", units)) ppx <- grams * 1e3
      if(grepl("ppb", units)) ppx <- grams * 1e6
      if(is.null(ppx)) stop(paste("units", units, "not available for conversion"))
      # Use the logarithm if specified
      if(grepl("log", units)) ppx <- log10(ppx)
      ppx
    }
    # Do the conversion for the conserved basis species, then each species
    value$loga.balance <- ppfun(value$loga.balance, units, grams.per.mole)
    value$loga.equil <- lapply(value$loga.equil, ppfun, units = units, grams.per.mole = grams.per.mole)
    # Identify the units in the function text
    value$fun <- paste(value$fun, units, sep = ".")
    # Return the updated object
    return(value)
  }

  ### Argument handling for non-list value
  if(is.null(value)) return(NULL)
  if(!is.character(units)) stop(paste('convert: please specify',
    'a character argument for the destination units.\n',
    'possibilities include (G or logK) (C or K) (J or cal) (cm3bar or calories) (Eh or pe)\n',
    'or their lowercase equivalents.\n'),call.=FALSE)
  Units <- units # For the possible message to user
  units <- tolower(units)

  # Tests and calculations for the specified units
  if(units %in% c('c','k')) {
    CK <- 273.15
    if(units=='k') value <- value + CK
    if(units=='c') value <- value - CK 
  }
  else if(units[1] %in% c('j','cal')) {
    Jcal <- 4.184
    if(units=='j') value <- value * Jcal
    if(units=='cal') value <- value / Jcal
  }
  else if(units %in% c('g','logk')) {
    #R <- 1.9872  # Gas constant, cal K^-1 mol^-1
    R <- 8.314445  # Gas constant, J K^-1 mol^-1  20220325
    if(units=='logk') value <- value / (-log(10) * R * T)
    if(units=='g') value <- value * (-log(10) * R * T)
  }
  else if(units %in% c('cm3bar','joules')) {
    if(units=='cm3bar') value <- value * 10
    if(units=='joules') value <- value / 10
  }
  else if(units %in% c('eh','pe')) {
    R <- 0.00831470
    F <- 96.4935
    if(units=='pe') value <- value * F / ( log(10) * R * T )
    if(units=='eh') value <- value * ( log(10) * R * T ) / F
  }
  else if(units %in% c('bar','mpa')) {
    barmpa <- 10
    if(units=='mpa') value <- value / barmpa
    if(units=='bar') value <- value * barmpa
  }
  else if(units %in% c('e0','logfo2')) {
    # Convert between Eh and logfO2
    supcrt.out <- suppressMessages(subcrt(c("H2O", "oxygen", "H+", "e-"), c(-1, 0.5, 2, 2), T=T, P=P, convert=FALSE))
    if(units=='logfo2') value <- 2*(supcrt.out$out$logK + logaH2O + 2*pH + 2*(convert(value,'pe',T=T)))
    if(units=='e0') value <- convert(( -supcrt.out$out$logK - 2*pH + value/2 - logaH2O )/2, 'Eh',T=T)
  }
  else cat(paste('convert: no conversion to ',Units,' found.\n',sep=''))
  return(value)
}

### Unexported functions ###

outvert <- function(value, units) {
  # Converts the given value from the given units to those specified in thermo()$opt
  units <- tolower(units)
  opt <- get("thermo", CHNOSZ)$opt
  if(units %in% c("c", "k")) {
    if(units == "c" & opt$T.units == "K") return(convert(value, "k"))
    if(units == "k" & opt$T.units == "C") return(convert(value, "c"))
  }
  if(units %in% c("j", "cal")) {
    if(units == "j" & opt$E.units == "cal") return(convert(value, "cal"))
    if(units == "cal" & opt$E.units == "J") return(convert(value, "j"))
  }
  if(units %in% c("bar", "mpa")) {
    if(units == "mpa" & opt$P.units == "bar") return(convert(value, "bar"))
    if(units == "bar" & opt$P.units == "MPa") return(convert(value, "mpa"))
  }
  return(value)
}

envert <- function(value,units) {
  # Convert values to the specified units
  # from those given in thermo()$opt
  if(!is.numeric(value[1])) return(value)
  units <- tolower(units)
  opt <- get("thermo", CHNOSZ)$opt
  if(units %in% c('c','k','t.units')) {
    if(units=='c' & opt$T.units=='K') return(convert(value,'c'))
    if(units=='k' & opt$T.units=='C') return(convert(value,'k'))
  }
  if(units %in% c('j','cal','e.units')) {
    if(units=='j' & opt$T.units=='Cal') return(convert(value,'j'))
    if(units=='cal' & opt$T.units=='J') return(convert(value,'cal'))
  }
  if(units %in% c('bar','mpa','p.units')) {
    if(units=='mpa' & opt$P.units=='bar') return(convert(value,'mpa'))
    if(units=='bar' & opt$P.units=='MPa') return(convert(value,'bar'))
  }
  return(value)
}

