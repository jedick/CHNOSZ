# CHNOSZ/util.args.R
# functions to create argument lists and get name of calling function

### unexported functions ###

# force T and P to equal length
# also looks for the keyword Psat in the value of P and substitutes calculated values of the saturation vapor pressure
TP.args <- function(T=NULL, P=NULL) {
  # keep the [1] here because some functions (e.g. subcrt) will repeat "Psat"
  if(identical(P[1], "Psat")) {
    P <- water("Psat", T, P="Psat")[, 1]
    # water.SUPCRT92 issues its own warnings about 
    # exceeding Psat's temperature limit
    if(get("thermo", CHNOSZ)$opt$water == "IAPWS95")
      if(sum(is.na(P))>0) 
        warning('TP.args: NAs in Psat (likely T > Tc where Tc = 647.096 K)',call.=FALSE)
  }
  if(length(P) < length(T) & !is.null(P)) P <- rep(P, length.out=length(T))
  else if(length(T) < length(P) & !is.null(T)) T <- rep(T, length.out=length(P))
  # something we do here so the SUPCRT water calculations work
  T[T==273.15] <- 273.16
  return(list(T=T, P=P))
}

# make the argument lowercase, then transform a, c, g, and l to aq, gas, cr, and liq
state.args <- function(state=NULL) {
  if(is.null(state) | is.numeric(state[1])) return(state)
  # normalize state arguments
  for(i in 1:length(state)) {
    if(tolower(state[i])=='a') state[i] <- 'aq'
    if(tolower(state[i])=='c') state[i] <- 'cr'
    if(tolower(state[i])=='g') state[i] <- 'gas'
    if(tolower(state[i])=='l') state[i] <- 'liq'
  }
  return(state)
}

caller.name <- function(n=2) {
  # returns the name of the calling function n frames up
  # (n=2: the caller of the function that calls this one)
  # or character() if called interactively
  if(sys.nframe() < n) name <- character()
  else {
    sc <- sys.call(-n)[[1]]
    name <- tryCatch(as.character(sc), error = function(e) character())
    # also return character() if the value from sys.call is
    # the function itself (why does this sometimes happen,
    # e.g. when called from affinity()?)
  }
  return(name)
}
