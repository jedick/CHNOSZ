# CHNOSZ/util.character.R
# functions to work with character objects

### unexported functions ###

# join the elements of a character object into a character object of length 1 (a string)
c2s <- function(x, sep=' ') {
  # make a string out of a character vector
  if(length(x) %in% c(0,1)) return(x)
  s <- paste(x,collapse=sep)
  return(s)
}

# split a string into elements of a character object of length n+1, where n is the number of separators in the string
# default sep=NULL indicates a separator at every position of x
# keep.sep is used to keep the separators in the output
s2c <- function(x,sep=NULL,keep.sep=TRUE) {
  # recursively split 'x' according to separation strings in 'sep'
  do.split <- function(x,sep,keep.sep=TRUE) {
    # split the elements of x according to sep
    # output is a list the length of x
    if(is.list(x)) stop("x is a list; it must be a character object (can have length > 1)")
    x <- as.list(x)
    for(i in 1:length(x)) {
      # do the splitting
      xi <- strsplit(x[[i]],sep,fixed=TRUE)[[1]]
      # paste the separation term term back in
      if(keep.sep & !is.null(sep)) {
        xhead <- character()
        xtail <- xi
        if(length(xi) > 1) {
          xhead <- head(xi,1)
          xtail <- tail(xi,-1)
          # in-between matches
          xtail <- paste("",xtail,sep=sep)
        } 
        # a match at the end ... grep here causes problems
        # when sep contains control characters (e.g. protein.refseq)
        #if(length(grep(paste(sep,"$",sep=""),x[[i]]) > 0)) xtail <- c(xtail,sep)
        # use substr instead
        nx <- nchar(x[[i]])
        ns <- nchar(sep)
        if(substr(x[[i]],nx-ns+1,nx) == sep) xtail <- c(xtail,sep)
        xi <- c(xhead,xtail)
      }
      x[[i]] <- xi
    }
    return(x)
  }
  # now do it!
  for(i in 1:length(sep)) x <- unlist(do.split(x,sep[i],keep.sep=keep.sep))
  return(x)
}

# return a value of TRUE or FALSE for each element of x
can.be.numeric <- function(x) {
  # return FALSE if length of argument is zero
  if(length(x) == 0) FALSE else
  if(length(x) > 1) as.logical(sapply(x, can.be.numeric)) else {
    if(is.numeric(x)) TRUE else
    if(!is.na(as.numeric.nowarn(x))) TRUE else
    if(x %in% c('.','+','-')) TRUE else FALSE
  }
}

# something like R's as.numeric(), but without the "NAs introduced by coercion" warnings
# (needed because testthat somehow detects the warnings suppressed by suppressWarnings) 20170427
as.numeric.nowarn <- function(x) {
  if(length(x) == 0) numeric() else
  if(length(x) > 1) sapply(x, as.numeric.nowarn) else
  # http://stackoverflow.com/questions/12643009/regular-expression-for-floating-point-numbers
  if(grepl("^[+-]?([0-9]*[.])?[0-9]+$", x)) as.numeric(x) else NA_real_
}

# convert to integer without NA coercion warnings
as.integer.nowarn <- function(x) {
  if(length(x) == 0) integer() else
  if(length(x) > 1) sapply(x, as.integer.nowarn) else
  if(grepl("[^0-9]", x)) NA_integer_ else as.integer(x)
}
