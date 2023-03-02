# CHNOSZ/util.character.R
# Functions to work with character objects

### Unexported functions ###

# Join the elements of a character object into a character object of length 1 (a string)
c2s <- function(x, sep=' ') {
  # Make a string out of a character vector
  if(length(x) %in% c(0,1)) return(x)
  s <- paste(x,collapse=sep)
  return(s)
}

# Split a string into elements of a character object of length n+1, where n is the number of separators in the string
# Default sep=NULL indicates a separator at every position of x
# keep.sep is used to keep the separators in the output
s2c <- function(x,sep=NULL,keep.sep=TRUE) {
  # Recursively split 'x' according to separation strings in 'sep'
  do.split <- function(x,sep,keep.sep=TRUE) {
    # Split the elements of x according to sep
    # Output is a list the length of x
    if(is.list(x)) stop("x is a list; it must be a character object (can have length > 1)")
    x <- as.list(x)
    for(i in 1:length(x)) {
      # Do the splitting
      xi <- strsplit(x[[i]],sep,fixed=TRUE)[[1]]
      # Paste the separation term term back in
      if(keep.sep & !is.null(sep)) {
        xhead <- character()
        xtail <- xi
        if(length(xi) > 1) {
          xhead <- head(xi,1)
          xtail <- tail(xi,-1)
          # In-between matches
          xtail <- paste("",xtail,sep=sep)
        } 
        # A match at the end ... grep here causes problems
        # when sep contains control characters
        #if(length(grep(paste(sep,"$",sep=""),x[[i]]) > 0)) xtail <- c(xtail,sep)
        # Use substr instead
        nx <- nchar(x[[i]])
        ns <- nchar(sep)
        if(substr(x[[i]],nx-ns+1,nx) == sep) xtail <- c(xtail,sep)
        xi <- c(xhead,xtail)
      }
      x[[i]] <- xi
    }
    return(x)
  }
  # Now do it!
  for(i in 1:length(sep)) x <- unlist(do.split(x,sep[i],keep.sep=keep.sep))
  return(x)
}

# Return a value of TRUE or FALSE for each element of x
can.be.numeric <- function(x) {
  # Return FALSE if length of argument is zero
  if(length(x) == 0) FALSE else
  if(length(x) > 1) as.logical(sapply(x, can.be.numeric)) else {
    if(is.numeric(x)) TRUE else
    if(!is.na(suppressWarnings(as.numeric(x)))) TRUE else
    if(x %in% c('.','+','-')) TRUE else FALSE
  }
}
