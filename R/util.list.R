# CHNOSZ/util.list.R
# functions to work with lists

# which list elements have the maximum (or minimum) values
# revised for speed 20200725
which.pmax <- function(x, maximum = TRUE) {
  # start with NA indices, -Inf (or Inf) working values, and a record of NA values
  iiNA <- tmp <- imax <- x[[1]]
  imax[] <- NA
  if(maximum) tmp[] <- -Inf else tmp[] <- Inf
  iiNA[] <- 0
  # loop over elements of x
  for(i in seq_along(x)) {
    # find values that are greater (or lesser) than working values
    if(maximum) iimax <- x[[i]] > tmp
    else iimax <- x[[i]] < tmp
    # keep NAs out
    iNA <- is.na(iimax)
    iiNA[iNA] <- 1
    iimax[iNA] <- FALSE
    # save the indices and update working values
    imax[iimax] <- i
    tmp[iimax] <- x[[i]][iimax]
  }
  imax[iiNA == 1] <- NA
  # keep attributes from x
  mostattributes(imax) <- attributes(x[[1]])
  imax
}

### unexported functions ###

lsum <- function(x,y) {
  # sum up the respective elements of lists
  # x and y to give list z of the same length
  z <- x
  for(i in 1:length(x)) z[[i]] <- x[[i]] + y[[i]]
  return(z)
}

pprod <- function(x,y) {
  # multiply each element of vector y
  # by corresponding value in list x
  pfun <- function(i) x[[i]]*y[i]
  lapply(1:length(y),pfun)
}
