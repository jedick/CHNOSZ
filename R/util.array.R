# CHNOSZ/util.array.R
# Functions to work on multidimensional arrays
# 20100314 jmd

list2array <- function(l) {
  # Turn a list of arrays (each with the same dimension) into an array,
  # which has an additional dimension with size equal to the number of initial arrays
  d <- dim(l[[1]])
  if(is.null(d)) d <- length(l[[1]])
  n <- length(l)
  # The size of each of the input arrays
  ni <- prod(d)
  # The total size of the output array
  nt <- prod(c(ni, n))
  # lapply is fast, but it just returns another list
  # and sapply (and unlist) is slow
  #arr <- sapply(l, as.numeric)
  arr <- numeric(nt)
  # Enter values into the vector
  for(i in 1:n) arr[(1:ni)+(i-1)*ni] <- l[[i]]
  # Turn the vector into an array with the required dimensions
  dim(arr) <- c(d, n)
  return(arr)
}

slice <- function(arr, d = NULL, i = 1, value = NULL) {
  # Extract/assign values from/to the ith slice(s) 
  # in the dth dimension of an array
  mydim <- dim(arr)
  nd <- length(mydim) 
  # Build an expression used to index the array 20101031
  prefix <- paste(rep(",", d-1), collapse = "")
  suffix <- paste(rep(",", nd-d), collapse = "")
  # The ith slices of the dth dimension
  expr <- paste(prefix, deparse(i), suffix, sep = "")
  # The old (maybe slower) way
  #expr <- rep("", nd)
  #expr[d] <- "i"
  #expr <- c2s(expr, sep = ",")
  if(is.null(value)) {
    expr <- paste("arr[", expr, "]", sep = "")
    return(eval(parse(text = expr)))
  } else {
    expr <- paste("arr[", expr, "] <- value", sep = "")
    # The following performs the given operation on arr
    eval(parse(text = expr))
    return(arr)
  }
}

dimSums <- function(arr, d = 1, i = NULL) {
  # Sum an n-dimensional array along the dth dimension
  # using only the ith slices in that dimension
  # Merged from mj/plot.R 20100314 jmd
  # e.g., if 'arr' is a matrix,  the following are TRUE:
  # identical(dimSums(arr, 1), colSums(arr))
  # identical(dimSums(arr, 2), rowSums(arr))
  # If i is NULL we use all slices in the dth dimension
  if(is.null(i)) i <- 1:dim(arr)[d]
  # Now take the sum
  for(j in 1:length(i)) {
    s <- slice(arr, d = d, i = i[j])
    if(j == 1) out <- s else out <- out + s
  }
  return(out)
}


