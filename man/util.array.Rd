\encoding{UTF-8}
\name{util.array}
\alias{util.array}
\alias{list2array}
\alias{slice}
\alias{dimSums}
\alias{slice.affinity}
\title{Functions to work with multidimensional arrays}
\description{
  These functions can be used to turn a list into an array and extract or replace values or take the sum along a certain dimension of an array.
}

\usage{
  list2array(l)
  slice(arr, d = NULL, i = 1, value = NULL)
  dimSums(arr, d = 1, i = NULL)
  slice.affinity(affinity, d = 1, i = 1)
}

\arguments{
  \item{l}{a list.}
  \item{arr}{an array.}
  \item{d}{numeric, what dimension to use.}
  \item{i}{numeric, what slice to use.}
  \item{value}{values to assign to the portion of an array specified by \code{d} and \code{i}.}
  \item{affinity}{list, output from \code{\link{affinity}} function.}
}

\details{
  \code{list2array} turns a list of \code{\link{array}}s, each with the same dimensions, into a new array having one more dimension whose size is equal to the number of initial arrays.

  \code{slice} extracts or assigns values from/to the \code{i}th slice(s) in the \code{d}th dimension of an array. Values are assigned to an array if \code{value} is not NULL. This function works by building an expression containing the extraction operator (\code{\link{[}}).

  \code{slice.affinity} performs a slice operation on the \samp{values} element of the \samp{affinity} variable (which should be the output of \code{\link{affinity}}).

  \code{dimSums} sums an array along the \code{d}th dimension using only the \code{i}th slices in that dimension. If \code{i} is NULL, all slices in that dimension are summed together. For matrices, \code{dimSums(x,1)} has the same result as \code{\link{colSums}(x)} and \code{dimSums(x,2)} has the same result as \code{\link{rowSums}(x)}.
}

\examples{
# Start with a matrix
mat <- matrix(1:12, ncol = 3)
# Pay attention to the following when
# writing examples that test for identity!
identical(1 * mat, mat)   # FALSE
# Create two matrices that are multiples of the first
a <- 1 * mat
b <- 2 * mat
# These both have two dimensions of lengths 4 and 3
dim(a)  # 4 3
# Combine them to make an array with three dimensions
x <- list2array(list(a, b))
# The third dimension has length 2
dim(x)  # 4 3 2
# The first slice of the third dimension
slice(x, 3)    # a
# The second slice of the third dimension
slice(x, 3, 2) # b
# 'slice' works just like the bracket operator
slice(x, 1)        # x[1, , ]
slice(x, 1, 2)     # x[2, , ]
slice(x, 2, 1)     # x[, 1, ]
slice(x, 2, 1:2)   # x[, 1:2, ]

# Replace part of the array
y <- slice(x, 3, 2, value = a)
# Now the second slice of the third dimension == a
slice(y, 3, 2)  # a
# and the sum across the third dimension == b
dimSums(y, 3)   # b
# Taking the sum removes that dimension
dim(y)              # 4 3 2
dim(dimSums(y, 1))  # 3 2
dim(dimSums(y, 2))  # 4 2
dim(dimSums(y, 3))  # 4 3

# Working with an 'affinity' object
\dontshow{reset()}
basis("CHNOS+")
species("alanine")
a1 <- affinity(O2 = c(-80, -60)) # i.e. pH = 7
a2 <- affinity(O2 = c(-80, -60), pH = c(0, 14, 7))
# In the 2nd dimension (pH) get the 4th slice (pH = 7)
a3 <- slice.affinity(a2, 2, 4)
all.equal(a1$values, a3$values)  # TRUE
}

\concept{Utility functions}
