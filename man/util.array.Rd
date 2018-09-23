\encoding{UTF-8}
\name{util.array}
\alias{util.array}
\alias{list2array}
\alias{slice}
\alias{dimSums}
\alias{slice.affinity}
\title{Functions to Work with Multidimensional Arrays}
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
# start with a matrix
x <- matrix(1:12,ncol=3)
# pay attention to the following when
# writing examples that test for identity!
identical(1*x,x)   # FALSE
# create two matrices that are multiples of the first
a <- 1*x
b <- 2*a
# these both have two dimensions of lengths 4 and 3
dim(a)  # 4 3
# combine them to make an array with three dimensions
c <- list2array(list(a,b))
# the third dimension has length 2
dim(c)  # 4 3 2
# the first slice of the third dimension == a
stopifnot(identical( slice(c,3), a ))
# the second slice of the third dimension == b
stopifnot(identical( slice(c,3,2), b ))
# 'slice' works just like the bracket operator
c11 <- slice(c,1)
c12 <- slice(c,1,2)
c21 <- slice(c,2,1)
c212 <- slice(c,2,1:2)
stopifnot(identical( c11, c[1,,] ))
stopifnot(identical( c12, c[2,,] ))
stopifnot(identical( c21, c[,1,] ))
stopifnot(identical( c212, c[,1:2,] ))
# let us replace part of the array
d <- slice(c,3,2,value=a)
# now the second slice of the third dimension == a
stopifnot(identical( slice(d,3,2), a ))
# and the sum across the third dimension == b
stopifnot(identical( dimSums(d,3), b ))
# taking the sum removes that dimension
dim(d)             # 4 3 2
dim(dimSums(d,1))  # 3 2
dim(dimSums(d,2))  # 4 2
dim(dimSums(d,3))  # 4 3

# working with an 'affinity' object
\dontshow{data(thermo)}
basis("CHNOS+")
species("alanine")
a1 <- affinity(O2=c(-80,-60))   # at pH=7
a2 <- affinity(O2=c(-80,-60),pH=c(0,14,7))
# in the 2nd dimension (pH) get the 4th slice (pH=7)
a3 <- slice.affinity(a2,2,4)
stopifnot(all.equal(a1$values,a3$values))
}

\concept{Utility functions}
