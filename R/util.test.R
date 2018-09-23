# CHNOSZ/util.test.R
# functions for writing tests 20171005

# the maximum (absolute) pairwise difference between x and y
maxdiff <- function(x, y) max(abs(y - x))

# modelled after the "expect_" functions in testthat ...
# is the maximum of the pairwise differences between two objects less than some value?
expect_maxdiff <- function(object, expected, maxdiff = 0) {
  if(!"testthat" %in% row.names(installed.packages())) {
    stop("please install the 'testthat' package to use this function")
  } else {
    # get the names for the object and expected values
    # we use a double substitute here because
    # as.character(substitute(x$a$b)) ==  c("$", "x$a", "b")
    # as.character(substitute(substitute(x$a$b)))[2] == "x$a$b" --> better for printing
    lab_act <- as.character(substitute(substitute(object)))[2]
    lab_exp <- as.character(substitute(substitute(expected)))[2]
    truemd <- maxdiff(object, expected)
    failmsg <- sprintf("maxdiff(%s, %s) not less than %s", lab_act, lab_exp, maxdiff)
    truemsg <- sprintf("actual value: %s", truemd)
    testthat::expect(maxdiff(object, expected) <= maxdiff, sprintf("%s\n%s", failmsg, truemsg))
    invisible(object)
  }
}
