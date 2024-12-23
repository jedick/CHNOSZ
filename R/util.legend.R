# CHNOSZ/util.legend.R
# Functions for making legend text
# 20190530 jmd first version

lNaCl <- function(x, digits = 2) substitute(NaCl == x~mol~kg^-1, list(x = round(x, digits)))

lS <- function(x, digits = 3) substitute(sum(S) == x~mol~kg^-1, list(x = round(x, digits)))

lT <- function(x, digits = 0) substitute(x~degree*C, list(x = round(x, digits)))

lP <- function(x, digits = 0) if(identical(x, "Psat")) quote(italic(P)[sat]) else substitute(x~bar, list(x = round(x, digits)))

lTP <- function(x, y, digits = 0) substitute(list(x, y), list(x = lT(x, digits), y = lP(y, digits)))

lex <- function(...) as.expression(c(...))
