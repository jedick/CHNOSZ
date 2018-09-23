# integrate Berman's equations using sympy 20170930
# add k4, k5, k6 20180328
library(rSymPy)

# create SymPy variables called T and P
sympy("var('T')")
sympy("var('P')")

sympy("k0 = Symbol('k0', real=True)")
sympy("k1 = Symbol('k1', real=True)")
sympy("k2 = Symbol('k2', real=True)")
sympy("k3 = Symbol('k3', real=True)")
sympy("k4 = Symbol('k4', real=True)")
sympy("k5 = Symbol('k5', real=True)")
sympy("k6 = Symbol('k6', real=True)")
sympy("v1 = Symbol('v1', real=True)")
sympy("v2 = Symbol('v2', real=True)")
sympy("v3 = Symbol('v3', real=True)")
sympy("v4 = Symbol('v4', real=True)")
sympy("Tr = Symbol('Tr', real=True)")
sympy("Pr = Symbol('Pr', real=True)")
sympy("VPrTr = Symbol('VPrTr', real=True)")
sympy("HPrTr = Symbol('HPrTr', real=True)")
sympy("SPrTr = Symbol('SPrTr', real=True)")

# Cp and its integrals
Cp <- "k0 + k1*T**-0.5 + k2*T**-2 + k3*T**-3 + k4*T**-1 + k5*T + k6*T**2"
intCp <- sympy(paste("integrate(", Cp,", (T, Tr, T))"))
message("intCp = ", intCp)

Cp_T <- "k0*T**-1 + k1*T**-1.5 + k2*T**-3 + k3*T**-4 + k4*T**-2 + k5 + k6*T"
intCp_T <- sympy(paste("integrate(", Cp_T,", (T, Tr, T))"))
message("intCp_T = ", intCp_T)

# V and its integrals
# as written in Berman, 1988
#V <- "VPrTr * (1 + v1*(P - Pr) + v2*(P - Pr)**2 + v3*(T- Tr) + v4*(T - Tr)**2)"
# simplified version with Pr==1 (Anderson and Crerar, 1993, p. 175)
V <- "VPrTr * (1 + v1*(P - 1) + v2*(P - 1)**2 + v3*(T - Tr) + v4*(T - Tr)**2)"
intV <- sympy(paste("integrate(", V,", (P, 1, P))"))
# check that intV is equal to the expression in Berman, 1988 and Anderson, 2005 eq. 5.36
#refintV <- "VPrTr * ( (v1/2 - v2) * (P**2 - Pr**2) + v2/3 * (P**3 - Pr**3) + (1 - v1 + v2 + v3*(T-Tr) + v4*(T-Tr)**2) * (P - Pr) )"
# simplified version with Pr==1
refintV <- "VPrTr * ( (v1 / 2 - v2) * (P**2 - 1) + v2 / 3 * (P**3 - 1) + (1 - v1 + v2 + v3*(T-Tr) + v4*(T-Tr)**2) * (P - 1) )"
# this doesn't collect and cancel terms, so we have to expand the expressions
#sympy(paste(intV, "- (", refintV, ")"))
expintV <- sympy(paste("expand (", intV, ")"))
exprefintV <- sympy(paste("expand (", refintV, ")"))
# the difference is zero!
diffintV <- sympy(paste(expintV, "- (", exprefintV, ")"))
message("diffintV = ", diffintV, " (should be 0)")

# another way to check that the expressions are equal
# a function to separate the terms of an expression
septerms <- function(x) {
  x <- paste("+", x)
  x <- gsub("+ ", "+", x, fixed=TRUE)
  x <- gsub("- ", "-", x, fixed=TRUE)
  strsplit(x, " ")[[1]]
}
intVterms <- septerms(expintV)
refintVterms <- septerms(exprefintV)
equalintVterms <- setequal(intVterms, refintVterms)
message("equalintVterms is ", equalintVterms, " (should be TRUE)")

# continuing with V integrals for S and H
dVdT <- sympy(paste("diff(", V, ", T, 1)"))
# again, integrate using Pr==1
intdVdT <- sympy(paste("integrate(", dVdT,", (P, 1, P))"))         # for S
message("intdVdT = ", intdVdT)

minusintdVdT <- sympy(paste("- (", intdVdT, ")"))                  # apply minus sign!
V_TdVdT <- sympy(paste(V, " - T * (", dVdT, ")"))
intV_TdVdT <- sympy(paste("integrate(", V_TdVdT,", (P, 1, P))"))   # for H
message("intV_TdVdT = ", intV_TdVdT)

# check that intV_TdVdT - T * minusintdVdT is equal to intV (for G)
intV2 <- sympy(paste(intV_TdVdT, "- T * (", minusintdVdT, ")") )
expintV2 <- sympy(paste("expand (", intV2, ")"))
## maybe we need to expand the terms first...
#expintV_TdVdT <- sympy(paste("expand (", intV_TdVdT, ")"))
#TintdVdT <- sympy(paste("T * (", intdVdT, ")") )
#expTintdVdT <- sympy(paste("expand (", TintdVdT, ")"))
#expintV3 <- sympy(paste(expintV_TdVdT, "- (", expTintdVdT, ")"))
# this is zero!
diffexpintV <- sympy(paste(expintV2, "- (", expintV, ")"))
message("diffexpintV = ", diffexpintV, " (should be 0)")

# H, S and G
#HCp <- sympy(paste("HPrTr +", intCp, "+", intV_TdVdT))
#SCp <- sympy(paste("SPrTr +", intCp_T, "- (", intdVdT, ")"))
#GCp <- sympy(paste(HCp, "- T * (", SCp, ")"))
#G <- sympy(paste(GCp, "+", intV))

# References
# Anderson, G. M. (2005) \emph{Thermodynamics of Natural Systems}, 2nd ed., Cambridge University Press, 648 p. \url{http://www.worldcat.org/oclc/474880901}
# Anderson, G. M. and Crerar, D. A. (1993) \emph{Thermodynamics in Geochemistry: The Equilibrium Model}, Oxford University Press. \url{http://www.worldcat.org/oclc/803272549}
# Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals in the system Na{\s2}O-K{\s2}O-CaO-MgO-FeO-Fe{\s2}O{\s3}-Al{\s2}O{\s3}-SiO{\s2}-TiO{\s2}-H{\s2}O-CO{\s2}. \emph{J. Petrol.} \bold{29}, 445-522. \url{https://doi.org/10.1093/petrology/29.2.445}

