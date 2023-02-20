# Load default settings for CHNOSZ
reset()

info <- "Chemical formulas with unknown elements cause a warning"
expect_warning(makeup("X"), "element\\(s\\) not in thermo\\(\\)\\$element", info = info)

info <- "Unparseable chemical formulas cause an error"
expect_error(makeup("h"), "is not a simple chemical formula", info = info)
expect_error(makeup("1H"), "is not a simple chemical formula", info = info)
expect_error(makeup("(H2O"), "has unpaired parentheses", info = info)
expect_error(makeup("H)2O"), "has unpaired parentheses", info = info)

info <- "Numeric species indices, and coefficients indicating charge can be parsed"
# these are all equivalent formulas for the electron
expect_equal(makeup("-1"), makeup("Z0-1"), info = info)
expect_equal(makeup("-1"), makeup("(Z-1)"), info = info)
expect_equal(makeup("-1"), makeup("Z-1+0"), info = info)
# The species index of the electron in thermo()$OBIGT
ie <- info("e-")
expect_equal(makeup("-1"), makeup(ie), info = info)

info <- "Signed and fractional coefficients can be parsed"
expect_equal(10*makeup("C10.6N1.6P0.1"), makeup("C106N16P1"), info = info)
expect_equal(as.numeric(makeup("C-1.567")), c(1, -1.567), info = info)
expect_equal(as.numeric(makeup("C-1.567+0")), -1.567, info = info)

info <- "Parenthetical and suffixed subformulas can be parsed"
expect_equal(makeup("(H)2O"), makeup("H2O"), info = info)

info <- "Summing and multiply formulas works as expected"
ispecies <- info(c("B(OH)3", "gypsum", "lysine:HCl"))
# The elemental composition all of those, added together
msaved <- as.array(c(B = 1, C = 6, Ca = 1, Cl = 1, H = 22, N = 2, O = 11, S = 1))
expect_equal(makeup(ispecies, sum = TRUE), msaved, info = info)
# The elemental composition in a 1:2:-1 ratio
msaved121 <- as.array(c(B = 1, C = -6, Ca = 2, Cl = -1, H = -4, N = -2, O = 13, S = 2))
expect_equal(makeup(ispecies, c(1,2,-1), sum = TRUE), msaved121, info = info)
expect_error(makeup(ispecies, c(1,2), sum = TRUE), "multiplier does not have .* length = number of formulas", info = info)

info <- "makeup() has a fall-through mechanism for matrices and named objects"
# This series of tests mimics situations encountered in residue.info() via species.basis()
# Ultimately we want to be able to use species.basis() for species indices, formulas or makeups
expect_equal(makeup(makeup("CH4")), makeup("CH4"), info = info)
expect_equal(makeup(list(makeup("CH4"))), list(makeup("CH4")), info = info)
basis("CHNOS")
# We turn the result into a vector using [1, ] so as to drop row names conditionally placed by species.basis
expect_equal(CHNOSZ:::species.basis("CH4")[1, ], CHNOSZ:::species.basis(info("CH4"))[1, ], info = info)
expect_equal(CHNOSZ:::species.basis(makeup("CH4"))[1, ], CHNOSZ:::species.basis("CH4")[1, ], info = info)
# A matrix should be turned into a list
protein <- c("LYSC_CHICK", "RNAS1_BOVIN", "CYC_BOVIN", "MYG_PHYCA", "MYG_HORSE")
pf <- protein.formula(protein)  # a matrix with elements on the columns
basis(protein)          # yup, a basis set made of proteins, just for fun
bmat <- basis.elements()  # a matrix with elements on the columns
expect_equal(makeup(pf)[[1]], makeup(as.chemical.formula(pf)[1]), info = info)
expect_equal(makeup(pf), makeup(bmat), info = info)

info <- "as.chemical.formula() moves charge to the end"
mkp <- makeup("Z-1HCO3")
expect_equal(as.chemical.formula(mkp), "HCO3-1", info = info)  # i.e. not -1HCO3

# TODO: Commented this test because detaching the package isn't always possible 20220131
#info <- "makeup() can process formulas if the package is not attached"
## test added 20200727
#CHNOSZattached <- "CHNOSZ" %in% (.packages())
#if(CHNOSZattached) detach("package:CHNOSZ", unload = TRUE)
#mH2O <- CHNOSZ::makeup("H2O")
#expect_identical(mH2O, c(H = 2, O = 1), info = info)
#if(CHNOSZattached) library(CHNOSZ)
