## Steps in calculation of chemical activities of two proteins
## in metastable equilibrium, after Dick and Shock, 2011
library(CHNOSZ)

protein <- pinfo(c("CSG_METVO", "CSG_METJA"))
# Use superseded properties of [Met], [Gly], and [UPBB] (Dick et al., 2006)
mod.OBIGT("[Met]", G = -35245, H = -59310, S = 40.38)
mod.OBIGT("[Gly]", G = -6075, H = -5570, S = 17.31)
mod.OBIGT("[UPBB]", G = -21436, H = -45220, S = 1.62)
# Set up the basis species to those used in DS11
basis("CHNOS+")
# Note this yields logaH2 = -4.657486
swap.basis("O2", "H2")
# Demonstrate the steps of the equilibrium calculation
protein.equil(protein, loga.protein = -3)
## We can also look at the affinities
# (Reaction 7, Dick and Shock, 2011)
# A/2.303RT for protein at unit activity (A-star for the protein)
a <- affinity(iprotein = protein[1], loga.protein = 0)
Astar.protein <- a$values[[1]]
# Divide affinity by protein length (A-star for the residue)
pl <- protein.length(protein[1])
Astar.residue <- a$values[[1]]/pl  # 0.1893, Eq. 11
# A/2.303RT per residue corresponding to protein activity of 10^-3
loga.residue <- log10(pl*10^-3)
Aref.residue <- Astar.residue - loga.residue  # 0.446, after Eq. 16
# A-star of the residue in natural log units (A/RT)
log(10) * Astar.residue  # 0.4359, after Eq. 23

# Forget about the superseded group properties for whatever comes next
reset()
