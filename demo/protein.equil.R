## steps in calculation of chemical activities of two proteins
## in metastable equilibrium, after Dick and Shock, 2011
protein <- pinfo(c("CSG_METVO", "CSG_METJA"))
# clear out amino acid residues loaded by the example above
# ( in affinity(iprotein=ip) )
data(thermo)
# use properties of the "old" [Met] sidechain group (Dick et al., 2006)
mod.obigt("[Met]", G=-35245, H=-59310)
# set up the basis species to those used in DS11
basis("CHNOS+")
# note this yields logaH2 = -4.657486
swap.basis("O2", "H2")
# demonstrate the steps of the equilibrium calculation
protein.equil(protein, loga.protein=-3)
## we can also look at the affinities
# (Reaction 7, Dick and Shock, 2011)
# A/2.303RT for protein at unit activity (A-star for the protein)
a <- affinity(iprotein=protein[1], loga.protein=0)
Astar.protein <- a$values[[1]]
# divide affinity by protein length (A-star for the residue)
pl <- protein.length(protein[1])
Astar.residue <- a$values[[1]]/pl  # 0.1893, Eq. 11
# A/2.303RT per residue corresponding to protein activity of 10^-3
loga.residue <- log10(pl*10^-3)
Aref.residue <- Astar.residue - loga.residue  # 0.446, after Eq. 16
# A-star of the residue in natural log units (A/RT)
log(10) * Astar.residue  # 0.4359, after Eq. 23
# forget about the old [Met] group for whatever comes next
data(thermo)
