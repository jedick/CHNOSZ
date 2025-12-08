# CHNOSZ/demo/phosphorylate.R

# Make T-pH and P-pH diagrams for Gibbs energy of phosphorylation reactions,
# using mosaic() to account for pH-dependent speciation 20251207

# Figure based on LaRowe and Dick (2025) doi:10.1029/2025JG009095

# adenosine_for_RNA uses only protonated species (i.e., no AMP-2 or PO4-3)
# to model RNA formation from monophosphate nucleotides

phospho_plot("adenosine_for_RNA", "P")
title(sub = "Only protonated species to model RNA formation (LaRowe and Dick, 2025)", line = 2.2, xpd = NA)
