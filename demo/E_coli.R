# CHNOSZ/demo/E_coli.R
# Calculate Gibbs energy of biomass synthesis in E. coli
# 20210316 jmd version 1

# After LaRowe and Amend (2016): https://doi.org/10.1038/ismej.2015.227
# Polymerization scheme from Amend et al. (2013): https://doi.org/10.1098/rstb.2012.0255

library(CHNOSZ)

# Concentrations of biomolecules (mol (g cell)-1)
# (from Table 1 of LaRowe and Amend, 2016)
concentrations <- c(
  # Amino acids
  alanine = 5.43e-04, arginine = 2.81e-04, asparagine = 2.29e-04, aspartate = 2.29e-04,
  cysteine = 8.70e-05, glutamate = 2.78e-04, glutamine = 2.50e-04, glycine = 5.82e-04,
  histidine = 9.00e-05, isoleucine = 2.76e-04, leucine = 4.28e-04, lysine = 3.26e-04,
  methionine = 1.46e-04, phenylalanine = 1.76e-04, proline = 2.10e-04, serine = 2.05e-04,
  threonine = 2.41e-04, tryptophan = 5.40e-05, tyrosine = 1.31e-04, valine = 4.02e-04,
  # Amines
  ethanolamine = 1.31e-04, `diaminopimelic acid` = 2.79e-05, putrescine = 3.40e-05, spermidine = 6.88e-06,
  # Nucleotides
  `AMP2-` = 1.65E-04, `GMP2-` = 1.26E-04, `CMP2-` = 2.03E-04, `UMP2-` = 1.36E-04,
  `dAMP2-` = 2.46E-05, `dGMP2-` = 2.54E-05, `dCMP2-` = 2.54E-05, `dTMP2-` = 2.46E-05,
  # Fatty acids
  palmitate = 1.12e-04, oleate = 6.22e-05, palmitoleate = 8.56e-05, myristate = 1.67e-05, `beta-hydroxymyristate` = 3.37e-05,
  # Saccharides and more
  glycerol = 1.61e-04, glucose = 2.50e-05, heptose = 2.52e-05, galactose = 8.33e-06, rhamnose = 8.53e-06,
  glucoseamine = 1.67e-05, `N-acetylglucosamine` = 3.62e-05, `N-acetylmuramic acid` = 2.76e-05
)
# Keep the names of the biomolecules here
biomolecules <- names(concentrations)

# Set temperature values
T <- 0:125
T.K <- convert(T, "K")
# Convert Eh to pe for oxidizing and reducing conditions
pe_ox <- convert(0.858, "pe", T = T.K)
pe_red <- convert(-0.384, "pe", T = T.K)
pe <- list(pe_ox, pe_red)

# Parameters for [UPBB] in OBIGT are from Kitadai (2014)
#   (https://doi.org/10.1007/s00239-014-9616-1)
# This command loads "old" parameters for [UPBB]
#   (Dick et al., 2006; https://doi.org/10.5194/bg-3-311-2006)
# - increases G.P278 by ca. 35-40% (closer to Figure 5 of Amend et al., 2013)
add.OBIGT("OldAA")

# Calculate polymerization contribution
# Standard Gibbs energy (J / mol) for AABB -> PBB + H2O
# (Figure 4 of Amend et al., 2013)
E.units("J")
G0.AABB_to_PBB_plus_H2O <- subcrt(c("[AABB]", "[UPBB]", "H2O"), c(-1, 1, 1), T = T)$out$G
# Standard Gibbs energy for 278 AA -> P[278] + 277H2O
G0.P278 <- 277 * G0.AABB_to_PBB_plus_H2O
# logQ for this reaction (decimal logarithm)
logQ <- log10(8.7e-6) - 278 * log10(6.5e-3)
# Gibbs energy for this reaction
# G = G0 + 2.303*RT*logQ
# (cf. Figure 5 of Amend et al., 2013)
G.P278 <- G0.P278 + log(10) * 8.3145 * T.K * logQ
# Gibbs energy (J / peptide bond)
G.P278_per_bond <- G.P278 / 277 / 6.02e23
# Gibbs energy of protein polymerization (J / g cell)
bonds_per_g_cell <- 2.82e21  # Table 2 of Amend et al., 2013
Gpoly_protein_per_g_cell <- G.P278_per_bond * bonds_per_g_cell
# The value calculated at 25 degrees C is equal to that given by Amend et al., 2013
stopifnot(round(Gpoly_protein_per_g_cell[26]) == 191)
# Calculate energy for non-protein polymerization (J / g cell)
Gpoly_nonprotein_per_g_cell <- 45 / 55 * Gpoly_protein_per_g_cell
Gpoly_per_g_cell <- Gpoly_protein_per_g_cell + Gpoly_nonprotein_per_g_cell

# Function to plot Gibbs energy of biomolecule synthesis
# for a given combination of C-, N- and S-bearing basis species
plot_G <- function(C, N, S) {

  # Retrieve logarithm of activity for given basis species
  # (from Table 2 of LaRowe and Amend, 2016)
  loga_C <- switch(C, "CO2" = -3, "CH3COO-" = -5, "CH4" = -6)
  loga_N <- switch(N, "NO3-" = -5, "NH4+" = -6)
  loga_S <- switch(S, "SO4-2" = -3, "HS-" = -6)
  # Set basis species
  # (Note: we set activity of e- in affinity())
  basis(c(C, N, S, "HPO4-2", "H2O", "H+", "e-"), c(loga_C, loga_N, loga_S, -5, 0, -7, 0))
  # Load formed species
  species(biomolecules, -9)

  # Start plot
  ylab <- quote(list(Delta*italic(G[synth]), kJ*(dry~g~cells)^-1))
  plot(c(0, 125), c(-15, 30), xlab = axis.label("T"), ylab = ylab, type = "n", xaxs = "i", yaxs = "i")
  axis(3, labels = FALSE)
  axis(4, labels = FALSE)
  # Loop over oxidizing/reducing conditions
  for(ipe in 1:2) {
    # Calculate dimensionless affinity (A/2.303RT) from 0 to 125 degC at 1 bar
    a <- affinity(T = T, `e-` = -pe[[ipe]])
    # Convert affinity to Gibbs energy (kJ/mol)
    G.cal <- lapply(a$values, convert, "G", T = T.K)
    G.J <- lapply(G.cal, convert, "J")
    G.kJ <- lapply(G.J, "*", 1e-3)
    # Calculate Gibbs energy (kJ (g cell)-1) for each biomolecule
    G.kJ.g_cell <- Map("*", G.kJ, concentrations)
    # Sum Gibbs energy for all biomolecules
    sum.G.kJ.g_cell <- Reduce("+", G.kJ.g_cell)
    # Add polymerization contribution
    total.G.kJ.g_cell <- sum.G.kJ.g_cell + Gpoly_per_g_cell / 1000
    # Add line to plot
    # (Note: ipe * 2 = 2 (red) or 4 (blue))
    lines(T, total.G.kJ.g_cell, col = ipe * 2, lwd = 2)
    # Add label
    dy_ox <- 3
    dy_red <- ifelse(C == "CH4", 3, -3)
    if(ipe == 1) text(T[25], total.G.kJ.g_cell[25] + dy_ox, "Oxidizing")
    if(ipe == 2) text(T[25], total.G.kJ.g_cell[25] + dy_red, "Reducing")
  }
  # Add legend
  Cexpr <- expr.species(C)
  Nexpr <- expr.species(N)
  Sexpr <- expr.species(S)
  legend <- bquote(list(.(Cexpr), .(Nexpr), .(Sexpr)))
  x <- ifelse(C == "CO2", "bottomright", "topright")
  legend(x, legend = legend, bty = "n")

}

# Make plots with different combinations of basis species
par(mfrow = c(3, 2))
par(mar = c(2.5, 3, 1.5, 1), mgp = c(1.5, 0.3, 0))
par(tcl = 0.25)
par(cex = 1)
plot_G("CO2", "NO3-", "SO4-2")
plot_G("CO2", "NH4+", "HS-")
plot_G("CH3COO-", "NO3-", "SO4-2")
plot_G("CH3COO-", "NH4+", "HS-")
plot_G("CH4", "NO3-", "SO4-2")
plot_G("CH4", "NH4+", "HS-")

# Reset CHNOSZ settings (units and OBIGT database)
reset()
