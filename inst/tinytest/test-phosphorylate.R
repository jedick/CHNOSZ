# Test added on 20251207

# First add thermodynamic data for sugar phosphates (from Table 3.2 of Alberty, 2003)
mod.OBIGT("glucose-6-phosphate-2", formula = "C6H11O9P-2", G = -1763940)
mod.OBIGT("glucose-6-phosphate-1", formula = "C6H12O9P-", G = -1800590)
# Alberty (2003) doesn't have ΔG° for neutral glucose-6-phosphate,
# so we calculate it from pKa1 = 1.5 (Degani and Halmann, 1966)
DG0_G6P <- -1800590 + convert(1.5, "G")
mod.OBIGT("glucose-6-phosphate", formula = "C6H13O9P", G = DG0_G6P)

# Calculate affinity at pH 5-9 with unit activities (loga = 0)
pH <- 5:9
info <- "phosphorylate() runs without error"
expect_silent(result <- phosphorylate("glucose", "ATP", pH = pH), info = info)

# Extract the overall affinity
A <- result$a12

# Convert affinity to Gibbs energy in kJ/mol
TK <- convert(25, "K")
deltaG_J <- convert(A, "G", T = TK)
deltaG_kJ <- deltaG_J / 1000

info <- "ΔG°' becomes more negative at higher pH"
expect_true(all(diff(deltaG_kJ) < 0), info = info)

info <- "ΔG°' has expected value at pH 7"
expect_equal(round(deltaG_kJ[3], 2), -32.15, info = info)

info <- "Using const_pH argument gives expected value"
result <- phosphorylate("glucose", "ATP", const_pH = 7)
A <- result$a12
TK <- convert(25, "K")
deltaG_J <- convert(A, "G", T = TK)
deltaG_kJ <- deltaG_J / 1000
expect_equal(round(as.numeric(deltaG_kJ), 2), -32.15, info = info)

# Tests added on 20260107

info <- "phospho_plot() doesn't error"
# Save plot to a temporary PNG file
pngfile <- tempfile(fileext = ".png")
png(pngfile)
expect_silent(phospho_plot("adenosine_to_AMP", "PP", loga_Mg = -5, res = 50), info = info)
# Close the graphics device and remove the temporary PNG file
dev.off()
file.remove(pngfile)

## Tests for more species

info <- "Acetate phosphorylation"
# Delta G0 from Table 3.2 of Alberty (2003)
mod.OBIGT("acetylphosphate0", formula = "C2H5O5P", G = -1298260)
mod.OBIGT("acetylphosphate-1", formula = "C2H4O5P-", G = -1268080)
mod.OBIGT("acetylphosphate-2", formula = "C2H3O5P-2", G = -1219390)
# Alberty doesn't have acetylphosphate-3, so we just use an arbitrarily high number
mod.OBIGT("acetylphosphate-3", formula = "C2H2O5P-3", G = 0)
# This is the dimensionless affinity (A/2.303RT) of the reaction
affinity_nodim <- as.numeric(phosphorylate("acetic acid", "P")$a12)
# This is a circular reference value (calculated with CHNOSZ, not from an outside source)
expect_equal(affinity_nodim, -6.20, tolerance = 0.01, scale = 1, info = info)

info <- "Glycerol phosphorylation"
mod.OBIGT("3-phosphoglycerol", formula = "C3H9O6P", G = 0)
mod.OBIGT("3-phosphoglycerol-1", formula = "C3H8O6P-", G = -1397040)
mod.OBIGT("3-phosphoglycerol-2", formula = "C3H7O6P-2", G = -1358960)
affinity_nodim <- as.numeric(phosphorylate("glycerol", "P")$a12)
# Reference value from p. 295 of Alberty (2003)
expect_equal(convert(affinity_nodim, "G"), -1780.75, tolerance = 1200, scale = 1, info = info)

info <- "Adenosine phosphorylation (to AMP)"
affinity_nodim <- as.numeric(phosphorylate("adenosine_to_AMP", "P")$a12)
# Reference value from p. 295 of Alberty (2003)
expect_equal(convert(affinity_nodim, "G"), 12969.6, tolerance = 500, scale = 1, info = info)

info <- "Adenosine phosphorylation (model 2)"
# adenosine_for_RNA excludes PO4-3 and AMP-2 (polymerization model from LaRowe and Dick, 2025)
affinity_nodim <- as.numeric(phosphorylate("adenosine_for_RNA", "P")$a12)
# For this let's use the same reference value as above, but with a greater expected difference
expect_equal(convert(affinity_nodim, "G"), 12969.6, tolerance = 7000, scale = 1, info = info)

info <- "Adenosine phosphorylation (to cAMP)"
# Delta G0 are placeholder values
mod.OBIGT("cyclic-HAMP", formula = "C10H12N5O6P", G = 0)
mod.OBIGT("cyclic-AMP-1", formula = "C10H11N5O6P-", G = 0)
affinity_nodim <- as.numeric(phosphorylate("adenosine_to_cAMP", "P")$a12)
# Circular reference
expect_equal(affinity_nodim, -149.20, tolerance = 0.01, scale = 1, info = info)

info <- "Ribose phosphorylation"
affinity_nodim <- as.numeric(phosphorylate("ribose", "P")$a12)
# Circular reference
expect_equal(affinity_nodim, -3.75, tolerance = 0.01, scale = 1, info = info)

info <- "Uridine phosphorylation"
affinity_nodim <- as.numeric(phosphorylate("uridine", "P")$a12)
# Circular reference
expect_equal(affinity_nodim, -2.24, tolerance = 0.01, scale = 1, info = info)

info <- "AMP phosphorylation"
affinity_nodim <- as.numeric(phosphorylate("AMP", "P")$a12)
# Circular reference
expect_equal(affinity_nodim, -5.83, tolerance = 0.01, scale = 1, info = info)

info <- "ADP phosphorylation"
affinity_nodim <- as.numeric(phosphorylate("ADP", "P")$a12)
# Circular reference
expect_equal(affinity_nodim, -6.52, tolerance = 0.01, scale = 1, info = info)

info <- "Pyruvate phosphorylation with ATP"
mod.OBIGT("phosphoenolpyruvate", formula = "C3H5O6P", G = 0)
mod.OBIGT("phosphoenolpyruvate-1", formula = "C3H4O6P-", G = 0)
mod.OBIGT("phosphoenolpyruvate-2", formula = "C3H3O6P-2", G = -1303610)
mod.OBIGT("phosphoenolpyruvate-3", formula = "C3H2O6P-3", G = -1263650)
affinity_nodim <- as.numeric(phosphorylate("pyruvic acid", "ATP")$a12)
# Reference value from p. 224 of Alberty (2003)
expect_equal(convert(affinity_nodim, "G"), 30622.4, tolerance = 4000, scale = 1, info = info)

info <- "Pyruvate phosphorylation with acetylphosphate"
affinity_nodim <- as.numeric(phosphorylate("pyruvic acid", "acetylphosphate")$a12)
# Circular reference
expect_equal(affinity_nodim, -5.06, tolerance = 0.01, scale = 1, info = info)

info <- "Error produced with unknown reactant or P_source"
# reactant should be pyruvic acid, not pyruvate
expect_error(phosphorylate("pyruvate", "acetylphosphate"), "unrecognized reactant", info = info)
# acetyphosphate is missing an "l"
expect_error(phosphorylate("pyruvic acid", "acetyphosphate"), "unrecognized P_source", info = info)
