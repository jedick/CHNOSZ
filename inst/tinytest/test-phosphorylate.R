# Tests added on 20251207

# Calculate affinity for glucose phosphorylation at pH 5-9 with unit activities (loga = 0)
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
expect_equal(round(deltaG_kJ[3], 2), -19.75, info = info)

info <- "Using const_pH argument gives expected value"
result <- phosphorylate("glucose", "ATP", const_pH = 7)
A <- result$a12
TK <- convert(25, "K")
deltaG_J <- convert(A, "G", T = TK)
deltaG_kJ <- deltaG_J / 1000
expect_equal(round(as.numeric(deltaG_kJ), 2), -19.75, info = info)

# Tests added on 20260107

info <- "phospho.plot() doesn't error"
# Save plot to a temporary PNG file
pngfile <- tempfile(fileext = ".png")
png(pngfile)
expect_silent(phospho.plot("adenosine_to_AMP", "PP", loga_Mg = -5, res = 50), info = info)
# Close the graphics device and remove the temporary PNG file
dev.off()
file.remove(pngfile)

## Tests for more species

info <- "Acetate phosphorylation"
# This is the dimensionless affinity (A/2.303RT) of the reaction
affinity_nodim <- as.numeric(phosphorylate("acetic acid", "P")$a12)
# This is a circular reference value (calculated with CHNOSZ, not from an outside source)
expect_equal(affinity_nodim, -2.35, tolerance = 0.01, scale = 1, info = info)

info <- "Adenosine phosphorylation (to AMP)"
affinity_nodim <- as.numeric(phosphorylate("adenosine_to_AMP", "P")$a12)
# Reference value from p. 295 of Alberty (2003)
expect_equal(convert(affinity_nodim, "G"), 12969.6, tolerance = 500, scale = 1, info = info)

info <- "Glycerol phosphorylation"
# Use Delta G0 from Table 3.2 of Alberty (2003)
mod.OBIGT("3-phosphoglycerol", E_units = "J", G = 0)
mod.OBIGT("3-phosphoglycerol-1", E_units = "J", G = -1397040)
mod.OBIGT("3-phosphoglycerol-2", E_units = "J", G = -1358960)
affinity_nodim <- as.numeric(phosphorylate("glycerol", "P")$a12)
# Reference value from p. 295 of Alberty (2003)
expect_equal(convert(affinity_nodim, "G"), -1780.75, tolerance = 1200, scale = 1, info = info)

info <- "Pyruvate phosphorylation with ATP"
mod.OBIGT("phosphoenolpyruvate", E_units = "J", G = 0)
mod.OBIGT("phosphoenolpyruvate-1", E_units = "J", G = 0)
mod.OBIGT("phosphoenolpyruvate-2", E_units = "J", G = -1303610)
mod.OBIGT("phosphoenolpyruvate-3", E_units = "J", G = -1263650)
affinity_nodim <- as.numeric(phosphorylate("pyruvic acid", "ATP")$a12)
# Reference value from p. 224 of Alberty (2003)
expect_equal(convert(affinity_nodim, "G"), 30622.4, tolerance = 4000, scale = 1, info = info)

info <- "Error produced with unknown reactant or P_source"
# reactant should be pyruvic acid, not pyruvate
expect_error(phosphorylate("pyruvate", "acetylphosphate"), "unrecognized reactant", info = info)
# acetyphosphate is missing an "l"
expect_error(phosphorylate("pyruvic acid", "acetyphosphate"), "unrecognized P_source", info = info)
