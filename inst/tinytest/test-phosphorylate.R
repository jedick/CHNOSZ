# 20251207

# First add thermodynamic data for sugar phosphates (from Table 3.2 of Alberty, 2003)
mod.OBIGT("glucose-6-phosphate-2", formula = "C6H11O9P-2", G = -1763940)
mod.OBIGT("glucose-6-phosphate-1", formula = "C6H12O9P-", G = -1800590)
# Alberty (2003) doesn't have ΔG° for neutral glucose-6-phosphate,
# so we calculate it from pKa1 = 1.5 (Degani and Halmann, 1966)
DG0_G6P <- -1800590 + convert(1.5, "G")
mod.OBIGT("glucose-6-phosphate", formula = "C6H13O9P", G = DG0_G6P)

# Calculate affinity at pH 5-9 with unit activities (loga = 0)
pH <- 5:9
result <- phosphorylate("glucose", "ATP", pH = pH)

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
