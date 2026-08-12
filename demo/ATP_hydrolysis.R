# CHNOSZ/demo/ATP_hydrolysis.R
# Calculate Gibbs energy of ATP hydrolysis in the biochemical standard state
# (as a function of pH) and compare with values in Table 4.1 of Alberty (2003)
# 20260810 jmd first version

library(CHNOSZ)

# Define T and pH grid
T <- c(10, 40, 3)
pH <- c(5, 9, 41)
# Colors for T: red, black, blue
col <- c(4, 1, 2)

# Calculate affinity of ADP phosphorylation
results <- phosphorylate("ADP", "P", T = T, pH = pH)

# Convert dimensionless affinities (A/2.303RT) to Gibbs energy (kJ/mol)
TK <- convert(results$m1$A.species$T_grid, "K")
DGr_ADP_P <- convert(results$a12, "G", T = TK) / 1000
# Take the opposite of ADP phosphorylation to get ATP hydrolysis
DGr_ATP_H2O <- - DGr_ADP_P

# Make plot with one line for each T value
ylab <- quote(Delta * italic(G[r]) * degree * "'" ~ "(kJ/mol)")
pH_vals <- results$m1$A.species$vals$pH
matplot(pH_vals, t(DGr_ATP_H2O), type = "l", col = col, lwd = 2, lty = 1, xlab = "pH", ylab = ylab)

# Add points from Table 4.1 of Alberty (2003)
DGr_Alb03_vals <- c(-34.73, -35.38, -36.98, -41.50, -46.93,
                    -35.34, -35.95, -37.64, -42.53, -48.32,
                    -35.95, -36.53, -38.31, -43.56, -49.71)

DGr_Alb03 <- matrix(DGr_Alb03_vals, nrow = 3, byrow = TRUE)
Alb03_pH_vals <- 5:9
matplot(Alb03_pH_vals, t(DGr_Alb03), type = "p", add = TRUE, pch = 19, col = col)

# Calculate mean absolute difference (MAD)
# Extract values of DGr calculated at whole pH values
DGr_pH <- DGr_ATP_H2O[, pH_vals %% 1 == 0]
MAD <- mean(abs(DGr_pH - DGr_Alb03))
text(5, -40, paste("Mean absolute difference: ", round(MAD, 2), "kJ/mol"), adj = 0)

# Add legends and title
T_vals <- results$m1$A.species$vals$T
T_text <- as.expression(sapply(T_vals, function(T) bquote(.(T) ~ degree * C)))
legend("topright", legend = T_text, title = "Alberty (2003)", pch = 19, col = col)
legend("bottomleft", legend = T_text, title = "CHNOSZ::mosaic()", lwd = 2, col = col)
title("ATP hydrolysis: Biochemical standard state", font.main = 1)
