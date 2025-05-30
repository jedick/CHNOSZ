# CHNOSZ/demo/comproportionation.R
# Gibbs energy of sulfur comproportionation,
# after Fig. 1 of Amend et al., 2020 (doi:10.1111/1462-2920.14982)
# 20191112 jmd first version

library(CHNOSZ)

# Set basis species and activities
basis(c("H2S", "SO4-2", "H2O", "H+"))
basis("H2S", -3)
basis("SO4-2", -2)
# Form native sulfur from sulfide and sulfate
species("S")

# If we calculate the affinity like this, we're stuck with H2S and SO4-2
#a <- affinity(T = c(0, 100), pH = c(0, 7))
# Instead, use mosaic() to speciate H2S/HS- and SO4-2/HSO4-
bases <- list(c("H2S", "HS-"), c("SO4-2", "HSO4-"))
m <- mosaic(bases, T = c(0, 100), pH = c(0, 7))
a <- m$A.species

# Get grid values for variables
T <- a$vals[[1]]
pH <- a$vals[[2]]
# The affinity as a function of T (rows) and pH (columns)
A <- a$values[[1]]
# Get values of Kelvin along the temperature scale
TK <- convert(T, "K")
## NOTE: If T was in the columns of A, we would need to transpose to get T into the rows of A
## (i.e., the first indexed dimension), then transpose again to get back to the original dimensions
#G.J <- t(convert(t(A), "G", T = TK))
# Since T is in the rows of A (first indexed dimension), no transposition is needed here
G.J <- convert(A, "G", T = TK)
G.kJ <- G.J / 1000
# Multiply by 4
# (formation reaction in CHNOSZ is for 1 S; reaction in paper has 4 S)
G.kJ.4 <- G.kJ * 4

# Use subcrt() to write the balanced reaction (shown on the plot)
rxn <- subcrt("S", 1)$reaction
rxn$coeff <- rxn$coeff * 4
rxntext <- describe.reaction(rxn)
# Get label for Delta G (kJ / mol)
DGlab <- axis.label("DGr", prefix = "k")

# Calculate pK of H2S and HSO4-
pK_H2S <- subcrt(c("HS-", "H+", "H2S"), c(-1, -1, 1), T = T)$out$logK
pK_HSO4 <- subcrt(c("SO4-2", "H+", "HSO4-"), c(-1, -1, 1), T = T)$out$logK

# Make contour plot
filled.contour(T, pH, G.kJ.4, xlab = axis.label("T"), ylab = axis.label("pH"),
  levels = -55:0,
  color.palette = ifelse(getRversion() >= "3.6.0", function(n) hcl.colors(n), topo.colors),
  # use plot.axes to label the contour plot (see ?filled.contour)
  plot.axes = {
    contour(T, pH, G.kJ.4, levels = c(-10, -30, -50), add = TRUE, col = "white", lwd = 2, labcex = 0.8)
    legend("topleft", legend = rxntext, bty = "n", inset = c(0, 0.03))
    legend("topleft", describe.basis(1:2), bty = "n", inset = c(0, 0.08))
    lines(T, pK_H2S, lty = 2)
    text(85, 6.7, expr.species("HS-"))
    text(85, 6.3, expr.species("H2S"))
    lines(T, pK_HSO4, lty = 2)
    text(85, 3.0, expr.species("SO4-2"))
    text(85, 2.5, expr.species("HSO4-"))
    axis(1)
    axis(2)
    title("Sulfur comproportionation, after Amend et al., 2020", font.main = 1)
  }
)

# Add legend text
par(xpd = NA)
text(87, 7.3, DGlab)
par(xpd = FALSE)
