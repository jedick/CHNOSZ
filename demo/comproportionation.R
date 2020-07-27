# CHNOSZ/demo/comproportionation.R
# Gibbs energy of sulfur comproportionation,
# after Fig. 1 of Amend et al., 2020 (doi:10.1111/1462-2920.14982)
# 20191112 jmd first version
library(CHNOSZ)

# set basis species and activities
basis(c("H2S", "SO4-2", "H2O", "H+"))
basis("H2S", -3)
basis("SO4-2", -2)
# form native sulfur from sulfide and sulfate
species("S")

# if we calculate the affinity like this, we're stuck with H2S and SO4-2
#a <- affinity(T = c(0, 100), pH = c(0, 7))
# instead, use mosaic() to speciate H2S/HS- and SO4-2/HSO4-
bases <- list(c("H2S", "HS-"), c("SO4-2", "HSO4-"))
m <- mosaic(bases, T = c(0, 100), pH = c(0, 7))
a <- m$A.species

# get plot values
T <- a$vals[[1]]
pH <- a$vals[[2]]
# the affinity as a function of T (rows) and pH (columns)
A <- a$values[[1]]
# convert dimensionless affinity (A/2.303RT) to delta G (cal)
TK <- convert(T, "K")
G.cal <- convert(A, "G", T = TK)
# convert cal to kJ
G.J <- convert(G.cal, "J")
G.kJ <- G.J / 1000
# multiply by 4
# (formation reaction in CHNOSZ is for 1 S; reaction in paper has 4 S)
G.kJ.4 <- G.kJ * 4

# use subcrt() to write the balanced reaction (shown on the plot)
rxn <- subcrt("S", 1)$reaction
rxn$coeff <- rxn$coeff * 4
rxntext <- describe.reaction(rxn)
# set units to get label for Delta G (kJ / mol)
E.units("J")
DGlab <- axis.label("DGr", prefix = "k")

# calculate pK of H2S and HSO4-
pK_H2S <- subcrt(c("HS-", "H+", "H2S"), c(-1, -1, 1), T = T)$out$logK
pK_HSO4 <- subcrt(c("SO4-2", "H+", "HSO4-"), c(-1, -1, 1), T = T)$out$logK

# make contour plot
filled.contour(T, pH, G.kJ.4, xlab = axis.label("T"), ylab = axis.label("pH"),
  levels = -55:0,
  color.palette = function(n) hcl.colors(n),
  # use plot.axes to label the contour plot (see ?filled.contour)
  plot.axes = {
    contour(T, pH, G.kJ.4, levels = c(-10, -30, -50), add = TRUE, col = "white", lwd = 2, labcex = 0.8)
    legend("topleft", legend = rxntext, bty = "n", inset = c(0, 0.03))
    legend("topleft", describe.basis(ibasis = 1:2), bty = "n", inset = c(0, 0.08))
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

# add legend text
par(xpd = NA)
text(87, 7.3, DGlab)
par(xpd = FALSE)
