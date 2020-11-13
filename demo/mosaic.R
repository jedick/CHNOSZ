# CHNOSZ/demo/mosaic.R
# 20141221 first version jmd
# 20200819 revision:
#   - comment out GC65 thermodynamic data
#   - use mash() to show S and C species together
#   - reorder plot layers to make better use of transparency
#   - add 'dy' argument to adjust positions of labels
#   - add legend to show activity of aqueous Fe species

library(CHNOSZ)

# Fe-bearing minerals and aqueous species in the Fe-S-O-H-C system
# after Figure 7.21 of Garrels and Christ, 1965

## Uncomment these lines to use thermodynamic data from Appendix 2 of GC65 (at 25 degrees C only)
#mod.OBIGT(c("Fe+2", "Fe+3"), G = c(-20300, -2520))
#mod.OBIGT(c("hematite", "magnetite", "pyrrhotite", "pyrite", "siderite"), G = c(-177100, -242400, -23320, -36000, -161060))
#mod.OBIGT(c("SO4-2", "HS-", "H2S", "HSO4-"), G = c(-177340, 3010, -6540, -179940))
#mod.OBIGT(c("CO2", "HCO3-", "CO3-2"), G = c(-92310, -140310, -126220))

# Define conditions
res <- 500
pH <- c(0, 14, res)
Eh <- c(-1, 1, res)
T <- 25
P <- 1
loga_S <- -6
loga_C <- 0
# Define chemical system
basis(c("FeO", "SO4-2", "H2O", "H+", "e-", "CO3-2"))
basis("SO4-2", loga_S)
basis("CO3-2", loga_C)
species(c("Fe+2", "Fe+3"), -4)
species(c("pyrrhotite", "pyrite", "hematite", "magnetite", "siderite"), add = TRUE)

# Use two sets of changing basis species:
#   speciate SO4-2, HSO4-, HS-, H2S as a function of Eh and pH
#   speciate CO3-2, HCO3-, CO2 as a function of pH
bases <- c("SO4-2", "HSO4-", "HS-", "H2S")
bases2 <- c("CO3-2", "HCO3-", "CO2")

# Make a diagram with log(activity of aqueous Fe species) = -4
m4 <- mosaic(bases, bases2, pH = pH, Eh = Eh, T = T, P = P)
diagram(m4$A.species, lty = 2, names = FALSE)

## Show the predominance fields for the sulfur and carbonate basis species
dS <- diagram(m4$A.bases, italic = TRUE, plot.it = FALSE)
dC <- diagram(m4$A.bases2, italic = TRUE, plot.it = FALSE)
dSC <- mash(dS, dC)
diagram(dSC, lty = 3, col = 4, col.names = 4, add = TRUE)

# Show lines for log(activity of aqueous Fe species) = -6
s6 <- species(c("Fe+2", "Fe+3"), -6)
m6 <- mosaic(bases, bases2, pH = pH, Eh = Eh, T = T, P = P)
srt <- dy <- numeric(nrow(s6))
dy[s6$name == "Fe+2"] <- 0.15
dy[s6$name == "hematite"] <- -0.4
dy[s6$name == "magnetite"] <- 0.07
dy[s6$name == "siderite"] <- 0.2
srt[s6$name == "pyrite"] <- -25
bold <- s6$state == "cr"
diagram(m6$A.species, add = TRUE, dy = dy, srt = srt, bold = bold)

# Add legend and title
TP <- describe.property(c("T", "P"), c(T, P))
SC <- c(
  bquote(log * italic(a)["S(aq)"] == .(loga_S)),
  bquote(log * italic(a)["C(aq)"] == .(loga_C))
)
Fe <- c(
  bquote(log * italic(a)["Fe(aq)"] == .(-4)),
  bquote(log * italic(a)["Fe(aq)"] == .(-6))
)
legend <- lex(TP, SC, Fe)
legend("topright", legend, lty = c(NA, NA, NA, NA, 1, 2), bty = "n")
title(main = paste("Iron minerals, sulfur, and carbonate in water,",
  "after Garrels and Christ, 1965, Figure 7.21", sep = "\n"), font.main = 1)

## Reset the database if we changed it using mod.OBIGT()
#reset()
