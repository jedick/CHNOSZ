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
basis(c("FeO", "SO4-2", "CO3-2", "H2O", "H+", "e-"))
basis("SO4-2", loga_S)
basis("CO3-2", loga_C)
# Start with log(activity of aqueous Fe species) = -4
species(c("Fe+2", "Fe+3", "HFeO2-"), -4)
species(c("pyrrhotite", "pyrite", "hematite", "magnetite", "siderite"), add = TRUE)

# Use two sets of changing basis species:
#   speciate SO4-2, HSO4-, HS-, H2S as a function of Eh and pH
#   speciate CO3-2, HCO3-, CO2 as a function of pH
bases <- list(
  c("SO4-2", "HSO4-", "HS-", "H2S"),
  c("CO3-2", "HCO3-", "CO2")
)

# Make a diagram with dotted lines and aqueous species labels for log(activity of aqueous Fe species) = -4
m4 <- mosaic(bases, pH = pH, Eh = Eh, T = T, P = P)
names.aq <- species()$name; names.aq[4:8] <- ""
diagram(m4$A.species, lty = 3, names = names.aq, col.names = 4)

# Overlay solid lines and mineral labels for log(activity of aqueous Fe species) = -6
species(c("Fe+2", "Fe+3", "HFeO2-"), -6)
m6 <- mosaic(bases, pH = pH, Eh = Eh, T = T, P = P)
names <- species()$name; names[1:3] <- ""
# Adjust labels
srt <- dy <- numeric(length(names))
dy[names == "hematite"] <- -0.4
dy[names == "siderite"] <- 0.2
srt[names == "pyrite"] <- -25
diagram(m6$A.species, add = TRUE, names = names, dy = dy, srt = srt, bold = TRUE)

# Replot the aqueous species labels for stronger contrast
diagram(m4$A.species, lty = 3, names = names.aq, col.names = 4, add = TRUE, fill = NA)

# Show the predominance fields for the sulfur and carbonate basis species
dS <- diagram(m4$A.bases[[1]], italic = TRUE, plot.it = FALSE)
dC <- diagram(m4$A.bases[[2]], italic = TRUE, plot.it = FALSE)
dSC <- mash(dS, dC)
diagram(dSC, lty = 2, col = 8, col.names = 8, add = TRUE, srt = 90)

# Add water lines
water.lines(dSC, lty = 5)

# Add legend and title
TP <- describe.property(c("T", "P"), c(T, P))
SC <- c(
  bquote(sum("S(aq)") == 10^.(loga_S)~m),
  bquote(sum("C(aq)") == 10^.(loga_C)~m)
)
legend1 <- lex(TP, SC)
legend("topright", legend1, bty = "n")

Fe <- c(
  bquote(10^-4~"m Fe"),
  bquote(10^-6~"m Fe")
)
legend2 <- lex(Fe)
legend("bottomleft", legend2, lty = c(3, 1), bty = "n")

title(main = paste("Iron oxides, sulfides, and carbonate in water,",
  "after Garrels and Christ, 1965, Figure 7.21", sep = "\n"), font.main = 1)

## Reset the database if we changed it using mod.OBIGT()
#reset()
