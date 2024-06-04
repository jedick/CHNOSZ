# CHNOSZ/demo/total_S.R
# 20240604 Use total S as a variable; modified from mosaic.R

library(CHNOSZ)

# Define conditions
res <- 500
pH <- c(0, 14, res)
loga_S <- c(-6, -2, res)
T <- 25
P <- 1
loga_C <- 0
Eh <- -0.5
# Define chemical system
basis(c("FeO", "SO4-2", "CO3-2", "H2O", "H+", "e-"))
basis("CO3-2", loga_C)
basis("Eh", Eh)
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
m4 <- mosaic(bases, pH = pH, "SO4-2" = loga_S, T = T, P = P)
names.aq <- species()$name; names.aq[4:8] <- ""
diagram(m4$A.species, lty = 3, names = names.aq, col.names = 4)

# Overlay solid lines and mineral labels for log(activity of aqueous Fe species) = -6
species(c("Fe+2", "Fe+3", "HFeO2-"), -6)
m6 <- mosaic(bases, pH = pH, "SO4-2" = loga_S, T = T, P = P)
names <- species()$name; names[1:3] <- ""
diagram(m6$A.species, add = TRUE, names = names, bold = TRUE)

# Replot the aqueous species labels for stronger contrast
diagram(m4$A.species, lty = 3, names = names.aq, col.names = 4, add = TRUE, fill = NA)

# Add legend and title
TP <- describe.property(c("T", "P"), c(T, P))
SC <- c(
  bquote(Eh == .(Eh) ~ V),
  bquote(sum("C(aq)") == 10^.(loga_C)~m)
)
legend1 <- lex(TP, SC)
legend("topleft", legend1, bty = "n")

Fe <- c(
  bquote(10^-4~"m Fe"),
  bquote(10^-6~"m Fe")
)
legend2 <- lex(Fe)
legend("bottomleft", legend2, lty = c(3, 1), bty = "n")

title(main = "Using total S as a variable; modified from mosaic.R", font.main = 1)

## Reset the database if we changed it using mod.OBIGT()
#reset()
