# CHNOSZ/demo/protbuff.R
# Buffer + ionization: Metastablilities of thiol peroxidases from model bactera
# (ECOLI, BACSU mesophile; AQUAE thermophile, THIDA acidophile, BACHD alkaliphile)

library(CHNOSZ)

basis("CHNOS+")
organisms <- c("ECOLI", "AQUAE", "BACSU", "BACHD", "THIDA")
species("TPX", organisms)
# Create a buffer with our proteins in it
mod.buffer("TPX", paste("TPX", organisms, sep = "_"))
# Set up the buffered activities
basis(c("CO2", "H2O", "NH3", "O2"), "TPX")
a <- affinity(return.buffer = TRUE, T = 50)
basis(c("CO2", "H2O", "NH3", "O2"), as.numeric(a[1:4]))
a <- affinity(pH = c(4, 10, 300), T = c(40, 60, 300))
e <- equilibrate(a, normalize = TRUE)
fill <- ZC.col(ZC(protein.formula(species()$name)))
diagram(e, fill = fill)
title(main = "Thiol peroxidases from bacteria")
legend("topleft", describe.basis(basis  =  thermo()$basis[-6,]), bg = "slategray1", box.lwd = 0)

## Buffer + ionization: relative stabilities
## of E. coli sigma factors on a T-pH diagram
# (sigma factors 24, 32, 38, 54, 70, i.e.
# RpoE, RpoH, RpoS, RpoN, RpoD)
proteins <- c("RPOE", "RP32", "RPOS", "RP54", "RPOD")
basis("CHNOS+")
basis("pH", 7.4)
# Define and set the buffer
mod.buffer("sigma", paste(proteins, "ECOLI", sep = "_"))
basis(c("CO2", "NH3", "H2S", "O2"), "sigma")
logact <- affinity(return.buffer = TRUE, T = 25)
# Set the activities of the basis species to constants 
# corresponding to the buffer, and diagram the relative
# stabilities as a function of T and pH
basis(c("CO2", "NH3", "H2S", "O2"), as.numeric(logact))
species(paste(proteins, "ECOLI", sep = "_"))
a <- affinity(pH = c(5, 10, 300), T = c(10, 40, 300))
fill <- ZC.col(ZC(protein.formula(species()$name)))
diagram(a, normalize = FALSE, fill = fill)
title(main = expression("Sigma factors in"~italic("E. coli")))
ptext <- c(describe.property("T", 25), 
  describe.basis(c(2, 6), oneline = TRUE))
legend("topleft", legend = c("preset conditions:", ptext), bg = "slategray1", box.lwd = 0)
btext <- describe.basis(c(1, 3, 4, 5), oneline = TRUE)
legend("bottomright", legend = c("buffered conditions:", btext), bg = "slategray1", box.lwd = 0)
