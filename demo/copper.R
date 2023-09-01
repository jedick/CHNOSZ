# CHNOSZ/demo/copper.R
## Eh-pH diagrams for copper-water-glycine
## After Fig. 2 of Aksu and Doyle, 2001
## (Aksu, S. and Doyle, F. M., 2001. Electrochemistry of copper in aqueous glycine 
## solutions. J. Electrochem. Soc., 148, B51-B57. doi:10.1149/1.1344532)
library(CHNOSZ)

# We need data for Cu-Gly complexes 20190206
add.OBIGT(system.file("extdata/adds/SK95.csv", package = "CHNOSZ"))
# Add some new species to thermo()$OBIGT
m1 <- makeup(info(c("Cu+", "glycinate", "glycinate")), sum = TRUE)
mod.OBIGT(name = "Cu(Gly)2-", formula = as.chemical.formula(m1))
m2 <- makeup(info(c("Cu+2", "glycinate", "H+")), sum = TRUE)
mod.OBIGT(name = "HCu(Gly)+2", formula = as.chemical.formula(m2))
# Names of species in AD01 Table 1 and Table II
Cu_s <- c("copper", "cuprite", "tenorite")
Gly <- c("glycinium", "glycine", "glycinate")
Cu_aq <- c("Cu+", "Cu+2", "CuO2-2", "HCuO2-")
CuGly <- c("Cu(Gly)+", "Cu(Gly)2", "Cu(Gly)2-", "HCu(Gly)+2")
names <- c(Cu_s, Gly, Cu_aq, CuGly)
G <- c(
  # Table I: Gibbs energies in kJ/mol
  c(0, -146, -129.7,
  -384.061, -370.647, -314.833,
  49.98, 65.49, -183.6, -258.5, -298.2)*1000,
  # Table II: Association constants, converted to Gibbs energy
  convert(c(15.64, 10.1, 2.92), "G")
)

# Run updates in order so later species take account of prev. species' values
getG <- function(x) info(info(x))$G
for(i in 1:length(G)) {
  myG <- G[i]
  if(i == 12) myG <- myG + getG("Cu+2") + 2*getG("glycinate")
  if(i == 13) myG <- myG + getG("Cu+") + 2*getG("glycinate")
  if(i == 14) myG <- myG + getG("Cu(Gly)+")
  # Energies are in Joules, so we have to change units of species in default OBIGT 20220325
  mod.OBIGT(names[i], G = myG, E_units = "J")
}  

# In Fig. 2b, total log activities of Cu (Cu_T) and glycine (L_T) are -4 and -1
basis(c("Cu+2", "H2O", "H+", "e-", "glycinium", "CO2"), c(999, 0, 999, 999, -1, 999))
# Add solids and aqueous species
species(Cu_s)
species(c(Cu_aq, CuGly), -4, add = TRUE)
names <- c(Cu_s, Cu_aq, CuGly)
# Mosaic diagram with glycine speciation as a function of pH
m <- mosaic(bases = Gly, pH = c(0, 16, 500), Eh = c(-0.6, 1.0, 500))
fill <- c(rep("lightgrey", 3), rep("white", 4), rep("lightblue", 4))
d <- diagram(m$A.species, fill = fill, names = FALSE, xaxs = "i", yaxs = "i", fill.NA = "pink2", limit.water = TRUE)
# Adjustments for labels
names <- names[sort(unique(as.numeric(d$predominant)))]
for(i in 1:length(names)) {
  if(i %in% 1:3) lab <- names[i] else lab <- expr.species(names[i])
  # Some manual adjustment so labels don't collide
  srt <- dy <- dx <- 0
  if(names[i] == "tenorite") dy <- -0.1
  if(names[i] == "CuO2-2") dy <- -0.1
  if(names[i] == "HCu(Gly)+2") srt <- 90
  if(names[i] == "HCu(Gly)+2") dx <- -0.2
  if(names[i] == "Cu(Gly)+") srt <- 90
  text(na.omit(d$namesx)[i]+dx, na.omit(d$namesy)[i]+dy, lab, srt = srt)
}

# Add glycine ionization lines
d <- diagram(m$A.bases, add = TRUE, col = "darkblue", lty = 3, names = FALSE)
text(d$namesx, -0.5, Gly, col = "darkblue")

# Add water lines and title
water.lines(d)
mtitle(expression("Copper-water-glycine at 25"~degree*"C and 1 bar",
  "After Aksu and Doyle, 2001 (Fig. 2b)"))

# Done!
reset()
