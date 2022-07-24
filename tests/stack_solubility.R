# Make stacked mosaic diagrams and solubility contours
# 20220723 jmd

# - Is solubility() function basically correct?
#   Check: Solubility contour (right column for each stack) coincides with mineral-aqueous stability boundary
#   NOTE: Contour and aqueous species are both set to logarithm of activity = -6

# - Are solubilities of bimetallic minerals (i.e. chalcopyrite) correct?
#   Check: ???
#   TODO: Is there a reference for the correct result?

# - Is an Fe(first)-Cu stack the same as a Cu(first)-Fe stack?
#   Check: Chalcopyrite field has consistent extents in Stack 2 (Diagram 4) and Stack 4 (Diagram 8)
#   TODO: Chalcopyrite field doesn't change, but assemblages with other minerals are different

# - Are solubilities correct for different stoichiometries?
#   Check: Change chalcopyrite (Cu:Fe = 1:1) to bornite (Cu:Fe = 5:1)
#   TODO: Solubility contours for bornite are incorrect

# NOTE: Among these diagrams, I think that Diagrams 4 and 8 are likely to be
# the most accurate depictions of the solubilities of Cu and Fe, respectively
# (Bornite and chalcopyrite should both be present, after the solubility() calculations for bornite are fixed.)

# This is a long test: don't run it for routine package checking
if(FALSE) {

library(CHNOSZ)
reset()

pdf("stack_solubility.pdf", width = 6, height = 9)
par(mfrow = c(4, 2))

# Define system
res <- 200
pH <- c(0, 14, res)
O2 <- c(-48, -33, res)
T <- 200
logmS <- -2
m_NaCl <- 0.1
logm_aq <- -6 # for both Fe- and Cu-bearing aq species
# Define basis species
S.aq <- c("H2S", "HS-", "HSO4-", "SO4-2")
# Define minerals
Fe.cr <- c("pyrite", "pyrrhotite", "magnetite", "hematite")
FeCu.cr <- c("chalcopyrite", "bornite")
Cu.cr <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
# Define aqueous species
iFe.aq <- retrieve("Fe", c("S", "O", "H", "Cl"), "aq")
Fe.aq <- info(iFe.aq)$name
iCu.aq <- retrieve("Cu", c("S", "O", "H", "Cl"), "aq")
Cu.aq <- info(iCu.aq)$name
# Apply NaCl concentration
nacl <- NaCl(T = T, P = "Psat", m_tot = m_NaCl)

### Setup basis species for Fe-Cu stacks
reset()
basis(c("Cu+", "pyrite", "H2S", "oxygen", "H2O", "H+", "Cl-"))
basis("H2S", logmS)
basis("Cl-", log10(nacl$m_Cl))

## Diagram 1: Only Fe-bearing minerals
species(Fe.cr)
# Mosaic with S species as basis species
mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T, IS = nacl$IS)
dFe <- diagram(mFe$A.species, col = 2)

## Diagram 2: Cu- (and FeCu-) bearing minerals and aqueous species
species(c("chalcopyrite", Cu.cr))
species(Cu.aq, -6, add = TRUE)
mFeCu <- mosaic(list(S.aq, Fe.cr), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dFe$predominant))
diagram(mFeCu$A.species)

## Diagram 2a: Overlay single solubility contour
species(c("chalcopyrite", Cu.cr))
sout1 <- solubility(iCu.aq, bases = list(S.aq, Fe.cr), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dFe$predominant))
diagram(sout1, add = TRUE, col = 4, lwd = 2, levels = -6)

mtext("Stack 1: Fe minerals only -> Cu(Fe) minerals and Cu aqueous", adj = 1.1, line = 1.1)

## Diagram 3: Fe-bearing minerals and aqueous species
species(Fe.cr)
species(Fe.aq, -6, add = TRUE)
mFe <- mosaic(S.aq, pH = pH, O2 = O2, T = T, IS = nacl$IS)
dFe <- diagram(mFe$A.species, col = 2)

## Diagram 4: Cu- (and FeCu-) bearing minerals and aqueous species
species(c("chalcopyrite", Cu.cr))
species(Cu.aq, -6, add = TRUE)
mFeCu <- mosaic(list(S.aq, c(Fe.cr, Fe.aq)), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dFe$predominant), loga_aq = c(NA, logm_aq))
diagram(mFeCu$A.species)

## Diagram 4a: Overlay single solubility contour
species(c("chalcopyrite", Cu.cr))
sout2 <- solubility(iCu.aq, bases = list(S.aq, c(Fe.cr, Fe.aq)), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dFe$predominant), loga_aq = c(NA, logm_aq))
diagram(sout2, add = TRUE, col = 4, lwd = 2, levels = -6)

mtext("Stack 2: Fe minerals and aqueous -> Cu(Fe) minerals and Cu aqueous", adj = 1, line = 1.1)

### Setup basis species for Cu-Fe stacks
reset()
basis(c("Fe+2", "copper", "H2S", "oxygen", "H2O", "H+", "Cl-"))
basis("H2S", logmS)
basis("Cl-", log10(nacl$m_Cl))

## Diagram 5: Only Cu-bearing minerals
species(Cu.cr)
mCu <- mosaic(S.aq, pH = pH, O2 = O2, T = T, IS = nacl$IS)
dCu <- diagram(mCu$A.species, col = 4)

## Diagram 6: Fe- (and CuFe-) bearing minerals and aqueous species
species(c("chalcopyrite", Fe.cr))
species(Fe.aq, -6, add = TRUE)
mCuFe <- mosaic(list(S.aq, Cu.cr), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dCu$predominant))
diagram(mCuFe$A.species)

## Diagram 6a: Overlay single solubility contour
species(c("chalcopyrite", Fe.cr))
sout3 <- solubility(iFe.aq, bases = list(S.aq, Cu.cr), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dCu$predominant))
diagram(sout3, add = TRUE, col = 2, lwd = 2, levels = -6)

mtext("Stack 3: Cu minerals only -> Fe(Cu) minerals and Fe aqueous", adj = 1.1, line = 1.1)

## Diagram 7: Cu-bearing minerals and aqueous species
species(Cu.cr)
species(Cu.aq, -6, add = TRUE)
mCu <- mosaic(S.aq, pH = pH, O2 = O2, T = T, IS = nacl$IS)
dCu <- diagram(mCu$A.species, col = 4)

## Diagram 8: Fe- (and CuFe-) bearing minerals and aqueous species
species(c("chalcopyrite", Fe.cr))
species(Fe.aq, -6, add = TRUE)
mCuFe <- mosaic(list(S.aq, c(Cu.cr, Cu.aq)), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dCu$predominant), loga_aq = c(NA, logm_aq))
diagram(mCuFe$A.species)

## Diagram 8a: Overlay single solubility contour
species(c("chalcopyrite", Fe.cr))
sout4 <- solubility(iFe.aq, bases = list(S.aq, c(Cu.cr, Cu.aq)), pH = pH, O2 = O2, T = T, IS = nacl$IS, stable = list(NULL, dCu$predominant), loga_aq = c(NA, logm_aq))
diagram(sout4, add = TRUE, col = 2, lwd = 2, levels = -6)

mtext("Stack 4: Cu minerals and aqueous -> Fe(Cu) minerals and Fe aqueous", adj = 1, line = 1.1)

dev.off()

}
