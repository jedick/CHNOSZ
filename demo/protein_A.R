# CHNOSZ/demo/protein_A.R
# Calculate affinity of ionization of proteins
# 20250510 demo extracted from anintro.Rmd

library(CHNOSZ)

ip <- pinfo(c("CYC_BOVIN", "LYSC_CHICK", "MYG_PHYCA", "RNAS1_BOVIN"))
basis("CHNOS+")
a_ion <- affinity(pH = c(0, 14), iprotein = ip)
basis("CHNOS")
a_nonion <- affinity(iprotein = ip)
plot(c(0, 14), c(50, 300), xlab = "pH", ylab = quote(italic(A/2.303*RT)), type = "n")
for(i in 1:4) {
  A_ion <- as.numeric(a_ion$values[[i]])
  A_nonion <- as.numeric(a_nonion$values[[i]])
  lines(a_ion$vals[[1]], A_ion - A_nonion, col=i)
}
legend("topright", legend = a_ion$species$name,
       col = 1:4, lty = 1, bty = "n", cex = 0.9)
title("Affinity of ionization of proteins", font.main = 1)

# The affinity is always positive, showing the strong energetic drive for ionization of proteins in aqueous solution.
# The degrees of ionization of amino and carboxyl groups increase at low and high pH, respectively, giving rise to the U-shaped lines.
