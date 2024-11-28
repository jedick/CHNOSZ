# CHNOSZ/demo/sum_S.R
# 20240604 Fe-S-C-O-H log a (sum S) - pH      modified from mosaic.R
# 20241127 Fe-S-O-H   logfO2 - log a (sum S)  added IS

library(CHNOSZ)

for(IS in c(0, 0.5)) {

  # Define conditions
  res <- 500
  loga_S <- c(-6, 0, res)
  logf_O2 <- c(-45, -20, res)
  pH <- 5
  T <- 300
  P <- "Psat"
  #IS <- 0
  # Setup chemical system
  # Use "oxygen" instead of "O2" to get the gas
  basis(c("FeO", "SO4-2", "H2O", "H+", "oxygen"))
  basis("pH", pH)
  # Add aqueous species followed by minerals
  species(c("Fe+2", "Fe+3", "HFeO2-"))
  species(c("pyrrhotite", "pyrite", "hematite", "magnetite"), add = TRUE)

  # List basis species to swap through
  bases <- list( c("SO4-2", "HSO4-", "HS-", "H2S") )

  # Start diagram for loga_Fe = -6
  species(1:3, -6)
  m1 <- mosaic(bases, "SO4-2" = loga_S, O2 = logf_O2, T = T, P = P, IS = IS)
  diagram(m1$A.species, bold = TRUE)

  # Overlay lines and labels for loga_Fe = -5
  species(1:3, -5)
  m2 <- mosaic(bases, "SO4-2" = loga_S, O2 = logf_O2, T = T, P = P, IS = IS)
  d <- diagram(m2$A.species, add = TRUE, names = NA, fill = NA, lty = 3)

  # Add dashed line for water stability limti
  water.lines(d)

  # Add legends and title
  TP <- describe.property(c("T", "P"), c(T, P))
  PI <- c(bquote(pH == .(pH)), bquote(italic(I) == .(IS)~mol~kg^-1))
  legend1 <- lex(TP, PI)
  legend("topright", legend1, bty = "n")

  Fe <- c(
    bquote(10^-5~italic(m)~Fe^"+2"),
    bquote(10^-6~italic(m)~Fe^"+2")
  )
  legend2 <- lex(Fe)
  legend("topleft", legend2, lty = c(3, 1), bty = "n")

  title(main = "Fe-S-O-H with summed S and ionic strength", font.main = 1)

}
