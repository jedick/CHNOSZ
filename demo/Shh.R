# CHNOSZ/demo/Shh.R
# Compare affinities of Sonic hedgehog and transcription factors involved in dorsal-ventral patterning
# (Dick, 2015. Chemical integration of proteins in signaling and development. https://doi.org/10.1101/015826)

library(CHNOSZ)

# To reproduce the calculations in the paper, use superseded data for [Gly] and [UPBB] 20190206
mod.OBIGT("[Gly]", G = -6075, H = -5570, S = 17.31)
mod.OBIGT("[UPBB]", G = -21436, H = -45220, S = 1.62)

# UniProt names of the proteins
pname <- c("SHH", "OLIG2", "NKX22", "FOXA2", "IRX3", "PAX6", "NKX62", "DBX1",
  "DBX2", "NKX61", "PAX7", "GLI1", "GLI2", "GLI3", "PTC1", "SMO", "GLI3R")[1:11]

# Colors modified from Dessaud et al., 2008 and Hui and Angers, 2011
fill <- c(SHH = "#c8c7c8", OLIG2 = "#f9a330", NKX22 = "#ef2b2c", FOXA2 = "#6ab0b9",
  NKX61 = "#b76775", DBX2 = "#35bcba", PAX7 = "#f6ef42",
  PAX6 = "#4d2a59", IRX3 = "#63c54e", NKX62 = "#f24e33", DBX1 = "#d4e94e",
  PTC1 = "#c7b0ee", SMO = "#8fd4ef", GLI1 = "#fcdc7e", GLI2 = "#c7e3b0", GLI3 = "#fcdc7e",
  GLI3R = "#f1b1ae")

# Names for plotting
names <- c(SHH = "Shh", OLIG2 = "Olig2", NKX22 = "Nkx2.2", FOXA2 = "Foxa2",
  NKX61 = "Nkx6.1", DBX2 = "Dbx2", PAX7 = "Pax7",
  PAX6 = "Pax6", IRX3 = "Irx3", NKX62 = "Nkx6.2", DBX1 = "Dbx1",
  PTC1 = "Ptch1", SMO = "Smo", GLI1 = "Gli1A", GLI2 = "Gli2A", GLI3 = "Gli3A",
  GLI3R = "Gli3R")

# Protein indices of Shh and the transcription factors
ip <- match(pname, thermo()$protein$protein)
aa <- thermo()$protein[ip, ]

# Set up basis species
basis("CHNOS")
basis("NH3", -7)

# Save as PDF file?
pdf <- FALSE
# Draw interpretive legend?
interp <- TRUE

# Set up ranges of logfO2 and logaH2O
O2 <- seq(-70, -100, length.out = 500)
H2O <- seq(0.5, -4.5, length.out = 500)
A <- affinity(H2O = H2O, O2 = O2, iprotein = ip)
# Plot affinities per residue, compared to SHH
pl <- protein.length(ip)
names(A$values) <- pname
for(i in 1:length(A$values)) A$values[[i]] <- A$values[[i]] / pl[i]
A.SSH <- A$values$SHH
for(i in 1:length(A$values)) A$values[[i]] <- A$values[[i]] - A.SSH
ylab <- expression(bold(A)/2.303*italic(RT)*" vs Shh")
xlab <- expression(log*italic(a)[H[2]][O])
# Set up normal plot, or plot with interpretive drawings
opar <- par(no.readonly = TRUE)
par(mar = c(5.1, 4.1, 4.1, 2.1))
if(interp) {
  if(pdf) pdf("tfactor_interp.pdf", width = 6, height = 6)
  plot.new()
  plot.window(rev(range(H2O)), c(-0.7, 4.7), xaxs = "i", yaxs = "i")
  par(xpd=FALSE)
  clip(range(H2O)[1], range(H2O)[2], -0.3, 4.2)
} else {
  if(pdf) pdf("tfactor_affinity.pdf", width = 6, height = 6)
  plot(range(H2O), c(-0.5, 4.5), type = "n", xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, xlim = rev(range(H2O)), mgp = c(2.5, 1, 0))
}
for(i in 1:length(pl)) {
  lty <- 3
  lwd <- 1
  # Highlight SHH with solid line
  if(pname[i] == "SHH") {
    lty <- 1
    lwd <- 3
  }
  lines(H2O, A$values[[i]], lty = lty, lwd = lwd)
}
# Highlight lines for reaction sequence from OLIG2 to SHH
A <- A$values
names(A) <- pname
# Olig2
iOLIG2 <- A$OLIG2 > A$SHH
lines(H2O[iOLIG2], A$OLIG2[iOLIG2], col = fill["OLIG2"], lwd = 3)
# Foxa2 with offset to distinguish it from Nkx6.1/Dbx2
iFOXA2 <- A$FOXA2 > A$OLIG2 & A$FOXA2 > A$NKX22
lines(H2O[iFOXA2], A$FOXA2[iFOXA2], col = fill["FOXA2"], lwd = 3)
# Nkx2.2
iNKX22 <- A$NKX22 > A$FOXA2 & A$NKX22 > A$SHH
lines(H2O[iNKX22], A$NKX22[iNKX22], col = fill["NKX22"], lwd = 3)
# Pax6
iPAX6 <- A$PAX6 > A$SHH
lines(H2O[iPAX6], A$PAX6[iPAX6], col = fill["PAX6"], lwd = 3)
# Nkx6.1 with overstepping
iNKX61 <- A$NKX61 > A$DBX2
imax <- max(which(iNKX61))
iNKX61[1: (imax-20)] <- FALSE
lines(H2O[iNKX61], A$NKX61[iNKX61], col = fill["NKX61"], lwd = 3)
# Dbx2
iDBX2 <- A$DBX2 > A$NKX61 & A$DBX2 < A$IRX3
lines(H2O[iDBX2], A$DBX2[iDBX2], col = fill["DBX2"], lwd = 3)
# Irx3
iIRX3 <- A$IRX3 > A$NKX62 & A$IRX3 > A$OLIG2
lines(H2O[iIRX3], A$IRX3[iIRX3], col = fill["IRX3"], lwd = 3)
# Nkx6.2
iNKX62 <- A$NKX62 > A$IRX3 & A$NKX62 > A$DBX1
lines(H2O[iNKX62], A$NKX62[iNKX62], col = fill["NKX62"], lwd = 3)
# Dbx1
iDBX1 <- A$DBX1 > A$NKX62 & A$DBX1 > A$SHH
lines(H2O[iDBX1], A$DBX1[iDBX1], col = fill["DBX1"], lwd = 3)
# Shh
#iSHH <- A$SHH > A$DBX1
#lines(H2O[iSHH], A$SHH[iSHH], col = fill["SHH"], lwd = 3)
# The lines need names
if(interp) {
  # Remove plot clip region
  par(xpd = NA)
  text(-2.12, -0.48, "Nkx2.2", srt = 90)
  text(-0.87, -0.65, "Olig2", srt = 90)
  text(0, -0.72, "Pax6", srt = 90)
  text(0.06, 4.3, "Olig2", srt = 90)
  text(-0.13, 4.3, "Irx3", srt = 90)
  text(-0.77, 4.3, "Nkx6.1", srt = 90)
  text(-.97, 4.3, "Dbx2", srt = 90)
  text(-3.45, 4.3, "Nkx6.2", srt = 90)
  text(-3.65, 4.3, "Dbx1", srt = 90)
} else {
  text(-1.5, 0.5, "Nkx2.2")
  text(-0.1, 0.3, "Pax6")
  text(-0.2, 3.8, "Olig2")
  text(-0.75, 2.45, "Irx3")
  text(-1.13, 1.3, "Nkx6.1")
  text(-1.5, 1, "Dbx2")
  text(-2.4, 1.3, "Nkx6.2")
  text(-3.8, 0.3, "Dbx1")
}
text(0.3, -0.15, "Shh")
text(-4.25, 0.15, "Shh")
text(-0.47, 0.5, "Pax7", srt = -35)
text(-0.22, 2, "Foxa2", srt = -61)
# Are we making an interpretive plot?
if(interp) {
  # The left,bottom x,y-position and horizontal width of the bottom gradient wedge
  xbot <- H2O[1]
  ybot <- -1.4
  wbot <- 3
  # The height of the bottom gradient wedge as a function of the x position
  hbot <- function(x) 0.3 + 0.08*(xbot - x)
  lines(c(xbot, xbot-wbot), ybot+hbot(c(xbot, xbot-wbot)))
  # Draw the base and sides of the bottom gradient wedge
  lines(c(xbot, xbot-wbot), c(ybot, ybot))
  lines(rep(xbot, 2), c(ybot, ybot+hbot(xbot)))
  lines(rep(xbot, 2)-wbot, c(ybot, ybot+hbot(xbot-wbot)))
  # Draw drop lines from reactions between TFs and Shh
  xNKX22 <- H2O[max(which(A$NKX22 > A$SHH))]
  lines(rep(xNKX22, 2), c(ybot, 0), lty = 2)
  xOLIG2 <- H2O[max(which(A$OLIG2 > A$SHH))]
  lines(rep(xOLIG2, 2), c(ybot, 0), lty = 2)
  xPAX6 <- H2O[max(which(A$PAX6 > A$SHH))]
  lines(rep(xPAX6, 2), c(ybot, 0), lty = 2)
  # The left,bottom x,y-position and horizontal width of the top gradient wedge
  xtop <- H2O[2]
  ytop <- 4.7
  wtop <- 5
  # The height of the top gradient wedge as a function of the x position
  htop <- function(x) 0.4 - 0.08*(xtop - x)
  lines(c(xtop, xtop-wtop), ytop+htop(c(xtop, xtop-wtop)))
  # Draw the base and sides of the top gradient wedge
  lines(c(xtop, xtop-wtop), c(ytop, ytop))
  lines(rep(xtop, 2), c(ytop, ytop+htop(xtop)))
  lines(rep(xtop, 2)-wtop, c(ytop, ytop+htop(xtop-wtop)))
  # Draw drop lines to reactions between TFs
  iIRX3 <- min(which(A$IRX3 > A$OLIG2))
  lines(rep(H2O[iIRX3], 2), c(A$IRX3[iIRX3], ytop+htop(H2O[iIRX3])), lty = 2)
  iDBX2 <- min(which(A$DBX2 > A$NKX61))
  lines(rep(H2O[iDBX2], 2), c(A$DBX2[iDBX2], ytop+htop(H2O[iDBX2])), lty = 2)
  iDBX1 <- min(which(A$DBX1 > A$NKX62))
  lines(rep(H2O[iDBX1], 2), c(A$DBX1[iDBX1], ytop+htop(H2O[iDBX1])), lty = 2)
  # Indicate plot variables
  arrows(-2.7, 2, -2.7, 2.5, 0.1)
  text(-2.8, 2.3, "affinity\nvs. Shh", adj = 0)
  arrows(-2.7, 2, -2.27, 2, 0.1)
  text(-2.3, 1.8, expression(list(log*italic(f)[O[2]], )), adj = 0)
  text(-2.3, 1.6, expression(list(log*italic(a)[H[2]*O])), adj = 0)
  # Label neural progenitors
  text(-2.37, -1.55, "FP")
  text(-1.6, -1.55, "p3")
  text(-0.53, -1.55, "pMN")
  text(0.2, -1.55, "p2")
  text(0.2, 5.25, "pMN")
  text(-0.47, 5.25, "p2")
  text(-2.2, 5.25, "p1")
  text(-4, 5.25, "p0")
  # Label Shh gradient and stages
  text(1.2, -1.2, "Shh\ngradient", adj = 0)
  text(1.2, 4.9, "Shh\ngradient", adj = 0)
  text(-2.6, -0.6, expression(italic("Stage 1: loading")), adj = 0)
  text(-1.5, 4.3, expression(italic("Stage 2: unloading")), adj = 0)
} else {
  # Add second axis: logfO2
  pu <- par("usr")
  pu[1:2] <- rev(range(O2))
  par(usr = pu)
  axis(3, at = seq(-75, -105, by = -5))
  mtext(expression(log*italic(f)[O[2]]), line = 2)
}
# All done!
par(opar)
if(pdf) dev.off()
reset()
