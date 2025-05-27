# CHNOSZ/demo/rubisco_Zc.R
# Plot Zc of Rubisco vs optimal growth temperature
# 20250510 demo extracted from anintro.Rmd

library(CHNOSZ)

datfile <- system.file("extdata/protein/rubisco.csv", package = "CHNOSZ")
aafile <- system.file("extdata/protein/rubisco_aa.csv", package = "CHNOSZ")
dat <- read.csv(datfile)
aa <- read.csv(aafile)
Topt <- (dat$T1 + dat$T2) / 2
idat <- match(dat$ID, substr(aa$protein, 4, 9))
aa <- aa[idat, ]
ZC <- ZC(protein.formula(aa))
pch <- match(dat$domain, c("E", "B", "A")) - 1
col <- match(dat$domain, c("A", "B", "E")) + 1
plot(Topt, ZC, pch = pch, cex = 2, col = col,
     xlab = expression(list(italic(T)[opt], degree*C)),
     ylab = expression(italic(Z)[C]))
text(Topt, ZC, rep(1:9, 3), cex = 0.8)
abline(v = c(36, 63), lty = 2, col = "grey")
legend("topright", legend = c("Archaea", "Bacteria", "Eukaryota"),
       pch = c(2, 1, 0), col = 2:4, pt.cex = 2)
title("Average oxidation state of carbon in Rubisco vs.\noptimal growth temperature of organisms", font.main = 1)

# The old code for making the SVG image
# - RSVGTipsDevice is no longer on CRAN: https://cran.r-project.org/web/packages/RSVGTipsDevice/index.html
# - I got a compilation error trying to install it with remotes::install_github("cran/RSVGTipsDevice") 20250510

if(require(RSVGTipsDevice, quietly = TRUE)) {
  datfile <- system.file("extdata/protein/rubisco.csv", package = "CHNOSZ")
  aafile <- system.file("extdata/protein/rubisco_aa.csv", package = "CHNOSZ")
  dat <- read.csv(datfile)
  aa <- read.csv(aafile)
  Topt <- (dat$T1 + dat$T2) / 2
  idat <- match(dat$ID, substr(aa$protein, 4, 9))
  aa <- aa[idat, ]
  ZC <- ZC(protein.formula(aa))
  pch <- match(dat$domain, c("E", "B", "A")) - 1
  col <- match(dat$domain, c("A", "B", "E")) + 1
  # Because the tooltip titles in the SVG file are shown by recent browsers,
  # we do not need to draw the tooltips explicitly, so set toolTipMode=0
  devSVGTips("rubisco_Zc.svg", toolTipMode = 0, title = "Rubisco")
  par(cex=1.4)
  # Unfortunately, plotmath can't be used with devSVGTips,
  # so axis labels here don't contain italics.
  plot(Topt, ZC, type = "n", xlab = "T, &#176;C", ylab = "ZC")
  n <- rep(1:9, 3)
  for(i in seq_along(Topt)) {
    # adjust cex to make the symbols look the same size
    cex <- ifelse(pch[i] == 1, 2.5, 3.5)
    points(Topt[i], ZC[i], pch = pch[i], cex = cex, col = col[i])
    URL <- dat$URL[i]
    setSVGShapeURL(URL, target = "_blank")
    setSVGShapeContents(paste0("<title>", dat$species[i], "</title>"))
    text(Topt[i], ZC[i], n[i], cex = 1.2)
  }
  abline(v = c(36, 63), lty = 2, col = "grey")
  legend("topright", legend = c("Archaea", "Bacteria", "Eukaryota"),
         pch = c(2, 1, 0), col = 2:4, cex = 1.5, pt.cex = c(3, 2.3, 3), bty = "n")
  dev.off()
}
