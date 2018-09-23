# plot temperature dependence of log K for some dehydration reactions

# the RSVGTipsDevice package allows us to create an SVG file with
# tooltips and hyperlinks
if(require("RSVGTipsDevice")) {

# because the tooltip titles in the SVG file are shown by recent browsers,
# we do not need to draw the tooltips explicitly, so set toolTipMode=0
devSVGTips("dehydration.svg", toolTipMode=0, title="Dehydration reactions")

# unfortunately, plotmath can't be used with devSVGTips,
# so axis labels here don't contain italics.
T <- seq(1, 175)
plot(range(T), c(-2, 1), type="n", xlab="T, &#176;C", ylab="log K")
title(main="Dehydration reactions")

reactants <- c("[AABB]", "[AABB]", "malate-2", "goethite", "gypsum", "epsomite", "ethanol")
products <- c("[UPBB]", "[PBB]", "fumarate-2", "hematite", "anhydrite", "hexahydrite", "ethylene")
rstate <- c("aq", "cr", "aq", "cr", "cr", "cr", "aq")
pstate <- c("aq", "cr", "aq", "cr", "cr", "cr", "gas")
rcoeff <- c(-1, -1, -1, -2, -0.5, -1, -1)
pcoeff <- c(1, 1, 1, 1, 0.5, 1, 1)
# position and rotation of the names
ilab <- c(140, 120, 60, 60, 20, 120, 120)
srt <- c(10, 29, 25, 12, 13, 20, 35)
# reference and temperature for examples of similar calculations
ex.T <- c(NA, NA, NA, NA, 40, NA, 170)
ex.txt <- c(NA, NA, NA, NA, "cf. Mercury et al., 2001", NA, "Shock, 1993")
ex.doi <- c(NA, NA, NA, NA, "10.1016/S0883-2927(00)00025-1", NA, "10.1016/0016-7037(93)90542-5")

for(i in 1:length(reactants)) {

  # lines
  s <- subcrt(c(reactants[i], products[i], "H2O"),
              c(rstate[i], pstate[i], "liq"),
              c(rcoeff[i], pcoeff[i], 1), T=T)
  lines(T, s$out$logK)

  # points
  if(!is.na(ex.T[i])) {
    URL <- paste0("https://doi.org/", ex.doi[i])
    setSVGShapeURL(URL, target="_blank")
    setSVGShapeContents(paste0("<title>", ex.txt[i], "</title>"))
    # we would use this instead with toolTipMode=1 :
    #setSVGShapeToolTip(title=ex.txt[i])
    points(ex.T[i], s$out$logK[ex.T[i]])
  }

  # names
  for(j in 1:2) {
    formula <- thermo$obigt$formula[s$reaction$ispecies[j]]
    key1 <- thermo$obigt$ref1[s$reaction$ispecies[j]]
    # remove suffix from the key (e.g. "DLH06 [S15]" --> "DLH06")
    key1 <- strsplit(key1, " ")[[1]][1]
    ikey1 <- which(thermo$refs$key==key1)
    URL1 <- thermo$refs$URL[ikey1]
    setSVGShapeURL(URL1, target="_blank")
    setSVGShapeContents(paste0("<title>", paste(formula, s$reaction$state[j]), "</title>"))
    if(j==1) dy <- 0.08 else dy <- -0.03
    if(j==1) dx <- 0 else dx <- 5
    # strip charge from names
    name <- gsub("-.*", "", s$reaction$name[j])
    text(T[ilab[i]] + dx, s$out$logK[ilab[i]] + dy, name, adj=1, srt=srt[i])
    # add a second reference link if needed
    key2 <- thermo$obigt$ref2[s$reaction$ispecies[j]]
    if(!is.na(key2)) {
      ikey2 <- which(thermo$refs$key==key2)
      URL2 <- thermo$refs$URL[ikey2]
      setSVGShapeURL(URL2, target="_blank")
      setSVGShapeContents("<title>2nd reference</title>")
      text(T[ilab[i]] + dx, s$out$logK[ilab[i]] + dy, "(*)", adj=0)
    }
  }

}

# dotted line for logaH2O = 0
abline(h=0, lty=3)
# done!
dev.off()

}
