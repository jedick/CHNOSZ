# formation potential of microbial proteins in colorectal cancer
# based on "bugstab" function in Supporting Information of Dick, 2016
# (https://doi.org/10.7717/peerj.2238)

# to reproduce the calculations in the paper, use superseded data for [Gly] 20190206
add.obigt("OldAA")
# set up graphics device
layout(cbind(matrix(sapply(list(c(1, 2), c(3, 4)), function(x) rep(rep(x, each=3), 3)), nrow=6, byrow=TRUE),
             matrix(rep(c(0, 5, 5, 5, 5, 0), each=4), nrow=6, byrow=TRUE)))
opar <- par(mar=c(3.3, 3.3, 1.5, 1.5), mgp=c(2.1, 0.7, 0), xaxs="i", yaxs="i", las=1, cex=0.9)
# resolution for plots
res <- 500
# basis can be "QEC" or "CHNOS"
basis <- "QEC"
# read bioproject ids, species names
mfile <- system.file("extdata/abundance/microbes.csv", package="CHNOSZ")
bugs <- read.csv(mfile, as.is=TRUE)
# where to keep the locations of healthy zones
healthbugs <- list()
if(basis=="CHNOS") {
  O2 <- c(-85, -65, res)
  H2O <- c(-15, 5, res)
} else {
  O2 <- c(-75, -55, res)
  H2O <- c(-10, 10, res)
}
# the datasets that are considered
datasets <- c("WCQ+12", "ZTV+14", "CTB+14", "FLJ+15")
# colors
blue <- "#4A6FE3"
lightblue <- "#9DA8E2"
neutral <- "#E2E2E2"
lightred <- "#E495A5"
red <- "#D33F6A"
# labels
logfO2lab <- expression(log*italic("f")[O[2]*group("(", italic("g"), ")")])
logaH2Olab <- expression(log*italic("a")[H[2]*O*group("(", italic("liq"), ")")])
for(i in 1:4) {
  ibug <- bugs$study==datasets[i] & (is.na(bugs$logodds) | abs(bugs$logodds) > 0.15)
  # get amino acid compositions
  aafile <- system.file("extdata/protein/microbial.aa.csv", package="CHNOSZ")
  microbial.aa <- read.csv(aafile, stringsAsFactors=FALSE)
  aa <- microbial.aa[match(bugs$bioproject[ibug], microbial.aa$organism), ]
  ip <- add.protein(aa)
  # set up system, calculate relative stabilities
  col <- ifelse(bugs$upcan[ibug], lightred, lightblue)
  basis(basis)
  a <- affinity(O2=O2, H2O=H2O, iprotein=ip, T=37)
  names <- bugs$abbrv[ibug]
  d <- diagram(a, names=names, fill=col, as.residue=TRUE, tplot=FALSE, xlab=logfO2lab, ylab=logaH2Olab, format.names=FALSE)
  if(i==1) title(main="Fecal 16S rRNA", cex.main=1)
  if(i==2) title(main="Fecal metagenome (ZTV+14)", cex.main=1)
  if(i==3) title(main="Co-abundance groups", cex.main=1)
  if(i==4) title(main="Fecal metagenome (FLJ+15)", cex.main=1)
  box()
  label.figure(LETTERS[i], yfrac=0.96, paren=FALSE, font=2, cex=1)
  # store locations of healthy bug zones
  p <- d$predominant
  p[p %in% which(bugs$upcan[ibug])] <- 0
  p[p != 0] <- 1
  healthbugs[[i]] <- p
}
# now show the healthy zones
xs <- seq(O2[1], O2[2], length.out=O2[3])
ys <- seq(H2O[1], H2O[2], length.out=H2O[3])
hhh <- healthbugs[[1]] + healthbugs[[2]] + healthbugs[[3]] + healthbugs[[4]]
image(xs, ys, hhh, col=c(red, lightred, neutral, lightblue, blue), useRaster=TRUE, xlab=logfO2lab, ylab=logaH2Olab)
title(main="Cumulative stability count")
box()
label.figure("E", yfrac=0.96, paren=FALSE, font=2, cex=1)

# reset graphics device to default
par(opar)
layout(1)
# reset thermodynamic database
data(thermo)
