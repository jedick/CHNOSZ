# CHNOSZ/demo/bison.R
# average oxidation state of carbon (ZC) of metagenome-derived
# proteins in different microbial phyla at Bison Pool  20171217

# set default xaxs, yaxs, and tcl here becuase they are affected by previous demos
par(mar=c(3.5, 4, 2.5, 2), las=1, mgp=c(2.2, 0.8, 0), xaxs="r", yaxs="r", tcl=-0.5)
# read the amino acid compositions
aa.annot <- read.csv(system.file("extdata/protein/DS11.csv", package="CHNOSZ"), as.is=TRUE)
aa.phyla <- read.csv(system.file("extdata/protein/DS13.csv", package="CHNOSZ"), as.is=TRUE)
sites <- c("N", "S", "R", "Q", "P")
sitenames <- paste("bison", sites, sep="")
# the names of the phyla in alphabetical order (except Deinococcus-Thermus at end)
phyla.abc <- sort(unique(aa.phyla$organism))[c(1:7,9:11,8)]
# an abbreviation for Dein.-Thermus
phyla.abbrv <- phyla.abc
phyla.abbrv[[11]] <- "Dein.-Thermus"
phyla.cols <- c("#f48ba5", "#f2692f", "#cfdd2a",
  "#962272", "#87c540", "#66c3a2", "#12a64a", "#f58656",
  "#ee3237", "#25b7d5", "#3953a4")
# set Chlorobi color to NA because no line is needed (it's identifed only at one location)
phyla.cols[4] <- NA
# chemical formula and ZC of proteins
pf.phyla <- protein.formula(aa.phyla)
ZC.phyla <- ZC(pf.phyla)
# set up plot
plot(0, 0, xlim=c(1, 5), ylim=c(-0.23, -0.14), xlab="distance from source, m", ylab=NA, xaxt="n")
distance <- c(0, 6, 11, 14, 22)
axis(1, at=1:5, labels=distance)
mtext(axis.label("ZC"), side=2, line=3, las=0)
for(i in 1:length(phyla.abc)) {
  # skip Euryarchaeota because it occurs at one location, on top of Dein.-Thermus and Firmicutes
  if(phyla.abc[i]=="Euryarchaeota") next
  # which of the model proteins correspond to this phylum
  iphy <- which(aa.phyla$organism==phyla.abc[i])
  # the locations (of 1, 2, 3, 4, 5) where this phylum is found
  ilocs <- match(aa.phyla$protein[iphy], sitenames)
  # the plotting symbol: determined by alphabetical position of the phylum
  points(ilocs, ZC.phyla[iphy], pch=i-1, cex=1.2)
  # lines to connect the phyla
  lines(ilocs, ZC.phyla[iphy], type="c", col=phyla.cols[i], lwd=2)
}
text(c(4.75, 2.0, 4.0, 4.0, 4.0, 2.0, 3.0, NA, 2.9, 1.3, 3.0),
     c(-0.146, -0.224, -0.161, -0.184, -0.145, -0.201, -0.144, NA, -0.176, -0.158, -0.192),
     phyla.abbrv, cex=0.9)
title(main="Bison Pool hot spring, Yellowstone National Park")
