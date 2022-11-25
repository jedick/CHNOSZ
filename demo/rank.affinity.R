## Affinity ranking for proteins coded by differentially expressed
## genes in response to carbon limitation in yeast
# Demo written on 20220620
# Adapted from an example previously in read.expr.Rd (CHNOSZ 1.1.0) and yeast.Rd (CHNOSZ 1.2.0-1.3.0)
library(CHNOSZ)

# Experimental data are from Tai et al. (2005)
# https://doi.org/10.1074/jbc.M410573200
file <- system.file("extdata/protein/TBD+05.csv", package = "CHNOSZ")
dat <- read.csv(file, row.names = 1, check.names = FALSE)

# The activities of ammonium and sulfate are similar to the
# non-growth-limiting concentrations used by Boer et al. (2003)
# https://doi.org/10.1074/jbc.M209759200 
basis(c("glucose", "H2O", "NH4+", "oxygen", "SO4-2", "H+"),
  c(-1, 0, -1.3, 999, -1.4, -7))

aafile <- system.file("extdata/protein/TBD+05_aa.csv", package = "CHNOSZ")
aa <- read.csv(aafile)
iprotein <- add.protein(aa, as.residue = TRUE)

res <- 200
aout <- affinity(C6H12O6 = c(-60, -20, res), O2 = c(-72, -60, res), iprotein = iprotein)
groups <- apply(dat, 2, which)
names(groups) <- paste0(names(groups), "\n(", colSums(dat), ")")
arank <- rank.affinity(aout, groups)
fill <- c("#d2b48c", "#b0e0e6", "#d3d3d3", "#d8bfd8")
diagram(arank, format.names = FALSE, fill = fill)

title(main = paste("Affinity ranking for proteins coded by differentially expressed\n",
  "genes under carbon limitation in yeast (data from Tai et al., 2005)"), font.main = 1)
