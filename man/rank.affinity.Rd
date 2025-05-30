\encoding{UTF-8}
\name{rank.affinity}
\alias{rank.affinity}
\title{Average ranks of chemical affinities}
\description{
Affinity rankings for groups of species.
}

\usage{
  rank.affinity(aout, groups, rescale = TRUE, percent = FALSE)
}

\arguments{
  \item{aout}{list, output of \code{\link{affinity}}}
  \item{groups}{named list of indices (integer or logical) for species in each group}
  \item{rescale}{logical, rescale average ranks for each group to have the same bounds?}
  \item{percent}{logical, take percentages of average ranks for each group (after rescaling)?}
}

\details{
The affinities for all species are \code{\link{rank}}ed, then the average rank for the species in each group is calculated.
The calculations are applied to each set of conditions individually (i.e., each grid point in the affinity \code{\link{affinity}} calculation).

Unless \code{rescale} is FALSE, the average rank of each group is rescaled to have the same bounds (from 1 to the total number of species).
If \code{percent} is TRUE, the average ranks of groups (after rescaling) are converted to percentages.
}

\value{
The average rankings are inserted into the \code{values} element of \code{aout}, and the names of the groups are inserted into the \code{species} element.
The result can be used by \code{\link{diagram}} to make line plots or predominance diagrams (the predominance fields correspond to the groups with highest average ranking of affinity).
}

\note{
The reaction coefficients in the \code{species} element of the returned value of \code{aout} are not valid.
Because balancing on a basis species (i.e., dividing by its reaction coefficient) would be incorrect, \code{diagram} enforces \code{balance = 1} so that that average rankings are used without further modification.
}

\seealso{
\code{demo("rank.affinity")}
}

\examples{
\dontshow{reset()}# Compare Rubisco proteins from three domains
datfile <- system.file("extdata/protein/rubisco.csv", package = "CHNOSZ")
aafile <- system.file("extdata/protein/rubisco_aa.csv", package = "CHNOSZ")
dat <- read.csv(datfile)
aa <- read.csv(aafile)
groups <- sapply(c("A", "B", "E"), "==", dat$domain, simplify = FALSE)
names(groups) <- c("Archaea", "Bacteria", "Eukaryota")
ip <- add.protein(aa, as.residue = TRUE)
basis("QEC")
aout <- affinity(O2 = c(-74, -66, 100), H2O = c(-4, 4, 100), iprotein = ip)
arank <- rank.affinity(aout, groups = groups)
nspecies <- sapply(groups, sum)
names <- paste0(names(groups), " (", nspecies, ")")
diagram(arank, fill = "terrain", font = 2, names = names, format.names = FALSE)
title("Average rank of affinity for Rubisco in different domains")
}

\concept{Extended workflow}
