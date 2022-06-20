\encoding{UTF-8}
\name{rank_affinity}
\alias{rank_affinity}
\title{Normalized Sums of Ranks of Chemical Affinities}
\description{
Affinity ranking for groups of species.
}

\usage{
  rank_affinity(aout, groups)
}

\arguments{
  \item{aout}{list, output of \code{\link{affinity}}}
  \item{groups}{named list of indices (integer or numeric) for species in each group}
}

\details{
The following calculations are applied to each set of conditions (i.e., grid point if \code{\link{affinity}} was called with two variables).
The \code{\link{rank}}s of affinities for all species are first computed.
Then, the ranks for the species in each group are summed and divided by the number of species in that group (this is the normalization step).
}

\value{
The normalized sum of ranks are inserted into the \code{values} element of \code{aout}, and the names of the groups are inserted into the \code{species} element.
The result can be used by \code{\link{diagram}} to show the groups with the highest normalized sum of ranks.
}

\note{
The reaction coefficients in the \code{species} element of the returned value of \code{aout} are not valid.
Because balancing on a basis species (i.e., dividing by its reaction coefficient) would be incorrect, \code{diagram} enforces \code{balance = 1} so that that the normalized sums of ranks are used as-is.
}

\seealso{
\code{demo("rank_affinity")}
}

\examples{
\dontshow{reset()}# Compare Rubisco proteins from three domains
datfile <- system.file("extdata/cpetc/rubisco.csv", package = "CHNOSZ")
fastafile <- system.file("extdata/protein/rubisco.fasta", package = "CHNOSZ")
dat <- read.csv(datfile)
aa <- read.fasta(fastafile)
groups <- sapply(c("A", "B", "E"), "==", dat$domain, simplify = FALSE)
names(groups) <- c("Archaea", "Bacteria", "Eukaryota")
ip <- add.protein(aa, as.residue = TRUE)
basis("QEC")
aout <- affinity(O2 = c(-74, -66, 100), H2O = c(-4, 4, 100), iprotein = ip)
arank <- rank_affinity(aout, groups = groups)
nspecies <- sapply(groups, sum)
names <- paste0(names(groups), " (", nspecies, ")")
diagram(arank, fill = "terrain", font = 2, names = names, format.names = FALSE)
title("Average affinity ranking of Rubisco proteins")
}

\concept{Extended workflow}
