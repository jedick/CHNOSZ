\encoding{UTF-8}
\name{util.blast}
\alias{util.blast}
\alias{read.blast}
\alias{id.blast}
\alias{write.blast}
\alias{def2gi}
\title{Functions to Work with BLAST Output Files}

\description{
  Read and filter BLAST tabular output files, make taxonomic identifications of the BLAST hits using gi numbers, write trimmed-down BLAST files. 
}

\usage{
  read.blast(file, similarity = 30, evalue = 1e-5, max.hits = 1, 
    min.length = NA, quiet = FALSE)
  id.blast(blast, gi.taxid, taxid.names, min.taxon = 0, 
    min.query = 0, min.phylum = 0, take.first = TRUE)
  write.blast(blast, outfile)
  def2gi(def)
}

\arguments{
  \item{file}{character, name of BLAST tabular output file}
  \item{similarity}{numeric, hits above this similarity score are kept}
  \item{evalue}{character, hits below this E value are kept}
  \item{max.hits}{numeric, up to this many hits are kept for each query sequence}
  \item{min.length}{numeric, hits with at least this alignment length are kept}
  \item{quiet}{logical, produce fewer messages?}
  \item{blast}{dataframe, BLAST table}
  \item{gi.taxid}{list, first component is sequence identifiers (gi numbers), second is taxon ids (taxids)}
  \item{taxid.names}{dataframe, with at least columns \samp{taxid} (taxon id), \samp{phylum} (name of phylum), \samp{species} (name of species)}
  \item{min.taxon}{numeric, this taxon is kept if it makes up at least this fraction of total}
  \item{min.query}{numeric, query sequence is counted if a single phylum makes up this fraction of its hits}
  \item{min.phylum}{numeric, this phylum is kept if it makes up at least this fraction of total}
  \item{take.first}{logical, keep only first hit after all other filtering steps?}
  \item{outfile}{character, name of output file}
  \item{def}{character, FASTA defline(s)}
}

\details{

  \code{read.blast} reads a BLAST (Altschul et al., 1997) tabular output \code{file} (such as generated using the -m 8 switch to the \sQuote{blastall} command), keeping only those hits with greater than or equal to \code{similarity} and less than or equal to \code{evalue} (expectation value). Furthermore, for each query sequence, only the top number of hits specified by \code{max.hits} are kept, and only hits with an alignment length of at least \code{min.length} are kept. One or more of these filters can be disabled by setting \code{similarity}, \code{evalue} and/or \code{max.hits} to NA. 

  \code{id.blast} takes a BLAST table (i.e., the output of \code{read.blast}) and finds the taxonomic ID, phylum and species name for each hit (subject sequence). The BLAST results are tied to taxids using \code{gi.taxid}, which is a list consisting of \samp{gi} and \samp{taxid} numeric vectors. Any subject sequence identifiers appearing in the BLAST file that do not match gi numbers in the \code{gi.taxid} list are dropped. The \code{taxid.names} dataframe lists the phylum and species names for each taxid.

  \code{id.blast} furthermore performs three possible filtering steps, which are all disabled by default. If one or more of the arguments is set to a non-zero value, its operation is performed, in this order. Any taxon that does not initially make up at least the fraction of total hits given by \code{min.taxon} is removed. Any query sequence that does not have a single phylum making up at least the fraction of hits (for each query sequence) given by \code{min.query} is removed. Finally, any phylum that does not make up at least the fraction of total hits given by \code{min.phylum} is removed.

  By default, for \code{take.first} equal to TRUE, \code{id.blast} performs a final filtering step (but \code{min.query} must be disabled). Only the first hit for each query sequence is kept.

\code{write.blast} takes a BLAST table (the output of \code{read.blast}) and writes to \code{outfile} a stripped-down BLAST file with empty values in the columns except for columns 1 (query sequence ID), 2 (hit sequence ID), 3 (similarity), 11 (E value).
In the process, \code{\link{def2gi}} is used to extract the GI numbers for the hit sequences that are then kept in the second column.
This function is used to reduce the size of the example BLAST files that are packaged with CHNOSZ (see the \sQuote{bison} section in \code{\link{extdata}}).

  \code{def2gi} extracts the GI number from a FASTA defline. 

}

\value{
  \code{read.blast} returns a dataframe with as many columns (12) as the BLAST file. \code{id.blast} returns a dataframe with columns \code{query}, \code{subject} (i.e., sequence id or gi number), \code{similarity}, \code{evalue}, \code{taxid}, \code{phylum} and \code{species}. \code{write.blast} \code{\link{invisible}}-y returns the results (that are also written to \code{outfile}).
}

\examples{
## using def2gi
def <- "gi|218295810|ref|ZP_03496590.1|"
stopifnot(all.equal(def2gi(def), "218295810"))

## process some of the BLAST output for proteins
## from Bison Pool metagenome (JGI, 2007)
# read the file that connects taxids with the sequence identifier
tfile <- system.file("extdata/bison/gi.taxid.txt.xz", package="CHNOSZ")
gi.taxid <- scan(tfile, what=as.list(character(2)), flush=TRUE)
# read the file that connects names with the taxids
nfile <- system.file("extdata/bison/taxid_names.csv.xz", package="CHNOSZ")
taxid.names <- read.csv(nfile)
# the BLAST files
sites <- c("N","S","R","Q","P")
bfile <- paste("extdata/bison/bison", sites, "_vs_refseq57.blastp.xz", sep="")
for(i in 1:5) {
  file <- system.file(bfile[i], package="CHNOSZ")
  # read the blast file, with default filtering settings
  bl <- read.blast(file)
  # process the blast file -- get taxon names
  ib <- id.blast(bl, gi.taxid, taxid.names, min.taxon=2^-7)
  # count each of the phyla
  bd <- as.matrix(sapply(unique(ib$phylum), function(x) (sum(x==ib$phylum))))
  colnames(bd) <- sites[i]
  # make a matrix -- each column for a different file
  if(i==1) bardata <- bd else {
    bardata <- merge(bardata, bd, all=TRUE, by="row.names")
    rownames(bardata) <- bardata$Row.names
    bardata <- bardata[,-1]
  }
}
# normalize the counts
bardata[is.na(bardata)] <- 0
bardata <- t(t(bardata)/colSums(bardata))
# make a bar chart
bp <- barplot(as.matrix(bardata), col=rainbow(nrow(bardata)),
  xlab="location", ylab="fractional abundance")
# add labels to the bars
names <- substr(row.names(bardata), 1, 3)
for(i in 1:5) {
  bd <- bardata[,i]
  ib <- bd!=0
  y <- (cumsum(bd) - bd/2)[ib]
  text(bp[i], y, names[ib])
}
title(main=paste("Phylum Classification of Protein Sequences",
  "in Part of the Bison Pool Metagenome", sep="\n"))
}

\references{
  Altschul, S. F., Madden, T. L., Schaffer, A. A., Zhang, J. H., Zhang, Z., Miller, W. and Lipman, D. J. (1997) Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. \emph{Nucleic Acids Res.} \bold{25}, 3389--3402. \url{https://doi.org/doi:10.1093/nar/25.17.3389} 
}

\concept{Other tools}
