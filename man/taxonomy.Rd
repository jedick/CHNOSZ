\encoding{UTF-8}
\name{taxonomy}
\alias{taxonomy}
\alias{getnodes}
\alias{getrank}
\alias{parent}
\alias{allparents}
\alias{getnames}
\alias{sciname}
\title{Extract Data from NCBI Taxonomy Files}

\description{
  Read data from NCBI taxonomy files, traverse taxonomic ranks, get scientific names of taxonomic nodes.
}

\usage{
  getnodes(taxdir)
  getrank(id, taxdir, nodes=NULL)
  parent(id, taxdir, rank=NULL, nodes=NULL)
  allparents(id, taxdir, nodes=NULL)
  getnames(taxdir)
  sciname(id, taxdir, names=NULL)
}

\arguments{
  \item{taxdir}{character, directory where the taxonomy files are kept.}
  \item{id}{numeric, taxonomic ID(s) of the nodes of interest.}
  \item{nodes}{dataframe, output from \code{getnodes} (optional).}
  \item{rank}{character, name of the taxonomic rank of interest.}
  \item{names}{dataframe, output from \code{getnames} (optional).}
}

\details{

  These functions provide a convenient way to read data from NCBI taxonomy files (i.e., the contents of \code{taxdump.tar.gz}, which is available from \url{https://ftp.ncbi.nih.gov/pub/taxonomy/}). 

  The \code{taxdir} argument is used to specify the directory where the \code{nodes.dmp} and \code{names.dmp} files are located. \code{getnodes} and \code{getnames} read these files into data frames. \code{getrank} returns the rank (species, genus, etc) of the node with the given taxonomic \code{id}. \code{parent} returns the taxonomic ID of the next-lowest node below that specified by the \code{id} in the argument, unless \code{rank} is supplied, in which case the function descends the tree until a node with that rank is found. \code{allparents} returns all the taxonomic IDs of all nodes between that specified by \code{id} and the root of the tree, inclusive. \code{sciname} returns the scientific name of the node with the given \code{id}.
  
  The \code{id} argument can be of length greater than 1 except for \code{allparents}. If \code{getrank}, \code{parent}, \code{allparents} or \code{sciname} need to be called repeatedly, the operation can be hastened by supplying the output of \code{getnodes} in the \code{nodes} argument and/or the output of \code{getnames} in the \code{names} argument.

}

\examples{
## get information about Homo sapiens from the
## packaged taxonomy files
taxdir <- system.file("extdata/taxonomy",package="CHNOSZ")
# H. sapiens' taxonomic id
id1 <- 9606
# that is a species
getrank(id1,taxdir)
# the next step up the taxonomy
id2 <- parent(id1,taxdir)
print(id2)
# that is a genus
getrank(id2,taxdir)
# that genus is "Homo"
sciname(id2,taxdir)
# we can ask what phylum is it part of?
id3 <- parent(id1,taxdir,"phylum")
# answer: "Chordata"
sciname(id3,taxdir)
# H. sapiens' complete taxonomy
id4 <- allparents(id1,taxdir)
sciname(id4,taxdir)

## the names of the organisms in the supplied taxonomy files
taxdir <- system.file("extdata/taxonomy",package="CHNOSZ")
id5 <- c(83333,4932,9606,186497,243232)
sciname(id5,taxdir)
# these are not all species, though
# (those with "no rank" are something like strains, 
# e.g. Escherichia coli K-12)
getrank(id5,taxdir)
# find the species for each of these
id6 <- sapply(id5,function(x) parent(x,taxdir=taxdir,rank="species"))
unique(getrank(id6, taxdir))  # "species"
# note that the K-12 is dropped
sciname(id6,taxdir)

## the complete nodes.dmp and names.dmp files are quite large,
## so it helps to store them in memory when performing multiple queries
## (this doesn't have a noticeable speed-up for the small files
## we use in this example)
taxdir <- system.file("extdata/taxonomy",package="CHNOSZ")
nodes <- getnodes(taxdir=taxdir)
# all of the node ids in this file
id7 <- nodes$id
# all of the non-leaf nodes
id8 <- unique(parent(id7,nodes=nodes))
names <- getnames(taxdir=taxdir)
sciname(id8,names=names)
}

\concept{Other tools}
