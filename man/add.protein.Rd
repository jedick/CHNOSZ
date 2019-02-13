\encoding{UTF-8}
\name{add.protein}
\alias{add.protein}
\alias{seq2aa}
\alias{aasum}
\title{Amino Acid Compositions of Proteins}
\description{
  Functions to get amino acid compositions and add them to protein list for use by other functions.
}

\usage{
  add.protein(aa)
  seq2aa(protein, sequence)
  aasum(aa, abundance = 1, average = FALSE, protein = NULL, organism = NULL)
}

\arguments{
  \item{aa}{data frame, amino acid composition in the format of \code{thermo$protein}}
  \item{protein}{character, name of protein; numeric, indices of proteins (rownumbers of \code{\link{thermo}$protein})}
  \item{sequence}{character, protein sequence}
  \item{...}{additional arguments passed to \code{\link{read.csv}}}
  \item{abundance}{numeric, abundances of proteins}
  \item{average}{logical, return the weighted average of amino acid counts?}
  \item{organism}{character, name of organism}
}

\details{
A \samp{protein} in CHNOSZ is defined by its identifying information and the amino acid composition, stored in \code{\link{thermo}$protein}.
The names of proteins in CHNOSZ are distinguished from those of other chemical species by having an underscore character ("_") that separates two identifiers, referred to as the \code{protein} and \code{organism}.
An example is \samp{LYSC_CHICK}. 
The purpose of the functions described here is to identify proteins and work with their amino acid compositions.
From the amino acid compositions, the thermodynamic properties of the proteins can be estimated by group additivity.

\code{seq2aa} returns a data frame of amino acid composition, in the format of \code{thermo$protein}, corresponding to the provided \code{sequence}.
Here, the \code{protein} argument indicates the name of the protein with an underscore (e.g. \samp{LYSC_CHICK}).

\code{aasum} returns a data frame representing the sum of amino acid compositions in the rows of the input \code{aa} data frame.
The amino acid compositions are multiplied by the indicated \code{abundance}; that argument is recycled to match the number of rows of \code{aa}.
If \code{average} is TRUE the final sum is divided by the number of input compositions.
The name used in the output is taken from the first row of \code{aa} or from \code{protein} and \code{organism} if they are specified.

Given amino acid compositions returned by the \code{*aa} functions described above, \code{add.protein} adds them to \code{thermo$protein} for use by other functions in CHNOSZ.
The amino acid compositions of proteins in \code{aa} with the same name as one in \code{thermo$protein} are replaced.
The value returned by this function is the rownumbers of \code{thermo$protein} that are added and/or replaced.
}

\examples{
\dontshow{reset()}
# manually adding a new protein
# Human Gastric juice peptide 1
aa <- seq2aa("GAJU_HUMAN", "LAAGKVEDSD")
ip <- add.protein(aa)
stopifnot(protein.length(ip)==10)
# the chemical formula of this peptide
as.chemical.formula(protein.formula(ip)) # "C41H69N11O18"
# we can also calculate a formula without using add.protein
aa <- seq2aa("pentapeptide_test", "ANLSG")
as.chemical.formula(protein.formula(aa))
}

\seealso{
\code{\link{read.fasta}}, \code{\link{uniprot.aa}}, \code{\link{yeast.aa}} for other ways of getting amino acid compositions.

\code{\link{pinfo}} for protein-level functions (length, chemical formulas, reaction coefficients of basis species).

\code{\link{protein}} for examples of affinity calculations and diagrams.
}

\concept{Protein properties}
