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
  add.protein(aa, as.residue = FALSE)
  seq2aa(sequence, protein = NA)
  aasum(aa, abundance = 1, average = FALSE, protein = NULL, organism = NULL)
}

\arguments{
  \item{aa}{data frame, amino acid composition in the format of \code{thermo()$protein}}
  \item{as.residue}{logical, normalize by protein length?}
  \item{sequence}{character, protein sequence}
  \item{protein}{character, name of protein; numeric, indices of proteins (rownumbers of \code{\link{thermo}()$protein})}
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

\code{seq2aa} returns a data frame of amino acid composition for the provided \code{sequence}, in the format of \code{thermo()$protein}.
In this function, the value of the \code{protein} argument is put into the \code{protein} column of the result.
If there is an underscore (e.g. \samp{LYSC_CHICK}), it is used to split the text, and the two parts are put into the \code{protein} and \code{organism} columns.

Given amino acid compositions returned by \code{seq2aa}, \code{add.protein} adds them to \code{thermo()$protein} for use by other functions in CHNOSZ.
The amino acid compositions of proteins in \code{aa} with the same name as one in \code{thermo()$protein} are replaced.
Set \code{as.residue} to TRUE to normalize by protein length; each input amino acid composition is divided by the corresponding number of residues, with the result that the sum of amino acid frequencies for each protein is 1.

\code{aasum} returns a data frame representing the sum of amino acid compositions in the rows of the input \code{aa} data frame.
The amino acid compositions are multiplied by the indicated \code{abundance}; that argument is recycled to match the number of rows of \code{aa}.
If \code{average} is TRUE the final sum is divided by the number of input compositions.
The name used in the output is taken from the first row of \code{aa} or from \code{protein} and \code{organism} if they are specified.
}

\value{
For \code{seq2aa}, a data frame of amino acid composition and identifying information for proteins.
For \code{add.protein}, the rownumbers of \code{thermo()$protein} that are added and/or replaced.
For \code{aasum}, a one-row data frame of amino acid composition and identifying information.
}

\examples{
\dontshow{reset()}
# Get the amino acid composition of a protein sequence
# (Human Gastric juice peptide 1)
aa <- seq2aa("LAAGKVEDSD", "GAJU_HUMAN")
# Add the protein to CHNOSZ
ip <- add.protein(aa)
# Calculate the protein length and chemical formula
protein.length(ip) # 10
as.chemical.formula(protein.formula(ip)) # "C41H69N11O18"

# Calculate a formula without using add.protein
aa <- seq2aa("ANLSG", "pentapeptide_test")
as.chemical.formula(protein.formula(aa))

# Sum the amino acid compositions of several poliovirus protein subunits
file <- system.file("extdata/protein/POLG.csv", package = "CHNOSZ")
aa <- read.csv(file, as.is = TRUE)
aasum(aa, protein = "POLG_sum")
}

\seealso{
\code{\link{read.fasta}} and \code{\link{uniprot.aa}} for other ways of getting amino acid compositions that can be used with \code{add.protein}.

\code{\link{pinfo}} for protein-level functions (length, chemical formulas, reaction coefficients of basis species).
}

\concept{Protein properties}
