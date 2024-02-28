\encoding{UTF-8}
\name{add.protein}
\alias{add.protein}
\alias{aasum}
\title{Amino Acid Compositions of Proteins}
\description{
  Functions to get amino acid compositions and add them to protein list for use by other functions.
}

\usage{
  add.protein(aa, as.residue = FALSE)
  aasum(aa, abundance = 1, average = FALSE, protein = NULL, organism = NULL)
}

\arguments{
  \item{aa}{data frame, amino acid composition in the format of \code{thermo()$protein}}
  \item{as.residue}{logical, normalize by protein length?}
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

Given a data frame of amino acid compositions in the format of \code{thermo()$protein}, \code{add.protein} adds them to \code{thermo()$protein} for use by other functions in CHNOSZ.
The amino acid compositions of proteins in \code{aa} with the same name as one in \code{thermo()$protein} are replaced.
Set \code{as.residue} to TRUE to normalize by protein length; each input amino acid composition is divided by the corresponding number of residues, with the result that the sum of amino acid frequencies for each protein is 1.

\code{aasum} returns a data frame representing the sum of amino acid compositions in the rows of the input \code{aa} data frame.
The amino acid compositions are multiplied by the indicated \code{abundance}; that argument is recycled to match the number of rows of \code{aa}.
If \code{average} is TRUE the final sum is divided by the number of input compositions.
The name used in the output is taken from the first row of \code{aa} or from \code{protein} and \code{organism} if they are specified.
}

\value{
For \code{add.protein}, the rownumbers of \code{thermo()$protein} that are added and/or replaced.
For \code{aasum}, a one-row data frame of amino acid composition and identifying information.
}

\examples{
\dontshow{reset()}
# Read a file with the amino acid compositions of several poliovirus protein subunits
file <- system.file("extdata/protein/POLG.csv", package = "CHNOSZ")
aa <- read.csv(file)

# Add the proteins to CHNOSZ
iprotein <- add.protein(aa)
# Calculate length and elemental formula
protein.length(iprotein)
protein.formula(iprotein)

# Sum the amino acid compositions and assign a new protein name
aasum(aa, protein = "POLG_sum")
}

\seealso{
\code{\link[canprot]{read.fasta}} for reading amino acid compositions from FASTA files.

\code{\link{pinfo}} for protein-level functions (length, chemical formulas, reaction coefficients of basis species).
}

\concept{Protein properties}
