\encoding{UTF-8}
\name{util.seq}
\alias{util.seq}
\alias{aminoacids}
\alias{nucleic.formula}
\alias{nucleic.complement}

\title{Functions to Work with Sequence Data}

\description{
  Return one- or three-letter abbreviations of amino acids; count nucleotides in nucleic acid sequences, calculate DNA and RNA complements of nucleic acid sequences.
}

\usage{
  aminoacids(nchar=1, which=NULL)
  nucleic.formula(nucleic = NULL)
  nucleic.complement(nucleic = NULL, type="DNA")
}

\arguments{
  \item{nchar}{numeric, \eqn{1} to return one-letter, \eqn{3} to return three-letter abbreviations for amino acids}
  \item{which}{character, which amino acids to name}
  \item{nucleic}{data frame, counts of nucleic-acid bases}
  \item{type}{character, target type of nucleic acid (DNA or RNA)}
}

\details{

\code{aminoacids} returns the one-letter abbreviations (\code{nchar}=\samp{1}) or the three-letter abbreviations (\code{nchar}=\samp{3}) or the names of the neutral amino acids (\code{nchar}=\samp{""}) or the names of the amino acids with ionized side chains (\code{nchar}=\samp{"Z"}).
The output includes 20 amino acids in alphabetic order by 1-letter abbreviation (the order used in \code{thermo()$protein}), unless \code{which} is provided, indicating the desired amino acids (either as 1- or 3-letter abbreviations or names of the neutral amino acids).


\code{nucleic.formula} returns a string representation of the chemical formula for each nucleic-acid composition contained in \code{nucleic}.
The names of the bases are indicated by the column names of \code{nucleic}.
At present, the formula is computed as the sum of the chemical formulas of the bases themselves, with no contribution from polymerization (dehydration) or phosphorylation.

\code{nucleic.complement} calculates the complement of the base composition given in \code{nucleic}. 
\code{type} specifies the type of nucleic acid of the complement - \samp{DNA} (A, G, C, T) or \samp{RNA} (A, G, C, U).
}

\seealso{\code{\link{count.aa}} for counting amino acids or nucleic-acid bases in a sequence; \code{\link{protein.formula}} for calculating the chemical formulas of proteins.}

\examples{\dontshow{reset()}
## count nucleobases in a sequence
bases <- count.aa("ACCGGGTTT", type="DNA")
# the DNA complement of that sequence
DNA.comp <- nucleic.complement(bases)
# the RNA complement of the DNA complement
RNA.comp <- nucleic.complement(DNA.comp, type="RNA")
# the formula of the RNA complement (bases only)
nucleic.formula(RNA.comp)  # C40H42N32O11
}

\concept{Protein properties}
