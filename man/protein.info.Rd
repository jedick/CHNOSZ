\encoding{UTF-8}
\name{protein.info}
\alias{protein.info}
\alias{pinfo}
\alias{protein.length}
\alias{protein.formula}
\alias{protein.OBIGT}
\alias{protein.basis}
\title{Summaries of thermodynamic properties of proteins}

\description{
  Protein information, length, chemical formula, thermodynamic properties by group additivity, reaction coefficients of basis species, and metastable equilibrium example calculation.
}

\usage{
  pinfo(protein, organism=NULL, residue=FALSE, regexp=FALSE)
  protein.length(protein, organism = NULL)
  protein.formula(protein, organism = NULL, residue = FALSE)
  protein.OBIGT(protein, organism = NULL, state=thermo()$opt$state)
  protein.basis(protein, T = 25, normalize = FALSE)
}

\arguments{
  \item{protein}{character, names of proteins; numeric, species index of proteins; data frame; amino acid composition of proteins}
  \item{organism}{character, names of organisms}
  \item{residue}{logical, return per-residue values (those of the proteins divided by their lengths)?}
  \item{regexp}{logical, find matches using regular expressions?}
  \item{normalize}{logical, return per-residue values (those of the proteins divided by their lengths)?}
  \item{state}{character, physical state}
  \item{T}{numeric, temperature in \degC}
}

\details{
For character \code{protein}, \code{pinfo} returns the rownumber(s) of \code{thermo()$protein} that match the protein names.
The names can be supplied in the single \code{protein} argument (with an underscore, denoting protein_organism) or as pairs of \code{protein}s and \code{organism}s.
NA is returned for any unmatched proteins, including those for which no \code{organism} is given or that do not have an underscore in \code{protein}.

Alternatively, if \code{regexp} is TRUE, the \code{protein} argument is used as a pattern (regular expression); rownumbers of all matches of \code{thermo()$protein$protein} to this pattern are returned.
When using \code{regexp}, the \code{organism} can optionally be provided to return only those entries that also match \code{thermo()$protein$organism}.

For numeric \code{protein}, \code{pinfo} returns the corresponding row(s) of \code{thermo()$protein}.
Set \code{residue} to TRUE to return the per-residue composition (i.e. amino acid composition of the protein divided by total number of residues).

For dataframe \code{protein}, \code{pinfo} returns it unchanged, except for possibly the per-residue calculation.

The following functions accept any specification of protein(s) described above for \code{pinfo}:

\code{protein.length} returns the lengths (number of amino acids) of the proteins.

\code{protein.formula} returns a stoichiometrix matrix representing the chemical formulas of the proteins that can be pased to e.g. \code{\link{mass}} or \code{\link{ZC}}.
The amino acid compositions are multiplied by the output of \code{\link{group.formulas}} to generate the result. 

\code{protein.OBIGT} calculates the thermodynamic properties and equations-of-state parameters for the completely nonionized proteins using group additivity with parameters taken from Dick et al., 2006 (aqueous proteins) and LaRowe and Dick, 2012 (crystalline proteins and revised aqueous methionine sidechain group).
The return value is a data frame in the same format as \code{thermo()$OBIGT}.
\code{state} indicates the physical state for the parameters used in the calculation (\samp{aq} or \samp{cr}).

The following functions also depend on an existing definition of the basis species:

\code{protein.basis} calculates the numbers of the basis species (i.e. opposite of the coefficients in the formation reactions) that can be combined to form the composition of each of the proteins.
The basis species must be present in \code{thermo()$basis}, and if \samp{H+} is among the basis species, the ionization states of the proteins are included.
The ionization state of the protein is calculated using \code{\link{ionize.aa}} at the pH defined in \code{thermo()$basis} and at the temperature specified by the \code{T} argument.
If \code{normalize} is TRUE, the coefficients on the basis species are divided by the lengths of the proteins. 

}

\seealso{
  \code{\link{ionize.aa}}
}

\examples{\dontshow{reset()}
# Search by name in thermo()$protein
# These are the same: ip1 == ip2
ip1 <- pinfo("LYSC_CHICK")
ip2 <- pinfo("LYSC", "CHICK")
# Two organisms with the same protein name
ip3 <- pinfo("MYG", c("HORSE", "PHYCA"))
# Their amino acid compositions
pinfo(ip3)
# Their thermodynamic properties by group additivity
protein.OBIGT(ip3)

# An unknown protein name gives NA
ip4 <- pinfo("MYGPHYCA")

## Example for chicken lysozyme C
# Index in thermo()$protein
ip <- pinfo("LYSC_CHICK")
# Amino acid composition
pinfo(ip)
# Protein length and chemical formula
protein.length(ip)
protein.formula(ip)
# Group additivity for thermodynamic properties and HKF equation-of-state
# parameters of non-ionized protein
protein.OBIGT(ip)
# Calculation of standard thermodynamic properties
# (subcrt uses the species name, not ip)
subcrt("LYSC_CHICK")
# NOTE: subcrt() only shows the properties of the non-ionized
# protein, but affinity() uses the properties of the ionized
# protein if the basis species have H+

## These are all the same
protein.formula("P53_PIG")
protein.formula(pinfo("P53_PIG"))
protein.formula(pinfo(pinfo("P53_PIG")))
}

\references{
  Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \doi{10.5194/bg-3-311-2006}

  LaRowe, D. E. and Dick, J. M. (2012) Calculation of the standard molal thermodynamic properties of crystalline peptides. \emph{Geochim. Cosmochim. Acta} \bold{80}, 70--91. \doi{10.1016/j.gca.2011.11.041}
}

\concept{Protein properties}
