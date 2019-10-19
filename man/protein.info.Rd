\encoding{UTF-8}
\name{protein.info}
\alias{protein.info}
\alias{pinfo}
\alias{protein.length}
\alias{protein.formula}
\alias{protein.obigt}
\alias{protein.basis}
\alias{protein.equil}
\title{Summaries of Thermodynamic Properties of Proteins}

\description{
  Protein information, length, chemical formula, thermodynamic properties by group additivity, reaction coefficients of basis species, and metastable equilibrium example calculation.
}

\usage{
  pinfo(protein, organism=NULL, residue=FALSE, regexp=FALSE)
  protein.length(protein, organism = NULL)
  protein.formula(protein, organism = NULL, residue = FALSE)
  protein.obigt(protein, organism = NULL, state=thermo()$opt$state)
  protein.basis(protein, T = 25, normalize = FALSE)
  protein.equil(protein, T=25, loga.protein = 0, digits = 4)
}

\arguments{
  \item{protein}{character, names of proteins; numeric, species index of proteins; data frame; amino acid composition of proteins}
  \item{organism}{character, names of organisms}
  \item{residue}{logical, return per-residue values (those of the proteins divided by their lengths)?}
  \item{regexp}{logical, find matches using regular expressions?}
  \item{normalize}{logical, return per-residue values (those of the proteins divided by their lengths)?}
  \item{state}{character, physical state}
  \item{T}{numeric, temperature in \degC}
  \item{loga.protein}{numeric, decimal logarithms of reference activities of proteins}
  \item{digits}{integer, number of significant digits (see \code{\link{signif}})}
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

\code{protein.obigt} calculates the thermodynamic properties and equations-of-state parameters for the completely nonionized proteins using group additivity with parameters taken from Dick et al., 2006 (aqueous proteins) and LaRowe and Dick, 2012 (crystalline proteins and revised aqueous methionine sidechain group).
The return value is a data frame in the same format as \code{thermo()$obigt}.
\code{state} indicates the physical state for the parameters used in the calculation (\samp{aq} or \samp{cr}).

The following functions also depend on an existing definition of the basis species:

\code{protein.basis} calculates the numbers of the basis species (i.e. opposite of the coefficients in the formation reactions) that can be combined to form the composition of each of the proteins.
The basis species must be present in \code{thermo()$basis}, and if \samp{H+} is among the basis species, the ionization states of the proteins are included.
The ionization state of the protein is calculated at the pH defined in \code{thermo()$basis} and at the temperature specified by the \code{T} argument.
If \code{normalize} is TRUE, the coefficients on the basis species are divided by the lengths of the proteins. 

  \code{protein.equil} produces a series of messages showing step-by-step a calculation of the chemical activities of proteins in metastable equilibrium. For the first protein, it shows the standard Gibbs energies of the reaction to form the nonionized protein from the basis species and of the ionization reaction of the protein (if \samp{H+} is in the basis), then the standard Gibbs energy/RT of the reaction to form the (possibly ionized) protein per residue. The per-residue values of \samp{logQstar} and \samp{Astar/RT} are also shown for the first protein. Equilibrium calculations are then performed, only if more than one protein is specified. This calculation applies the Boltzmann distribution to the calculation of the equilibrium degrees of formation of the residue equivalents of the proteins, then converts them to activities of proteins taking account of \code{loga.protein} and protein length. If the \code{protein} argument is numeric (indicating rownumbers in \code{thermo()$protein}), the values of \samp{Astar/RT} are compared with the output of \code{\link{affinity}}, and those of the equilibrium degrees of formation of the residues and the chemical activities of the proteins with the output of \code{\link{diagram}}. If the values in any of these tests are are not \code{\link{all.equal}} an error is produced indicating a bug. 
}

\examples{\dontshow{reset()}
# search by name in thermo()$protein
ip1 <- pinfo("LYSC_CHICK")
ip2 <- pinfo("LYSC", "CHICK")
# these are the same
stopifnot(all.equal(ip1, ip2))
# two organisms with the same protein name
ip3 <- pinfo("MYG", c("HORSE", "PHYCA"))
# their amino acid compositions
pinfo(ip3)
# their thermodynamic properties by group additivity
protein.obigt(ip3)

# an example of an unrecognized protein name
ip4 <- pinfo("MYGPHYCA")
stopifnot(is.na(ip4))

## example for chicken lysozyme C
# index in thermo()$protein
ip <- pinfo("LYSC_CHICK")
# amino acid composition
pinfo(ip)
# length and chemical formula
protein.length(ip)
protein.formula(ip)
# group additivity for thermodynamic properties and HKF equation-of-state
# parameters of non-ionized protein
protein.obigt(ip)
# calculation of standard thermodynamic properties
# (subcrt uses the species name, not ip)
subcrt("LYSC_CHICK")
# affinity calculation, protein identified by ip
basis("CHNOS+")
affinity(iprotein=ip)
# affinity calculation, protein loaded as a species
species("LYSC_CHICK")
affinity()
# NB: subcrt() only shows the properties of the non-ionized
# protein, but affinity() uses the properties of the ionized
# protein if the basis species have H+

## these are all the same
protein.formula("P53_PIG")
protein.formula(pinfo("P53_PIG"))
protein.formula(pinfo(pinfo("P53_PIG")))

\dontshow{opar <- par(no.readonly=TRUE)}
# using pinfo() with regexp=TRUE:
# plot ZC and nH2O/residue of HOX proteins
# basis species: glutamine-glutamic acid-cysteine-O2-H2O
basis("QEC")
# device setup
par(mfrow=c(2, 2))
# a red-blue scale from 1-13
col <- ZC.col(1:13)
# axis labels
ZClab <- axis.label("ZC")
nH2Olab <- expression(bar(italic(n))[H[2]*O])
# loop over HOX gene clusters
for(cluster in c("A", "B", "C", "D")) {
  # get protein indices
  pattern <- paste0("^HX", cluster)
  ip <- pinfo(pattern, "HUMAN", regexp=TRUE)
  # calculate ZC and nH2O/residue
  thisZC <- ZC(protein.formula(ip))
  thisH2O <- protein.basis(ip)[, "H2O"] / protein.length(ip)
  # plot lines
  plot(thisZC, thisH2O, type="l", xlab=ZClab, ylab=nH2Olab)
  # the number of the HOX gene
  pname <- pinfo(ip)$protein
  nHOX <- as.numeric(gsub("[A-Za-z]*", "", pname))
  # plot colored points
  points(thisZC, thisH2O, pch=19, col=col[nHOX], cex=3.5)
  points(thisZC, thisH2O, pch=19, col="white", cex=2.5)
  # plot the number of the HOX gene
  text(thisZC, thisH2O, nHOX)
  # add title
  title(main=paste0("HOX", cluster))
}
\dontshow{par(opar)}
}

\references{
  Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \url{https://doi.org/10.5194/bg-3-311-2006}

  LaRowe, D. E. and Dick, J. M. (2012) Calculation of the standard molal thermodynamic properties of crystalline peptides. \emph{Geochim. Cosmochim. Acta} \bold{80}, 70--91. \url{https://doi.org/10.1016/j.gca.2011.11.041}

  Dick, J. M. (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \url{https://doi.org/10.1098/rsif.2013.1095}
}

\concept{Protein properties}
