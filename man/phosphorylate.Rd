\name{phosphorylate}
\alias{phosphorylate}
\title{Calculate affinity of phosphorylation reactions with pH-dependent speciation}
\description{
  Calculate the affinity and Gibbs energy of phosphorylation reactions accounting for pH-dependent speciation of reactants, products, and phosphate sources using the mosaic approach.
}

\usage{
phosphorylate(reactant, P_source, loga_reactant = 0, loga_product = 0, 
                loga_P_source = 0, loga_P_remainder = 0, const_pH = 7, ...)
}
\arguments{
  \item{reactant}{
    character, the compound to be phosphorylated.
    Options include \code{"acetic acid"}, \code{"glycerol"}, \code{"adenosine_to_AMP"}, \code{"adenosine_to_cAMP"},
    \code{"ribose"}, \code{"uridine"}, \code{"AMP"}, \code{"ADP"}, \code{"glucose"}, or \code{"pyruvic acid"}.}
  \item{P_source}{
    character, the source of phosphate.
    Use \code{"P"} for basic reactions with inorganic phosphate,
    or \code{"PP"}, \code{"acetylphosphate"}, or \code{"ATP"} for extended reactions with other phosphate donors.}
  \item{loga_reactant}{numeric, logarithm of activity of the reactant. Default: 0}
  \item{loga_product}{numeric, logarithm of activity of the product. Default: 0}
  \item{loga_P_source}{numeric, logarithm of activity of the phosphate source. Default: 0}
  \item{loga_P_remainder}{numeric, logarithm of activity of the phosphate remainder (e.g., Pi when PP is the source). Default: 0}
  \item{const_pH}{numeric, pH value for the calculation. Default: 7}
  \item{\dots}{
    additional arguments passed to \code{\link{mosaic}},
    such as \code{pH}, \code{T} (temperature), or \code{P} (pressure) ranges for multi-dimensional calculations.}
}

\details{
  This function calculates the affinity of phosphorylation reactions by summing the affinities of component reactions
  while accounting for pH-dependent speciation of all complex species (reactants, products, and phosphate forms).
  Complex species are those that exist in multiple ionization states depending on pH.
  The \code{\link{mosaic}} function is used to calculate affinities accounting for the predominant species at each pH condition.

  The calculation method depends on whether the reaction is \dQuote{basic} or \dQuote{extended}:
  
  \strong{Basic reactions} (\code{P_source = "P"}): A single reaction of the form:
  \preformatted{  reactant + Pi = product + H2O}
  
  \strong{Extended reactions} (other \code{P_source} values): The sum of two reactions:
  \preformatted{  (1) reactant + Pi = product + H2O
  (2) P_source + H2O = P_remainder + Pi}
  
  The overall affinity is calculated as A_total = A1 - A2, where the negative sign accounts for the reverse direction of reaction 2.

  Use \code{const_pH} to perform calculations at constant pH (when \code{pH} is not provided as a range in \code{...}).

  If pressure (\code{P}) is one of the arguments in \code{...},
  then call this function with \code{P_source} as a named argument to avoid conflicting assignments due to partial argument matching.
  
  The function returns affinity values that can be converted to standard transformed Gibbs energy using \code{\link{convert}}.
}
\section{Warning}{
  This function depends on thermodynamic data for phosphorylated species that have not yet been added to OBIGT (2025-12-07).
  The data for glucose-6-phosphate used in the examples are provisional and are known to produce inaccurate results for the Gibbs energy of phosphorylation
  when used with the data in OBIGT for glucose, ATP, and ADP.
}

\value{
  A list with the following components:
  \item{m1}{the \code{\link{mosaic}} output for reaction 1 (reactant + Pi = product + H2O)}
  \item{m2}{the \code{\link{mosaic}} output for reaction 2 (NULL for basic reactions with \code{P_source = "P"})}
  \item{a12}{
    numeric, the overall affinity (dimensionless, A/2.303RT).
    For basic reactions this equals the affinity from m1; for extended reactions this is the sum of affinities from m1 and m2.}
  \item{P_reaction}{the formation reaction definition for the P_source (NULL for basic reactions)}
}
\seealso{
  \code{\link{mosaic}}, \code{\link{affinity}}, \code{\link{convert}}
}

\examples{
\dontshow{data(thermo)}
# Calculate standard transformed Gibbs energy (DeltaG°') for 
# glucose + ATP = glucose-6-phosphate + ADP at pH 7, 25 °C

# First add thermodynamic data for sugar phosphates (from Table 3.2 of Alberty, 2003)
mod.OBIGT("glucose-6-phosphate-2", formula = "C6H11O9P-2", G = -1763940)
mod.OBIGT("glucose-6-phosphate-1", formula = "C6H12O9P-", G = -1800590)
# Alberty (2003) doesn't have DeltaG° for neutral glucose-6-phosphate,
# so we calculate it from pKa1 = 1.5 (Degani and Halmann, 1966)
DG0_G6P <- -1800590 + convert(1.5, "G")
mod.OBIGT("glucose-6-phosphate", formula = "C6H13O9P", G = DG0_G6P)

# Calculate affinity at pH 7 with unit activities (loga = 0)
result <- phosphorylate("glucose", "ATP", const_pH = 7)

# Extract the overall affinity
A <- result$a12

# Convert affinity to Gibbs energy in kJ/mol
TK <- convert(25, "K")
deltaG_J <- convert(A, "G", T = TK)
deltaG_kJ <- deltaG_J / 1000

print(paste("DeltaG°' at pH 7:", round(deltaG_kJ, 2), "kJ/mol"))

# Example with non-standard activities
# log(a) = -3 for glucose, ATP, and products
result2 <- phosphorylate("glucose", "ATP", 
                         loga_reactant = -5,
                         loga_product = -4, 
                         loga_P_source = -3,
                         loga_P_remainder = -2,
                         const_pH = 7)
A2 <- result2$a12
deltaG2_kJ <- convert(A2, "G", T = convert(25, "K")) / 1000
print(paste("DeltaG at pH 7 with non-unit activities:", 
            round(deltaG2_kJ, 2), "kJ/mol"))
}

\references{
  Alberty RA (2003) \emph{Thermodynamics of Biochemical Reactions}. John Wiley & Sons, Hoboken, New Jersey, 397 p. \url{https://www.worldcat.org/oclc/51242181}

  Degani C and Halmann M. (1966) Solvolysis of phosphoric acid esters. Hydrolysis of glucose 6-phosphate. Kinetic and tracer studies.
  \emph{J Am Chem Soc} \bold{88}:4075-4082. \doi{10.1021/ja00969a033}
}

\concept{Thermodynamic calculations}
