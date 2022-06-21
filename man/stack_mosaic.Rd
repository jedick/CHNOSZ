\encoding{UTF-8}
\name{stack_mosaic}
\alias{stack_mosaic}
\alias{add_alpha}
\title{Stacked mosaic diagram}
\description{
Create a stacked mosaic diagram, where the species formed in the first layer become the basis species for the species formed in the second layer.
The species in each layer are usually minerals with different metals; any bimetallic species are added to the second layer.
}

\usage{
  stack_mosaic(bases, species1, species2, species12, names = NULL,
    col = list(4, 3, 6), col.names = list(4, 3, 6), fill = NULL,
    dx = list(0, 0, 0), dy = list(0, 0, 0), srt = list(0, 0, 0),
    lwd = list(1, 1, 1), lty = list(1, 1, 1), ...)
  add_alpha(col, alpha)
}

\arguments{
  \item{bases}{basis species to be changed for each layer (commonly S-bearing aqueous species)}
  \item{species1}{species (minerals and/or aqueous species) with metal 1}
  \item{species2}{species with metal 2}
  \item{species12}{bimetallic species}
  \item{names}{character, species names (or chemical formulas) for labeling fields}
  \item{col}{line color}
  \item{col.names}{text color}
  \item{fill}{field color}
  \item{dx}{label x-offset}
  \item{dy}{label y-offset}
  \item{srt}{label rotation}
  \item{lwd}{line width}
  \item{lty}{line type}
  \item{...}{arguments for \code{\link{mosaic}} and \code{\link{affinity}}}
  \item{alpha}{character, hexadecimal value of color transparency (alpha)}
}

\details{
\code{stack_mosaic} creates a stacked mosaic diagram following steps that are described in detail in the vignette \viglink{multi-metal}.
Briefly, the first layer of the diagram is made by speciating the species in \code{bases} across the diagram to form the first set of species in \code{species1}.
Then, both \code{bases} \emph{and} \code{species1} (the stable species at each point on the diagram) are used to form the second set of species, consisting of both \code{species2} \emph{and} \code{species12}.

The plot parameters \code{col}, \code{col.names}, \code{fill}, \code{dx}, \code{dy}, \code{srt}, \code{lwd}, and \code{lty} should be length-3 lists (not vectors).
The values of elements 1--3 of the list are recycled to the number of species in \code{species1}, \code{species2}, and \code{species12}, respectively.

For \code{fill}, the default is to use no fill except for \code{species12}, where the fill color is taken from \code{col.names} with added transparency.
The default definition of \code{fill} is \code{list(NA, NA, add_alpha(col.names[3], "50"))}.
\code{add_alpha} adds transparency to a color by appending the value of \code{alpha} to the hexadecimal representation of the color given in \code{col}.
}

\section{Warning}{
The bimetallic species in \code{species12} are shown as part of the second layer, although their formation is sensitive to the presence of stable species in the first layer.
It follows that changing the order of layers (i.e., swapping \code{species1} and \code{species2}) can affect the depiction of mineral assemblages that have \code{species12}.
Only one of the alternatives is thermodynamically correct, but currently there is no check to determine which one it is.
}

\value{
A list of length two containing the output of each of the \code{\link{diagram}} calls use to make the diagram.
}

\examples{
# Define temperature (degrees C), pressure (bar), pH and logfO2 ranges
T <- 200
P <- "Psat"
res <- 200
pH <- c(0, 14, res)
O2 <- c(-48, -33, res)

# Define system: Fe-Cu-O-S-Cl
# NOTE: the basis species must include the first species listed
# in each of bases, species1, and species2 below
basis(c("pyrite", "Cu", "Cl-", "H2S", "H2O", "oxygen", "H+"))
basis("H2S", -2)
# Calculate solution composition for 1 mol/kg NaCl
NaCl <- NaCl(T = T, P = P, m_tot = 1)
basis("Cl-", log10(NaCl$m_Cl))

# Define arguments for stack_mosaic: Speciate aqueous sulfur
bases <- c("H2S", "HS-", "HSO4-", "SO4-2")
# Calculate stabilities of Fe-bearing minerals first
species1 <- c("pyrite", "pyrrhotite", "magnetite", "hematite")
# Calculate stabilities of Cu-bearing and FeCu-bearing minerals second
species2 <- c("copper", "cuprite", "tenorite", "chalcocite", "covellite")
species12 <- c("chalcopyrite", "bornite")
# Use abbreviations for Fe-bearing minerals and formulas for Cu-bearing minerals
names1 <- c("Py", "Po", "Mag", "Hem")
names2 <- info(info(species2))$formula
names12 <- info(info(species12))$formula
names <- list(names1, names2, names12)
# Adjust x-position for one species (chalcocite, Cu2S)
dx <- list(c(0, 0, 0, 0), c(0, 0, 0, 1, 0), c(0, 0))
# Use thick dashed lines for the bimetallic species
lwd <- list(1, 1, 2)
lty <- list(1, 1, 2)

# Make the diagram
stack_mosaic(bases, species1, species2, species12, names = names,
  dx = dx, lwd = lwd, lty = lty,
  pH = pH, O2 = O2, T = T, P = P, IS = NaCl$IS)
# Add legend and title
lTP <- lex(lTP(T, P))
db <- describe.basis(ibasis = c(3:4))
legend("topright", c(lTP, db), bg = "white")
title("Fe-Cu-S-O-H-Cl", font.main = 1)
}

\concept{Extended workflow}
