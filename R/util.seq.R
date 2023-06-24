# CHNOSZ/util.seq.R
# Functions to work with sequence data

aminoacids <- function(nchar = 1, which = NULL) {
  # Return the abbreviations or names of the amino acids
  # The following are all in the same order as thermo()$protein
  # The single-letter codes
  aa1 <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
           "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  # The 3-letter codes
  aa3 <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
            "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  # Names of the neutral amino acids
  aaneutral <- c("alanine", "cysteine", "aspartic acid", "glutamic acid", "phenylalanine", 
    "glycine", "histidine", "isoleucine", "lysine", "leucine", "methionine", "asparagine", 
    "proline", "glutamine", "arginine", "serine", "threonine", "valine", "tryptophan", "tyrosine")
  # Names of the amino acids with ionized side chains
  aacharged <- c("alanine", "cysteinate", "aspartate", "glutamate", "phenylalanine", 
    "glycine", "histidinium", "isoleucine", "lysinium", "leucine", "methionine", "asparagine", 
    "proline", "glutamine", "argininium", "serine", "threonine", "valine", "tryptophan", "tyrosinate")
  # Defaults are in the same order as in thermo()$protein
  if(is.null(which)) which <- aa1
  # Figure out which amino acids are wanted (can use 1- or 3-letter codes, or neutral names)
  if(all(nchar(which) == 1)) iaa <- match(which, aa1)
  else if(all(nchar(which) == 3)) iaa <- match(which, aa3)
  else iaa <- match(which, aaneutral)
  # Return the desired abbreviations or names
  if(nchar == 1) return(aa1[iaa])
  else if(nchar == 3) return(aa3[iaa])
  else if(nchar == "") return(aaneutral[iaa])
  else if(nchar == "Z") return(aacharged[iaa])
}

