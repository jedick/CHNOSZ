# CHNOSZ/add.protein.R
# calculate properties of proteins 20061109 jmd
# reorganize protein functions 20120513

# add.protein - add amino acid counts to thermo()$protein (returns iprotein)
# seq2aa - calculate amino acid counts from a sequence
# aasum - combine amino acid counts (sum, average, or weighted sum by abundance)

seq2aa <- function(protein, sequence) {
  # remove newlines and whitespace
  sequence <- gsub("\\s", "", gsub("[\r\n]", "", sequence))
  # make a data frame from counting the amino acids in the sequence
  caa <- count.aa(sequence)
  colnames(caa) <- aminoacids(3)
  # a protein with no amino acids is sort of boring
  if(all(caa==0)) stop("no characters match an amino acid")
  ip <- pinfo(protein)
  # now make the data frame
  po <- strsplit(protein, "_")[[1]]
  aa <- data.frame(protein=po[1], organism=po[2], ref=NA, abbrv=NA, stringsAsFactors=FALSE)
  aa <- cbind(aa, chains=1, caa)
  return(aa)
}

aasum <- function(aa, abundance=1, average=FALSE, protein=NULL, organism=NULL) {
  # returns the sum of the amino acid counts in aa,
  # multiplied by the abundances of the proteins
  abundance <- rep(abundance, length.out=nrow(aa))
  # drop any NA rows or abundances
  ina.aa <- is.na(aa$chains)
  ina.ab <- is.na(abundance)
  ina <- ina.aa | ina.ab
  if(any(ina)) {
    aa <- aa[!ina, ]
    abundance <- abundance[!ina]
    message("aasum: dropped ", sum(ina), " proteins with NA composition and/or abundance")
  }
  # we don't know how to deal with different numbers of polypeptide chains
  if(!all(aa$chains==aa$chains[1])) stop("different numbers of polypeptide chains")
  # multiply
  aa[, 6:25] <- aa[, 6:25] * abundance
  # sum
  out <- aa[1, ]
  out[, 5:25] <- colSums(aa[, 5:25])
  # average if told to do so
  if(average) {
    # polypeptide chains by number of proteins, residues by frequence
    out[, 5] <- out[, 5]/nrow(aa)
    out[, 6:25] <- out[, 6:25]/sum(abundance)
  }
  # add protein and organism names if given
  if(!is.null(protein)) out$protein <- protein
  if(!is.null(organism)) out$organism <- organism
  return(out)
}

add.protein <- function(aa, as.residue = FALSE) {
  # Add a properly constructed data frame of 
  # amino acid counts to thermo()$protein
  thermo <- get("thermo", CHNOSZ)
  if(!identical(colnames(aa), colnames(thermo$protein)))
    stop("'aa' does not have the same columns as thermo()$protein")
  # Check that new protein IDs are unique 20220418
  po <- paste(aa$protein, aa$organism, sep = "_")
  idup <- duplicated(po)
  if(any(idup)) stop(paste("some protein IDs are duplicated:", paste(unique(po[idup]), collapse = " ")))
  # Normalize by protein length if as.residue = TRUE 20220416
  if(as.residue) {
    pl <- protein.length(aa)
    aa[, 5:25] <- aa[, 5:25] / pl
  }
  # Find any protein IDs that are already present
  ip <- pinfo(po)
  ip.present <- !is.na(ip)
  # Now we're ready to go
  tp.new <- thermo$protein
  if(!all(ip.present)) tp.new <- rbind(tp.new, aa[!ip.present, ])
  if(any(ip.present)) tp.new[ip[ip.present], ] <- aa[ip.present, ]
  rownames(tp.new) <- NULL
  thermo$protein <- tp.new
  assign("thermo", thermo, CHNOSZ)
  # Return the new rownumbers
  ip <- pinfo(po)
  # Make some noise
  if(!all(ip.present)) message("add.protein: added ", nrow(aa)-sum(ip.present), " new protein(s) to thermo()$protein")
  if(any(ip.present)) message("add.protein: replaced ", sum(ip.present), " existing protein(s) in thermo()$protein")
  return(ip)
}
