# CHNOSZ/add.protein.R
# Calculate properties of proteins 20061109 jmd
# Reorganize protein functions 20120513

# Add amino acid counts to thermo()$protein (returns iprotein)
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

