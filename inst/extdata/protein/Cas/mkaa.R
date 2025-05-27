# Cas/mkaa.R
# Generate amino acid compositions from protein sequences
# 20250522

dat <- read.csv("Cas_uniprot.csv")
# Use UniProt ID as the file name
ID <- dat$UniProt
# In case UniProt ID is missing, use alternate ID
ID[ID == ""] <- dat$Protein[ID == ""]
# Store ID in data frame
dat$ID <- ID
# Remove missing IDs
dat <- subset(dat, ID != "")

# Get amino acid composition for each protein
aalist <- lapply(1:nrow(dat), function(iID) {
  file <- file.path("fasta", paste0(dat$ID[iID], ".fasta"))
  aa <- canprot::read_fasta(file)
  # Store systematic name and ID
  aa$protein <- dat$Systematic[iID]
  aa$ref <- dat$ID[iID]
  aa
})

# Convert list to data frame
aa <- do.call(rbind, aalist)
# Capitalize protein names
aa$protein <- gsub("cas", "Cas", aa$protein)
aa$protein <- gsub("csx", "Csx", aa$protein)
aa$protein <- gsub("c2c", "C2c", aa$protein)
aa$protein <- gsub("din", "Din", aa$protein)
aa$protein <- gsub("tns", "Tns", aa$protein)
# Save results
write.csv(aa, "Cas_aa.csv", row.names = FALSE, quote = 1)
