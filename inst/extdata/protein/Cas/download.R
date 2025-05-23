# Cas/download.R
# Download Cas sequences from UniProt
# 20250522

dat <- read.csv("Cas_uniprot.csv")

# Loop over UniProt entries
for(UniProt in dat$UniProt) {
  if(UniProt == "") next
  file <- paste0(UniProt, ".fasta")
  # Skip if already downloaded
  if(file.exists(file.path("fasta/", file))) next
  URL <- paste0("https://rest.uniprot.org/uniprotkb/", file)
  cmd <- paste("wget", URL)
  print(cmd)
  system(cmd)
  # Move downloaded file to fasta directory
  file.rename(file, file.path("fasta/", file))
}

# Loop over Protein to get UniParc sequences
for(Protein in dat$Protein) {
  if(Protein == "") next
  if(!grepl("^UPI", Protein)) next
  file <- paste0(Protein, ".fasta")
  # Skip if already downloaded
  if(file.exists(file.path("fasta/", file))) next
  URL <- paste0("https://rest.uniprot.org/uniparc/", file)
  cmd <- paste("wget", URL)
  print(cmd)
  system(cmd)
  # Move downloaded file to fasta directory
  file.rename(file, file.path("fasta/", file))
}
