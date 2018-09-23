# the following data files support calculations using the 
# RefSeq database (release 61, 2013-09-09)
protein_refseq.csv: overall (average) amino acid composition of all proteins for each
  microbial genome in the RefSeq collection (n=6758)
taxid_names.csv: taxid, phylum name and species name for 788 microbial taxa

# these functions/scripts have the following purpose (output files listed in parentheses):
gencat.sh - extract gi number, taxid, sequence length from RefSeq release catalog (gi.taxid.txt)
protein.refseq.R - get average amino acid composition for each taxid from gzipped sequence files (protein_refseq.csv)
taxid.names.R - get taxonomic names for each taxid represented (taxid_names.csv)
mkfaa.sh - combine gzipped sequence files into one big FASTA file (refseq61.faa)

# bash scripts assume a GNU/Linux-like operating system
# timings were made for processing RefSeq 61 on a recent (2009) intel laptop

# download stuff
1. download 'release61.files.installed', 'RefSeq-release61.catalog.gz',
   and 'release61.multispecies_WP_accession_to_taxname.gz' from NCBI
   (ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog)
2. list URLS for the microbial protein sequence files:
     grep microbial.*.protein.faa* release61.files.installed | \
       sed -e "s/^/ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/microbial\//g" > urllist
3. download the files using 'wget -i urllist' [232 files, 5.5 GB]
4. move the .gz files to a directory named 'protein'

# protein stuff
5. gzip -d RefSeq-release61.catalog.gz [3.5 GB]
6. use 'gencat.sh' to generate gi.taxid.txt for microbial proteins in the catalog [3 minutes]
   for RefSeq61, 'cat gi.taxid.txt | wc -l' is 28548891
7. generate protein_refseq_complete.csv in R:  [~14 hours]
   > source("protein.refseq.R")
   > protein.refseq()
   note that this depends on gi.taxid.txt and the .faa.gz files in the 'protein' directory
8. trim entries to produce protein_refseq.csv (smaller size, better for package distribution)
   > source("trim_refseq.R")

# taxonomy stuff
9. edit 'taxid.names.R' so that 'taxdir' points to the directory where the files
    'names.dmp' and 'nodes.dmp' are present. these files can be downloaded from
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz (accessed on 2013-09-18)
10. source 'taxid.names.R' to generate the file 'taxid_names.csv' [~5.5 hours]
10a. 20170926: To save space for the package, the file has been trimmed to
     hold only those taxids listed in extdata/bison/gi.taxid.txt.

# BLAST stuff (optional)
11. run ls protein/*.gz > filelist
12. use 'mkfaa.sh' to combine the sequences into a single file 'refseq61.faa' [11 GB, 11 minutes]
    for RefSeq61, 'grep "^>" refseq61.faa | wc -l' is 28548891 (same as catalog)
13. make a BLAST database, e.g. formatdb -t refseq61 -i refseq61.faa -p T
