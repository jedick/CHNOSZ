#/bin/sh
# extract microbial, genomic records from the RefSeq catalog
RELEASE=61
INFILE=RefSeq-release$RELEASE.catalog 
OUTFILE=RefSeq-release$RELEASE.catalog.microbial.protein

# extract the microbial protein records
grep "[[:space:]]NP_.*microbial" $INFILE  > $OUTFILE  # 449849 records
grep "[[:space:]]YP_.*microbial" $INFILE >> $OUTFILE  # 7898382 records
grep "[[:space:]]WP_.*microbial" $INFILE >> $OUTFILE  # 20200660 records

# to save only the gi and taxid columns
# the field separator (tab) is defined in command line, not in awk program,
#   otherwise the first line gets processed incorrectly
cat $OUTFILE | awk -F\\t '{print $4,$1}' > gi.taxid.unsrt

# sort the file on gi so that it can be used with e.g. the unix 'join' command
# (note: for join to work, -n for numeric sort is *not* used here)
cat gi.taxid.unsrt | sort > gi.taxid.txt
