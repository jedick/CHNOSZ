# send the contents of all the .faa.gz files to a single file ("OUTFILE")

OUTFILE="refseq61.faa"
FILELIST="filelist"

# start with an empty file
rm $OUTFILE
touch $OUTFILE

while read -r LINE
do
  echo $LINE
  gzip -dc $LINE >> $OUTFILE
done < <( cat ${FILELIST} )
