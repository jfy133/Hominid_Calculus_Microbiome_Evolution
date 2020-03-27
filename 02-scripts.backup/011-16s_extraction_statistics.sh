#!/usr/bin/env bash

## Give folder with BWA aln '.stats' files to make a list of number of reads 
## each 'stats' file should be in a unique directory with the directory name 
## as the name of the sample.

INDIR=$(readlink -f "$1")

echo "File,Count" >> "$INDIR"/temp.txt

find "$INDIR" -maxdepth 2 -name *.stats -type f -exec grep 'mapped (' {} /dev/null \; >> \
"$INDIR"/temp.txt

sed -i -e "s|$INDIR||g" -e 's/:/,/g' -e 's/ + 0 mapped\(.*\)$//g' \
"$INDIR"/temp.txt >> "$INDIR"/temp.txt

cut -d/ -f3-99 "$INDIR"/temp.txt >> "$INDIR"/16s_extraction_statistics_"$(date +"%Y%m%d")".csv

rm "$INDIR"/temp.txt

exit 0
