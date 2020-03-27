#!/usr/bin/env bash

## To do: finish collating all the interesting information
##  - fragment lenths
## Convert to loop, with header and remove header (leaving only the first)

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INDIR=$(readlink -f "$1")

echo "Sample_Name,Raw_Reads,Post_Clip_and_Merge,Percent_Merged_Content_Post_Clip_and_Merge,Mapped_Reads,Percent_Endogenous_DNA,Post_Duplicate_Removal_Reads,Cluster_Factor,X_Mean_Coverage,mtDNA_Reads,Median_Fragment_Length,GC_Content_Percent,Non_Human_Reads" >> human_filtering_statistics_"$(date +"%Y%m%d")".csv

for TMPDIR in "$INDIR"/*/; do
	
	LIBDIR="$(echo $TMPDIR | rev | cut -d/ -f2-999 | rev)"
	SAMPLE="$(echo "$LIBDIR" | rev | cut -d/ -f1 | rev)"
	
	echo "Reading $SAMPLE"
	
	LIBREADS=(0)

	for ZIP in "$LIBDIR"/fastqc/*.zip; do
	  RAWREADS=$(unzip -p "$ZIP" */fastqc_data.txt | grep 'Total Sequences' | cut -d$'\t' -f2)
	  LIBREADS+="$(echo ' '+ $RAWREADS)"
	done

	TOTALREADS="$(echo "$LIBREADS" | bc)"
	POSTADAPTERREMOVAL="$(cat "$LIBDIR"/*settings | grep "Number of retained reads:" | rev | cut -d" " -f1 | rev)"
	
	
	if grep -q 'Trimming of single-end reads' "$LIBDIR"/*settings; then
		MERGEDPC="NA"
	elif grep -q 'Trimming of paired-end reads' "$LIBDIR"/*settings; then
		FULLMERGED="$(cat "$LIBDIR"/*settings | grep "Number of full-length collapsed pairs:" | rev | cut -d" " -f1 | rev)"
		TRUNMERGED="$(cat "$LIBDIR"/*settings | grep "Number of truncated collapsed pairs:" | rev | cut -d" " -f1 | rev)"
		MERGEDPC=$(echo "(( $FULLMERGED + $TRUNMERGED ) / $POSTADAPTERREMOVAL) * 100" | bc -l)
	else
		"FAIL"
	fi
	
	MAPPEDREADS=$(cat "$LIBDIR"/*stats | grep "mapped (" | cut -d" " -f1)
	ENDODNA=$(echo "( $MAPPEDREADS / $POSTADAPTERREMOVAL ) * 100" | bc -l)
	POSTDEDUP=$(($(echo MAPPEDREADS) - $(cat "$LIBDIR"/*.log | grep "Total removed:" | cut -d" " -f3)))
	CLUSTERFACTOR=$(echo "$MAPPEDREADS / $POSTDEDUP" | bc -l)
	MEANCOV="$(grep 'mean coverageData' "$LIBDIR"/qualimap/genome_results.txt | cut -d= -f 2 | cut -d" " -f 2 | cut -dX -f 1)"
	MTREADS="$(grep 'chrMT' "$LIBDIR"/*idxstats | cut -d$'\t' -f 3)"
	FRAGLEN="$($SCRIPTDIR/006_2-calculate_frag_lengths.R $LIBDIR | cut -d" " -f 2)"
	GCPC="$(cat "$LIBDIR"/qualimap/genome_results.txt | grep GC | cut -d= -f 2 | cut -d" " -f 2 | cut -d% -f1)"
	UNMAPREADS=$(echo "$(cat "$LIBDIR"/*stats | grep "in total" | cut -d" " -f1) - $(cat "$LIBDIR"/*stats | grep "mapped (" | cut -d" " -f1)" | bc -l )
	
	echo "$SAMPLE,$TOTALREADS,$POSTADAPTERREMOVAL,$MERGEDPC,$MAPPEDREADS,$ENDODNA,$POSTDEDUP,$CLUSTERFACTOR,$MEANCOV,$MTREADS,$FRAGLEN,$GCPC,$UNMAPREADS" >> human_filtering_statistics_"$(date +"%Y%m%d")".csv

done