#!/usr/bin/env bash

if [ "$#" -ne 3 ]; then
  echo "Usage: 009-preprocessing_16s_mapping.sh </path/to/sample_output_directory> </path/to/sample_input_directory> <name_of_library>"
  exit 1
else
  ## To load specific software versions defined for this analysis
  source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/000-analysis_profile
  CPU=4 ## number of cores
  MEM=32 ## amount of memory (in GB)

  LIBOUTDIR="$1"
  INDIR="$2"
  LIBNAME="$3"

  INFILE="$(find -L $INDIR -name '*.fq.gz' -type f)"

  printf "+++ Starting +++ \n"
  printf "Outdir: $LIBOUTDIR \n"
  printf "Indir: $INDIR \n"
  printf "Lib name: $LIBNAME \n"
  printf "In file: $INFILE \n"
  printf "BWA: $BWA \n"
  printf "CPU: $CPU \n"
  printf "SILVADB: $SILVADB \n" 

  ## To map reads in library that can match positions on human reference genome,
  ## convert to standard BAM file and get statistics
  printf "\n +++ BWA aln +++ \n"
  "$BWA" aln \
  -t "$CPU" \
  "$SILVADB" \
  "$INFILE" \
  -n 0.01 \
  -l 32 \
  -f "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.sai

  printf "\n +++ BWA samse +++ \n"
  "$BWA" samse \
  -r "@RG\tID:ILLUMINA-"$LIBNAME"\tSM:"$LIBNAME"\tPL:illumina" \
  "$SILVADB" \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.sai \
  "$INFILE" \
  -f "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.sam

  printf "\n +++ SAMTOOLS view (sam to bam) +++ \n"
  "$SAMTOOLS" view -@ "$CPU" -bS \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.sam \
  -o "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.bam

  printf "\n +++ SAMTOOLS flagstat +++ \n"
  "$SAMTOOLS" flagstat \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.bam > \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.stats

  ## To extract reads with 'mapped' flag respectively
  printf "\n +++ SAMTOOLS view (extract mapped) +++ \n"
  "$SAMTOOLS" view -@ "$CPU" -F4 -q 0 -b \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.bam \
  -o "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.bam

  ## To sort mapped reads prior conversion to FASTA
  printf "\n +++ SAMTOOLS sort (mapped) +++ \n"
  "$SAMTOOLS" sort -@ "$CPU" -m "$MEM"G \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.bam \
  -o "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.sorted.bam

  ## To convert unmapped reads back to FASTA for downstream
  printf "\n +++ SAMTOOLS fasta (mapped) +++ \n"
  "$SAMTOOLS" fasta \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.sorted.bam > \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.fa

  ## To convert FASTA headers to make compatible with QIIME
  	HEADER="$(echo "$LIBNAME" | sed 's/_/./g')"
	
	awk -v name="$HEADER" '/^>/{print ">"name "_" ++i; next}{print}' \
	< "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.fa > \
	"$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped_renamed.fa

	gzip "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped_renamed.fa

   ## To clean up unecessary files
  printf "\n +++ CLEANUP +++ \n"
  rm \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.fq.gz \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.sai \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.sam \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.bam \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.bam \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.fa \
  "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.fastq.merged.prefixed.16smapped.sorted.bam

  exit 0

fi
