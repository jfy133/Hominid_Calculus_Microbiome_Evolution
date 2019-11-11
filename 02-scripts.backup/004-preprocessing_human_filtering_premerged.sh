#!/usr/bin/env bash

##
## THIS ONLY WORKS FOR SHORT READ DATA THAT CAN BE MERGED PROPERLY!
## DO NOT USE FOR LONG READ ILLUMINA DATA (WHERE MATE PAIRS DO NOT OVERLAP)
##

if [ "$#" -ne 3 ]; then
  echo "Usage: 01-preprocessing_human_filtering.sh </path/to/sample_output_directory> </path/to/sample_input_directory> <name_of_library>"
  exit 1
else
  ## To load specific software versions defined for this analysis
  source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/00-analysis_profile
  CPU=4 ## number of cores
  MEM=32 ## amount of memory (in GB)

  ## To assign directories
  printf "\n +++ PREPARATION +++ \n"
  LIBOUTDIR="$1"

  ## Find input file
  FASTQ=($(find -L $2 -name '*merged.fq.gz' -type f | sort))
  ln -s $(readlink -f $FASTQ) $(readlink -f "$LIBOUTDIR") 

  ## Prefix modification

  "$ADAPTERREMOVALFIXPREFIX" \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.fq.gz \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.fq.gz

  ## To map reads in library that can match positions on human reference genome,
  ## convert to standard BAM file and get statistics
  printf "\n +++ BWA aln +++ \n"
  "$BWA" aln \
  -t "$CPU" \
  "$HG19REF" \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.fq.gz \
  -n 0.01 \
  -l 32 \
  -f "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.sai

  printf "\n +++ BWA samse +++ \n"
  "$BWA" samse \
  -r "@RG\tID:ILLUMINA-"$3"\tSM:"$3"\tPL:illumina" \
  "$HG19REF" \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.sai \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.fq.gz \
  -f "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.sam

  printf "\n +++ SAMTOOLS view (sam to bam) +++ \n"
  "$SAMTOOLS" view -@ "$CPU" -bS \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.sam \
  -o "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.bam

  printf "\n +++ SAMTOOLS flagstat +++ \n"
  "$SAMTOOLS" flagstat \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.bam > \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.stats

  ## To extract reads with 'mapped' and 'unmapped' flags respectively
  printf "\n +++ SAMTOOLS view (extract mapped) +++ \n"
  "$SAMTOOLS" view -@ "$CPU" -F4 -q 0 -b \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.sam \
  -o "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.bam

  printf "\n +++ SAMTOOLS view (extract unmapped) +++ \n"
  "$SAMTOOLS" view -@ "$CPU" -f4 -q 0 -b \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.sam -o \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.bam

  ## To convert unmapped reads back to FASTQ for downstream
  printf "\n +++ SAMTOOLS fastq (unmapped) +++ \n"
  "$SAMTOOLS" fastq \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.bam | \
  gzip > "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz

  ## To get remaining HG19 statistics, clean up BAM file for PCR duplicate removal
  printf "\n +++ SAMTOOLS sort (mapped) +++ \n"
  "$SAMTOOLS" sort -@ "$CPU" -m "$MEM"G \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.bam \
  -o "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.bam

  printf "\n +++ SAMTOOLS index (mapped) +++ \n"
  "$SAMTOOLS" index \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.bam

  printf "\n +++SAMTOOLS idxstats (mapped) +++ \n"
  "$SAMTOOLS" idxstats \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.bam \
  >> "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.bam.idxstats

  printf "\n +++ CLEANUP (post-samtools) +++ \n"
  rm \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.fq.gz \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.sai \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.sam \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.bam \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.bam \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.bam


  printf "\n +++ PICARDTOOLS CleanSam +++ \n"
  "$PICARDTOOLS" CleanSam \
  INPUT="$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.bam \
  OUTPUT="$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned.bam \
  VALIDATION_STRINGENCY=SILENT

  rm \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.bam

  ## To remove PCR dupliate removal and clean up
  printf "\n +++ DEDUP +++ \n"
  "$DEDUP" \
  -i "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned.bam \
  -m -o "$LIBOUTDIR"

  printf "\n +++ SAMTOOLS sort (dedupped mapped) +++ \n"
  samtools sort -@ "$CPU" -m "$MEM"G \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned_rmdup.bam \
  -o "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned_rmdup.sorted.bam

  printf "\n +++ SAMTOOLS index (dedupped mapped) +++ \n"
  "$SAMTOOLS" index \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned_rmdup.sorted.bam

  ## To generate final mapping statistics and damage
  printf "\n +++ QUALIMAP (dedupped mapped) +++ \n"
  "$QUALIMAP" bamqc -bam \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned_rmdup.sorted.bam \
  -nt 4 -outdir "$LIBOUTDIR"/qualimap -outformat HTML --java-mem-size=32G

  printf "\n +++ DAMAGEPROFILER (dedupped mapped) +++ \n"
  mkdir "$LIBOUTDIR"/damageprofiler
  "$DAMAGEPROFILER" \
  -i "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned_rmdup.sorted.bam \
  -r "$HG19REF" \
  -l 100 \
  -o "$LIBOUTDIR"/damageprofiler \
  -t 25

  ## To clean up unecessary files
  printf "\n +++ CLEANUP (final) +++ \n"
  rm -r \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.bam.bai \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned.bam \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned_rmdup.bam \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned_rmdup.sorted.bam \
  "$LIBOUTDIR"/"$3"_S0_L000_R1_000.fastq.merged.prefixed.hg19mapped.sorted.cleaned_rmdup.sorted.bam.bai \
  "$LIBOUTDIR"/qualimap/*/ \
  "$LIBOUTDIR"/qualimap/*.html

  exit 0

fi
