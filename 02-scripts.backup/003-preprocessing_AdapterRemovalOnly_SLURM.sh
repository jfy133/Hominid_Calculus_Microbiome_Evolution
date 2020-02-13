#!/bin/bash

#SBATCH -n 4                      					# number of CPUs (here 4)
#SBATCH --mem 32000                					# memory pool for all cores (here 8GB)
#SBATCH -t 0-23:00            						# walltime (D-HH:MM) (here 0 days, 8hours, 0 minutes)
#SBATCH --partition=medium						# partition (queue) to submit to
#SBATCH -o /projects1/clusterhomes/fellows/slurm_logs/slurm.%j.out      # STDOUT (the standard output stream)
#SBATCH -e /projects1/clusterhomes/fellows/slurm_logs/slurm.%j.err      # STDERR (the output stream for errors)
#SBATCH --array=0-4%4
#SBATCH --mail-type=FAIL,ARRAY_TASKS          					# notifications for job to abort
#SBATCH --mail-type=time_limit          				# notifications for job to abort
#SBATCH --mail-use=fellows@shh.mpg.de 					# these notifications will be sent to your shh.mpg.de email-address.
#SBATCH -J "AdapterRemovalTrimOnly_SlurmArray"						# name of job

INDIR=/projects1/microbiome_calculus/evolution/10-publication/ENA_upload/clipped/input/paired_fastq3_missing
OUTDIR=/projects1/microbiome_calculus/evolution/10-publication/ENA_upload/clipped/output/paired_fastq3_missing

DIRS=($(find -L "$INDIR" -maxdepth 1 -mindepth 1 -name '*' -type d))
DIRNAME=${DIRS[$SLURM_ARRAY_TASK_ID]}

LIBNAME=$(echo "$DIRNAME" | rev | cut -d/ -f1 | rev)

CPU=4 ## number of cores
MEM=32 ## amount of memory (in GB)

  ## To assign directories and find input Forward and Reverse input files
LIBOUTDIR="$OUTDIR"/"$LIBNAME"
FORWARD=($(find -L "$INDIR"/"$LIBNAME" -name '*_R1_*.gz' -type f | sort))
REVERSE=($(find -L "$INDIR"/"$LIBNAME" -name '*_R2_*.gz' -type f | sort))

if [ -z "$REVERSE" ]; then
  echo "SINGLE END"
  /projects1/users/fellows/bin.backup/miniconda3/envs/adapterremoval/bin/AdapterRemoval \
  --file1 "${FORWARD[@]}" \
  --basename "$LIBOUTDIR"/"$LIBNAME"_S0_L000_ \
  --gzip \
  --threads "$CPU" \
  --trimns \
  --trimqualities \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
  --minlength 0 \
  --minquality 20 \
  --minadapteroverlap 1    
  rename 's/_.truncated.gz/_R1_000.trimmed.fastq.gz/' "$LIBOUTDIR"/"$LIBNAME"_S0_L000_.truncated.gz
  MACHINE=$(zcat "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.trimmed.fastq.gz | head -n 1 | cut -f1 -d ':')
  for i in "$LIBOUTDIR"/*trimmed.fastq.gz; do
    zcat "$i" | sed -e '/^+$/,+1s/^$/!/' -e "/^$MACHINE/,+1s/^$/N/" | pigz -p 4 > "${i%%.fastq.gz}".len.fq.gz
  done
else
  echo "$PAIRED END"
  AdapterRemoval \
  --file1 "${FORWARD[@]}" \
  --file2 "${REVERSE[@]}" \
  --basename "$LIBOUTDIR"/"$LIBNAME"_S0_L000_ \
  --gzip \
  --threads "$CPU" \
  --trimns \
  --trimqualities \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
  --minlength 0 \
  --minquality 20 \
  --minadapteroverlap 1
  AdapterRemoval \
  --file1 "${FORWARD[@]}" \
  --file2 "${REVERSE[@]}" \
  --basename "$LIBOUTDIR"/"$LIBNAME"_S0_L000_ \
  --gzip \
  --threads "$CPU" \
  --trimns \
  --trimqualities \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
  --minlength 0 \
  --minquality 20 \
  --minadapteroverlap 1
  rename s/_.pair1.truncated.gz/_R1_000.trimmed.fastq.gz/ "$LIBOUTDIR"/"$LIBNAME"_S0_L000_.pair1.truncated.gz
  rename s/_.pair2.truncated.gz/_R2_000.trimmed.fastq.gz/ "$LIBOUTDIR"/"$LIBNAME"_S0_L000_.pair2.truncated.gz
  MACHINE=$(zcat "$LIBOUTDIR"/"$LIBNAME"_S0_L000_R1_000.trimmed.fastq.gz | head -n 1 | cut -f1 -d ':')
  for i in "$LIBOUTDIR"/*trimmed.fastq.gz; do
    zcat "$i" | sed -e '/^+$/,+1s/^$/!/' -e "/^$MACHINE/,+1s/^$/N/" | pigz -p 4 > "${i%%.fastq.gz}".len.fq.gz
  done
fi


