#!/bin/bash

#SBATCH -n 1                      					# number of CPUs (here 4)
#SBATCH --mem 8000                					# memory pool for all cores (here 8GB)
#SBATCH -t 0-23:00            						# walltime (D-HH:MM) (here 0 days, 8hours, 0 minutes)
#SBATCH --partition=medium						# partition (queue) to submit to
#SBATCH --array=0-187%8
#SBATCH --mail-type=FAIL,ARRAY_TASKS          					# notifications for job to abort
#SBATCH --mail-type=time_limit          				# notifications for job to abort
#SBATCH -J "AdapterRemovalTrimOnly_SlurmArray"						# name of job

INDIR=04-analysis/screening/eager/polyGremoval_input/input
OUTDIR=04-analysis/screening/eager/polyGremoval_input/output

DIRS=($(find -L "$INDIR" -maxdepth 2 -mindepth 2 -name '*' -type d))
DIRNAME=${DIRS[$SLURM_ARRAY_TASK_ID]}

LIBNAME=$(echo "$DIRNAME" | rev | cut -d/ -f1 | rev)

CPU=1 ## number of cores
MEM=8 ## amount of memory (in GB)

  ## To assign directories and find input Forward and Reverse input files
LIBOUTDIR="${DIRNAME/\/input/\/output}"
FORWARD=($(find -L "$DIRNAME" -name '*_R1_*.gz' -type f | sort))
REVERSE=($(find -L "$DIRNAME" -name '*_R2_*.gz' -type f | sort))

if [ -z "$REVERSE" ]; then
  echo "SINGLE END"
  AdapterRemoval \
  --file1 "${FORWARD[@]}" \
  --basename "$LIBOUTDIR"/"$LIBNAME"_S0_L000_ \
  --gzip \
  --threads "$CPU" \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
  --minlength 0 \
  --minquality 20 \
  --minadapteroverlap 1    
  rename 's/_.truncated.gz/_R1_000.trimmed.fastq.gz/' "$LIBOUTDIR"/"$LIBNAME"_S0_L000_.truncated.gz
else
  echo "$PAIRED END"
  AdapterRemoval \
  --file1 "${FORWARD[@]}" \
  --file2 "${REVERSE[@]}" \
  --basename "$LIBOUTDIR"/"$LIBNAME"_S0_L000_ \
  --gzip \
  --threads "$CPU" \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
  --minlength 0 \
  --minquality 20 \
  --minadapteroverlap 1
  rename s/_.pair1.truncated.gz/_R1_000.trimmed.fastq.gz/ "$LIBOUTDIR"/"$LIBNAME"_S0_L000_.pair1.truncated.gz
  rename s/_.pair2.truncated.gz/_R2_000.trimmed.fastq.gz/ "$LIBOUTDIR"/"$LIBNAME"_S0_L000_.pair2.truncated.gz
fi


