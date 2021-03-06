#!/bin/bash

#SBATCH -n 1                      					# number of CPUs (here 4)
#SBATCH --mem 4000                					# memory pool for all cores (here 8GB)
#SBATCH -t 0-02:00            						# walltime (D-HH:MM) (here 0 days, 8hours, 0 minutes)
#SBATCH --partition=short						# partition (queue) to submit to 
#SBATCH --array=0-353%8
#SBATCH --mail-type=fail          					# notifications for job to abort
#SBATCH --mail-type=time_limit          				# notifications for job to abort
#SBATCH -J "fastp_polyGtrim_SlurmArray"					# name of job


FILES=($(find -L 04-analysis/screening/eager/polyGremoval_output/input -name  '*.gz.pG.fq.gz' -type f ))
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

MACHINE=$(zcat "$FILENAME" | head -n 1 | cut -f1 -d ':')

zcat "$FILENAME" | sed -e '/^+$/,+1s/^$/!/' -e "/^$MACHINE/,+1s/^$/N/" -e '/^$/d' | pigz -p 4 > "$FILENAME".fixed.fq.gz