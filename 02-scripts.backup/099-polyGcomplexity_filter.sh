#!/bin/bash
 
#SBATCH -n 1                      					# number of CPUs (here 4)
#SBATCH --mem 4000                					# memory pool for all cores (here 8GB)
#SBATCH -t 0-02:00            						# walltime (D-HH:MM) (here 0 days, 8hours, 0 minutes)
#SBATCH --partition=short						# partition (queue) to submit to 
#SBATCH --array=0-183%8
#SBATCH -o slurm.%j.out      # STDOUT (the standard output stream)
#SBATCH -e slurm.%j.err      # STDERR (the output stream for errors)
#SBATCH --mail-type=fail          					# notifications for job to abort
#SBATCH --mail-type=time_limit          				# notifications for job to abort
#SBATCH --mail-use=fellows@shh.mpg.de 					# these notifications will be sent to your shh.mpg.de email-address.
#SBATCH -J "fasp_polyGtrim_SlurmArray"					# name of job


FILES=($(find -L 04-analysis/screening/eager/polyGremoval_input/ -name '*mappedonly.sorted.cleaned_rmdup.sorted.bam.fq.gz' -type f))
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

fastp -i "$FILENAME" -o ${FILENAME/.fq.gz/.polyGtrimmed.fq.gz} --trim_poly_g --disable_quality_filtering --disable_length_filtering --disable_adapter_trimming
