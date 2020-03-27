#!/bin/bash
 
#SBATCH -n 4                      					# number of CPUs (here 4)
#SBATCH --mem 32000                					# memory pool for all cores (here 8GB)
#SBATCH -t 0-47:00            						# walltime (D-HH:MM) (here 0 days, 8hours, 0 minutes)
#SBATCH --partition=medium						# partition (queue) to submit to 
#SBATCH --array=0-199%8
#SBATCH -o slurm.%j.out      # STDOUT (the standard output stream)
#SBATCH -e slurm.%j.err      # STDERR (the output stream for errors)
#SBATCH --mail-type=FAIL,ARRAY_TASKS          					# notifications for job to abort
#SBATCH --mail-type=time_limit          				# notifications for job to abort
#SBATCH -J "EAGER_SlurmArray"						# name of job

FILES=($(find -L 04-analysis/screening/EMN_Neanderthal_phylogeny_check/eager/output -name '*OFN*' -type d))
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

unset DISPLAY && eagercli ${FILENAME}
