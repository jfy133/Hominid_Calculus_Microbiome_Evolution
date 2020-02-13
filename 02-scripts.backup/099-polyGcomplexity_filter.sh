#!/bin/bash
 
#SBATCH -n 1                      					# number of CPUs (here 4)
#SBATCH --mem 4000                					# memory pool for all cores (here 8GB)
#SBATCH -t 0-02:00            						# walltime (D-HH:MM) (here 0 days, 8hours, 0 minutes)
#SBATCH --partition=short						# partition (queue) to submit to 
#SBATCH --array=0-187%8
#SBATCH -o /projects1/clusterhomes/fellows/slurm_logs/slurm.%j.out      # STDOUT (the standard output stream)
#SBATCH -e /projects1/clusterhomes/fellows/slurm_logs/slurm.%j.err      # STDERR (the output stream for errors)
#SBATCH --mail-type=fail          					# notifications for job to abort
#SBATCH --mail-type=time_limit          				# notifications for job to abort
#SBATCH --mail-use=fellows@shh.mpg.de 					# these notifications will be sent to your shh.mpg.de email-address.
#SBATCH -J "fastp_polyGtrim_SlurmArray"					# name of job


DIRS=($(find -L /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/output -maxdepth 2 -mindepth 2 -type d ))
DIRNAME=${DIRS[$SLURM_ARRAY_TASK_ID]}
OUTDIR="${DIRNAME/_input\/output/_output/input}"

FORWARD=($(find -L $DIRNAME -name '*_R1_*.gz' -type f | sort))
REVERSE=($(find -L $DIRNAME -name '*_R2_*.gz' -type f | sort))

mkdir -p $OUTDIR

if [ -z "$REVERSE" ]; then
  #SINGLE
  B_NAME=$(basename "${FORWARD[0]}")
  FORFILE="${B_NAME/_L00[1-9]_/_L000_}"

  /projects1/users/fellows/bin.backup/fastp/fastp \
  -i "$FORWARD" \
  -o "$OUTDIR"/"${FORFILE}".pG.fq.gz \
  --trim_poly_g \
  --disable_quality_filtering \
  --disable_length_filtering \
  --disable_adapter_trimming \
  -j "$OUTDIR"/"${FORFILE}".json \
  -h "$OUTDIR"/"${FORFILE}".html
 
else
   #PAIRED
  FOR_B_NAME=$(basename "${FORWARD[0]}")
  REV_B_NAME=$(basename "${REVERSE[0]}")
  FORFILE="${FOR_B_NAME/_L00[1-9]_/_L000_}"
  REVFILE="${REV_B_NAME/_L00[1-9]_/_L000_}"


  /projects1/users/fellows/bin.backup/fastp/fastp \
  -i "$FORWARD" \
  -I "$FORWARD" \
  -o "$OUTDIR"/"${FORFILE}".pG.fq.gz \
  -O "$OUTDIR"/"${REVFILE}".pG.fq.gz \
  --trim_poly_g \
  --disable_quality_filtering \
  --disable_length_filtering \
  --disable_adapter_trimming \
  -j "$OUTDIR"/fastp.json \
  -h "$OUTDIR"/fastp.html
 
fi



