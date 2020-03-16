#! /usr/bin/env bash

#SBATCH -c 4                      # number of CPUs (here 1)
#SBATCH --mem 1000                # memory pool for all cores (here 1GB)
#SBATCH -t 0-24:00                 # walltime (D-HH:MM) (here 0 days, 24hours, 0 minutes)
#SBATCH -o slurm.%j.out        # STDOUT (the standard output stream)
#SBATCH -e slurm.%j.err        # STDERR (the output stream for errors)
#SBATCH --mail-type=begin         # notifications for job to begin 
#SBATCH --mail-type=end           # notifications for job to end
#SBATCH --mail-type=fail          # notifications for job to abort

FILE="$(readlink -f $1)"
OUTDIR="$(readlink -f $2)"
CORES="$3"

## TURNED OFF AS WGET FROM ENA WAS SLOW, FASTER VIA DOWNLOADING FROM SRA
# Give output directory and ERR code
download_err () {
	local DIR="$1"
	local LINE="$2"
	mkdir "$DIR"/$("echo $LINE")
	
	if [[ "$(echo "$LINE" | wc -m)" == 9 ]]; then
		wget -P "$DIR"/"$LINE"/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"$(echo $LINE | cut -c1-6)"/"$LINE"/"$LINE"*.fastq.gz

	elif [[ "$(echo "$LINE" | wc -m)" == 10 ]]; then
		wget -P "$DIR"/"$LINE"/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"$(echo $LINE | cut -c1-6)"/00"$(echo $LINE | cut -c1-6 | rev | cut -c1 |rev)"/"$LINE"/"$LINE"*.fastq.gz

	elif [[ "$(echo "$LINE" | wc -m)" == 11 ]]; then
		wget -P "$DIR"/"$LINE"/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"$(echo $LINE | cut -c1-6)"/0"$(echo $LINE | cut -c1-6 | rev | cut -c1-2 |rev)"/"$LINE"/"$LINE"*.fastq.gz

	elif [[ "$(echo "$LINE" | wc -m)" == 12 ]]; then
		wget -P "$DIR"/"$LINE"/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"$(echo $LINE | cut -c1-6)"/"$(echo $LINE | cut -c1-6 | rev | cut -c1-3 |rev)"/"$LINE"/"$LINE"*.fastq.gz

	else
		printf "\n $SAMPLE is not compatible with this script \n"
	fi
}

export -f download_err

if [[ "$#" -ne 3 ]]; then
	echo "Usage: SRR_ERR_download_script.sh <file_of_S/ERR codes>.txt /<out>/<dir> <no. avail cores>"
else 
	cat "$FILE" | parallel -j "$CORES" "download_err $OUTDIR $SAMPLE"
fi