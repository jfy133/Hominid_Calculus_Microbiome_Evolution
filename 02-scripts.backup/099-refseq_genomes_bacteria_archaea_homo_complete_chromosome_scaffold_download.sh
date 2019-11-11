#!/bin/bash
 
#SBATCH -n 16                                                           # number of CPUs (here 4)
#SBATCH --mem 16000                                                     # memory pool for all cores (here 8GB)
#SBATCH --partition=long
#SBATCH -o slurm.%j.out      # STDOUT (the standard output stream)
#SBATCH -e slurm.%j.err      # STDERR (the output stream for errors)
#SBATCH --mail-type=fail                                                # notifications for job to abort
#SBATCH --mail-type=time_limit                                          # notifications for job to abort
#SBATCH -J "genomes_RefSeq_Download"                                           # name of job

tail -n+2 01-data/reference_databases/refseq/genomes/bacteria_archea_homo_20181122/docs/refseq_genomes_bacteria_archaea_homo_complete_chromosome_scaffold_downloadtable_20181122.txt | awk -F '\t' '{ print $20 }' | parallel \
--jobs 16 \
'mkdir 01-data/reference_databases/refseq/genomes/bacteria_archea_homo_20181122/raw/$(basename {}) && wget {}/*fna.gz \
-P 01-data/reference_databases/refseq/genomes/bacteria_archea_homo_20181122/raw/$(basename {}) && \
wget {}/*gff.gz \
-P 01-data/reference_databases/refseq/genomes/bacteria_archea_homo_20181122/raw/$(basename {}) && \
rm 01-data/reference_databases/refseq/genomes/bacteria_archea_homo_20181122/raw/$(basename {})/*_rna* \
01-data/reference_databases/refseq/genomes/bacteria_archea_homo_20181122/raw/$(basename {})/*_cds_*'

