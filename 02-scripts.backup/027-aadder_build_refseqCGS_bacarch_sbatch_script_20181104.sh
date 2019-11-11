#!/usr/bin/env bash
 
#SBATCH -n 112                      # number of CPUs (here 4)
#SBATCH --mem 1950G                # memory pool for all cores (here 8GB)
#SBATCH --partition=supercruncher
#SBATCH -o slurm.%j.out        # STDOUT (the standard output stream)
#SBATCH -e slurm.%j.err        # STDERR (the output stream for errors)
#SBATCH -J "AADDER-build-RefSeq_CGS_bacarchomo_2018-11"


aadder-build \
-igff 01-data/reference_databases/refseq/genomes/bacteria_archea_homo_20181122/raw \
-d  01-data/databases/aadder \
-a2t 01-data/malt/databases/acc2tax/nucl_acc2tax-Nov2018.abin \
-ex \
-v
