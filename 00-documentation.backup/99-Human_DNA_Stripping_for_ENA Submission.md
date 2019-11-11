# Human DNA Stripping for ENA Submission

Map to human reference genome

```bash
INDIR=/projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered/input
OUTDIR=/projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered/output
SCRIPTS=/projects1/microbiome_calculus/evolution/02-scripts.backup/

for LIBDIR in "$INDIR"/*/; do
  LIBNAME=$(echo "$LIBDIR" | rev | cut -d/ -f2 | rev)
  if [ -d "$OUTDIR"/"$LIBNAME" ]; then
    printf "\n $LIBNAME already Processed \n\n"
  else
    mkdir "$OUTDIR"/"$LIBNAME"
    sbatch \
    -c 4 \
    --mem 32000 \
    -o ~/slurm_logs/slurm.%j.out \
    -e ~/slurm_logs/slurm.%j.err \
    -t 48:00:00 \
    --mail-type=fail \
    --mail-type=time_limit \
    --mail-user=fellows@shh.mpg.de \
    --export=ALL \
    -J "PREPRO_$LIBNAME" \
    --wrap="$SCRIPTS/003-preprocessing_human_filtering_forHumanStripping.sh $OUTDIR/$LIBNAME $LIBDIR $LIBNAME"
    sleep 1
  fi
done
```

Extract read IDs from mapped BAM files, using GNU parallel [Tange 2011](http://dx.doi.org/10.5281/zenodo.16303)

```bash
find -name '*hg19mapped.bam' -type f | parallel -j 4 "samtools view {} | cut -d' ' -f 1 > {}.readids.txt"
find -name '*readids.txt' -type f | parallel -j 4 "sed -i 's/M_//g; s/F_//g; s/R_//g' {}"

```

Use `seqkit` v0.10.1 ([Shen et al. 2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163962)) 
to remove any read with read ID from list above from the original FASTQs

```bash
cd /projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered/input

for i in */*.gz; do
	echo "$i"
	zcat "$i" | seqkit grep -j 4 -f /projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered/output/"$(dirname $i)"/*readids.txt -v | gzip > /projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered/output/"$i".stripped.fastq.gz
done
```

To validate the stripping worked correctly, I then re-run the mapping script on 
the stripped files.

```bash
INDIR=/projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered_validation/input
OUTDIR=/projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered_validation/output
SCRIPTS=/projects1/microbiome_calculus/evolution/02-scripts.backup/

for LIBDIR in "$INDIR"/*/; do
  LIBNAME=$(echo "$LIBDIR" | rev | cut -d/ -f2 | rev)
  if [ -d "$OUTDIR"/"$LIBNAME" ]; then
    printf "\n $LIBNAME already Processed \n\n"
  else
    mkdir "$OUTDIR"/"$LIBNAME"
    sbatch \
    -c 4 \
    --mem 32000 \
    -o ~/slurm_logs/slurm.%j.out \
    -e ~/slurm_logs/slurm.%j.err \
    -t 48:00:00 \
    --mail-type=fail \
    --mail-type=time_limit \
    --mail-user=fellows@shh.mpg.de \
    --export=ALL \
    -J "PREPRO_$LIBNAME" \
    --wrap="$SCRIPTS/003-preprocessing_human_filtering.sh $OUTDIR/$LIBNAME $LIBDIR $LIBNAME"
    sleep 1
  fi
done
```