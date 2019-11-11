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

> IMPORTATN MAKE SURE THE -d IS ACTUALLY A TAB AND NOT SPACE!

```bash
find -name '*hg19mapped.bam' -type f | parallel -j 4 "samtools view {} | cut -d$'\t' -f 1 > {}.readids.txt"
find -name '*readids.txt' -type f | parallel -j 4 "sed -i 's/M_//g; s/F_//g; s/R_//g' {}"

```

Use `seqkit` v0.10.1 ([Shen et al. 2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163962)) 
to remove any read with read ID from list above from the original FASTQs

```bash
cd /projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered/input

for i in */*.gz; do
	echo "$i"
	zcat "$i" | seqkit grep -v -j 4 -f /projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered/output/"$(dirname $i)"/*readids.txt -v | gzip > /projects1/microbiome_calculus/evolution/01-data/ENA_upload/human_filtered/output/"$i".stripped.fastq.gz
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

## Adapter Clipping of Whole dataset

Can run the following SLURM array, on the normal raw FASTQs for most libaries,
but the stripped ones for the VLC modern day humans.

`003-preprocessing_AdapterRemovalOnly_SLURM.sh` (switch single and paired_fastq accordingly)

> Note that AdapterRemoval will leave an 'empty' read (entry exists,
> but sequence length of 0) if the entire read is adapter.
> The script aboves fixes this issue by replacing any 'blank' read with a
> single nucleotide of N and phred-score of ! 

> Make sure to switch `single_fastq` to `paired_fastq` accordingly

To check there are no 'empty' reads, you can run the following.

To identify e.g.:

```bash
zgrep -B 1 '^$' *_R1_*trimmed.fastq.gz -n > blank_reads_R1.txt
```

To check the differences between the trimmed and 'fixed' files (`len.fq.gz`)
you can run 

```bash
zdiff *trimmed.fastq.gz *trimmed.len.fq.gz
```

Next check that there are no singletons or discarded reads in those files (file
size should be 20bytes)

e.g.

```bash
cd output
ls -lh */* | less
```

And if so, remove the singleton and discard files

```bash
rm */*discarded* */*singleton* */*settings
```

> Make sure to do this also on single end files

Check for a couple of examples that same amount of lines in input
and output FASTQ files using `zcat <FILE> | wc -l`.

You can then rename these with

```bash
rename 's/.len.fq.gz/.fastq.gz/' */*
```



## MD5 calculations

For ENA upload we also need to create md5 hashs for validating the file
was not corrupted during upload.

```bash
touch paired_fastqs.md5

for i in */*.fastq.gz; do
  echo "$i"
  md5sum "$i" >>  paired_fastqs.md5
done
```

```bash

touch single_fastqs.md5

for i in */*.fastq.gz; do
  echo "$i"
  md5sum "$i" >>  single_fastqs.md5
done
```

