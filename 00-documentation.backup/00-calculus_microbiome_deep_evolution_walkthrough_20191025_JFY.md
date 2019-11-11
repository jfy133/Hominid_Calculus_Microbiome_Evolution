# Evolution of the Plaque Oral Microbiome Walkthrough

By James Fellows Yates

## Introduction

This walkthrough goes through each step of the analysis carried out for the 
**ARTICLE TITLE** paper.

This analysis was performed on Ubuntu 14.04 LTS, and a SLURM submission
scheduler. Note that the SLURM submission commands are directly related to
our system, so you will need to modify these commands accordingly throughout 
so it works with your cluster.

### Resources

Here is a list of programs and databases that will be used in this analysis
and that you should have ready installed/downloaded prior carrying out the
analysis. Note that the download and set up of the databases are described 
below.

**Software**:

Name                    | Version                 | Citation
------------------------|-------------------------|------------------------------------------------------------
GNU bash                | 4.3.11                  | NA
sratoolkit              | 2.8.0                   | https://www.ncbi.nlm.nih.gov/books/NBK158900/
EAGER                   | 1.92.55                 | Peltzer et al. 2016
FastQC                  | 0.11.4                  | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
AdapterRemoval2         | 2.2.0                   | Schubert et al. 2016
AdapterRemovalFixPrefix | 0.0.3                   | https://github.com/apeltzer/AdapterRemovalFixPrefix/
bwa                     | 0.7.12                  | Li and Durbin 2009
samtools                | 1.3                     | Li et al. 2009
PicardTools             | 1.140                   | https://github.com/broadinstitute/picard
DeDup                   | 0.12.2                  | Peltzer et al. 2016
QualiMap                | 2.2.1                   | Okonechnikov et al. 2016
DamageProfiler          | 0.3.8                   | Judith Neukamm (Unpublished)
MALT                    | 0.4.0                   | Herbig et al. 2016 bioRxiv, Vagene et al. 2018
MEGANCE                 | 6.X.X                   | Huson et al. 2016
QIIME                   | 1.9.1                   | Caporaso et al. 2010
R                       | 3.4.4                   | https://www.R-project.org/
MetaPhlAn2              | 2.7.1                   | Truong et al. 2015
Inkscape                | 0.92                    | www.inkscape.org
MultiVCFAnalyzer        | v0.91                   | https://github.com/alexherbig/MultiVCFAnalyzer
blastn                  | 2.7.1+                  | Package: blast 2.7.1, build Oct 18 2017 19:57:24
megax                   | 10.0.4                  | Tamura et al. 2018 
Entrez Direct           | Jun 2016                | http://www.ncbi.nlm.nih.gov/books/NBK179288
pdftotext               | 0.41.0                  | http://poppler.freedesktop.org
seqtk                   | 1.2-r95-dirty           | https://github.com/lh3/seqtk
pigz                    | 2.3                     | https://zlib.net/pigz/
mapDamage               | 2.0.6                   | Jónsson et al. ‎2013 Bioinformatics
bedtools                | 2.25.0                  | Quinlan et al. 2010 Bioinformatics

**R Packages**:

Name             | Version |  Citation
----------------|---------|-------------------------------------
tidyverse       | 1.2.1   | Wickham 2017 https://CRAN.R-project.org/package=tidyverse
ggplot2         | 2.2.1   | Wickham 2009
tibble          | 1.3.3   | Müller and Wickham 2017
tidyr           | 0.6.3   | Wickham 2017
readr           | 1.1.1   | Wickham, Hester and Francois 2017
dplyr           | 0.7.1   | Wickham et al. 2017
stringr         | 1.3.1   | https://github.com/tidyverse/stringr
zCompositions   | 1.1.1   | Palarea-Albaladejo et al. ‎2015
decontam        | 1.1.0   | Davis et al. 2018 bioRxiv, https://github.com/benjjneb/decontam
plotly          | 4.7.1   | https://help.plot.ly/citations/
philr           | 1.2.0   | Silverman et al. 2017
SourceTracker   | 1.0.1   | Knights et al. 2011
patchwork       | 0.0.1   | https://github.com/thomasp85/patchwork
phyloseq        | 1.20.0  | McMurdie and Holmes 2013
gridExtra       | 2.3     | Auguie 2017 https://cran.r-project.org/web/packages/gridExtra/index.html
vegan           | 2.5.2   | Oksanen et al. 2018 https://CRAN.R-project.org/package=vegan
viridis         | 0.5.1   | Garnier 2018 https://CRAN.R-project.org/package=viridis
amap            | 0.8-16  | Lucas 2018 https://CRAN.R-project.org/package=amap
ggtree          | 1.8.2   | Yu et al. 2017
treeio          | 1.0.2   | Yu et al. 2017 https://guangchuangyu.github.io/treeio 
ape             | 5.1     | Paradis et al. 2004
adegenet        | 2.1.1   | Jombart and Ahmed 2011
indicspecies    | 1.7.6   | De Caceres and Legendre 2009

**Databases**:

Name            | Version | Date  | Download Location
----------------|-----------------|---------------|-------------------------------------------------------
NCBI Nucleotide | nt.gz   | Dec. 2016     | ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/
SILVA           | 128_SSURef_Nr99 | Mar. 2017     | http://ftp.arb-silva.de/release_128/Exports/
GreenGenes      | gg_13_8_otus    | Mar. 2016     | http://qiime.org/
MetaPhlAn 2     | mpa_v20_m200    | Nov. 2017     | https://bitbucket.org/biobakery/metaphlan2/downloads/
lac2strep       | custom  | 2017  | NA  | blastn database generated as described in the lac2 operon below section from lac2 sequences from: Richards, V. P., Choi, S. C., Pavinski Bitar, P. D., Gurjar, A. A., & Stanhope, M. J. (2013). Transcriptomic and genomic evidence for Streptococcus agalactiae adaptation to the bovine environment. BMC Genomics, 14(1), 920. http://doi.org/10.1186/1471-2164-14-920

**Reference Genomes**:

Species                                       | Strain      | Date           | Completeness | Type           | Source
----------------------------------------------|-------------|----------------|--------------|----------------|-----------------------------------------------------
_Homo sapiens_                                | HG19        | 2016-01-14     | Complete     | Reference      | http://hgdownload.cse.ucsc.edu/downloads.html#human
_Aggregatibacter aphrophilus_                 | W10433      | 2017-06-14     | Complete     | Representative | https://www.ncbi.nlm.nih.gov/genome
_Desulfobulbus_ sp. oral taxon 041            | Dsb1-5      | 2018-06-06     | Contigs      | Assembly       | https://www.ncbi.nlm.nih.gov/genome
_Fusobacterium nucleatum_ subsp. _nucleatum_  | ATCC 25586  | 2017-06-14     | Complete     | Representative | https://www.ncbi.nlm.nih.gov/genome
_Pseudopropionibacterium propionicum_         | F0230a      | 2017-06-14     | Complete     | Representative | https://www.ncbi.nlm.nih.gov/genome
_Rothia dentocariosa_                         | ATCC 17831  | 2017-06-14     | Complete     | Representative | https://www.ncbi.nlm.nih.gov/genome
_Streptococcus gordonii_                      | CH1         | 2017-06-14     | Complete     | Representative | https://www.ncbi.nlm.nih.gov/genome
_Treponema socranskii_                        | ATCC 35535  | 2018-05-31     | Scaffold     | Representative | https://www.ncbi.nlm.nih.gov/genome

Note: for species with plasmids or multiple chromosomes, these were concatenated
into a single multi-FASTA file.


### General Directory Structure
The general struture of this project is typically as follows (although
variants will occur):

```
project/
  00-documentation/
    00-document_1.txt
    01-document_2.txt
  01-data/
    raw_data/
      screening/
      deep/
    databases/
      <database_1>/
  02-scripts/
    00-ANALYSIS_CONFIG
    01-script.sh
    02-notebook.Rmd
  03-preprocessing
    screening/
      human_filtering/
        input/
        output/
      library_merging/
        input/
        output/
    deep/
      human_filtering
        input/
        output/
      library_merging/
        input/
        output/
  04-analysis/
    screening/
      analysis_1/
        input
        output
      analysis_2/
        input
        output
    deep/
      analysis_1/
        input
        output
      analysis_2/
        input
        output
```

### Variable Assignment

We assume `bash` and `R` are set as default paths.

Once you've downloaded and set up all the software listed above, we need to
prepare a profile that stores the paths to some of the programs that 
are required in some of the longer scripts we use (namely the human 
preprocessing, 16s preprocessing and MALT scripts).

Make a new file under `02-scripts` named `00-analysis_profile`. In here, we
want to now set program variables to each path. You will need to change the
path to your own installation path for each bit of software before saving!


```bash
## SOFTWARE PATHS
FASTQC=/projects1/tools/fastqc_0.11.4/fastqc
ADAPTERREMOVALFIXPREFIX=/projects1/tools/adapterremovalfixprefix/0.0.3/bin/AdapterRemovalFixPrefix
ADAPTERREMOVAL=/projects1/tools/adapterremoval/2.2.0-multiple_input_files/bin/AdapterRemoval
BWA=/projects1/tools/bwa-0.7.12/bwa
SAMTOOLS=/projects1/tools/samtools_1.3/bin/samtools
PICARDTOOLS=/projects1/tools/picardtools/1.140/bin/picard
DEDUP=/projects1/tools/dedup/0.12.2/bin/dedup
QUALIMAP=/projects1/tools/qualimap/2.2.1/bin/qualimap
DAMAGEPROFILER=/projects1/tools/damageprofiler/0.3.8/bin/damageprofiler

METAPHLAN2=/projects1/tools/metaphlan2/biobakery-metaphlan2-7898bf1/metaphlan2.py

MULTIVCFANALYZER=/projects1/tools/multivcfanalyzer/0.0.87/bin/multivcfanalyzer

## For project only scripts (give directory only)
SCRIPTS=/projects1/microbiome_calculus/evolution/02-scripts.backup

## For QIIME we will use multiple tools so we will just set the bin directory
QIIME=/projects1/tools/qiime-environment/1.9.1/bin
```

**For all other commands, we manually define the path to the software directly**
**in the code block as we go through the walkthrough.** So please ensure you
check each code block before submitting.



### Database Setup

To build the MALT nt database, once you have downloaded the NCBI Nucleotide (NT)
database in compressed FASTA format, we can build the corresponding indexed
database. The 'acc2tax' file can be downloaded from the MEGAN website
[here](http://ab.inf.uni-tuebingen.de/data/software/megan6/download/welcome.html)

For the SILVA database and HG19, we just need to download and run 
`bwa index` on the FASTA file.

```bash

DBDIR=/projects1/microbiome_calculus/evolution/01-data/databases

## MALT indexed NT database
mkdir -p $DBDIR/malt/raw $DBDIR/malt/indexed
cd $DBDIR/malt

### Download nucleotide database fasta and md5sum file into a database directory
wget ftp://ftp-trace.ncbi.nih.gov/blast/db/FASTA/nt.gz .
wget ftp://ftp-trace.ncbi.nih.gov/blast/db/FASTA/nt.gz.md5 .

### Generate the md5sum of the downloaded file, and comapre with contents of
### the NCBI .md5 version
md5sum nt.gz
cat nt.gz.md5

### Download into a different directory the accession to taxonomy mapping file
### as provided on the MEGAN6 website, and unzip
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/nucl_acc2tax-May2017.abin.zip
unzip nucl_acc2tax-May2017.abin.zip

### Now to build...
BWA=/projects1/tools/bwa-0.7.12/bwa
SAMTOOLS=/projects1/tools/samtools_1.3/bin/samtools
PICARDTOOLS=/projects1/tools/picardtools/1.140/bin/picard


sbatch \
-c 112 \
--mem 1950000 \
--partition=supercruncher \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=fail \
--mail-user=fellows@shh.mpg.de \
-J "MALT-Build-Nt_2017-10" \
--wrap="/projects1/malt/versions/malt040/malt-build \
--step 2 \
-i $DBDIR/malt/raw/nt.gz \
-s DNA \
-d $DBDIR/malt/indexed \
-t 112 -a2taxonomy $DBDIR/malt/raw/nucl_acc2tax-May2017.abin"

## BWA indexed SILVA database
mkdir $DBDIR/silva
cd !$

wget http://ftp.arb-silva.de/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz
"$BWA" index SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta
"$SAMTOOLS" faidx hg19_complete.fasta
"$PICARDTOOLS" CreateSequenceDictionary \
REFERENCE=SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta \
O=SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.dict

## BWA index HG19 - note I'm not 100% sure how this was downloaded as this was
## done centrally in our department, so I'm making a semi-guess here.
mkdir $DBDIR/hg19
cd !$

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*.fa.gz
cat *fa.gz >> hg19_complete.fasta
"$BWA" index hg19_complete.fasta
"$SAMTOOLS" faidx hg19_complete.fasta
"$PICARDTOOLS" CreateSequenceDictionary \
REFERENCE=hg19_complete.fasta \
O=hg19_complete.fasta.dict
```

For the custom NCBI Genome RefSeq database containing bacterial and archaea
assemblies at scaffold, chromosome and complete levels - we follow the the
notebook here: `/projects1/microbiome_sciences/reference_databases/refseq/genomes/bacteria_archea_homo_20181029/docs/refseq_genomes_bacteria_archaea_homo_complete_chromosome_scaffold_walkthrough_20181029.r`

To build the `aadder` database for functional analysis, we run the command
`/projects1/microbiome_calculus/evolution/02-scripts.backup/027-aadder_build_refseqCGS_bacarch_sbatch_script_20181104.sh`
Note we have to change the `MEGAN.vmoptions` to have a large enough memory
allocation.


Now, in our analysis_profile file, lets add the the database and reference
genome (DB) paths. Make sure again to change each path to the one your system
before saving.

```bash
## MALT DB Directory containing all database files
MALTDB=/projects1/malt/databases/indexed/index038/full-nt_2017-10
## SILVA DB directory containing the converted U to T FASTA file and associated bwa indexed files
SILVADB=/projects1/microbiome_sciences/reference_databases/silva/release_128_DNA/
## GreenGenes DB directory, as provided in QIIME
GREENGENESDB=/projects1/tools/qiime-environment/1.9.1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus

HG19REF=/projects1/Reference_Genomes/Human/HG19/hg19_complete.fasta
```

## Preprocessing

### Raw Data Structure
When downloading the raw sequencing data, the following structure is required 
where each sequenced library has a unique directory and all fastq files 
associated with that library from a given sequencing run is together. For 
example, a library sequenced twice on two different runs should be considered
separate (at this stage).

The data files should be in unmerged gziped FASTQ format.

```
data/
  raw/
    screening/
      library_1/ # Single End HiSeq
        library_1_L001_R1.fq.gz
      library_2/ # Paired End HiSeq
        library_2_L001_R1.fq.gz
        library_2_L001_R2.fq.gz
      library_3/ # Single End NextSeq
        library_3_L001_R1.fq.gz
        library_3_L002_R1.fq.gz
        library_3_L003_R1.fq.gz
        library_3_L004_R1.fq.g
      library_4/ # Paired End NextSeq
        library_4_L001_R1.fq.gz
        library_4_L001_R2.fq.gz
        library_4_L002_R1.fq.gz
        library_4_L002_R2.fq.gz
        library_4_L003_R1.fq.gz
        library_4_L003_R2.fq.gz
        library_4_L004_R1.fq.gz
        library_4_L004_R2.fq.gz
      library_1_20170131/# Resequenced Library 1
        library_1_L001_R1.fq.gz
    deep/
      library_1_deep/ # Single End HiSeq
      library_1_deep_L001_R1.fq.gz
      library_2_deep/ # Paired End HiSeq
        library_2_deep_L001_R1.fq.gz
        library_2_deep_L001_R2.fq.gz
```

Note that when a library has been sequenced twice, in some cases these FASTQ 
files will be in a unique directory with additional date tag in the name OR will
have a 'SG1.1' appended to it's name, where SG1 can alternatively be SG2 to 
indicate the second sequencing of the same library.

These resequenced libraries will be merged with the original sequencing below. 
It is also important to keep the screening and deep sequenced libraries separate
as they will be processed independently.

Furthermore, in some cases libraries were sequenced over multiple lanes
of a HiSeq. These are considered independent until FastQC verification has 
occured. This directories are appended with the lane number e.g. `_L6`.

Two additional directories were created here, one for all laboratory library 
and extraction controls, and another for FASTQ files needed for source tracker.

### Public Data

All raw FASTQ files should be downloaded sample specific directories in
 `01-data/public_data/raw`.

**Sourcetracker**

Sourcetracker source files were selected based on the following criteria:

  * Had to be shotgun metagenomes
  * Have had no modification or treatment made to DNA selection or host (e.g. no pesticide)
  * Must have been generated on the Illumina platform
  * Must have more than 10 million reads
  * Must have contain than 1000 16s rRNA reads in the shotgun data (detected 
  during analysis)
  
In addition the human microbiome project samples had the additional criteria of:

  * Must be from unique individuals (checked using the biospecimen column) 
  * Aim for approximately 50/50 male and female where possible

The sourcetracker files were downloaded from the NCBI SRA and EBI ENA 
databases. A list of these libraries can be seen in the 'mapping' file 
`02-microbiome_calculus-deep_evolution-individualscontrolssources_metadata.tsv`. Metadata on the HMP samples can be 
seen in `00-documentation/99-sourcemetadata-hmp_SraRunTable_allRuns.tsv, with 
samples selected for each source type coming from unique individuals, as 
inferred by the hmp_subject_id column. Metadata for the other datasets can be 
seen in the `00-documentation/99-sourcemetadata-*` files

A file containing the ERR and SRR numbers of each library on each line was 
given to the script `03-SRR_ERR_download_script.sh`. This script includes 
sbatch information and downloads each file sequentially producing the forward 
and reverse files. This script was then parallelised in 'v_0_2' of this script, 
working in the same way but requiring an additional number of cores parameter.

Some HMP project samples were not avaliable directly from the FTP server,
in which case these were directly pulled using either `prefetch -v`, which is 
a similar command as in the `03-SRR_ERR_download_script.sh` file, but without 
the `wget` step.

```bash
SRATOOLKIT=/projects1/tools/sratoolkit/2.8.0/sratoolkit.2.8.0-ubuntu64/bin

"$SRATOOKIT"/prefetch -v SRR514306
"$SRATOOKIT"/fastq-dump -F --split-files --readids /projects1/clusterhomes/fellows/ncbi/public/sra/SRR514306.sra --gzip --outdir .
```

These were then renamed using

```bash
cd /projects1/microbiome_calculus/evolution/01-data/public_data/raw/
rename s/_1.fastq.gz/_S0_L001_R1_000.fastq.gz/ */*.fastq.gz
rename s/_2.fastq.gz/_S0_L001_R2_000.fastq.gz/ */*.fastq.gz

```

The final fastq files were then placed symlinked to individual directories in
`/projects1/microbiome_calculus/evolution/01-data/public_data/prepped`

The same thing happened with the Slon data, so had to resort to 
downloading the FASTQ files directly. However, unfortunately, the uploaded
data was actually not 'raw' but the already merged data from the Slon 2017 paper. 
We will download the FASTQ data anyway and do a modified pre-processing.

I did this with the following command, and utilising the file generated from
 the parsing script of the two metadata files which is stored here
 `02-scripts/99-Slon2017_DataFinderScript.R`.

```bash
cat /projects1/microbiome_calculus/evolution/00-documentation.backup/99-Slon2017_AccessionsToDownload_2.tsv | while read LINE; do
  sbatch \
  -c 1 \
  --mem 1000 \
  -t 48:00:00 \
  -o ~/slurm_logs/slurm.%j.out \
  -e ~/slurm_logs/slurm.%j.err \
  -J "ENA_Download_$LINE" \
  --wrap="mkdir /projects1/microbiome_calculus/evolution/01-data/raw_data/screening/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev) && \
  wget $(echo $LINE | cut -d ';' -f1) -P /projects1/microbiome_calculus/evolution/01-data/raw_data/screening/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)/ && \
  wget $(echo $LINE | cut -d ';' -f2) -P /projects1/microbiome_calculus/evolution/01-data/raw_data/screening/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)/ && \
  wget $(echo $LINE | cut -d ';' -f3) -P /projects1/microbiome_calculus/evolution/01-data/raw_data/screening/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)/ && \
  rename s/.fastq.gz/_S0_L001_R1_000.merged.fastq.gz/ \
  /projects1/microbiome_calculus/evolution/01-data/raw_data/screening/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev).fastq.gz && \
  rename s/_1.fastq.gz/_S0_L001_R1_000.fastq.gz/ /projects1/microbiome_calculus/evolution/01-data/raw_data/screening/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)/*.fastq.gz && \
  rename s/_2.fastq.gz/_S0_L001_R2_000.fastq.gz/ /projects1/microbiome_calculus/evolution/01-data/raw_data/screening/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)/*.fastq.gz"
done
```

With these downloaded data, we need to do one more step to make them 
ready for preprocessing. As the data is already merged with singletons
we want to have put them in a single file.

```r
cat /projects1/microbiome_calculus/evolution/00-documentation.backup/99-Slon2017_AccessionsToDownload_2.tsv | while read LINE; do
  sbatch \
  -c 1 \
  --mem 1000 \
  -t 48:00:00 \
  -o ~/slurm_logs/slurm.%j.out \
  -e ~/slurm_logs/slurm.%j.err \
  -J "ENA_Download_$LINE" \
  --wrap="mkdir /projects1/microbiome_calculus/evolution/03-preprocessing/screening/human_filtering/input/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)/ && \
  cat /projects1/microbiome_calculus/evolution/01-data/raw_data/screening/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)/*.fastq.gz >> \
  /projects1/microbiome_calculus/evolution/03-preprocessing/screening/human_filtering/input/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)/$(echo $LINE | cut -d ';' -f1 | rev | cut -d/ -f 2 | rev)_S0_L001_R1_000.merged.fq.gz"
done
```

The final merged individual fastq files were moved to individual directories in 
`/projects1/microbiome_calculus/evolution/01-data/public_data/prepped`

**Additional Individuals**

Finally, we want to download the shotgun sequenced individuals from Weyrich et
al. 2017 (Nature).

For this, I downloaded the 'processed' files from the OAGR database, as the
original data had barcodes, and they used AdapterRemoval anyway.

The list of files and hard links can be seen in the documentation under
`99-publicdata-Weyrich_Neanderthals.txt`.

I downloaded each file with `wget`, renaming with the file as listed in the 
OAGR website, concatenated multi-file samples if required (i.e. ElSidron1) then 
renamed to the EAGER standard.

An example:

```bash
cd /projects1/microbiome_calculus/evolution/01-data/public_data/
mkdir ElSidron1

wget https://www.oagr.org.au/api/v1/dataset_file/3086/download
rename s/download/2NoAdapt_ELSIDRON1L7_lTACTG_rCTCGA_R1R2_Collapsed.fastq.gz/ download
wget https://www.oagr.org.au/api/v1/dataset_file/3109/download
rename s/download/2NoAdapt_ELSIDRON1L7_lTACTG_rCTCGA_R1R2_Collapsed_Truncated.fastq.gz/ download

cat *fastq.gz > ElSidron1_S0_L000_R1_000.fastq.merged.fq.gz

## If no concatenation required
# mv 2NoAdapt_CHIMP_150519_lACGTG_rATTGA_R1R2_Collapsed.fastq.gz Chimp.150519_S0_L000_R1_000.fastq.merged.fq.gz

```

The final merged individual fastq files were moved to individual directories in 
`/projects1/microbiome_calculus/evolution/01-data/public_data/prepped`. 

The just renamed files were just symlinked into the above.

### Data Preprocessing

#### Script Version

Now downloaded, we can begin preprocessing by performing a sequencing quality 
check, merge any of the PE data (as well as trimming low quality bases), remove 
any DNA that maps to the human genome (but preserving the statistics) and 
extract all-non human DNA for downstream processing.

We can do this with the script `01-preprocessing_human_filtering.sh`, or with 
on EAGER v1.92.55. 

**Note:** The reason why EAGER itself wasn't actually original for most
of the samples on this step was due to a few bugs in the pipeline when I 
needed to process them. The `01-preproce[...]` commands themselves are pulled
directly from a EAGER log, so it is indistinguishable from EAGER itself.

Below I provide a loop that makes sure to run the
script on libraries you have not already run that are present in the output
directory. All you need to do is change the paths in the variables `INDIR` and
`OUTDIR`. 

If you need to increase the number of cores of memory, you can
modify this in the `CPU` and `MEM` variables in the script itself.

You will need to check if this manual script is working correctly
manually, as I didn't have time to inbuild checks. You can do this by checking 
all the fields in the output of the next script (below) are filled with numbers.

```bash
INDIR=/projects1/microbiome_calculus/evolution/03-preprocessing/screening/human_filtering/input/hiseq_single
OUTDIR=/projects1/microbiome_calculus/evolution/03-preprocessing/screening/human_filtering/output
SCRIPTS=/projects1/microbiome_calculus/evolution/02-scripts.backup

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
    --wrap="$SCRIPTS/03-preprocessing_human_filtering_premerged.sh $OUTDIR/$LIBNAME $LIBDIR $LIBNAME"
    sleep 1
  fi
done

```

#### EAGER Version

For the Slon et al. 2017 and Weyrich et al. 2017 datasets, we ran EAGER without 
AdapterRemoval as the ENA uploaded data was already trimmed and merged. 

For the VLC/ARS/DLV/PLV/RIG/OFN samples, some the blanks and deep 
sequenced samples, I ran the same as below, but with FastQC and AdapterRemoval 
turned on with default settings. 

```
Organism: Human
Age: Ancient
Treated Data: non-UDG [Deep : UDG]
Pairment: Single End (VLC/EXB/LIB - Paired End)
Input is already concatenated (skip merging): Y (VLC/EXB/LIB - N)
Concatenate Lanewise together: N (VLC -Y, EXB/LIB - N)
Reference: HG19
Name of mitocondrial chromosome: chrMT
FastQC: Off (VLC/EXB/LIB/ARS - On)
AdapterRemoval: Off (VLC/EXB/LIB - On)
Mapping: BWA
  Seedlength: 32
  Max # diff: 0.01
  Qualityfilter: 0
  Filter unmapped: Off
  Extract Mapped/Unmapped: On
Remove Duplicates: DeDup
	Treat all reads as merged: On
Damage Calculation: Damageprofiler
CleanUp: On
Create Report: On

```

and submitted with

```bash
EAGER=/projects1/tools/eager/1.92.55/bin/eager
EAGERCLI=/projects1/tools/eager/1.92.55/bin/eagercli

for FILE in $(find /projects1/microbiome_calculus/evolution/03-preprocessing/screening/human_filtering/output/EXB* -name '*.xml'); do
  sbatch \
  -c 4 \
  --mem 32000 \
  -p medium \
  -o ~/slurm_logs/slurm.%j.out \
  -e ~/slurm_logs/slurm.%j.err \
  -t 48:00:00 \
  --mail-type=fail \
  --mail-user=fellows@shh.mpg.de \
  -J "EAGER_$(echo "$FILE" | rev | cut -d/ -f2 | rev)" \
  --wrap="unset DISPLAY && $EAGERCLI $(readlink -f $FILE)"
done

```

Deep Sequencing files were submitted using the array 
`021-eager_microbiota_slurm_array_deep.sh `

For the screening data to clean up the EAGER results directories so 
that they are in the same format as the others you can run the following couple 
of commands.

I didn't need to clean up the EAGER results for the deep sequenced data, so
skipped that.

```bash
cd /projects1/microbiome_calculus/evolution/03-preprocessing/screening/human_filtering/output

## Clean up
SAMPLE=LIB007.A0124
SAMTOOLS=/projects1/tools/samtools_1.3/bin/samtools

for DIR in $(readlink -f "$SAMPLE"); do
  cd "$DIR"
  mkdir fastqc
  mv 0-FastQC/*zip fastqc
  mv 1-AdapClip/*settings .
  mv 4-Samtools/*mapped.bam.stats .
  mv 4-Samtools/*extractunmapped.bam .
  mv 5-DeDup/*.mapped.sorted.cleaned.hist .
  mv 5-DeDup/*.log .
  mv 6-QualiMap/*/ qualimap
  mv 7-DnaDamage/*/ damageprofiler
done

cd ..

for DIR in $(readlink -f "$SAMPLE"); do
  cd "$DIR"
  "$SAMTOOLS" idxstats $(readlink -f 4-Samtools/*.mapped.sorted.bam) >> $(readlink -f 4-Samtools/*.mapped.sorted.bam).idxstats
  mv 4-Samtools/*idxstats .
  rm -r 0-FastQC
  rm -r 1-AdapClip
  rm -r 3-Mapper
  rm -r 4-Samtools
  rm -r 5-DeDup
  rm -r 6-QualiMap
  rm -r 7-DnaDamage
  rm DONE.CleanUpRedundantData
  rm EAGER.log
done


```

or

```bash
SAMPLES=($(find -name '2019*.xml' -type f -exec readlink -f {} \;))

SAMTOOLS=/projects1/tools/samtools_1.3/bin/samtools

for DIR in ${SAMPLES[@]}; do 
  cd $(dirname "$DIR")
  mkdir fastqc
  mv 0-FastQC/*zip fastqc
  mv 1-AdapClip/*settings .
  mv 4-Samtools/*mapped.bam.stats .
  mv 4-Samtools/*extractunmapped.bam .
  mv 5-DeDup/*.mapped.sorted.cleaned.hist .
  mv 5-DeDup/*.log .
  mv 6-QualiMap/*/ qualimap
  mv 7-DnaDamage/*/ damageprofiler
done

cd ..

for DIR in ${SAMPLES[@]}; do 
  cd $(dirname "$DIR")
  "$SAMTOOLS" idxstats $(readlink -f 4-Samtools/*.mapped.sorted.bam) >> $(readlink -f 4-Samtools/*.mapped.sorted.bam).idxstats
  mv 4-Samtools/*idxstats .
  rm -r 0-FastQC
  rm -r 1-AdapClip
  rm -r 3-Mapper
  rm -r 4-Samtools
  rm -r 5-DeDup
  rm -r 6-QualiMap
  rm -r 7-DnaDamage
  rm DONE.CleanUpRedundantData
  rm EAGER.log
done
```

#### Post Human Filtering Steps

##### BAM2FASTQ
For bam2fastq converstion of the deep sequenced data by running the script
`02-scripts.backup/029-bam2fastq_array.sh`


##### Statistics
We can then extract the statistics of this pre-processing with the script 
`05-statistics_human_filtering.sh`. Once checked, we move the resulting 
`human_filtering_statistics.csv` script to our `00-documents` folder and
rename as `03-human_filtering_statistics.csv`.

```bash
SCRIPTS=/projects1/microbiome_calculus/evolution/02-scripts.backup

cd /projects1/microbiome_calculus/evolution/03-preprocessing/screening/human_filtering

"$SCRIPTS"/005-statistics_human_filtering.sh output/

## Copy to documentation folder
mv human_filtering_statistics_"$(date +"%Y%m%d")".csv ..
cp ../human_filtering_statistics_"$(date +"%Y%m%d")".csv ../../../00-documentation.backup/03-human_filtering_statistics_"$(date +"%Y%m%d")".csv

```

NOTE: if you get a "Runtime error (func=(main), adr=6): Divide by zero" error,
don't worry, this is related to the cluster factor calculation for when there 
are no human reads after de-duplication. I just manually fill in with NA.

For the deep data, the EAGER table was used for report statistics (see below for
further information)

#### Library Merging

Next, we can merge together libraries that have been sequenced multiple times.

A list of libraries that have been merged together can be seen 
04-library_merging_information.csv, however in general this can be worked out 
by merging together any library that shares the first six character section of 
each library. For example:

```
CDC005.A0101
CDC005.A0101.170817

or

ECO002.B0101
ECO002.C0101
```

These represent distinct individuals (but different libraries of a single
sample or multiple sequencing of the same library), and for this study we are 
not considering differences in sampling source of the calculus.

Into our final output preprocessing directory (library_merging), we can quickly 
import the files not needing merging like so:

```bash
cd 03-preprocessing/screening/library_merging

## Make directory
for DIR in ../../screening/human_filtering/output/*/; do 
  mkdir "$(echo $DIR | rev | cut -d/ -f2 | rev)/"
done

## Make symlink of FastQ into new directory above.
for DIR in ../../screening/human_filtering/output/*/; do 
  ln -s $(readlink -f "$DIR")/*fq.gz "$(echo $DIR | rev | cut -d/ -f2 | rev)/"; 
done
```

Now remove those 'extra' libraries that will be collapsed into the 
first entry 

```bash
rm \
ABM007.A0101/*.fq.gz \
ABM008.A0101/*.fq.gz \
BIT001.A0101/*.fq.gz \
CDC005.A0101/*.fq.gz \
CDC011.A0101/*.fq.gz \
DJA006.A0101/*.fq.gz \
DJA007.A0101/*.fq.gz \
EBO008.A0101/*.fq.gz \
EBO009.A0101/*.fq.gz \
ECO002.B0101/*.fq.gz \
ECO004.B0101/*.fq.gz \
FDM001.A0101/*.fq.gz \
GDN001.A0101/*.fq.gz \
IBA001.A0101/*.fq.gz \
IBA002.A0101/*.fq.gz \
LOB001.A0101/*.fq.gz \
MOA001.A0101/*.fq.gz \
MTK001.A0101/*.fq.gz \
MTM003.A0101/*.fq.gz \
MTM010.A0101/*.fq.gz \
MTM011.A0101/*.fq.gz \
MTM012.A0101/*.fq.gz \
MTM013.A0101/*.fq.gz \
MTS001.A0101/*.fq.gz \
MTS002.A0101/*.fq.gz \
MTS003.A0101/*.fq.gz \
TAF017.A0101/*.fq.gz \
TAF018.A0101/*.fq.gz \
WAL001.A0101/*.fq.gz
```

These extras are merged into a single FASTQ file and placed in an 
independent file with the following example command:

```bash
INDIR=/projects1/microbiome_calculus/evolution/03-preprocessing/screening/human_filtering/output
OUTDIR=/projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging

cat "$INDIR"/ABM007*/*.fq.gz >> "$OUTDIR"/ABM007.A0101/ABM007_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/ABM008*/*.fq.gz >> "$OUTDIR"/ABM008.A0101/ABM008_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/BIT001*/*.fq.gz >> "$OUTDIR"/BIT001.A0101/BIT001_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/CDC005*/*.fq.gz >> "$OUTDIR"/CDC005.A0101/CDC005_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/CDC011*/*.fq.gz >> "$OUTDIR"/CDC011.A0101/CDC011_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/DJA006*/*.fq.gz >> "$OUTDIR"/DJA006.A0101/DJA006_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/DJA007*/*.fq.gz >> "$OUTDIR"/DJA007.A0101/DJA007_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/EBO008*/*.fq.gz >> "$OUTDIR"/EBO008.A0101/EBO008_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/EBO009*/*.fq.gz >> "$OUTDIR"/EBO009.A0101/EBO009_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/ECO002*/*.fq.gz >> "$OUTDIR"/ECO002.B0101/ECO002_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/ECO004*/*.fq.gz >> "$OUTDIR"/ECO004.B0101/ECO004_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/FDM001*/*.fq.gz >> "$OUTDIR"/FDM001.A0101/FDM001_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/GDN001*/*.fq.gz >> "$OUTDIR"/GDN001.A0101/GDN001_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/IBA001*/*.fq.gz >> "$OUTDIR"/IBA001.A0101/IBA001_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/IBA002*/*.fq.gz >> "$OUTDIR"/IBA002.A0101/IBA002_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/LOB001*/*.fq.gz >> "$OUTDIR"/LOB001.A0101/LOB001_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MOA001*/*.fq.gz >> "$OUTDIR"/MOA001.A0101/MOA001_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MTK001*/*.fq.gz >> "$OUTDIR"/MTK001.A0101/MTK001_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MTM003*/*.fq.gz >> "$OUTDIR"/MTM003.A0101/MTM003_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MTM010*/*.fq.gz >> "$OUTDIR"/MTM010.A0101/MTM010_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MTM011*/*.fq.gz >> "$OUTDIR"/MTM011.A0101/MTM011_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MTM012*/*.fq.gz >> "$OUTDIR"/MTM012.A0101/MTM012_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MTM013*/*.fq.gz >> "$OUTDIR"/MTM013.A0101/MTM013_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MTS001*/*.fq.gz >> "$OUTDIR"/MTS001.A0101/MTS001_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MTS002*/*.fq.gz >> "$OUTDIR"/MTS002.A0101/MTS002_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/MTS003*/*.fq.gz >> "$OUTDIR"/MTS003.A0101/MTS003_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/TAF017*/*.fq.gz >> "$OUTDIR"/TAF017.A0101/TAF017_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/TAF018*/*.fq.gz >> "$OUTDIR"/TAF018.A0101/TAF018_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz
cat "$INDIR"/WAL001*/*.fq.gz >> "$OUTDIR"/WAL001.A0101/WAL001_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz


```

For screening data, the final number of reads going downstream analysis is also 
recorded in the file `04-samples_library_merging_information.csv`. For libraries 
sequenced twice, I manually added the two values together.

For deep sequencing data, The same concatenating of multiple samples and/or 
lanes was done for the deep sequencing samples. However the statistics were 
summarised across replicates using the R notebook 
`099-eager_table_individual_summarised.Rmd`.


### PolyG Trimming Assessment

#### Human

The human DNA GC content could be a bit off in some of the new libraries 
generated in this study, as we are sequencing with Illumina NextSeqs, which 
have a 2 colour chemistry that considers no light emission as a 'G'. Thus, 
empty or finished clusters can be read as long poly G reads - which can still 
map to the human genome.

To get improved human DNA content calculations, we can run EAGER to get the 
mapped reads. Then run `fastp` on the bam2fastq mapped reads and use their 
`--trim_poly_g` function to remove complete G reads or remove from reads that 
has the last 10 reads as Gs this tail. This procedure is recorded in 
`04-analysis/screening/eager` under the `polyGremoval*` directories, but the 
unprocessed reads as initial input as I deleted the mapped files previously.

The polyG issue would not likely affect our reference microbial data because
these should be so short that each read would have both ends of the read
entirely sequenced (the whole read with both adapters or most of both adapters 
would fit in 75 cycles). Thus, adapter removal would remove anything that 
comes after the adapters which would be the poly Gs. Poly G tails would 
only affect single index reads, where the read itself is too short and doesn't
have an adapter to indicate the read has ended.

For the first round of mapping we use the input parameters per sequencing 
configuration, and the mapping pipeline as follows: 

```
Organism: Human
Age: Ancient
Treated Data: non-UDG
Pairment: Single End (R1) or Paired End (R1/R2)
Input is already concatenated (skip merging): N
Concatenate Lanewise together: N (HiSeq) or Y (NextSeq)
Reference: HG19
Name of mitocondrial chromosome: chrMT

AdapterRemoval: Y (New Data)
Mapping: BWA
  Seedlength: 32
  Max # diff: 0.01
  Qualityfilter: 0
  Filter unmapped: On
  Extract Mapped/Unmapped: Off
Remove Duplicates: DeDup
  Treat all reads as merged: On
CleanUp: On
Create Report: Off

```

We then update the `find` path and number of jobs in the array in the script 
`18-eager_microbiota_slurm_array.sh` to the new input directory and number
of libraries.

Once completed, we can perform a cleanup FASTQs to reduce our footprint 

```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input

rm -r /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/*.collapsed.gz \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/*.collapsed.truncated.gz \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/*.discarded.gz \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/*.pair1.truncated.gz \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/*.pair2.truncated.gz \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/*.settings \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/*.singleton.truncated.gz \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/DONE.AdapterRemovaldefault \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/DONE.AdapterRemovalFixReadPrefix \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/DONE.CombineFastQ \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/DONE.TrackFastQ \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/1-AdapClip/track_fastq.log \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/3-Mapper/ \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/4-Samtools/ \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/6-QualiMap \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/*.bam \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/*.bai \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/*.hist \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/*.log \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/*.bam.mtnucratio \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/DONE.CleanSam \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/DONE.DeDup \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/DONE.MTToNucRatioCalculator \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/DONE.SamtoolsIndexDedup \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/*/5-DeDup/DONE.SamtoolsSortDeDup
```

and now use samtools to convert our mapped Human DNA reads (after DeDup) back to 
FASTQ.

```bash
INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/

for LIBDIR in "$INDIR"/*/; do
  LIBNAME=$(echo "$LIBDIR" | rev | cut -d/ -f2 | rev)
  if [ -f "$INDIR"/"$LIBNAME"/5-DeDup/*mappedonly.sorted.cleaned_rmdup.sorted.bam.fq.gz ]; then
    printf "\n $LIBNAME already Processed \n\n"
  else
    echo "$LIBNAME"
    sbatch \
  -c 1 \
  --mem 4G \
  --partition=short \
  -t 0-02:00 \
  -o ~/slurm_logs/slurm.%j.out \
  -e ~/slurm_logs/slurm.%j.err \
  --mail-type=fail \
  --mail-type=time_limit \
  --mail-user=fellows@shh.mpg.de \
  -J "bam2fq_$LIBNAME" \
  --wrap="samtools fastq $(readlink -f $INDIR/$LIBNAME/5-DeDup/*mappedonly.sorted.cleaned_rmdup.sorted.bam) | gzip > $(readlink -f $INDIR/$LIBNAME/5-DeDup/*mappedonly.sorted.cleaned_rmdup.sorted.bam).fq.gz"
  fi
done

```

Then we run the complexity filter on these to remove the poly-Gs using 
[fastp](https://github.com/OpenGene/fastp). This will remove/trim any read
where the last 10 reads are all Gs.

```bash

## When running for the first time - REMEMBER TO UPDATE NUMBER OF ENTRIES IN ARRAY
sbatch /projects1/microbiome_calculus/iberian/02-scripts.backup/99-plotGcomplexity_filter.sh

## If running with extra samples
INDIR=/projects1/microbiome_calculus/iberian/03-preprocessing/screening/human_filtering/output

for LIBDIR in "$INDIR"/*/; do
  LIBNAME=$(echo "$LIBDIR" | rev | cut -d/ -f2 | rev)
  if [ -f "$INDIR"/"$LIBNAME"/4-Samtools/*.polyGtrimmed.fq.gz ]; then
    printf "\n $LIBNAME already Processed \n\n"
  else
    echo "$LIBNAME"
    sbatch \
  -c 1 \
  --mem 4G \
  --partition=short \
  -t 0-02:00 \
  -o ~/slurm_logs/slurm.%j.out \
  -e ~/slurm_logs/slurm.%j.err \
  --mail-type=fail \
  --mail-type=begin \
  --mail-type=time_limit \
  --mail-user=fellows@shh.mpg.de \
  -J "fastp-ployG_$LIBNAME" \
  --wrap="/projects1/users/fellows/bin.backup/fastp/fastp \
  -i $(readlink -f $INDIR/$LIBNAME/4-Samtools/*mapped.sorted.bam.fq.gz) \
  -o $(readlink -f $INDIR/$LIBNAME/4-Samtools/*mapped.sorted.bam.fq.gz).polyGtrimmed.fq.gz \
  --trim_poly_g \
  --disable_quality_filtering \
  --disable_length_filtering \
  --disable_adapter_trimming"
  fi
done
```

Once completed, we prepare the final round of EAGER mapping.

```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_output/
mkdir input output
cd input

find /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/polyGremoval_input/ -name '*.polyGtrimmed.fq.gz' -type f | while read LINE; do 
  mkdir $(echo $LINE | rev | cut -d/ -f3 | rev) && ln -s "$LINE" $(echo $LINE | rev | cut -d/ -f3 | rev); 
done

```

And set up EAGER with the following settings

```
Organism: Human
Age: Ancient
Treated Data: non-UDG
Pairment: Single End 
Input is already concatenated (skip merging): Y
Concatenate Lanewise together: N
Reference: HG19
Name of mitocondrial chromosome: chrMT

AdapterRemoval: N
Mapping: BWA
  Seedlength: 32
  Max # diff: 0.01
  Qualityfilter: 0
  Filter unmapped: On
  Extract Mapped/Unmapped: Off
Remove Duplicates: DeDup
  Treat all reads as merged: Y
Damage Calculation: Y
CleanUp: On
Create Report: Off

```

## PRESERVATIONAL ASSESSMENT

### 16s Extraction and QIIME Clustering

We firstly want to make an assessment of how well preserved the calculus
microbiome signal is in our archaeological and museum samples. The first step
here is to run the program Sourcetracker, however this typically uses an OTU 
table of 16s reads for the proportion estimation of the sources in your sample.

The script `05-preprocessing-16s_mapping` works very much the same way as with 
the human DNA preprocessing, but mapping to the SILVA database and converting 
the mapped only reads to FASTA format with a particular header format 
(>sample.name_1) that works with QIIME.


```bash
INDIR=/projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging
OUTDIR=/projects1/microbiome_calculus/evolution/03-preprocessing/screening/silva_16s_reads
SCRIPTS=/projects1/microbiome_calculus/evolution/02-scripts.backup


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
    -J "16SMAP_$LIBNAME" \
    --wrap="$SCRIPTS/009-preprocessing_16s_mapping.sh $OUTDIR/$LIBNAME $LIBDIR $LIBNAME"
    sleep 1
  fi
done
```

A SLURM array version for this is in the file 
`010-preprocessing_16s_mapping_slurmarray.sh`

To get the number of 16s rRNA reads that mapped, we can run the script 
'06-16s_extraction_statistics.sh', and copy the resulting file into our 
documents directory with

```bash
SCRIPTS=/projects1/microbiome_calculus/evolution/02-scripts.backup

cd /projects1/microbiome_calculus/evolution/03-preprocessing/screening

## Generate stats
"$SCRIPTS"/011-16s_extraction_statistics.sh silva_16s_reads/

## Copy into documentation folder
mv silva_16s_reads/16s_extraction_statistics_"$(date +"%Y%m%d")".csv .
cp 16s_extraction_statistics_"$(date +"%Y%m%d")".csv ../../00-documentation.backup/05-16s_extraction_statistics_"$(date +"%Y%m%d")".csv
```

While we have these reads, we still don't know what taxa they are derived from. 
For this we can use QIIME to cluster the reads by similarity then assign a
taxonomic 'name'. We will stick with OTU clustering rather than the more 
recent (and more powerful/statistically valid) ASV (amplicon sequence variant) 
methods. This is because we are not dealing directly with amplicon data, and we 
has also have damage which would affect reliable of assignment. Finally, the 
more recent version of QIIME, QIIME2, has been made much less flexible for 
non-amplicon data and I couldn't work out how to adapt the data. For our 
rough preservational screening purposes using QIIME 1 is sufficient.

Firstly, we have to make a combined FASTA file containing all the 16s reads 
from all the samples.

```bash
INDIR=/projects1/microbiome_calculus/evolution/03-preprocessing/screening/silva_16s_reads
OUTDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/input

rm "$OUTDIR"/silva_16s_reads_concatenated.fna.gz
touch "$OUTDIR"/silva_16s_reads_concatenated.fna.gz

for SAMPLE in "$INDIR"/*; do
  find "$SAMPLE"  -maxdepth 1 -name '*_renamed.fa.gz' -exec cat {} >> \
  "$OUTDIR"/silva_16s_reads_concatenated.fna.gz \;
done
```

For OTU clustering tself we need to define some parameter files that have 
been adapted for shotgun data at LMAMR in Oklahoma. So in a file named 
'10-params_CrefOTUpick.txt' that is in the '02-scripts' folder, add the 
following lines:

```
pick_otus:max_accepts    100
pick_otus:max_rejects    500
pick_otus:word_length    12
pick_otus:stepwords    20
pick_otus:enable_rev_strand_match    True
```

To actually run the clustering analysis and generate our OTU table we need to 
do the following. 

1. Install the QIIME conda environment (make sure to have set up conda previously, see: http://qiime.org/install/install.html)

```bash
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
```

2. activate the conda QIIME environment

```bash
conda activate qiime1
```

3. Run the OTU clustering!

```bash
INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/input
OUTDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/output
GREENGENESDB=/projects1/users/fellows/bin.backup/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/
SCRIPTS=/projects1/microbiome_calculus/evolution/02-scripts.backup

gunzip "$INDIR/silva_16s_reads_concatenated.fna.gz"

sbatch \
-c 16 \
--mem 128000 \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
-p medium \
-t 47:00:00 \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fellows@shh.mpg.de \
-J "QIIME_Pick_OTUs_Closed" \
--wrap="pick_closed_reference_otus.py \
-i $INDIR/silva_16s_reads_concatenated.fna \
-o $OUTDIR/otu_picking \
-a \
-O 16 \
-r $GREENGENESDB/rep_set/97_otus.fasta \
-t $GREENGENESDB/taxonomy/97_otu_taxonomy.txt \
-p $SCRIPTS/010-qiime_1_9_1_params_CrefOTUpick.txt"


```

To deactivate (note yous below need to reactive for all QIIME steps below)

```bash
source deactivate qiime1
```

Once finished we can then re-gzip the input fasta file with:

```bash
INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/input

pigz -p 4 "$INDIR/silva_16s_reads_concatenated.fna"
```

To get some basic statistics about the OTU picking we can use the BIOM package
that came with the QIIME environment.

```bash
conda activate qiime

INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/input
OUTDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/output

cd "$OUTDIR"

biom summarize-table -i "$OUTDIR"/otu_picking/otu_table.biom >> \
"$OUTDIR"/otu_picking/otu_table_summary.txt

## Check with: 
head "$OUTDIR"/otu_picking/otu_table_summary.txt -n 20
```

Sourcetracker (I realised later) doesn't do rarefaction properly as it 
allows sampling with replacement of the OTUS. Therefore, we need to remove
samples that have less than 1000 OTUs (which is the default rarefaction level
in sourcetracker - see below). So next we need to filter our OTU table
to remove samples that have less than that threshold. Looking at the table
summary shows that fortunately this will remove very few samples and 
will remove mostly blanks.

 ```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/output/otu_picking
INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/output/otu_picking

filter_samples_from_otu_table.py \
-i "$INDIR"/otu_table.biom \
-o "$INDIR"/otu_table_1000OTUsfiltered.biom \
-n 1000


biom summarize-table -i "$INDIR"/otu_table_1000OTUsfiltered.biom >> \
"$INDIR"/otu_table_1000OTUsfiltered_summary.txt

## Check with: 
head "$INDIR"/otu_table_1000OTUsfiltered_summary.txt -n 20
 ```

Next, to try and 'standardise' the files for SourceTracker, we can 
rarefy down the samples. Although statistically 'unforgiveable' (see McMurdie 
and Holmes 2016 2014 PLoS Comp. Bio), in this case we aren't using the data for
for any high-rigour analysis, and is an acceptable compromise for validation.

For filtering to Genus (L6) level:

```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/output/otu_picking
INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/output/otu_picking

## For Genus
summarize_taxa.py \
-i "$INDIR"/otu_table_1000OTUsfiltered.biom \
-o "$INDIR" \
-a \
-L 6

sed s#\'Solanum#Solanum#g otu_table_1000OTUsfiltered_L6.txt > otu_table_1000OTUsfiltered_L6_cleaned.tsv


```

With this final OTU table, we can now run Sourcetracker

### Sourcetracker

Sourecetracker requires a OTU table (generated above) and a metadata file
that tells the program what libraries in the OTU are a 'sink' or a 'source'.
This metadata file in this case is recorded here, 
`02-microbiome_calculus-deep_evolution-individualscontrolssources_metadata_<DATE>.tsv`,
which I try to use for all downstream analysis. In particular here we need to 
ensure we have a 'Env' and a 'SourceSink' column. 

```bash
## change to new directory for output directories
INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/qiime/output/otu_picking
OUTDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/sourcetracker.backup/

MAPPINGFILE=/projects1/microbiome_calculus/evolution/00-documentation.backup/02-calculus_microbiome-deep_evolution-individualscontrolssources_metadata_20190429.tsv

sbatch \
-c 2 \
--mem 8000 \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
-t 48:00:00 \
-p long \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fellows@shh.mpg.de \
--export=ALL \
-J "Sourcetracker" \
--wrap="Rscript \
/projects1/microbiome_calculus/evolution/04-analysis/screening/sourcetracker.backup/sourcetracker-1.0.1/sourcetracker_for_qiime.r \
-i "$INDIR"/otu_table_1000OTUsfiltered_L6_cleaned.tsv \
-m "$MAPPINGFILE" \
-o "$OUTDIR"/otu_table_L6_1000_"$(date +"%Y%m%d")" \
-r 1000 \
-v"

```

For plotting of these, I use the following R notebook to summarise the 
results in stacked plots: `12-sourcetracker_stackedbarplots_evo_<DATE>.Rmd`. 
Within this notebook, we also perform a hard threshold filtering to remove 
samples that potentially have the largest identified proportion as sediment or 
skin. The  list of samples to be removed from the analysis is saved in the file 
`06-sourcetracker_samples_to_filter_out.csv` in the `00-documentation` 
directory.

### MALT

#### MALT Running

However, as our libraries are shotgun sequenced, we are only utilising a very
small fraction of our samples that have 1000 16s reads. We can also perform
taxonomic binning to assign all the reads to different taxa, and use this for
preservational assessment and basic compositional analysis.

For this binning we will use MALT, with a wrapper script for efficient 
submission. The settings as set in 
`07-malt-genbank-nt_2017_2step_85pc_supp_0.01` are as follows: a read requires 
a minimum of 85% sequence identity (-id), it will only retain a second hit 
if the alignment is within 1% of the bit-score of the top scoring hit (-top), 
it will only keep a leaf-node on the tree if it has more than 0.01 of the hits
over all the hits in the sample (-supp). It will only retain 100 possible hits
per node. Note that you may have to change the `INSTALL4J_JAVA_HOME` path 
within the script, as I'm currently not sure what that does. This script aligns
to the NCBI Nucleotide database (nt).

We can then run our taxonomic binning job with the following command, providing
a input directory with wildcards to all the files and an output directory.

```bash
sbatch \
-c 112 \
--mem 1900000 \
--partition=supercruncher \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fellows@shh.mpg.de \
-J "MALT040-FullNT_2017" \
--wrap="/projects1/microbiome_calculus/evolution/02-scripts.backup/07-malt-genbank-nt_2017_2step_85pc_supp_0.01 /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/*/*.fq.gz /projects1/microbiome_calculus/evolution/04-analysis/screening/malt/nt"
```

For the RefSeq database, see AADDER.

#### Summary Statistics

To extract the summary statistics for all the MALT runs, we can run
the following on the MALT logs.

```bash
grep -e "Loading MEGAN File:" \
-e "Total reads:" \
-e "With hits:" \
-e "Alignments:" \
-e "Assig. Taxonomy" \
-e "Min-supp. changes" \
-e "Numb. Tax. classes:" \
-e "Class. Taxonomy:" \
-e "Num. of queries:" \
-e "Aligned queries:" \
-e "Num. alignments:" \
-e "MinSupport set to:" \
/projects1/microbiome_calculus/evolution/04-analysis/screening/malt/nt/*log | cut -d":" -f 2-99 > /projects1/microbiome_calculus/evolution/00-documentation.backup/99-maltAlignedReadsSummary_raw_nt_$(date "+%Y%m%d").txt
```

Then we do a few clean up and calculation steps as in the R script `99-MALT_Summary_statistics.R`,
the output of which is recorded in `00-documentation.backup/05-MALT_taxonomic_binning_summary_statistics_<DATE>.tsv`.

The MinSupport value column(s) was then manually added to the indiviudals column of
our main metadata file. 


#### Effects of Removal of Human Reads

At the beginning of the project I decided to filter out any reads mapping to 
human based on the assumption that 
  a) a lot of reference genomes in NCBI are full of human reads, 
  b) it would remove a lot of stretches of uninformative
reads such as poly-A reads and with the additional justification of 
  c) that it would make my FASTQ files smaller and faster to analyse (as we 
  aren't interested in the human content anyway).

However, I decided to see exactly how much of a different it makes. So I 
selected a bunch of samples in three approximate cateogies of amount of 
human DNA (high: >10%, medium: >1%, and low <1). For an approximate level
of consistency I tried to aim for 2 Gorillas, 2 Neanderthals and 2 Humans (
Chimps were generally overall low, so couldn't identify any with very high
endogenous DNA).

The selected samples can be found in the documentation file:
`99-SamplesSelectedHumanRemovalMALTTest.csv`.

I then used the clipped-and-merged FASTQ files from the output of the 
polyG removal test  step 1 (so finding the human DNA by mapping, _prior_ 
to the polyG removal itself.) and ran MALT, as above.

```bash
sbatch \
-c 112 \
--mem 1900000 \
--partition=supercruncher \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fellows@shh.mpg.de \
-J "MALT040-FullNT_2017" \
--wrap="/projects1/microbiome_calculus/evolution/02-scripts.backup/07-malt-genbank-nt_2017_2step_85pc_supp_0.01 \
/projects1/microbiome_calculus/evolution/04-analysis/screening/malt/nt_withhumanreads/input/*.gz \
/projects1/microbiome_calculus/evolution/04-analysis/screening/malt/nt_withhumanreads/output"
```

Analysis of the output of this test can be seen in the notebook 
`099-MALT_HumanReadsRemovalEffectsTest.Rmd`

### MEGAN

Once MALT has completed, we need to generate an OTU table of all the alignments.
 
For this we need to use the GUI tool MEGAN6 locally on your desktop. To 
generate the OTU table, open the MALT RMA6 files in MEGAN with absolute counts, 
and ignore all unassigned reads.

To uncollapse the tree: Tree > Rank > <Taxonomic level>, then select species 
nodes with: Select > Rank > <Taxonomic level>. 

Now go File > Export > Text (CSV) format, select taxonName_to-count, summarised 
and tab as a separator.

For a 'microbial-only' OTU table (basically not prokaryotes or synthetic DNA 
sequences), we can before uncollapsing the tree select 
'collapse non-prokaryotes' under the 'Tree' menu and then do Select > Rank > 
<taxonomic level> to select only Bacteria and Archaea. We also need to export 
the same data with as a tree from MEGAN with the option: file > export > tree, 
and save as newick. This should also be done for each non-prokaryote taxonomic 
level. 


### Cumulative Proportion Decay Plots

The OTU table itself does not give us much information about the oral signature.

Instead I came up with a simple visualisation to show how abundant the oral
signal in the samples are. This visualisation needs two things: an OTU table 
from MEGAN at species level and a database of taxa with their 'sources'. 

We have already generated the OTU above.

For the database, you can follow the steps as recorded in 
`09-Organism_Isolation_Source_Database_Generation.Rmd`. The database also 
requires manual curation over time. This database is stored under
`00-documentation` as `07-master_oralgenome_isolationsource_database_<DATE>.tsv`

With these two things, the R notebook 
`10-cumulative_proportion_decay_curves_<DATE>.Rmd` shows you how to generate
the visualisation.

### DECONTAM

Tp remove possible contaminating species from our samples that are derived from
the laboratory environment, we can use the R package `decontam`. The idea here 
is to use this to reduce the number of noisey taxa in the downstream 
compositional analysis, e.g. reducing the overplotting of the loadings of PCAs.

Firstly we need to obtain the summed quant values for all our individuals, as
we merge the samples/libraries during pre-processing. This script uses the 
labwork data which is stored externally from SDAG/CDAG in Dropbox. The script
is in the scripts folder under '14-decontam_quant_preparations.R' The output
file is in `99-decontam_quants.tsv`. We then manually add these values to our
main metadata file 
`02-microbiome_calculus-deep_evolution-individualscontrolssources_metadata_<MALT_DB>_<DATE>.tsv`.

We then run Decontam following the Decontam tutorial in the script 
`14-decontam.Rmd`


## Compositional Analysis

### MetaPhlAn2

To begin, based on Velsko et al. (2018, mSystems), in addition to MALT 
MetaPhLan2 performs pretty well in recalling a synthetic oral microbiome
dataset from damaged shotgun data, so we will use this for our dataset.

For the moment we will keep all viruses, eukaryotes, bacteria and archaea. 
In the future we may remove e.g. eukaryotes.

```bash
METAPHLAN2=/projects1/tools/metaphlan2/biobakery-metaphlan2-7898bf1/metaphlan2.py

INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/input
OUTDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/output

for LIBFILE in "$INDIR"/*/*.fq.gz; do
  LIBNAME=$(echo "$LIBFILE" | rev | cut -d/ -f2 | rev)
  if [ -d "$OUTDIR"/"$LIBNAME" ]; then
    printf "\n $LIBNAME already Processed \n\n"
  else
    mkdir "$OUTDIR"/"$LIBNAME"
    sbatch \
    -c 2 \
    --mem 16000 \
    -o ~/slurm_logs/slurm.%j.out \
    -e ~/slurm_logs/slurm.%j.err \
    -t 02:00:00 \
    -p short \
    --mail-type=fail \
    --mail-type=time_limit \
    --mail-user=fellows@shh.mpg.de \
    --export=ALL \
    -J "MP2_$LIBNAME" \
    --wrap="$METAPHLAN2 $LIBFILE \
    -o $OUTDIR/$LIBNAME/$LIBNAME.mp2profile.tsv \
    --input_type fastq \
    --nproc 2 \
    -t rel_ab"
  fi
done

## if re-running with different parameters, here changing fastq to bowtie2out

METAPHLAN2=/projects1/tools/metaphlan2/biobakery-metaphlan2-7898bf1/metaphlan2.py

INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/input
OUTDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/output

for LIBFILE in "$INDIR"/*/*bowtie2out.txt; do
  LIBNAME=$(echo "$LIBFILE" | rev | cut -d/ -f2 | rev)
  if [ -d "$OUTDIR"/"$LIBNAME" ]; then
    printf "\n $LIBNAME already Processed \n\n"
  else
    mkdir "$OUTDIR"/"$LIBNAME"
    sbatch \
    -c 2 \
    --mem 16000 \
    -o ~/slurm_logs/slurm.%j.out \
    -e ~/slurm_logs/slurm.%j.err \
    -t 02:00:00 \
    -p short \
    --mail-type=fail \
    --mail-type=time_limit \
    --mail-user=fellows@shh.mpg.de \
    --export=ALL \
    -J "MP2_$LIBNAME" \
    --wrap="$METAPHLAN2 $LIBFILE \
    -o $OUTDIR/$LIBNAME/$LIBNAME.mp2profile.tsv \
    --input_type bowtie2out \
    --nproc 2 \
    -t rel_ab"
  fi
done

## Re-running to get estimated read counts with rel_ab to rel_ab_w_read_stats

METAPHLAN2=/projects1/tools/metaphlan2/biobakery-metaphlan2-7898bf1/metaphlan2.py
INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/input
OUTDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/output_readcounts

for LIBFILE in "$INDIR"/*/*bowtie2out.txt; do
  LIBNAME=$(echo "$LIBFILE" | rev | cut -d/ -f2 | rev)
  if [ -d "$OUTDIR"/"$LIBNAME" ]; then
    printf "\n $LIBNAME already Processed \n\n"
  else
    mkdir "$OUTDIR"/"$LIBNAME"
    sbatch \
    -c 2 \
    --mem 4000 \
    -t 02:00:00 \
    -p short \
    -o ~/slurm_logs/slurm.%j.out \
    -e ~/slurm_logs/slurm.%j.err \
    --mail-type=fail \
    --mail-type=time_limit \
    --mail-user=fellows@shh.mpg.de \
    --export=ALL \
    -J "MP2_$LIBNAME" \
    --wrap="$METAPHLAN2 $LIBFILE \
    -o $OUTDIR/$LIBNAME/$LIBNAME.mp2profile_readstats.tsv \
    --input_type bowtie2out \
    --nproc 2 \
    -t rel_ab_w_read_stats"
  fi
done


## Re-running to get actual rad counts with rel_ab to rel_ab_w_read_stats

METAPHLAN2=/projects1/tools/metaphlan2/biobakery-metaphlan2-7898bf1/metaphlan2.py
INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/input
OUTDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/output_readmappedcounts

for LIBFILE in "$INDIR"/*/*bowtie2out.txt; do
  LIBNAME=$(echo "$LIBFILE" | rev | cut -d/ -f2 | rev)
  if [ -d "$OUTDIR"/"$LIBNAME" ]; then
    printf "\n $LIBNAME already Processed \n\n"
  else
    mkdir "$OUTDIR"/"$LIBNAME"
    sbatch \
    -c 2 \
    --mem 4G \
    -t 02:00:00 \
    -p short \
    -o ~/slurm_logs/slurm.%j.out \
    -e ~/slurm_logs/slurm.%j.err \
    --mail-type=fail \
    --mail-type=time_limit \
    --mail-user=fellows@shh.mpg.de \
    --export=ALL \
    -J "MP2_$LIBNAME" \
    --wrap="$METAPHLAN2 $LIBFILE \
    -o $OUTDIR/$LIBNAME/$LIBNAME.mp2profile_readsmappedstats.tsv \
    --input_type bowtie2out \
    --nproc 2 \
    -t reads_map"
  fi
done
```

To merging all of the MetaPhlAn2 percentage profiles of all the libraries

```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/output

## Note, remove metaphlan2.py script in the variable here as using util scripts
METAPHLAN2=/projects1/tools/metaphlan2/biobakery-metaphlan2-7898bf1/

INDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/input
OUTDIR=/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/output

"$METAPHLAN2"/utils/merge_metaphlan_tables.py "$OUTDIR"/*/*mp2profile.tsv > "$OUTDIR"/mp2_merged_abundance_table_all_"$(date +%Y%m%d)".txt
```

For merging of the estimated read count files, run the 
`016-metaphlan2_readcount_table_generator.Rmd` notebook.


For the mapped reads, method, the output actually the name of the read and 
a given taxonomic ID. For this we run the following:

```bash

for i in */*; do 
  echo $(echo "$i" | cut -d "/" -f 1) $(zcat "$i" | wc -l) ; 
done > mp2_merged_readsmapped_table_all_"$(date +%Y%m%d)".txt

```

Note that all of those files need to be -1 because the count includes a header.

Finally, some read statistics by applying the same 0.01% threshold used in MALT
are gained via the `99-MetaPhlan2_Summary_statistics.R` script. These were
then manually added to the metadata file.

### Raw OTU Tables

Raw MALT OTU tables with and without bad samples and at different min-support
values are generated by the Notebook `016-MALT_otutable_generation.Rmd`. These
tables are stored as `.tsv` files in the `04-analysis/screening/megan.backup` 
directory.

### PC(o)A

To explore if we have a structure in our data that can describe differences 
between each group we want to explore, we can perform a PCA or PCoAs to reduce
the variation between the samples to human-readable dimensions.

We perform a Phylogenetic Isometric-Log-Ratio transform (Silverman et al. 2017 
_Frontiers in XXX_), for 'proper' CoDa ordination as described in the notebook
`22-PhILR_PC0A.Rmd`.

Due to the notebook ending up having lots of options, I also used `knitr::purl()`
to create a Script version that accepts arguments.

To generate the script

```{r}
knitr::purl("/home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190130.Rmd", "/home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190130_script.R", documentation = 2)

```

Then cycle through all variants of the script I would like to run

```bash
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt genus withSources withControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt genus withSources withControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt genus noSources withControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt genus noSources withControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt genus withSources noControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt genus withSources noControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt genus noSources noControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt genus noSources noControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt species withSources withControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt species withSources withControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt species noSources withControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt species noSources withControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt species withSources noControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt species withSources noControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt species noSources noControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R nt species noSources noControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq genus withSources withControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq genus withSources withControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq genus noSources withControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq genus noSources withControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq genus withSources noControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq genus withSources noControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq genus noSources noControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq genus noSources noControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq species withSources withControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq species withSources withControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq species noSources withControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq species noSources withControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq species withSources noControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq species withSources noControls out 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq species noSources noControls in 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190320_script.R refseq species noSources noControls out 12




```


Again but changing the sample filtering method

```bash
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls in none 1
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out sourcetracker 1
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out onepcburnin 1
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out twopcburnin 1
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out fivepcburnin 1
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out tenpcburnin 1
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out withinvariation 1

Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls in none 5
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out sourcetracker 5
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out onepcburnin 5
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out twopcburnin 5
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out fivepcburnin 5
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out tenpcburnin 5
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out withinvariation 5

Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls in none 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out sourcetracker 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out onepcburnin 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out twopcburnin 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out fivepcburnin 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out tenpcburnin 12
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190408_script.R nt genus withSources withControls out withinvariation 12



```


Finally by cycling through the MinSupport values

```bash
for i in 1 2 3 4 5 6 7; do
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190415_script.R nt genus withSources withControls out withinvariation "$i"
done

for i in 1 2 3 4 5 6 7; do
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190415_script.R nt genus withSources withControls in withinvariation "$i"
done

```

Final one

```bash
for i in 7; do
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt genus withSources withControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt genus noSources withControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt genus withSources noControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt genus noSources noControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt genus withSources withControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt genus noSources withControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt genus withSources noControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt genus noSources noControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq genus withSources withControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq genus noSources withControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq genus withSources noControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq genus noSources noControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq genus withSources withControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq genus noSources withControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq genus withSources noControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq genus noSources noControls out withinvariation "$i"
done

for i in 4; do
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt species withSources withControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt species noSources withControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt species withSources noControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt species noSources noControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt species withSources withControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt species noSources withControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt species withSources noControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R nt species noSources noControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq species withSources withControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq species noSources withControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq species withSources noControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq species noSources noControls in withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq species withSources withControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq species noSources withControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq species withSources noControls out withinvariation "$i"
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/017-PhILR_PCoA_20190430_script.R refseq species noSources noControls out withinvariation "$i"
done

```


The R script for summarising the results across all runs is named
`017-PhILR_PCoA_summaries_<DATE>.R`, and the output files are in 
  `00-documentation` under `philr_permanova_*_summary_<DATE>.tsv`.


--- OLD ----

For this we can use our MALT OTU table (see below), perform zero-count 
correctation and centered-log-ratio normalisation, as recommend by Gloor et al. 
(2017, _Front. In Microbio._) as described in 
`16-compositional_pca_20180907.Rmd`.

This is alternatively implemented for MetaPhlan2 in the `CLR_MetaPhlan2.Rmd` 
notebook.

However, as there are mathematical issues with CLR, we additionally perform
an Phylogenetic Isometric-Log-Ratio transform (Silverman et al. 2017 
_Frontiers in XXX_), for 'proper' 'CoDa' ordination as described in the notebook
`22-PhILR_PC0A.Rmd`.

----

### Core Microbiome Selection

This is performed in the notebook `018-CoreMicrobiome.Rmd`.

A script version is also provided as `018-CoreMicrobiome_20190318_script.R`

This was then run with the following commands, switching database, minimum
support, fraction of individuals per population a taxon needs to be in to be 
core to that population, and the fraction of populations a taxon needs to be 
core in to be core to the genus.

```bash
## NT tests
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 1 0.8 0.66 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 2 0.8 0.66 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.8 0.66 F

## Refseq tests
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "refseq" 1 0.8 0.66 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "refseq" 2 0.8 0.66 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "refseq" 5 0.8 0.66 F

## core threshold tests, with 0.8 and 0.66 our 'standard' - mutate ind threshold
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.5 0.66 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.66 0.66 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.8 0.66 F

## core threshold tests, with 0.8 and 0.66 our 'standard' - mutate pops threshold
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.8 0.5 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.8 0.65 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.8 0.75 F

## core threshold tests, with 0.8 and 0.66 our 'standard' - mutate pops threshold
## But after selection of 0.5 as ind threshold
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.5 0.5 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.5 0.66 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190319_script.R "nt" 5 0.5 0.75 F

## Now tweaking minsupport again
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190320_script.R "nt" 5 0.5 0.66 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190320_script.R "nt" 6 0.5 0.66 F
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190320_script.R "nt" 7 0.5 0.66 F
```

```bash
## Final being minSupport of X, 0.5 inds of population, and 66 pops of genus, 
## minsupport 0.7 for genus and 0.4 for species, and testing whether dropping
## single individual populations makes a difference.

for i in 4 7; do
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190429_script.R "nt" "$i" 0.5 0.66 F
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190429_script.R "nt" "$i" 0.5 0.66 T
done

for i in 4 7; do
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190429_script.R "refseq" "$i" 0.5 0.66 F
  Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/018-CoreMicrobiome_20190429_script.R "refseq" "$i" 0.5 0.66 T
done


```

Due to the large number of variants of parameters prodocued in the above 
notebook, an R script was used to condense all values from tables to a single 
summary table. The R script is `018-CoreMicrobiome_summaries_<DATE>.R`, and the 
output files are in `00-documentation` under 
`23-intersection_proktaxapassingthresholdsstats_20190211_<DATE>.tsv`.

### MaltExtract

To further verify the authenticity the calculus samples - we can also run 
maltExtract on a subset of the core taxa to check for damage patterns and 
short fragment lengths. 

We took the species level 'core' of the Anthropoids, Hominids and Homininis 
based on the MALT nt database, with a minsupport value of 0.04, requiring 
a taxon being in 50% of each population having a taxon to be core to the 
population, and 66% of the populations to be core to the genus. We retain
single individual populations, and calculate intersections between each
host genus. 

This list is then given to maltExtract as follows.

```bash
sbatch \
-c 32 \
--mem 256G \
-p long \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
-t 0-12:00 \
--mail-type=all \
--mail-user=fellows@shh.mpg.de \
-J "MaltExtract_Evo" \
--wrap="which MaltExtract && MaltExtract \
--destackingOff \
-f def_anc \
-i /projects1/microbiome_calculus/evolution/04-analysis/screening/maltExtract/input/all/ECO*rma6 \
-o /projects1/microbiome_calculus/evolution/04-analysis/screening/maltExtract/output/test/ \
-p 32 \
-r /projects1/clusterhomes/huebler/RMASifter/RMA_Extractor_Resources \
-t /projects1/microbiome_calculus/evolution/04-analysis/screening/maltExtract/taxon_lists/08b-coremicrobiome_presenceabsence_upsettable_allsoftware_maltdbnt_0.04_fracinds0.5_fracpops0.66_singleindpopsdroppedF_hostgenus_species_AnthropoidsHominidaeHoiminini_20190509.tsv \
-v"

```

We then run the HOPS postprocessing script with

```bash
 Rscript /projects1/clusterhomes/huebler/RMASifter/AMPS/PostProcessing/amps-master/postprocessing.AMPS.r -r AnthropoidsHominidaeHoiminini_core_20190509/ -m def_anc -t 4 -n ../taxon_lists/08b-coremicrobiome_presenceabsence_upsettable_allsoftware_maltdbnt_0.04_fracinds0.5_fracpops0.66_singleindpopsdroppedF_hostgenus_species_AnthropoidsHominidaeHoiminini_20190509.tsv 
```

We also run maltExtract on just the samples and taxa we perform for the 
phylogenies.

```bash
sbatch \
-c 32 \
--mem 256000M \
-p long \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=all \
--mail-user=fellows@shh.mpg.de \
-J "MaltExtract_Evo" \
--wrap="/projects1/tools/java/jdk-9.0.4/bin/java -Xmx255G -jar /projects1/clusterhomes/huebler/RMASifter/RMAExtractor_jars/MaltExtract1.5.jar \
--destackingOff \
-f def_anc \
-i /projects1/microbiome_calculus/evolution/04-analysis/screening/maltExtract/input/all/*rma6 \
-o /projects1/microbiome_calculus/evolution/04-analysis/screening/maltExtract/output/PhylogenyTaxa_20190715/ \
-p 32 \
-r /projects1/clusterhomes/huebler/RMASifter/RMA_Extractor_Resources \
-t /projects1/microbiome_calculus/evolution/04-analysis/screening/maltExtract/taxon_lists/phylogeny_list.txt \
-v"

```

### Indicator Analysis

For compostional analysis we often have the issue we have way too
many taxa to easily visualise. As a way of improving this we can run 'Indicator 
Analysis' to select taxa that are indicators of a particular group or 
group-of-groups, based on the taxa being prevelent across the group and at high 
abundance compared to other groups. This is performed following the script 
`21-Indicator_analysis.R`

Submitted as follows:

```bash
sbatch \
-c 8 \
--mem 64G \
-t 48:00:00 \
-p medium \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=all \
--mail-user=fellows@shh.mpg.de \
--export=ALL \
-J "IndicSpecies" \
--wrap="Rscript /projects1/microbiome_calculus/evolution/02-scripts.backup/21-Indicator_analysis.R"

```

This saves a file in the documentation which has a list of the species, the 
group that the species has a high IndVal and significance value (to that group
combination) of <0.05. It also has a column whether the species is considered
an indicator species core (to a single group), This file is saved as
`15-Evolution IndicSpecies - MP2_Species_HostCommon_cleaned_20180817.tsv`. The
list of taxa in there can be then used as a filtering variable, when we run 
the analysis below.

### Heatmap

To visualise the whole MetaPhlAn2 dataset, we can either run the default 
MetaPhlAn2 heatmap utility script as follows.

```bash
## Note, in METAPHLAN2 variable remove metaphlan2.py script here as using util scripts
METAPHLAN2=~/.bin/biobakery-metaphlan2-5c606797370a

INDIR=/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/input
OUTDIR=/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/metaphlan2/output

python "$METAPHLAN2"/utils/metaphlan_hclust_heatmap.py \
-c bbcry \
-d euclidean \
-f euclidean \
--top 40 \
--minv 0.1 \
--tax_lev g \
-s log \
--in "$OUTDIR"/mp2_merged_abundance_table_nosource_"$(date +%Y%m%d)".txt \
--out "$OUTDIR"/mp2_merged_abundance_heatmap_nosource_top40_genus_"$(date +%Y%m%d)".svg
```

Note that here I have had to use local version of MetaPhlAn2, as something
screwy on our cluster that would cause a wierd problem due to one of the
libraries. Furthermore I had to switch from bray-curtis default to euclidean 
distances due to a 'ValueError: Linkage 'Z' contains negative indices.' error. 
This switch was recommended by the MetaPhlAn2 developer here -
 https://groups.google.com/forum/#!topic/metaphlan-users/h3nQqpnsz84

**Alternatively**, I reconstructed the heatmap script above but in R, giving us
more flexibility (grouping by our known host groups, removing lab contaminants 
and bad samples etc.). This can be run using the notebook here: 
`17-compositional_heatmaps_20180906.Rmd`


## Genome Reconstruciton

### Screening Data

#### Mapping to Reference

For intitial genome reconstruction assessment, EAGER was used to map
to the genomes of different oral microbiota of interest, based on
species that were found to be relatively consistant across two or more
Host species, as reported in the MetaPhlAn2 species heatmap.

The Intitial species that were selected are:

  * Pseudopropionibacterium propionicum
  * Treponema socranskii (<- only scaffolds)
  * Rothia dentocariosa 
  * Desulfobulbus sp. oral taxon 041 (<- only contigs)
    - For this selected the sp. oral taxon 041 assembly with the largest  
    assembly size, Desulfobulbus sp. oral taxon 041 str. Dsb1-5. 
    - Downloaded from the Genbank FTP server the contigs from here: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/403/865/GCA_000403865.1_Dsb1-5
  * Fusobacterium nucleatum

And for others of interest

  * Aggregatibacter aphrophilus
  * Streptococcus gordonii

In addition added red complex bacteria at request from Tina

 * Treponema denticola
 * Tannerella forsythia
 * Treponema denticola

Reference sequences were downloaded from NCBI Genome using the reference or
representative strain. The downloaded FASTAs were then run through the following
indexing commands

```bash
### Now to build...
BWA=/projects1/tools/bwa-0.7.12/bwa
SAMTOOLS=/projects1/tools/samtools_1.3/bin/samtools
PICARDTOOLS=/projects1/tools/picardtools/1.140/bin/picard



"$BWA" index <reference>.fa
"$SAMTOOLS" faidx <reference>.fa
"$PICARDTOOLS" CreateSequenceDictionary REFERENCE=<NAME>.fa O=<NAME>.fa.dict

```

For testing purposes, two sets of mapping parameters were used to each 
different reference genome. These two mapping parameters are "relaxed" (-n 0.1, 
seed length 32) and "strict" (-n 0.01, -l 16).

Otherwise, the config files are set up as follows:

```
Organism: Bacteria/Other
Age: Ancient
Treated Data: non-UDG
Pairment: Single End
Input is already concatenated (skip merging): Y

Reference: [listed above]
Name of mitocondrial chromosome: [none]

FastQC: Off
AdapterRemoval: Off
Mapping: BWA
 	Seedlength: 16 or 32
 	Max # diff: 0.01 or 0.1
 	Qualityfilter: 37
 	Filter unmapped: On
 	Extract Mapped/Unmapped: Off
Complexity Estimation: Off
Remove Duplicates: DeDup
	Treat all reads as merged: On
Damage Calculation: Damageprofiler
SNP Calling: Unified Genotyper
	Emit All Sites: On
CleanUp: On
Create Report: On
```

and submitted with the following array in the script `13-eager_microbiota_slurm_array.sh`

```bash
#!/bin/bash
 
#SBATCH -n 4                                                            # number of CPUs (here 4)
#SBATCH --mem 32000                                                     # memory pool for all cores (here 8GB)
#SBATCH -t 0-01:00                                                      # walltime (D-HH:MM) (here 0 days, 8hours, 0 minutes)
#SBATCH --partition=short                                              # partition (queue) to submit to 
#SBATCH --array=0-535%8
#SBATCH -o /projects1/clusterhomes/fellows/slurm_logs/slurm.%j.out      # STDOUT (the standard output stream)
#SBATCH -e /projects1/clusterhomes/fellows/slurm_logs/slurm.%j.err      # STDERR (the output stream for errors)
#SBATCH --mail-type=fail                                                # notifications for job to abort
#SBATCH --mail-type=time_limit                                          # notifications for job to abort
#SBATCH --mail-use=fellows@shh.mpg.de                                   # these notifications will be sent to your shh.mpg.de email-address.
#SBATCH -J "EAGER_relaxed_SlurmArray"                                           # name of job


FILES=($(find -L /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output/relaxed -name '*xml' -type f))
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

EAGERCLI=/projects1/tools/eager/1.92.55/bin/eagercli
$EAGERCLI ${FILENAME}
```

```bash
sbatch /projects1/microbiome_calculus/evolution/02-scripts.backup/13-eager_microbiota_slurm_array.sh
```

Note, that there is a maximum of 1001 jobs (by default) that can be submitted 
with the array.

Therefore I run the relaxed/strict directories separately.

After running, the script `99-eager_microbiotamapping_parameter_comparison.Rmd` 
is used to summarise the mapping results, in particular showing the difference
between the two mapping parameters. 

Summaries of DamageProfiler output can be seen in 24-DamageProfiler_Summary.Rmd

#### PreSeq

For deeper sequencing calculations, we wanted to do an estimation of what
our sequencing limits would be, where we start reaching extremely high 
cluster factors (i.e. duplication rate, or sequencing the same DNA fragment
again). To get this limit, run on each EAGER output folder the script in 
`99-sbatch_preseq array.sh`

Then download the PDFs in each species folder to your local computer 

```bash
scp sdag:/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output/strict/tannerella_*/*/*_preseq.pdf .
```

and extract the number of reads until the 2.5 and 4.0 cluster factors are reached

```bash
for i in */*pdf; do echo "$i" "$(pdftotext "$i" - | grep CF:2.5)" $(pdftotext "$i" - | grep CF:4.0); done
```

Copy and paste this output into a spreadsheet and convert to a table.

Merge this table into EAGER report Table using the following R commands, changging
the paths in `eager_data`. Where `cf_data` is the data pulled from the pdftotext 
command above. 

```r
library(tidyverse)
eager_data <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/00-documentation.backup/10-oral_microbiota_mapping_eager_results_20180814_strictonly.tsv")
cf_data <- read_tsv("~/Downloads/Evolution Deeper Sequencing Calculuation - Felix_Preseq_CF.tsv")


final_data <- left_join(eager_data, cf_data)

write_tsv(final_data, "~/Downloads/Evolution_Deeper_Sequencing_Calculus_Bacterial_CF_merged.tsv")
```

#### MultiVCFAnalyzer

Next we want to check if the EAGER mappings have resulted in cross-mapping 
from other strains or species, and also then perform some phylogenies. Our 
input data will be the SNP table from MultiVCFAnalyzer.

To prepare the input data in a format MultiVCFAnalyzer works well in, as it 
takes the sample name from the folder of the VCF file, not the sample name.

```bash
mkdir /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer
cd !$
mkdir input output

## Make input directories and symlink VCF files
for i in $(find /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output/ -name '*.vcf.gz' -type f); do
  PARAM="$(echo $i | rev | cut -d/ -f 5 | rev)"
  SPECIES="$(echo $i | rev | cut -d/ -f 4 | rev)"
  SAMPLE="$(echo $i | rev | cut -d/ -f 3 | rev)"
  mkdir -p /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/input/"$PARAM"/"$SPECIES"/"$SAMPLE"/
  ln -s "$i" /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/input/"$PARAM"/"$SPECIES"/"$SAMPLE"/ 
done



## Make output directories

cd output 
  mkdir -p multiallelic_allowed/relaxed multiallelic_allowed/strict
  mkdir -p multiallelic_notallowed/relaxed multiallelic_notallowed/strict
cd ..


## Make species output directories
for i in /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output/*/*; do
  for j in /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output/relaxed/*/; do 
    mkdir "$i"/"$(echo "$j" | rev | cut -d/ -f2 | rev)" 
  done
done

```

We want to export allele frequencies, set a genotyping quality of 30, a minimum 
of 3X coverage, homozygous calls require a frequency of 90% at a position to be 
reported and 'heterozygous' calls require 10% to be reported as in Vågene et 
al. 2018 (Nat. Eco. Evo.). This setting is recorded below as 
'multiallelic_allowed'

For the phylogenies (multiallelic_notallowed) with confident SNP calls we 
switch 'het' calls to 0.9 as well.

I couldn't quickly come up with a convienent way of automating submission
of these, so for the moment I just copy and past the script below, changing
the variables, and the het call flag (just before the second NA).

i.e. I change species for each species, then for each of these change between
the relaxed/strict EAGER mapping parameters. and also multallelic_allowed 
and *_notallowed (but for the latter two also changing the het call 
from 0.1 for allowed to 0.9 for notallowed).

```bash
MULTIVCFANALYZER=/projects1/tools/multivcfanalyzer/0.0.87/bin/multivcfanalyzer
ANALYSIS="multiallelic_allowed"
SPECIES="fusobacterium_nucleatum"
PARAM="relaxed"
FILES=($(find -L /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/input/$PARAM/$SPECIES -name '*.vcf.gz' -type f))

sbatch \
-c 4 \
--mem 32G \
--partition=short \
-t 02:00:00 \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=fail \
--mail-user=fellows@shh.mpg.de \
-J "MultiVCFAnalyzer_$PARAM_$SPECIES" \
--wrap="$MULTIVCFANALYZER \
NA \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/references/${SPECIES^}/*.fa \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/references/${SPECIES^}/*.gff \
/projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output/$ANALYSIS/$PARAM/$SPECIES \
T \
30 \
3 \
0.9 \
0.1 \
NA \
${FILES[*]}

find /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output/$ANALYSIS/$PARAM/$SPECIES -name '*fasta' -type f -exec pigz -p 4 {} \;
find /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output/$ANALYSIS/$PARAM/$SPECIES -name '*tsv' -type f -exec pigz -p 4 {} \;
gunzip /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output/$ANALYSIS/$PARAM/$SPECIES/snpStatistics.tsv.gz"

```

Unfortunately MultiVCFAnalyzer does not accept multi-fastas references currently,
so _Treponema socranskii_ and _Desulfobulbus sp. oral taxon 41_ cannot be
used here.

To assess for cross mapping, we can plot the SNP frequencies of 
'heterozygous plots', as described in `19-multi_allelic_SNP_plots.Rmd`. The 
plots themselves can be found in 
`/projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/plots`


#### Low Coverage Phylogenies

One problem we have is we get very very low coverage of most of the species.
It as been shown however, it is still possible to get topology information
it seems even with 0.2X coverage (see _Y. pestis_ section in Damgaard et al. 
2018, Nature)

The block below is how Marcel has done this for his data.

```
fellows [14:31]
Hey Marcel,
how did you do your phylogenies with really low coverage data?

mkeller [14:37]
i did it once now for a sample with 0.2 or so
did alexander tell you?
it was non-udg, so i ran it with non-udg parameters through eager to get max. coverage, not caring about damage after that, i used the bam file from the gatk folder and turned it into a fastq file with bedtools

fellows [14:39]
after RMDUP?

mkeller [14:40]
yes?geom
then i multiplied it 5x and concatenated the files
i ran the 5xfastq file again through eager, just mapping and snp calling, but without duplicate removal

fellows [14:41]
What do you mean by multiplied it 5x?

mkeller [14:41]
i literally multiplied it 5 times
copypaste
:smile:

fellows [14:42]
You mean concatenated the same FASTQ to itself 5 times?

mkeller [14:42]
yes

fellows [14:42]
Ok, understood

mkeller [14:42]
i’m sure you know more elegant ways :wink:

fellows [14:42]
Nope :wink:

mkeller [14:43]
for the phylogenetic tree, i ran multivcfanalyzer with min 5x coverage after copying the vcf-file to a folder with outgroup_ in the name, this way the unique snps (most of it damage) will be removed

mkeller [14:46]
initially it was meant for real outgroups like y. pseudotuberculosis, because the root would be unnecessarily long

but you can do it with multiple samples
then all singletons will be removed, only homoplastic sites will remain

```


To implement what is described:

```bash

SAMTOOLS=/projects1/tools/samtools_1.3/bin/samtools

mkdir /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/input_haplotypecaller_concat2x
cd !$
mkdir strict relaxed

## make species directories
    for i in ../output/relaxed/*/; do 
      mkdir relaxed/"$(basename $i)" strict/"$(basename $i)"; 
    done

mkdir /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output_haplotypecaller_concat2x
cd !$
mkdir strict relaxed

## make species directories
    for i in ../output/relaxed/*/; do 
      mkdir relaxed/"$(basename $i)" strict/"$(basename $i)"; 
    done

## Set species and mapping parameters! SWITCH HERE!
species="pseudopropionibacterium_propionicum"
mappingparam="relaxed"

## Find the post dedup-bams normally used for GATK and convert to FASTQ
sbatch \
-c 16 \
--mem 64000 \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fellows@shh.mpg.de \
-p short \
-t 02:00:00 \
-J "BAM2FASTQ_$species_$mappingparam" \
--wrap="find /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output/$mappingparam/$species/ -name '*.sorted.cleaned_rmdup.sorted.bam' -type f | parallel -j 16 '$SAMTOOLS fastq {} | gzip > {}.fq.gz'"

## Concatenate the newly created FASTQs onto themselves 5 times for 5x coverage (or 2 times for 2x coverage)
## but put in the directoeies we made above 
sbatch \
-c 16 \
--mem 64000 \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fellows@shh.mpg.de \
-p short \
-t 02:00:00 \
-J "FASTQ_CONCAT_$species_$mappingparam" \
--wrap="find /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output/$mappingparam/$species/ -name '*.sorted.cleaned_rmdup.sorted.bam.fq.gz' -type f | parallel -j 16 'cat {} {} > {}_concat2x.fq.gz'"

## Move these into a new EAGER input directory
find \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output/$mappingparam/$species/ \
-name '*_concat2x.fq.gz' | \
parallel \
-j 16 \
'mkdir /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/input_unifiedgenotyper_concat2x/$(echo {}| rev | cut -d/ -f5 | rev)/$(echo {}| rev | cut -d/ -f4 | rev)/$(echo {}| rev | cut -d/ -f3 | rev) && \
mv {} /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/input_unifiedgenotyper_concat2x/$(echo {}| rev | cut -d/ -f5 | rev)/$(echo {}| rev | cut -d/ -f4 | rev)/$(echo {}| rev | cut -d/ -f3 | rev)'


```

Then generated EAGER configs with the same parameters above, but into a new 
EAGER output directory (here: output_lowcoverage_mappings_concat5x) and with 
only modules Mapping (e.g. relaxed (bwa mapping -l 32, -n 0.01) and SNP Calling 
(EMIT_ALL_Sites) turned on.

We can submit as an array as before, by modifying the 
`13-eager_microbiota_slurm_array.sh` to the corresponding directory and 
the number of jobs in the array range to that of the number of config files


```bash
#!/bin/bash
 
#SBATCH -n 4                                                            # number of CPUs (here 4)
#SBATCH --mem 32000                                                     # memory pool for all cores (here 8GB)
#SBATCH -t 0-01:00                                                      # walltime (D-HH:MM) (here 0 days, 8hours, 0 minutes)
#SBATCH --partition=short                                              # partition (queue) to submit to 
#SBATCH --array=0-106%8
#SBATCH -o /projects1/clusterhomes/fellows/slurm_logs/slurm.%j.out      # STDOUT (the standard output stream)
#SBATCH -e /projects1/clusterhomes/fellows/slurm_logs/slurm.%j.err      # STDERR (the output stream for errors)
#SBATCH --mail-type=fail                                                # notifications for job to abort
#SBATCH --mail-type=time_limit                                          # notifications for job to abort
#SBATCH --mail-use=fellows@shh.mpg.de                                   # these notifications will be sent to your shh.mpg.de email-address.
#SBATCH -J "EAGER_relaxed_SlurmArray"                                           # name of job


FILES=($(find -L /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output_lowcoverage_mappings_concat5x/aggregatibacter_aphrophilus/relaxed -name '*xml' -type f))
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

EAGERCLI=/projects1/tools/eager/1.92.55/bin/eagercli
$EAGERCLI ${FILENAME}
```

After running, the script `99-eager_microbiotamapping_parameter_comparison.Rmd` 
is used again to summarise the mapping results, in particular showing the 
difference between the two mapping parameters. When not comparing parameters,
a reduced version of this `18-eager_reporttable_concatenator.R` can be used.

Then to prepare for MultiVCFAnalyzer


```{r}
cd /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer
mkdir input_unifiedgenotyper_concat2x output_unifiedgenotyper_concat2x

## Make input directories and symlink VCF files
for i in $(find /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output_unifiedgenotyper_concat2x/relaxed/ -name '*.vcf.gz' -type f); do
  PARAM="$(echo $i | rev | cut -d/ -f 5 | rev)"
  SPECIES="$(echo $i | rev | cut -d/ -f 4 | rev)"
  SAMPLE="$(echo $i | rev | cut -d/ -f 3 | rev)"
  mkdir -p /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/input_unifiedgenotyper_concat2x/"$PARAM"/"$SPECIES"/"$SAMPLE"/
  ln -s "$i" /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/input_unifiedgenotyper_concat2x/"$PARAM"/"$SPECIES"/"$SAMPLE"/ 
done



## Make output directories
cd output_unifiedgenotyper_concat2x
  mkdir -p multiallelic_notallowed/relaxed multiallelic_notallowed/strict
cd ..


## Make species output directories
for i in /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output_unifiedgenotyper_concat2x/*/*; do
  for j in /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output/relaxed/*/; do 
    mkdir "$i"/"$(echo "$j" | rev | cut -d/ -f2 | rev)" 
  done
done

```

Then run MultiVCFAnalyzer as follows, switching the variable as 
needed.

```bash
MULTIVCFANALYZER=/projects1/tools/multivcfanalyzer/0.0.87/bin/multivcfanalyzer
ANALYSIS="multiallelic_notallowed"
SPECIES="pseudopropionibacterium_propionicum"
PARAM="relaxed"
FILES=($(find -L /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/input_unifiedgenotyper_concat2x/$PARAM/$SPECIES/ -name '*.vcf.gz' -type f))

sbatch \
-c 4 \
--mem 32G \
--partition=short \
-t 02:00:00 \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=fail \
--mail-user=fellows@shh.mpg.de \
-J "MultiVCFAnalyzer_$PARAM_$SPECIES" \
--wrap="$MULTIVCFANALYZER \
NA \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/references/${SPECIES^}/*.fa \
/projects1/microbiome_calculus/evolution/04-analysis/screening/eager/references/${SPECIES^}/*.gff \
/projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output_unifiedgenotyper_concat2x/$ANALYSIS/$PARAM/$SPECIES \
T \
30 \
2 \
0.9 \
0.9 \
NA \
${FILES[*]}

find /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output_unifiedgenotyper_concat2x/$ANALYSIS/$PARAM/$SPECIES -name '*fasta' -type f -exec pigz -p 4 {} \;
find /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output_unifiedgenotyper_concat2x/$ANALYSIS/$PARAM/$SPECIES -name '*tsv' -type f -exec pigz -p 4 {} \;
gunzip /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output_unifiedgenotyper_concat2x/$ANALYSIS/$PARAM/$SPECIES/snpStatistics.tsv.gz"


```

To prepare for phylogenetic analysis, the Rmarkdown script
`20-eager_microbiota_phylogeny_alignerReducer_20180915_ISBA.Rmd` is used. This 
filters the alignments to those mappings that have at least 10% percent 
coverage at 1x.

#### MEGA Phylogenies

With the filtered alignments produced above, performed neighbour joining
phylogenies as in MEGAX with the following settings.

In MEGAPROTO set up a config file as follows:

1. Data input file type: Nucleotide
2. Phylogeny > Construct/Test Nighbout Joining-tree
3. Test of phylogeny: Bootstrap
4. No. of Bootstrap Replications: 100
5. Model/Method: Jukes-Cantor model
6. Gaps/Missing Data Treatment: Pairwise deletion
7. Select Codon Positions: 1st, 2nd 3rd Codon Positions and Noncoding sites
8. Number of Threads: 4
9. Missing Base: N

Save settings as: `infer_NJ_nucleotide_boots100_jukescantor_uniformrates_pairwisedeletion_noncodingssites_4threads.mao`
in the MEGAX folder. Note: have to select all codon positions even if don't
consider as protein-coding.

Now prepare the MEGAX folders

```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/screening/megax
mkdir input output
for i in /projects1/microbiome_calculus/evolution/04-analysis/screening/megax/*/; do
  mkdir -p \
  "$i"/multiallelic_notallowed/strict \
  "$i"/multiallelic_notallowed/relaxed \
  "$i"/multiallelic_allowed/strict \
  "$i"/multiallelic_allowed/relaxed
done

for i in /projects1/microbiome_calculus/evolution/04-analysis/screening/megax/*put/*/*; do
  for j in /projects1/microbiome_calculus/evolution/04-analysis/screening/eager/output/strict/*/; do
    SPECIES="$(echo $j | rev | cut -d/ -f 2 | rev)"
    mkdir "$i"/$SPECIES
  done
done

## Now import the reduced fasta files for running
for i in $(find /projects1/microbiome_calculus/evolution/04-analysis/screening/multivcfanalyzer/output_unifiedgenotyper_concat2x/multiallelic_notallowed/ -name '*_pseudo2x_coverage_10percent.fasta.gz' -type f); do
  PARAM=$(echo "$i" | rev | cut -d/ -f 4 | rev)
  MAPPING=$(echo "$i" | rev | cut -d/ -f 3 | rev)
  SPECIES=$(echo "$i" | rev | cut -d/ -f 2 | rev)
  ln -s "$i" /projects1/microbiome_calculus/evolution/04-analysis/screening/megax/input/"$PARAM"/"$MAPPING"/"$SPECIES"/
done

```

Now submit MEGACC runs as follows.

```bash
for i in $(find -L /projects1/microbiome_calculus/evolution/04-analysis/screening/megax/input/ -name '*_pseudo2x_coverage_10percent.fasta.gz' -type f); do
  PARAM=$(echo "$i" | rev | cut -d/ -f 4 | rev)
  MAPPING=$(echo "$i" | rev | cut -d/ -f 3 | rev)
  SPECIES=$(echo "$i" | rev | cut -d/ -f 2 | rev)
  
  sbatch \
  -c 4 \
  --mem 32G \
  --partition=short \
  -t 02:00:00 \
  -o ~/slurm_logs/slurm.%j.out \
  -e ~/slurm_logs/slurm.%j.err \
  --mail-type=fail \
  --mail-user=fellows@shh.mpg.de \
  -J "MEGAX-CC_$i" \
  --wrap="gunzip -f $i && megacc \
  --analysisOptions /projects1/microbiome_calculus/evolution/04-analysis/screening/megax/infer_NJ_nucleotide_boots100_jukescantor_uniformrates_pairwisedeletion_allcodonposnoncodingssites_4threads.mao \
  --data ${i%.gz} \
  -format Fasta \
  --outfile /projects1/microbiome_calculus/evolution/04-analysis/screening/megax/output/$PARAM/$MAPPING/$SPECIES/ \
  && gzip -f $i"
done

```

Once completed, we can open the R notebook `25-Tree_visualisation.Rmd` to 
make prettified trees.

### Deep Sequencing Data (Competitive Mapping)

This section describes an attempt to assess if competitive mapping across
a genus can help reduce the number of multi-allelic SNPs than reduce the
reliability of phylogenetic analysis.

The competitive mapping procedure consists of:

1. Map to a superreference of all species of a genus.
2. Select a single genome of the most likely 'most abundant' species
3. Map to the single genome
4. Extract multi-allelic SNP information by MultiVCFAnalyzer
5. Compare the multi-allelic SNP rates of the single genome species vs the same 
genome but in the superreference

#### Superreference Setup

To begin, we want to assess what taxa we do actually have in our sample.

For this we took the core genera selected by prevalence (see above), then
ran the notebook `99-phylogeny_reference_genome_collection_YYYMMDD` which was
used to generate 'superreference' genomes of one representative assembly per 
species of a genus.

Note that Psuedopropionibacterium was done separately, due to recent renaming
procedures but closely related skin taxa. 

With this we then map all of our data to this reference genome with EAGER in
the form of 'competitive mapping' to see which taxa we more likely have (and
if we reduce the amount of cross mapping)

Make our directories

```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/deep_temp/eager/superreference_mapping

## Input
cd input
ln -s /projects1/microbiome_calculus/evolution/03-preprocessing/deep/library_merging/output/* .
cd ..

## References
cd references
for i in /projects1/microbiome_calculus/evolution/01-data/genomes/*; do 
  mkdir "$(basename $i)"
  ln -s "$i"/*.* "$(pwd)"/"$(basename $i)"; 
done
cd ..

## Output
cd output
for i in /projects1/microbiome_calculus/evolution/01-data/genomes/*; do 
  mkdir "$(basename $i)"
done
cd ..
```

Then set up EAGER as the following


```
Organism: Bacteria/Other
Age: Ancient
Treated Data: UDG
Pairment: Single End
Input is already concatenated (skip merging): Y

Reference: [listed above]
Name of mitocondrial chromosome: [none]

FastQC: Off
AdapterRemoval: Off
Mapping: BWA
  Seedlength: 32 
  Max # diff: 0.1
  Qualityfilter: 0
  Filter unmapped: On
  Extract Mapped/Unmapped: Off
Complexity Estimation: Off
Remove Duplicates: DeDup
  Treat all reads as merged: On
Damage Calculation: On
SNP Calling: On
  Emit All Sites: Y
CleanUp: On
Create Report: On
```

The EAGER runs are stored in `04-analysis/deep/eager/superreference_mapping`

Once complete, we can then clean up the non-useful intermediate files:

```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/superreference_mapping/output
find . -name '*only.sorted.bam' -delete
find . -name '*only.sorted.bam.bai' -delete
find . -name '*sorted.cleaned.bam' -delete

```

#### Superreference Coverage Statistics

We need to generate coverage statistics with `bedtools coverage`. We can run 
this with the script `035-bedtools_stats_array.sh`. This makes two output 
files: one with `-mean` for depth of coverage, and another without for number
of reads and breadth of coverage. 

#### Species Selection

To check how well this competitive mapping works, we also need a single
genome mapping for each genus to compare against. We used the stats generated
above to visualise a variety of mapping metrics (e.g. depth, breadth, no. reads, 
as performed in `031-superreferencemapped_genotyping_stats_*.Rmd` (the plots are in 
`04-analysis/deep/competitive_mapping/species_selection/`).

A list of taxa were collected based on visual inspection of plots and an automated
system which selected the most prevalent taxa across all individuals, and that 
also were the had the highest metric that was above one mean + standard 
deviation of the metric across all samples.

This list is as follows:

Genus                    | Majority Vote                                
-------------------------|-----------------------------------------------
Actinomyces              | Actinomyces_dentalis_DSM_19115 
Campylobacter            | Campylobacter_gracilis
Capnocytophaga           | Capnocytophaga_gingivalis_ATCC_33624
Corynebacterium          | Corynebacterium_matruchotii_ATCC_14266
Fretibacterium           | Fretibacterium_fastidiosum
Fusobacterium            | Fusobacterium_hwasookii_ChDC_F206
Olsenella                | Olsenella_sp_oral_taxon_807
Ottowia                  | Ottowia_sp_oral_taxon_894
Porphyromonas            | Porphyromonas_gingivalis_ATCC_33277
Prevotella               | Prevotella_loescheii_DSM_19665_=_JCM_12249_=_ATCC_15930
Pseudopropionibacterium  | Pseudopropionibacterium_propionicum_F0230a
Selenomonas              | Selenomonas_sp_F0473
Streptococcus            | Streptococcus_sanguinis_SK36
Tannerella               | Tannerella_forsythia_92A2
Treponema                | Treponema_socranskii_subsp_paredis_ATCC_35535

These specific genomes are then copied from the `01-data/genomes` folder into
`04-analysis/deep/eager/initial_single_genome/references`. Those with 
multiple chromosomes are collapsed as by the 
script described in `99-phylogeny_reference_genome_collection_YYYMMDD`, as also 
done for the superreference mapping genome generation performed above.

We thus generate indices of the references, run EAGER (and clean up) with the 
same settings above, but with the above single references instead of the 
super-references, and GATK UnifiedGenotyper on, with 2 ploidy and emit_all_sites.
The output of these runs are in `04-analysis/deep/eager/initial_single_genome`.

##### Superreference

Make output directories and prepare input files.

```bash
mkdir -p /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/superreference_mapping
cd !$
mkdir output input

for i in /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/superreference_mapping/output/*/; do 
  mkdir input/"$(basename $i)"
  for j in /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/superreference_mapping/output/"$(basename $i)"/*/; do
    mkdir input/"$(basename $i)"/"$(basename $j)"
    ln -s /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/superreference_mapping/output/"$(basename $i)"/"$(basename $j)"/10-GATKGenotyper/*.vcf.gz $(readlink -f input/"$(basename $i)"/"$(basename $j)")
  done
done

for i in /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/superreference_mapping/output/*/; do 
  mkdir -p output/superreference_mapping_2X_0.7/"$(basename $i)"; 
done
```

To run MultiVCFAnalyzer, I here had to do this manually as the superreferences
are so large and had to keep increasing the java heap size by changing the SLURM
`--mem` and the java `-Xmx` flag (e.g. `-Xmx 64G`.

To run MultiVCFAnayzer itself run the following sbatch script:
`036-crossmappingassessment_multivcfanalyzer_superreference_2X_slurm_array.sh`.

Once completed, we need to make a subset of the resulting tables and 
statistics to just the entry in the superreference that matches our selected
single genome for each genus

This subsetting can be performed using the script: 
`036-generate_MultiVCFAnalyzer_subset.R`


```bash
for i in /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/*; do 
  species="$(basename $i)"
  genus="$(echo $species | cut -d_ -f 1)"
  Rscript /projects1/microbiome_calculus/evolution/02-scripts.backup/036-generate_multiVCFanalyzer_subset_statistics_20190322.R \
  /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/superreference_mapping/output/"$genus"/ \
  /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/superreference_mapping/references/"$genus"/collapsed_"$genus"_superreference.bed \
  "$species"
done
```

##### Single Genome 

Next we want to generate an estimation of the number of multi-allelic-sites.
For this we can run MultiVCFAnalyzer which summarises all of this information
from the VCF files.

```bash
mkdir -p /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome
cd !$
mkdir output input

for i in /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/output/*/; do 
  mkdir input/"$(basename $i)"
  for j in /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/output/"$(basename $i)"/*/; do
    mkdir input/"$(basename $i)"/"$(basename $j)"
    ln -s /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/output/"$(basename $i)"/"$(basename $j)"/10-GATKGenotyper/*.vcf.gz $(readlink -f input/"$(basename $i)"/"$(basename $j)")
  done
done

for i in /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/output/*/; do 
  mkdir output/"$(basename $i)"; 
done
```

We can then use a bash array to find the corresponding references and then
submit each MultiVCFAnalyzer run.

```bash
unset astr
## Make bash array, storing the species name as key and value as reference genome
declare -A astr
for i in /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/references/*; do
  species="$(basename $i)"
  astr["$species"]="$(readlink -f /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/references/$species/*fna)"
done

## Submit MultiVCFAnalyzer runs by looping through array keys and values
for species in "${!astr[@]}"; do
  MULTIVCFANALYZER=/projects1/tools/multivcfanalyzer/0.0.87/bin/multivcfanalyzer
  FILES=($(find -L /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/input/$species/ -name '*.vcf.gz' -type f))
  sbatch \
  -c 4 \
  --mem 32G \
  --partition=short \
  -t 02:00:00 \
  -o ~/slurm_logs/slurm.%j.out \
  -e ~/slurm_logs/slurm.%j.err \
  --mail-type=fail \
  --mail-user=fellows@shh.mpg.de \
  -J "MultiVCFAnalyzer_single_$species" \
  --wrap="$MULTIVCFANALYZER \
  NA \
  ${astr[$species]} \
  NA \
  /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/$species \
  T \
  30 \
  3 \
  0.9 \
  0.1 \
  NA \
  ${FILES[*]}
  pigz -p 4 /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/$species/*
  gunzip /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/$species/snpStatistics.tsv.gz"
done



```

For those reference genomes with superreferences or `.fa` ending.

```bash
declare -A astr
for i in /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/references/*; do
  species="$(basename $i)"
  astr["$species"]="$(readlink -f /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/references/$species/*fa)"
done

## Submit MultiVCFAnalyzer runs by looping through array keys and values
for species in "${!astr[@]}"; do
  MULTIVCFANALYZER=/projects1/tools/multivcfanalyzer/0.0.87/bin/multivcfanalyzer
  FILES=($(find -L /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/input/$species/ -name '*.vcf.gz' -type f))
  sbatch \
  -c 4 \
  --mem 32G \
  --partition=short \
  -t 02:00:00 \
  -o ~/slurm_logs/slurm.%j.out \
  -e ~/slurm_logs/slurm.%j.err \
  --mail-type=fail \
  --mail-user=fellows@shh.mpg.de \
  -J "MultiVCFAnalyzer_single_$species" \
  --wrap="$MULTIVCFANALYZER \
  NA \
  ${astr[$species]} \
  NA \
  /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/$species \
  T \
  30 \
  3 \
  0.9 \
  0.1 \
  NA \
  ${FILES[*]}
  pigz -p 4 /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/$species/*
  gunzip /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/$species/snpStatistics.tsv.gz"
done

```

#### Amylase Binding Protein Statistics

For the Amylase binding protein A/B stuff, we also want to look at the ratio of
streptococcus coverage vs abpA/abP specific coverage. Using the mapping to
_Streptococcus gordonii_, we extract the coverage of the following genes

  * abpA - NC_009785.1     Protein Homology        CDS     2,167,409 2167996 .       -       0       ID=cds2039;Parent=gene2087;Dbxref=Genbank:WP_008810020.1,GeneID:25052597;Name=WP_008810020.1;gbkey=CDS;inference=COORDINATES: similar to AA sequence:RefSeq:WP_008810020.1;product=hypothetical protein;protein_id=WP_008810020.1;transl_table=11
  * abpB -  NC_009785.1     Protein Homology        CDS     171195  173153  .       +       0       ID=cds160;Parent=gene161;Dbxref=Genbank:WP_011999704.1,GeneID:25051779;Name=WP_011999704.1;gbkey=CDS;inference=COORDINATES: similar to AA sequence:RefSeq:WP_006597304.1;product=dipeptidase;protein_id=WP_011999704.1;transl_table=11

We can extract the coverages with the following command

```bash

cd /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/output/Streptococcus_gordonii_str_Challis_substr_CH1

for BAM in */5-DeDup/*_rmdup.sorted.bam ; do
  bedtools coverage -b "$BAM" -a /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/references/Streptococcus_gordonii_str_Challis_substr_CH1/GCF_000017005.1_ASM1700v1_genomic.gff -mean | pigz -p 4 > "$BAM"_featurecoverage.txt.gz;
done

zgrep -e "WP_008810020.1" -e "WP_011999704.1" */5-DeDup/*featurecoverage.txt.gz > amylasebindingprotein_averagedepthcoverages.txt
```

And for read counts

```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/output/Streptococcus_gordonii_str_Challis_substr_CH1

for BAM in */5-DeDup/*_rmdup.sorted.bam ; do
  bedtools coverage -b "$BAM" -a /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/initial_single_genome/references/Streptococcus_gordonii_str_Challis_substr_CH1/GCF_000017005.1_ASM1700v1_genomic.gff | pigz -p 4 > "$BAM"_featurereadcount.txt.gz;
done

zgrep -e "WP_008810020.1" -e "WP_011999704.1" */5-DeDup/*_featurereadcount.txt.gz > amylasebindingprotein_averagereadcounts.txt

```

Then we run the R notebook `031-superreference_genotyping_statistics.Rmd` for analysis
of these files.

#### Multi Allellic SNP Assessment

The `032-multibasesites_mappingstrategycomparison.Rmd` notebook describes
the assessment of whether the superreference mapping reduces the amount
of multi-allelic SNPS. 

While in most cases this does occur, the number of SNPS left over is reduced
suggesting that while the quality increases of SNP calling, the resulting
data will result in less data that will lead to lower resolution phylogenetic 
analysis.

Instead, we will focus on the single genome mapping. To select which 
SNP calling frequency threshold, the notebook also assessment on what frequency
cut off to use - the most common of which appears to be 0.7.


#### Final SNP Alignments for Phylogenies

We clearly see that we do not have clean alignments, and likely will be
unable to get these 'clean' even with the competitive mapping approach.
Regardless, we will test both single genome and competitive mapping alignments 
but focusing on a single genome.

For this, we will need to run both the single genome and superreference 
MultiVCFAnalyzer runs three times at each threshold level. However we will
turn off the SNP frequencies, and set both 'heterogyzous' and 'homozygous' 
thresholds to 0.7, (i.e. not allowing heterozygous calls to be made).

The output of these runs are stored in `/projects1/microbiome_calculus/evolution/04-analysis/deep/snp_alignment`.
The runs use the same input files as above, however with tweaked MultiVCFAnalyzer
settings. These settings can be seen in the scripts:

```bash
sbatch /projects1/microbiome_calculus/evolution/02-scripts.backup/037-snpalignment_multivcfanalyzer_superreference_2X_slurm_array.sh
sbatch /projects1/microbiome_calculus/evolution/02-scripts.backup/037-snpalignment_multivcfanalyzer_initialsinglegenome_2X_slurm_array.sh

```

> We again retain 2X as a depth coverage, due to the low coverage across most
> of the species.

For the single genome mappings, we can immediately proceed with the output
'snpAlignment.fasta' file for phylogenies (see next section).

For the superreference mapping, we firstly need to extract the same genomes
as in the single genome mappings from the superreference mapping. This we
can use with the same 'subset' script as above: `036-generate_MultiVCFAnalyzer_subset.R`

For the 2X coverage parameters

```bash
for i in /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/initial_single_genome_2X_0.7_0.7/*; do 
  species="$(basename $i)"
  genus="$(echo $species | cut -f 1 -d_)"

  echo "$species" "$genus"
  
  Rscript /projects1/microbiome_calculus/evolution/02-scripts.backup/036-generate_multiVCFanalyzer_subset_statistics_20190404.R \
  /projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/superreference_mapping/output/superreference_mapping_2X_0.7_0.7/"$genus"/ \
  /projects1/microbiome_calculus/evolution/04-analysis/deep/eager/superreference_mapping/references/"$genus"/collapsed_"$genus"_superreference.bed \
  "$species"
done
```

For prevotella due to stupid strain name formatting

```bash
genus="Prevotella"
species='Prevotella_loescheii_DSM_19665_=_JCM_12249_=_ATCC_15930'

Rscript /projects1/microbiome_calculus/evolution/02-scripts.backup/036-generate_multiVCFanalyzer_subset_statistics_20190404.R \
/projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/superreference_mapping/output/superreference_mapping_2X_0.7_0.7/"$genus"/ \
/projects1/microbiome_calculus/evolution/04-analysis/deep/eager/superreference_mapping/references/"$genus"/collapsed_"$genus"_superreference.bed \
"$species"
```

#### Deep Phylogenies

I later experimented more with fold coverage SNP calling thresholds (see above)
and found that in some cases this provided enough data that the MEGAX 
phylogenies did not crash as much during bootstrapping. I also think this
was the cause of bootstrapping failing in R with ape/adegenet.

Also, after further exploration, I identified that I could do filtering
_within_ R of samples with low numbers of SNPs. Therefore, I wrote a script
which does both sample filtering (i.e. removing samples with less than X
number of SNPs), then does bootstrapping of pairwise neighbour-joining trees. 
The script produces a newick file we can use for visualisaion.

We can run this script  to generate the tree files.

```bash

PROJDIR=/home/fellows/projects1/microbiome_calculus/evolution

## Initial Single Genomes
for i in "$PROJDIR"/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/initial_single_genome_2X_0.7_0.7/*/snpAlignment.fasta.gz; do
  echo "$i"
  Rscript "$PROJDIR"/02-scripts.backup/042-generate_NJ_tree.R "$i" 1000 JC69 100 none
  echo ""
done

Rscript "$PROJDIR"/02-scripts.backup/042-generate_NJ_tree.R "$PROJDIR"/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/initial_single_genome_2X_0.7_0.7/Porphyromonas_gingivalis_ATCC_33277/snpAlignment.fasta.gz 2000 JC69 100 OME003

## Superreference 
for i in "$PROJDIR"/04-analysis/deep/multivcfanalyzer/superreference_mapping/output/superreference_mapping_2X_0.7_0.7/*/snpAlignment_subset.fasta.gz; do
  echo "$i"
  Rscript "$PROJDIR"/02-scripts.backup/042-generate_NJ_tree.R "$i" 1000 JC69 100 none
  echo ""
done


Rscript "$PROJDIR"/02-scripts.backup/042-generate_NJ_tree.R "$PROJDIR"/04-analysis/deep/multivcfanalyzer/superreference_mapping/output/superreference_mapping_2X_0.7_0.7/Tannerella/snpAlignment_subset.fasta.gz 2000 JC69 100 none
Rscript "$PROJDIR"/02-scripts.backup/042-generate_NJ_tree.R "$PROJDIR"/04-analysis/deep/multivcfanalyzer/superreference_mapping/output/superreference_mapping_2X_0.7_0.7/Streptococcus/snpAlignment_subset.fasta.gz 2500 JC69 100 none

Rscript "$PROJDIR"/02-scripts.backup/042-generate_NJ_tree.R "$PROJDIR"/04-analysis/deep/multivcfanalyzer/superreference_mapping/output/superreference_mapping_2X_0.7_0.7/Selenomonas/snpAlignment_subset.fasta.gz 1000 JC69 100 PES001,OME003,TAF008


```

* Initial genome: Porphyromonas_gingivalis_ATCC_33277 at bootstrapping 
  regardless of min SNP threshold
* Supperreferference: Selenomas fails at initial tree regardless of min SNP 
 threshold

_Porphyromonas gingivalis_ failed during bootstrapping because OME003 has 
an extremely high variability from the reference across the entire alignment. 
During some steps of the resampling for the alignment during bootstrapping, 
this variability would lead to the randomly sampled alignment to have a 
very high difference from the frequencies compared to the other samples (i.e.
most other samples having an 'A' basefrequency of 0.25, and for OME003 it would
be 0.19), violating an assumption of JC69 of equal base frequencies. This was 
noticed as re-running bootstrapping would fail at different bootstrapping 
interations. This can possibly be demonstrated by looking at the ratio of
differences/shared SNPS between each sample combination (see the `*_overlappingNucleaotidesReport.csv` in
the `04-analysis/deep/multivcfanalyzer/superreference_mapping/output/superreference_mapping_2X_0.7_0.7` folders), 
where OME003 dominates the highest ratio of different/shared SNPS having almost 
up to three times more differences than SNPs.

The superreference Selenomonas alignment also suffers the same problem but
to a greater extent (i.e. 5-6 times more differences than shared), which 
is why this fails with the initial distance matrix calculation.

In both cases this suggests that the present species that has mapped to
the reference is a very highly diverged relative, resulting in an inordinate 
number of SNPs, and those samples should be excluded. For Selenomonas, this
occurs across all samples so likely the diveristy of that genus is currently
too undersampled in terms of reference genomes to allow further investigation. 

##### El Mirón Observation Check

We observed that El Mirón always clustered with Neanderthals in all the deep
phylogenies. We wanted to check wether this is specific to upper palaeolithic
Eurpoean individuals, potentially indicating a proxy of population history
patterns.

We decided to try confirm this single by creating phylogenies including
another well-preserved individual from Upper Palaeolithic Europe (in particular
pre-LGM), which was not deep-sequenced. The only individual matching this
criteria is PLV001 and RIG001, based on RIG human DNA falling in the El Mirón 
cluster in Fu et al 2016, and an Dolni Vestonice individua lfrom 31ka in Posth
et al, and PLV falling in the earlier Vestonice cluster in Fu et al.

Therefore we will try to make phylogenies based on the screening data, with 
the caveats that low-coverage and damage may cause problems.

The analysis for this is here: 

```
04-analysis/screening/EMN_Neanderthal_phylogeny_check
```

I set up EAGER runs for the screening versions of the samples used in
deep sequencing, plus RIG001 and PLV001, with the same settings for the 
production phylogenies (reminder below), to the references of the better 
support phylogenies.

```
Organism: Bacteria/Other
Age: Ancient
Treated Data: UDG
Pairment: Single End
Input is already concatenated (skip merging): Y

Reference: [listed above]
Name of mitocondrial chromosome: [none]

FastQC: Off
AdapterRemoval: Off
Mapping: BWA
  Seedlength: 32 
  Max # diff: 0.1
  Qualityfilter: 0
  Filter unmapped: On
  Extract Mapped/Unmapped: Off
Complexity Estimation: Off
Remove Duplicates: DeDup
  Treat all reads as merged: On
Damage Calculation: On
SNP Calling: On
  Emit All Sites: Y
CleanUp: On
Create Report: On
```


Once completed, coverage statistics were assessed by concatenation of EAGER
tables as in `056-Phylogenies_Screening_EMNCheck_EAGERResults.Rmd`. Coverages
remain very low for most ancient samples, so we won't expect much from 
downstream analysis...

Next we run MultiVCFAnalyzer on the VCF files. To prepare

```bash
mkdir -p /projects1/microbiome_calculus/evolution/04-analysis/screening/EMN_Neanderthal_phylogeny_check/multivcfanalyzer/
cd !$
mkdir output input

for i in /projects1/microbiome_calculus/evolution/04-analysis/screening/EMN_Neanderthal_phylogeny_check/eager/output/*/; do 
  mkdir input/"$(basename $i)"
  for j in /projects1/microbiome_calculus/evolution/04-analysis/screening/EMN_Neanderthal_phylogeny_check/eager/output/"$(basename $i)"/*/; do
    mkdir input/"$(basename $i)"/"$(basename $j)"
    ln -s /projects1/microbiome_calculus/evolution/04-analysis/screening/EMN_Neanderthal_phylogeny_check/eager/output/"$(basename $i)"/"$(basename $j)"/10-GATKGenotyper/*.vcf.gz $(readlink -f input/"$(basename $i)"/"$(basename $j)")
  done
done

for i in /projects1/microbiome_calculus/evolution/04-analysis/screening/EMN_Neanderthal_phylogeny_check/eager/output/*/; do 
  mkdir output/"$(basename $i)"; 
done
```

Then run the SBATCH script `057-snpalignment_multivcfanalyzer_singlegenome_2X_screening_slurm_array.sh`

We can attempt to make some phylogenies with

```bash
PROJDIR=/home/fellows/projects1/microbiome_calculus/evolution

## Initial Single Genomes
for i in "$PROJDIR"/04-analysis/screening/EMN_Neanderthal_phylogeny_check/multivcfanalyzer/output/*/snpAlignment.fasta.gz; do
  echo "$i"
  Rscript "$PROJDIR"/02-scripts.backup/042-generate_NJ_tree.R "$i" 1000 JC69 100 none
  echo ""
done
```

With these we can load into the R markdown document `058-Tree_visualisation_20191025_screening.Rmd` to see the trees.

## Functional Analysis

### HuMANN2

#### Installation
To create conda environment and install software

```bash
conda create --name humann2 -c bioconda humann2
```

To load

```bash
conda activate humann2
```

From now on following dependency checks as described in:
http://huttenhower.sph.harvard.edu/humann2
https://bitbucket.org/biobakery/humann2/wiki/Home#markdown-header-configuration

Small bug is the MetaPhlAn2 database doesn't come with the bioconda humann2 install.

We can link that into the directory it was looking for with the following from our local install

```bash
/projects1/users/fellows/bin.backup/miniconda3/envs/humann2/bin
ln -s /projects1/tools/metaphlan2/biobakery-metaphlan2-27f7e0c86785/db_v20/ .
```

ChocoPhlan I had already installed and it somehow picked it up (wtf?!), but
uniref wasn't installed.

```bash
/projects1/microbiome_calculus/evolution/01-data/databases/uniref90
humann2_databases --download uniref uniref90_ec_filtered_diamond "$(pwd)"
```

the humann2 utility scripts came with the conda environment.

#### Running

Load conda environment

```bash
conda activate humann2
```

Then submit the batch array `02-scripts.backup/030-humann2_slurm_array.sh`

Note, the output temporary files (which don't appear to be removed) are HUGE.
We need to remove them after successful running by doing the following:

```bash
cd /projects1/microbiome_calculus/evolution/04-analysis/screening/humann2/output
rm -r */*_temp/
```

which reduces our footprint from 4.9 TERABYTES(!?!?!) down to 3.9 gigabytes.
This is mostly due to the uncompressed SAM files.

Next we need to normalise the data of each file. I'm not familiar which is 
better, so for the moment we will go for the default. We only need
to do this on the genefamilies and pathabundance data, as you don't have to
do this for the [pathwaycoverage](https://bitbucket.org/biobakery/biobakery/wiki/humann2#rst-header-manipulating-humann2-output-tables)

```bash

conda activate humann2

cd /projects1/microbiome_calculus/evolution/04-analysis/screening/humann2/output

for SAMPLE in */*_genefamilies.tsv; do
    humann2_renorm_table --input $SAMPLE --output ${SAMPLE%.tsv}_cpm.tsv --units cpm --update-snames
done

for SAMPLE in */*_pathabundance.tsv; do
    humann2_renorm_table --input $SAMPLE --output ${SAMPLE%.tsv}_cpm.tsv --units cpm --update-snames
done

rm */*_genefamilies.tsv
rm */*_pathabundance.tsv

# -s to search subdirectories of current directory to look for the files
humann2_join_tables --input . --output humann2_genefamilies.tsv --file_name genefamilies_cpm -s 
humann2_join_tables --input . --output humann2_pathcoverage.tsv --file_name pathcoverage -s 
humann2_join_tables --input . --output humann2_pathabundance.tsv --file_name pathabundance_cpm -s

```

Irina would like the output in KEGG format. To do this we firstly need to 
download the re-grouping database

```bash
humann2_databases --download utility_mapping full /projects1/microbiome_calculus/evolution/01-data/databases/humann2
humann2_regroup_table -i humann2_pathabundance.tsv -g uniref90_ko -o humann2_pathabundance_ko.tsv
```

Then we can regroup our table

```bash
humann2_regroup_table -i humann2_genefamilies.tsv -g uniref90_ko -o humann2_genefamilies_ko.tsv
```

Next need to add metadata to the file...

### AADDER Analysis

We want to experiment with the unpublished AADDER tool that comes with
MEGAN6. This uses `.gff` files to compare taxonomic assignments with
annotations. 

We will run this against the filtered RefSeq genomes we built above both
with MALT (fasta files) and aadder (gff files).

We then run MALT but instead of immediately producing RMA6 files, we 
generate SAM files. This can be run with the following:

```bash
sbatch \
-c 112 \
--mem 1900G \
--partition=supercruncher \
-o /projects1/microbiome_calculus/evolution/04-analysis/screening/malt/refseq_bacarchhomo_gcs_20181122/slurm.%j.out \
-e /projects1/microbiome_calculus/evolution/04-analysis/screening/malt/refseq_bacarchhomo_gcs_20181122/slurm.%j.err \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fellows@shh.mpg.de \
-J "MALT040-RefSeq_GCS_2018" \
--wrap="/projects1/microbiome_calculus/evolution/02-scripts.backup/008-malt-genbank-refseq_bacarchhomo_gcs_20181122_4step_85pc_supp_0.01 /projects1/microbiome_calculus/evolution/04-analysis/screening/malt/temp_input/*/*.fq.gz /projects1/microbiome_calculus/evolution/04-analysis/screening/malt/refseq_bacarchhomo_gcs_20181122/"
```

To extract the summary statistics for the Refseq runs, we can run
the following on the MALT logs.

```bash
grep -e "Loading MEGAN File:" \
-e "Total reads:" \
-e "With hits:" \
-e "Alignments:" \
-e "Assig. Taxonomy" \
-e "Min-supp. changes" \
-e "Numb. Tax. classes:" \
-e "Class. Taxonomy:" \
-e "Num. of queries:" \
-e "Aligned queries:" \
-e "Num. alignments:" \
-e "MinSupport set to:" \
/projects1/microbiome_calculus/evolution/04-analysis/screening/malt/refseq_bacarchhomo_gcs_20181122/*log | cut -d":" -f 2-99 > /projects1/microbiome_calculus/evolution/00-documentation.backup/99-maltAlignedReadsSummary_raw_refseq_bacarchhomo_gcs_20181122_$(date "+%Y%m%d").txt
```


Then we run AADDER

```bash

sbatch \
-c 112 \
--mem 1900G \
--partition=supercruncher \
-o /projects1/microbiome_calculus/evolution/04-analysis/screening/aadder/output/slurm.%j.out \
-e /projects1/microbiome_calculus/evolution/04-analysis/screening/aadder/output/slurm.%j.err \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fellows@shh.mpg.de \
-J "aadder-run" \
--wrap="/projects1/users/fellows/bin.backup/megan6/tools/aadder-run \
-i /projects1/microbiome_calculus/evolution/04-analysis/screening/aadder/input/*.sam.gz \
-d /projects1/microbiome_calculus/evolution/01-data/databases/aadder/ \
-o /projects1/microbiome_calculus/evolution/04-analysis/screening/aadder/output/ \
-v
pigz -p 112 /projects1/microbiome_calculus/evolution/04-analysis/screening/aadder/output/*"
```

Finally run blast2rma

```bash
sbatch \
-c 112 \
--mem 1900G \
--partition=supercruncher \
-o /projects1/microbiome_calculus/evolution/04-analysis/screening/aadder/output/slurm.%j.out \
-e /projects1/microbiome_calculus/evolution/04-analysis/screening/aadder/output/slurm.%j.err \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fellows@shh.mpg.de \
-J "blast2rma" \
--wrap="/projects1/users/fellows/bin.backup/megan6/tools/blast2rma \
--format SAM \
-i /projects1/microbiome_calculus/evolution/04-analysis/screening/aadder/output/*out.gz \
-o /projects1/microbiome_calculus/evolution/04-analysis/screening/aadder/rma6_2/ \
-a2seed /projects1/malt/databases/acc2seed/acc2seed-May2015XX.abin \
-a2t /projects1/malt/databases/acc2tax/nucl_acc2tax-Nov2018.abin \
-v"
```

Once completed, the resulting RMA6 files were opened in MEGAN6CE (v6.15.2) via 
MEGAN SERVER, opened in compare mode using absolute counts and 
'ignore all unassigned reads'. I then switched to the 'SEED' categories, 
'uncollapse-all', select everythign then exported the table as TSV with 
seedPath_to_count.

# Compositional Heatmaps

## Generation

For these we can run the script that performs CLR transform of the OTU matrix
and unsupervised clustering of host taxa and microbial taxa. Clustering
algorithm is selected automatically within the script.

The script allows filtering as above (database, taxonomic levels, with/without
sources, with/without controls, with/without bad samples (+ bad sample removal
method option)) and also additional taxon filtering, zero-replacement method,
additional min-support filtering and a prevalence filter (i.e. a taxon is
only kept if it is in _n_ number of indiviudals across the dataset)

We look for the parameters with the best overall bootstrap support 
in the deepest nodes (i.e. the ones we are most interested in - splits between
host genus level clades)

```bash

## No additional min support filtering, no genus filtering
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd nt species noSources noControls out withinvariation none pseudo 0 0
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd nt species noSources noControls out withinvariation none pseudo 0 5

## No additional min support filtering,  genus filtering to core only
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd nt species noSources noControls out withinvariation "Actinomyces|Campylobacter|Capnocytophaga|Corynebacterium|Desulfomicrobium|Fusobacterium|Fretibacterium|Mogibacterium|Mycobacterium|Olsenella|Ottowia|Parvimonas|Prevotella|Porphyromonas|Pseudopropionibacterium|Selenomonas|Streptococcus|Treponema|Tannerella" pseudo 0 0
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd nt species noSources noControls out withinvariation "Actinomyces|Campylobacter|Capnocytophaga|Corynebacterium|Desulfomicrobium|Fusobacterium|Fretibacterium|Mogibacterium|Mycobacterium|Olsenella|Ottowia|Parvimonas|Prevotella|Porphyromonas|Pseudopropionibacterium|Selenomonas|Streptococcus|Treponema|Tannerella" pseudo 0 5

## Additional min support filtering, no genus filtering
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd nt species noSources noControls out withinvariation none pseudo 4 0
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd nt species noSources noControls out withinvariation none pseudo 4 5

## Additional min support filtering, genus filtering to core only
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd nt species noSources noControls out withinvariation "Actinomyces|Campylobacter|Capnocytophaga|Corynebacterium|Desulfomicrobium|Fusobacterium|Fretibacterium|Mogibacterium|Mycobacterium|Olsenella|Ottowia|Parvimonas|Prevotella|Porphyromonas|Pseudopropionibacterium|Selenomonas|Streptococcus|Treponema|Tannerella" pseudo 4 0
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd nt species noSources noControls out withinvariation "Actinomyces|Campylobacter|Capnocytophaga|Corynebacterium|Desulfomicrobium|Fusobacterium|Fretibacterium|Mogibacterium|Mycobacterium|Olsenella|Ottowia|Parvimonas|Prevotella|Porphyromonas|Pseudopropionibacterium|Selenomonas|Streptococcus|Treponema|Tannerella" pseudo 4 5
```

Output PDFs were then visually checked for the bootstrap values, and the 
bootstraps of the first 6 major (approximate) bifurcations were averaged. 
Major bifurcations were those that described large clade splits (i.e. not
for single 'orphan' tips), and bootstraps on the branch leading to a major
clade.

A table of resuts of the visual checking can be seen here:

| Database | Taxon Filtering | Min Support | Sample Prevalence | Approximate Deep Bootstrap Support | Phylogeny_Follows_Host?  | Notes                                                                    |
|----------|-----------------|-------------|-------------------|------------------------------------|--------------------------|--------------------------------------------------------------------------|
| nt       | CoreGenusOnly   | 0.04        | 0                 | 75.83                              | F                        |                                                                          |
| nt       | CoreGenusOnly   | 0.04        | 5                 | 82.50                              | F                        |                                                                          |
| nt       | CoreGenusOnly   | 0           | 0                 | 78.67                              | Close                    |                                                                          |
| nt       | CoreGenusOnly   | 0           | 5                 | 79.33                              | V. close                 |                                                                          |
| nt       | none            | 0.04        | 0                 | 88.17                              | V. close                 |                                                                          |
| nt       | none            | 0.04        | 5                 | 82.67                              | T                        | FAVOURITE: less messy (less single taxa; clean clustering, good support) |
| nt       | none            | 0           | 0                 | 71.17                              | F                        |                                                                          |
| nt       | none            | 0           | 5                 | 81.33                              | F                        |                                                                          |

The sample clustering with no taxon filtering, additional minsupport of 0.04, 
and prevalence filtering set to 5 was selected. This was based on it having
generally the highest bootstrap values in the internal nodes, and the phylogeny
showing 'cleanest' clustering of individuals of the same host genus falling 
together.

## Zero replacement validation

Next we can check whether this result is affected by the zero replacement model,
and generate with the same settings heatmaps from the alternative databases
as checks for the pattern in the nt heatmap.

```bash
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd nt species noSources noControls out withinvariation none czm 4 5
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd refseq species noSources noControls out withinvariation none pseudo 4 5
Rscript /home/fellows/projects1/microbiome_calculus/evolution/02-scripts.backup/045-Compositional_Heatmaps_20190806_script.Rmd mp2 species noSources noControls out withinvariation none pseudo 4 5
```

Additional observations below.

| Database | Taxon Filtering | Min Support | Sample Prevalence | Approximate Deep Bootstrap Support | Phylogeny_Follows_Host?                                                                                           | Notes                                                       |   |
| refseq   | none            | 0.04        | 5                 | 88.83                              | Close(ish) (modern day humans; then Gorilla/Alouattas; the ancient human (a); then chimps from ancient human (b)) |                                                             |   |
| mp2      | none            | 0           | 5                 | 49.17                              | Close(ish) (Gorilla/alouatta split first; but chimps in single clade within humans)                               | Tail filtering not applicable                               |   |

Comparing the zero replacement methods shows no difference between clustering.
There is cosmetic tree topology changes but  only by clade flipping - no 
structural chananges.

## Taxon block selection

I manually took the 0.04 min support and 5 individual prevalence threshold
heatmaps of the nt, refseq and MetaPhlAn2 databases, and visually defined
blocks of interest that defined either single host genus or host genus 
combination clades.

### Lac2 operon

This section is to look for and reconstruct the lac2 operon from various 
Streptococcus species. This operon is used for the catabolism of lactose

Search steps (detailed below):

1. make lac2 operon BLAST database
2. run all cleaned and merged evolution project screening files against the 
lac2 operon with blastn
3. list out the reads that had hits in the BLAST search
4. pull out the reads that had hits from the fastq files as fasta files
5. align the reads against each of the individual operons
6. check the reads that aligned with bwa for appropriate damage profiles
7. profile the cleaned & merged evolution project screening files with Kaiju to 
get read assignents and estimated proportion of Strep in each sample
8. read stats from Kaiju and bwa
9. assemble the reads that aligned to each operon


#### lac2 operon BLAST database Generation

Firstly, within an interactive SLURM session, we make the BLAST database out of 
a fasta file that contains the full lac2 operon sequence from 21 different 
Strep species, plus the reversed sequence of 3 of those species 
(24 sequences total). This fasta file is stored here 
`04-analysis/screening/lac2operon/lac2blastdb`

```bash
srun -c 1 --mem 4000 -t 0-00:30 --pty -u bash -i
source ~/.bash_profile

cd /projecst1/microbiome_calculus/evolution/04-analysis/screening/lac2operon/lac2blastdb

makeblastdb -in all.fasta -dbtype nucl -title lac2strezp -out lac2strep
```


#### BLAST Reads Against Database 

Next we convert the `.fq.gz` files into fasta so BLAST can read them. 

To do this we make a file called samples.paths.list

```bash
awk '{print "/projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/"$1"/"$1"_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz"}' samples.list > samples.paths.list
```

Then run the script `101-qtoa.sh` to convert the FASTQs to FASTA using `seqtk`,
and then gzip the FASTAS with `pigz` as in `103-zipit.sh`.

We then run the script script `102-lac2evolblast.sh` to search the FASTA files 
against the lac2strep database with BLASTn, and again gzip with `103-zipit.sh`. 
The search script streams as standard input the gziped FASTA into BLAST, which
aligns to the lac2operon database, and exports the top 20 alignments.

#### Extract BLAST lac2operon Positive Reads 

Script `104-ncbiparser.pl` was downloaded from 
`https://github.com/mel-astar/mel-ngs/blob/master/utils/ncbiblastparser.pl` 
and the perl path was updated for our cluster. 

To run the read extraction using the script above, the script `105-gethits.sh` 
was used and submitted with.

```bash
sbatch ./105-gethits.sh
```

The parser script gives a lot of additional statistics regarding the BLAST hit 
itself. Instead, we just want to know the actual reads themselves. This we 
can get by taking the first column of each file and storing this in a new file, 
as here `106-hitlist.sh`

#### Map BLAST lac2operon Positive Reads

Next we need to find these reads in the original FASTQ file, which we can then
use for mapping with `bwa`. This we do by running `seqtk` on  each file, and 
using the `grep` function with the list of files given above. We can then store
these in new fastq files. This procedure is in performed in `107-lac2fastq.sh`.

With these reads, we can now align all the reads to the each of the original 
lac2operon reference file(s) separately using `bwa`, as in `108-alignlac2.sh`. 
When aligning, we use the same parameters as for human, 16s and genome mapping 
above to account for damage.

We can then convert the resulting samfiles to bam with samtools. We run this in 
an interactive SLURM session due to the 1000 job array limit, and we have over 
4000 sam files. This command is stored in `110-sam2bam.sh`

Again in an SLURM sesson, we remove PCR duplicate reads from bam files and get a 
list of the mapped reads with `samtools`. The command is also stored in 
`111-rmdupgetmapped.sh`

Alignment stats are generated from the mappings with the commands stored in 
`112-imv-alignmentstats.sh` and `113-imv-uniqrmdupcounts.sh`. These commands 
give you the following information: Number of overall reads aligned (i.e. 
including multiple alignments), the number unique reads of reads to each 
operon, and the number of unique non-duplicated reads to each operon.

To calculate alignment coverage and depth of coverage on each operon, we first 
index lac2 fasta files in the file path below

```bash
for f in /projects1/microbiome_calculus/evolution/04-analysis/screening/lac2operon/bwaindex/*.fasta; do
sbatch -c 2 --mem=4000 --wrap="samtools faidx $f"
done
```

Then we organise the bam files in each folder replacing the host source e.g.
`B-Gorilla` to each source.

```bash
for f in /projects1/microbiome_calculus/evolution/04-analysis/screening/lac2operon/bwaout/B-Gorilla/*.bam; do 
sbatch -c 2 --mem=4000 --partition=short --wrap="samtools sort -o `basename $f .bam`.sort.bam $f"
done
```

We then convert lac2 reference fasta files to bed format.

```bash
for f in /projects1/microbiome_calculus/evolution/04-analysis/screening/lac2operon/bwaindex/*.fai; do
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $f > `basename $f .fai`.bed
done
```

We can then calculate depth and breadth of coverage of each sample to each 
lac2operon reference with bedtools to get the # of bases with each fold 
coverage level using `114-imv-getcoverage.sh`

To summarise the lowest and highest coverage for each, we can run 
`115-imv-highlow.sh` and move the results to 
`/projects1/microbiome_calculus/evolution/04-analysis/screening/lac2operon/bwaout/**XX**/`

We submit generate the summaries with 

```bash
sbatch -c 2 --mem=4000 ./115-imv-highlow.sh
```

NOTE FOR IRINA: THE SCRIPT IN THE FOLDER DOESN'T HAVE THE RENAMING COMMANDS 
COULD YOU UPDATE THIS AND THEN REMOVE THE CODE BLOCK BELOW.

```bash
#!/bin/bash

for f in *.genomecov
do
grep genome $f | sed -n '1p;$p' | awk '{print $2}' > `basename $f .genomecov`.highlow
done

paste *ALQD_Sagallac2.aln.sort.highlow > all.ALQD.aln.sort.highlow.txt
paste *ALRK_Sagallac2.aln.sort.highlow > all.ALRK.aln.sort.highlow.txt
paste *ALSA_Sagallac2.aln.sort.highlow > all.ALSA.aln.sort.highlow.txt
paste *ALSH_Sagallac2.aln.sort.highlow > all.ALSH.aln.sort.highlow.txt
paste *ANDK_Sagallac2.aln.sort.highlow > all.ANDK.aln.sort.highlow.txt
paste *ANES_Sagallac2.aln.sort.highlow > all.ANES.aln.sort.highlow.txt
paste *ANET_Sagallac2.aln.sort.highlow > all.ANET.aln.sort.highlow.txt
paste *ANEU_Sagallac2.aln.sort.highlow > all.ANEU.aln.sort.highlow.txt
paste *dysgalactiaelac2.aln.sort.highlow > all.dysgalactiae.aln.sort.highlow.txt
paste *gordoniilac2.aln.sort.highlow > all.gordonii.aln.sort.highlow.txt
paste *infantariuslac2.aln.sort.highlow > all.infantarius.aln.sort.highlow.txt
paste *mutanslac2.aln.sort.highlow > all.mutans.aln.sort.highlow.txt
paste *oralislac2.aln.sort.highlow > all.oralis.aln.sort.highlow.txt
paste *parasanguinislac2.aln.sort.highlow > all.parasanguinis.aln.sort.highlow.txt
paste *pneumoniaelac2.aln.sort.highlow > all.pneumoniae.aln.sort.highlow.txt
paste *pyogeneslac2.aln.sort.highlow > all.pyogenes.aln.sort.highlow.txt
paste *.sanguinislac2.aln.sort.highlow > all.sanguinis.aln.sort.highlow.txt
paste *Sdd_Sdysgallac2.aln.sort.highlow > all.Sdd_Sdysgal.aln.sort.highlow.txt
paste *Sdd_Sdysgal-revlac2.aln.sort.highlow > all.Sdd_Sdysgal-rev.aln.sort.highlow.txt
paste *staphaureus-revlac2.aln.sort.highlow > all.staphaureus.aln.sort.highlow.txt
paste *suislac2.aln.sort.highlow > all.suis.aln.sort.highlow.txt
paste *uberislac2.aln.sort.highlow > all.uberis.aln.sort.highlow.txt
paste *uberis-revlac2.aln.sort.highlow > all.uberis-rev.aln.sort.highlow.txt
paste *ALQD_Sagallac2.aln_rmdup.sort.highlow > all.ALQD.aln_rmdup.sort.highlow.txt
paste *ALRK_Sagallac2.aln_rmdup.sort.highlow > all.ALRK.aln_rmdup.sort.highlow.txt
paste *ALSA_Sagallac2.aln_rmdup.sort.highlow > all.ALSA.aln_rmdup.sort.highlow.txt
paste *ALSH_Sagallac2.aln_rmdup.sort.highlow > all.ALSH.aln_rmdup.sort.highlow.txt
paste *ANDK_Sagallac2.aln_rmdup.sort.highlow > all.ANDK.aln_rmdup.sort.highlow.txt
paste *ANES_Sagallac2.aln_rmdup.sort.highlow > all.ANES.aln_rmdup.sort.highlow.txt
paste *ANET_Sagallac2.aln_rmdup.sort.highlow > all.ANET.aln_rmdup.sort.highlow.txt
paste *ANEU_Sagallac2.aln_rmdup.sort.highlow > all.ANEU.aln_rmdup.sort.highlow.txt
paste *dysgalactiaelac2.aln_rmdup.sort.highlow > all.dysgalactiae.aln_rmdup.sort.highlow.txt
paste *gordoniilac2.aln_rmdup.sort.highlow > all.gordonii.aln_rmdup.sort.highlow.txt
paste *infantariuslac2.aln_rmdup.sort.highlow > all.infantarius.aln_rmdup.sort.highlow.txt
paste *mutanslac2.aln_rmdup.sort.highlow > all.mutans.aln_rmdup.sort.highlow.txt
paste *oralislac2.aln_rmdup.sort.highlow > all.oralis.aln_rmdup.sort.highlow.txt
paste *parasanguinislac2.aln_rmdup.sort.highlow > all.parasanguinis.aln_rmdup.sort.highlow.txt
paste *pneumoniaelac2.aln_rmdup.sort.highlow > all.pneumoniae.aln_rmdup.sort.highlow.txt
paste *pyogeneslac2.aln_rmdup.sort.highlow > all.pyogenes.aln_rmdup.sort.highlow.txt
paste *.sanguinislac2.aln_rmdup.sort.highlow > all.sanguinis.aln_rmdup.sort.highlow.txt
paste *Sdd_Sdysgallac2.aln_rmdup.sort.highlow > all.Sdd_Sdysgal.aln_rmdup.sort.highlow.txt
paste *Sdd_Sdysgal-revlac2.aln_rmdup.sort.highlow > all.Sdd_Sdysgal-rev.aln_rmdup.sort.highlow.txt
paste *staphaureus-revlac2.aln_rmdup.sort.highlow > all.staphaureus.aln_rmdup.sort.highlow.txt
paste *suislac2.aln_rmdup.sort.highlow > all.suis.aln_rmdup.sort.highlow.txt
paste *uberislac2.aln_rmdup.sort.highlow > all.uberis.aln_rmdup.sort.highlow.txt
paste *uberis-revlac2.aln_rmdup.sort.highlow > all.uberis-rev.aln_rmdup.sort.highlow.txt

rm *.highlow

rename 's/all./D-all./' all.*
rename 's/all./E-all./' all.*
rename 's/all./F-all./' all.*
rename 's/all./G-all./' all.*
rename 's/all./Z-Plaque-all./' all.*
rename 's/all./Z-Skin-all./' all.*
rename 's/all./B-Gorilla-all./' all.*

```

#### lac2operon Damage Patterns

We also want to check whether the mapped reads also exhibit C to T fragmentation
patterns indicative of authentic ancient DNA (i.e. that these lac2operons come
from stuff originally in the calculus and not from modern environmental 
contmamination). For this run `mapDamage` as per `116-imv-mapdamageset.sh`, 
submitted with

```bash
sbatch -c 2 --mem=4000 ./116-imv-getcoverage.sh
```

#### Kaiju Streptoccus Proportion Estimation

#### Complete lac2operon read assembly

### Amylase-binding proteins apbA, abpB

in /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/
make samples.paths.list

> Note: Extra sequencing data was generated for FUM002, GOY006 GDN001 PES001 in
> the form of non-UDG treated single stranded libraries. These have been added
> to this analysis only and is not in the data described above.

```bash
awk '{print "/projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/"$1"/"$1"_S0_L000_R1_000.fastq.merged.prefixed.hg19unmapped.fq.gz"}' newforfasta.list > samples.paths.list
```

also run the Radcliffe samples through this
/projects1/microbiome_sciences/raw_data/external/velsko2017

in an interactive session run 

```bash
ls | while read line; do cat $line/$line.fastq.gz $line/$line.pair1.fastq.gz $line/$line.pair2.fastq.gz $line/$line.trunc.fastq.gz > $line.all.fastq.gz; done &
```

then run 

```bash
cat samples.paths.list | while read line; do seqtk seq -a $line > ./$(basename $line .fastq.gz).fasta; done &
```

```bash
cat samples.paths.list | while read line; do seqtk seq -a $line > ./$(basename $line .fastq).fasta; done
```

`117-imv-qtoa-extra.sh`

make the BLAST database

```bash
makeblastdb -in abpA.selected.fasta -dbtype nucl -title apbA -out abpA
makeblastdb -in abpB.selected.fasta -dbtype nucl -title apbB -out abpB
makeblastdb -in abpOther1.selected.fasta -dbtype nucl -title apbOther1 -out abpOther1
```

`118-imv-abpblast.sh`
``119-imv-abpBblast.sh`

```bash
rename 's/.fastq.combined.fq.prefixed.extractunmapped.bam.fasta.gz.abpA.out.gz/.abpA.out.gz/' *.fastq.combined.fq.prefixed.extractunmapped.bam.fasta.gz.abpA.out.gz
```

```bash
rename 's/.fastq.combined.fq.prefixed.extractunmapped.bam.fasta.gz.abpB.out.gz/.abpB.out.gz/' *.fastq.combined.fq.prefixed.extractunmapped.bam.fasta.gz.abpB.out.gz
```

`120-imv-abp-gethits.sh`
`125-imv-abp-hitlist.sh`
`125b-imv-abp-hitlist.sh`

remove the first line from all hits.list files (line is query_name)

```bash
for f in *.list; do tail -n +2 $f > $(basename $f .list).c.list; done
```

 pull out the reads that were identified as hits by BLAST

`126-abpfastq.sh`

and for the Radcliffe samples too
pull out the reads that were identified as hits by BLAST

`126c-abpfastq.sh`
`126cR-abpfastq.sh`


these ones were missed

```bash
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC001.A0101/VLC001.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpA-blastout/VLC001.A0101_S0_L001_R1_001.abpA.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAfastq/VLC001.A0101.abpA.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC002.A0101/VLC002.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpA-blastout/VLC002.A0101_S0_L001_R1_001.abpA.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAfastq/VLC002.A0101.abpA.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC003.A0101/VLC003.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpA-blastout/VLC003.A0101_S0_L001_R1_001.abpA.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAfastq/VLC003.A0101.abpA.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC004.A0101/VLC004.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpA-blastout/VLC004.A0101_S0_L001_R1_001.abpA.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAfastq/VLC004.A0101.abpA.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC005.A0101/VLC005.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpA-blastout/VLC005.A0101_S0_L001_R1_001.abpA.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAfastq/VLC005.A0101.abpA.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC008.A0101/VLC008.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpA-blastout/VLC008.A0101_S0_L001_R1_001.abpA.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAfastq/VLC008.A0101.abpA.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC009.A0101/VLC009.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpA-blastout/VLC009.A0101_S0_L001_R1_001.abpA.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAfastq/VLC009.A0101.abpA.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC010.A0101/VLC010.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpA-blastout/VLC010.A0101_S0_L001_R1_001.abpA.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAfastq/VLC010.A0101.abpA.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/JAE007.A0101/JAE007.A0101 _S0_L001_R1_001.fastq.merged.prefixed.hg19unmapped.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpA-blastout/JAE007.A0101_S0_L001_R1_001.abpA.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAfastq/JAE007.A0101.abpA.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC001.A0101/VLC001.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpB-blastout/VLC001.A0101_S0_L001_R1_001.abpB.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBfastq/VLC001.A0101.abpB.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC002.A0101/VLC002.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpB-blastout/VLC002.A0101_S0_L001_R1_001.abpB.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBfastq/VLC002.A0101.abpB.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC003.A0101/VLC003.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpB-blastout/VLC003.A0101_S0_L001_R1_001.abpB.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBfastq/VLC003.A0101.abpB.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC004.A0101/VLC004.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpB-blastout/VLC004.A0101_S0_L001_R1_001.abpB.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBfastq/VLC004.A0101.abpB.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC005.A0101/VLC005.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpB-blastout/VLC005.A0101_S0_L001_R1_001.abpB.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBfastq/VLC005.A0101.abpB.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC008.A0101/VLC008.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpB-blastout/VLC008.A0101_S0_L001_R1_001.abpB.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBfastq/VLC008.A0101.abpB.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC009.A0101/VLC009.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpB-blastout/VLC009.A0101_S0_L001_R1_001.abpB.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBfastq/VLC009.A0101.abpB.fastq
zcat /projects1/microbiome_calculus/evolution/03-preprocessing/screening/library_merging/VLC010.A0101/VLC010.A0101_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.fq.gz | seqkit grep -f /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpB-blastout/VLC010.A0101_S0_L001_R1_001.abpB.out.gz.hits.c.list > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBfastq/VLC010.A0101.abpB.fastq
```


now index the amylase binding proteins for bwa in an interactive session

in both 

/projects1/microbiome_calculus/evolution/01-data/databases/amylaseBPs/amylaseBPA/individual
/projects1/microbiome_calculus/evolution/01-data/databases/amylaseBPs/amylaseBPB/individual

```bash
for f in *.fasta; do bwa index -p `basename $f .fasta` $f; done
```

align the reads to each abpA sequence 
`128-imv-alignabpa.sh`
align the reads to each abpB sequence 
`128b-imv-alignapbb.sh`

run in an interactive session
the maximum # of arrays is 1000 and there re >4000 sam files

```bash
srun -c 1 --mem 2000 -t 0-01:00 --pty -u bash -i
source ~/.bash_profile
for f in *.sam; do samtools view -bS $f | samtools sort -o ./$(basename $f .sam).bam -T $(basename $f .sam); done &

for f in *.bam; do samtools index $f ./$(basename $f .bam).bai; done &
```

`129-imv-abprmdupgetmapped.sh`
run in an interactive session
`130-imv-abpalignmentstats.sh`
run in an interactive session
get a matrix of all reads that aligned to each amylase-binding protein (includes multiple alignments)

```bash
paste /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAbwaout/*.all.alignments > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAbwaout/all.abpA.alignments
paste /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBbwaout/*.all.alignments > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBbwaout/all.abpB.alignments
```

get a matrix of all unique reads that aligned to each amylase-binding protein

```bash
paste /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAbwaout/*.uniq.alignments > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAbwaout/uniq.abpA.alignments
paste /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBbwaout/*.uniq.alignments > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBbwaout/uniq.abpB.alignments
```

get a matrix of all unique non-duplicate reads that aligned to each amylase-binding protein

```bash
paste /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAbwaout/*.rmdup.alignments > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpAbwaout/rmdup.abpA.alignments
paste /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBbwaout/*.rmdup.alignments > /projects1/microbiome_calculus/evolution/04-analysis/screening/amylaseBPs/abpBbwaout/rmdup.abpB.alignments
```


`113-imv-uniqrmdupcounts.sh` 
run in an interactive session

map again with looser parameters to see if any additional samples pick up either abpA or abpB
`139a-imv-mapabpa-l16.sh`
`139b-imv-mapapbb-l16.sh`

run in an interactive session

the maximum # of arrays is 1000 and there re >4000 sam files

```bash
srun -c 1 --mem 2000 -t 0-01:00 --pty -u bash -i
source ~/.bash_profile
for f in *.sam; do samtools view -bS $f | samtools sort -o ./$(basename $f .sam).bam -T $(basename $f .sam); done &
```

get the stats `142-imv-getabpstats.sh`

`143-imv-abpA_beast.sh`

`143b-imv-abpB_beast.sh`


# Misc

To convert PDF to PNG on Ubuntu

```bash
for i in *.pdf; do
   pdftoppm -png -rx 300 -ry 300 $i ${i%.pdf*}
done
```