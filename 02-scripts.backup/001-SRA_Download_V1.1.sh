#!/usr/bin/env bash 

### Takes as input a text file containing a list of SRA Experiment (SRX), Sample (SRS) (SRR) and Run numbers []not study (SRP)], and the output directory!
### Downloads with wget a .sra file into a SRA specific directory and converts into .fastq.gz format with SRAtoolkit's fastq-dump
### Written by James. A Fellows Yates (jfy133@gmail.com) on 2016-08-03

INPUTLIST=$1
OUTPUTDIR=$(basename $2)

while read CODE; do
        if [ "$(echo $CODE | cut -c1-3)" == "SRR" ]; then
                mkdir $OUTPUTDIR/$("$CODE")
                wget -P $OUTPUTDIR/$CODE/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/"$(echo $CODE | cut -c1-6)"/"$CODE"/"$CODE".sra
                sratoolkit.2.8.0-ubuntu64/bin/fastq-dump -F --split-files --readids $OUTPUTDIR/"$CODE"/"$CODE".sra -O $OUTPUTDIR/"$CODE" --gzip
                rm $OUTPUTDIR/"$CODE"/"$CODE".sra

        elif [ "$(echo $CODE | cut -c1-3)" == "ERR" ]; then
                mkdir $OUTPUTDIR/$("$CODE")
                wget -P $OUTPUTDIR/$CODE/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/"$(echo $CODE | cut -c1-6)"/"$CODE"/"$CODE".sra
                sratoolkit.2.8.0-ubuntu64/bin/fastq-dump -F --split-files --readids $OUTPUTDIR/"$CODE"/"$CODE".sra -O $OUTPUTDIR/"$CODE" --gzip
                rm $OUTPUTDIR/"$CODE"/"$CODE".sra

        elif [ "$(echo $CODE | cut -c1-3)" == "SRX" ]; then
                mkdir $OUTPUTDIR/$("$CODE")
                wget -r -nH -nd -np -P $OUTPUTDIR/$CODE/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/"$(echo $CODE | cut -c1-6)"/"$CODE"/ -A "*.sra"
                #sratoolkit.2.8.0-ubuntu64/bin/fastq-dump $OUTPUTDIR/"$CODE"/"$CODE".sra -O $OUTPUTDIR/"$CODE" --gzip
                #rm $OUTPUTDIR/"$CODE"/"$CODE".sra
                
        elif [ "$(echo $CODE | cut -c1-3)" == "SRS" ]; then
                mkdir $OUTPUTDIR/$("$CODE")
                wget -r -nH -nd -np -P $OUTPUTDIR/$CODE/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/BySample/sra/SRS/"$(echo $CODE | cut -c1-6)"/"$CODE"/ -A "*.sra"
                #sratoolkit.2.8.0-ubuntu64/bin/fastq-dump $OUTPUTDIR/"$CODE"/"$CODE".sra -O $OUTPUTDIR/"$CODE" --gzip
                #rm $OUTPUTDIR/"$CODE"/"$CODE".sra

        else
                echo "$CODE: code not recognised"
        fi
done <$INPUTLIST
