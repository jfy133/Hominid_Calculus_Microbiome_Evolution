#!/usr/bin/env bash

path=$1
species="$(echo $1 | rev | cut -d/ -f2 | rev)"

gunzip "$path"
sed -i "s/ /_$species /" ${path%.gz}
#sed "s/    /_$species /" ../01-data/genomes/"$genus"/"$species"/*gff
gzip ${path%.gz}

exit 0
