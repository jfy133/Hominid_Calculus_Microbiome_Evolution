library(tidyverse)
library(data.table)

## Get my superreference DNA seuqence coordinates
superref_genome_coords <- fread("/home/fellows/projects1/microbiome_calculus/evolution/01-data/genomes/Streptococcus/collapsed_Streptococcus_superreference.coords", col.names = c("input_file", "fasta_entry_header", "fasta_entry_length", "start_coordinate", "end_coordinate")) %>% 
  as_tibble() %>%
  mutate(file_name = basename(input_file) %>% tools::file_path_sans_ext(.) %>% tools::file_path_sans_ext(.),
         seqname = map(fasta_entry_header, ~str_split(.x, "_")[[1]][1]) %>% unlist)

## Get locations
gff <- Sys.glob("/home/fellows/projects1/microbiome_calculus/evolution/01-data/genomes/Streptococcus/*/") %>%
  list.files(".gff.gz", full.names = T) %>%
  enframe(name = NULL, value = "file_path") %>%
  mutate(Taxon = map(file_path, ~.x %>% str_split(., "\\/") %>% unlist %>% .[10]) %>% unlist,
         file_contents = map(file_path, ~read_tsv(.x, 
                                                  comment = '#', 
                                                  col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
                                                  col_types = cols(seqname = "c", source = "c", feature = "c", start = "d", end = "d", score = "c", strand = "c", frame = "c", attribute = "c")) %>% 
                               filter(feature == "gene"))) %>%
  mutate(file_name = basename(file_path) %>% tools::file_path_sans_ext(.) %>% tools::file_path_sans_ext(.)) %>% 
  unnest() %>% 
  mutate(attribute = map(attribute, ~.x %>% str_split(";") %>% unlist %>% .[grepl("locus_tag", .)] %>% paste0(collapse = ";") %>% gsub("locus_tag=", "", .)) %>% unlist)


## Get selected abpA/b similar annotations
list_abpA <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/streptococcus_investigation.backup/panX/selected_panX_abpA_annotations.tsv") %>% separate(Annotation, into = c("Assembly", "locus_tag", "number", "description", "extra"), "-")

list_abpB <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/streptococcus_investigation.backup/panX/selected_panX_abpB_annotations.tsv") %>% separate(Annotation, into = c("Assembly", "locus_tag", "number", "description", "extra"), "-")

## Filter gffs to get locations on original genome

abpA_locs <- gff %>% filter(attribute %in% (list_abpA %>% pull(locus_tag))) 
abpB_locs <- gff %>% filter(attribute %in% (list_abpB %>% pull(locus_tag))) 

## note that D2908_07695 is not present in my annotations because it wasn't 
## annotated at the time - and has since been done. Will ignore for the 
## purposes of this project.

## Combine the reference coordinates and abp* data

abpA_combined <- left_join(abpA_locs, superref_genome_coords)

abpB_combined <- left_join(abpB_locs, superref_genome_coords)


## calculate where the gene is in the superreference - while have to check 
## which contig the gene falls in for non-complete assemblies

abpA_bed <- abpA_combined %>% mutate(gene_start_coordinate = (start_coordinate + start) - 1, 
                         gene_end_coordinate = (start_coordinate + end) - 1)

abpB_bed <- abpB_combined %>% mutate(gene_start_coordinate = (start_coordinate + start) - 1, 
                         gene_end_coordinate = (start_coordinate + end) - 1)

## Now create bed files
abpA_bed <- abpA_bed %>% 
  mutate(chrom = "Streptococcus_superreference") %>% 
  select(chrom, gene_start_coordinate, gene_end_coordinate, Taxon, file_name, attribute) %>%
  unite(col = "name", Taxon, file_name, attribute, sep = "-") %>%
  write_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/streptococcus_investigation.backup/Streptococcus_superreference_coordinates_abpAlike_genes.bed", col_names = F)

abpB_bed <- abpB_bed %>% 
  mutate(chrom = "Streptococcus_superreference") %>% 
  select(chrom, gene_start_coordinate, gene_end_coordinate, Taxon, file_name, attribute) %>%
  unite(col = "name", Taxon, file_name, attribute, sep = "-") %>%
  write_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/streptococcus_investigation.backup/Streptococcus_superreference_coordinates_abpBlike_genes.bed", col_names = F)

