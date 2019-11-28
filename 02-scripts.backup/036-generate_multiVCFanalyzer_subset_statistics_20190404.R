#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

if (is.na(args[1]) | args[1] == "" ) {
  stop("Usage: Rscript generate_multiVCFanalzer_subset_statistics.R <multiVCFanalyzer_outdir> <bed_file> <bed_feature_string> 

        Description: Given a MultiVCFAnalyzer 0.87 output directory, this script 
        will output a subsetted snpTable, genotype matrix, FASTA alignment and 
        snpStatistics file based on the features of interest as indicated by 
        your bed file. Output will be in the directory of the input table.

        Requires the following R packages: tidyverse, data.table

        TODO: Input validation, Error messages, accept non-GZIP.

        Written by: James A. Fellows Yates, 2019
       ")
} else if (args[1] == "-h" | args[1] == "--help") {
  cat("Usage: Rscript generate_multiVCFanalzer_subset_statistics.R <multiVCFanalyzer_outdir> <bed_file> <bed_feature_string> 

       Description: Given a MultiVCFAnalyzer 0.87 output directory, this script 
       will output a subsetted snpTable, genotype matrix, FASTA alignment and 
       snpStatistics file based on the features of interest as indicated by 
       your bed file. Output will be in the directory of the input table.

       Requires the following R packages: tidyverse, data.table

       TODO: Input validation, Error messages, accept non-GZIP.

       Written by: James A. Fellows Yates, 2019
      ")
  quit()
}

library(tidyverse)
library(data.table)

input_dir <- args[1]
input_bed <- args[2]
input_feature <- args[3]

## Testing ##
 #input_dir <- "../04-analysis/deep/multivcfanalyzer/superreference_mapping/output/Treponema/"
 #input_bed <- "../01-data/genomes/Treponema/collapsed_Treponema_superreference.bed"
 #input_feature <- "Treponema_socranskii"
 #out_dir <- "~/Downloads/"
#############

print("loading SNP table")

data_snptable <- fread(paste(input_dir, "/snpTable.tsv.gz", sep = "")) %>% 
  as_tibble()

print("loading Genotype Matrix")

data_genotypematrix <- fread(paste(input_dir, "/genotyeMatrix.tsv.gz", sep = "")) %>% 
  as_tibble()

print("loading bed file")

data_bedfile <- fread(input_bed) %>% 
  as_tibble() %>% 
  rename(contig = V1, start = V2, stop = V3, feature = V4) %>% 
  mutate(start = start + 1)

print("extracting your feature of interest")

data_feature_range <- c()

data_feature_range <- data_bedfile %>% 
  filter(grepl(input_feature, feature)) %>%
  head(n = 1) %>%
  select(start) %>%
  append(data_feature_range)

data_feature_range <- data_bedfile %>% 
  filter(grepl(input_feature, feature)) %>%
  tail(n = 1) %>%
  select(stop) %>%
  append(data_feature_range)

print("filter")

data_snptable_filtered <- data_snptable %>% 
  filter(Position >= data_feature_range$start, 
         Position <= data_feature_range$stop) %>%
  write_tsv(., paste(input_dir, "/snpTable_subset.tsv.gz", sep = ""))

data_genotypematrix_filtered <- data_genotypematrix %>% 
  filter(Position >= data_feature_range$start, 
         Position <= data_feature_range$stop) %>%
  write_tsv(., paste(input_dir, "/genotyeMatrix_subset.tsv.gz", sep = ""))

print("calculating statistics (Please ignore 'expected 2 pieces' warnings)")

data_snptable_filtered %>%
  gather(sample, Call, 3:ncol(.)) %>%
  mutate(Call = gsub("\\(|\\)", "", Call)) %>%
  separate(Call, c("Call", "Frequency"), sep = " ") %>%
  group_by(sample) %>%
  summarise(
    `SNP Calls (all)` = sum(Call != "N" & Call != "."),
    `SNP Calls (het)` = sum(!is.na(Frequency) & 
                              !Call %in% c("C", "G", "T", "A")),
    `coverage(fold)` = NA,
    `coverage(percent)` = NA,
    `refCall` = NA,
    `allPos` = NA,
    `noCall` = NA,
    `discardedRefCall` = NA,
    `discardedVarCall` = NA,
    `filteredVarCall` = NA,
    `unhandledGenotype` = NA
  ) %>%
    write_tsv(., paste(input_dir, "/snpStatistics_subset.tsv.gz", sep = ""))

print("creating FASTA (Please ignore 'expected 2 pieces' warnings)")

data_fasta_tab <- data_snptable_filtered %>%
  gather(sample, Call, 3:ncol(.)) %>%
  filter(Call != "N") %>%
  mutate(Call = if_else(Call == ".", Ref, Call)) %>%
  separate(Call, c("Call", "Frequency")) %>%
  select(sample, Position, Call) %>%
  spread(Position, Call, fill = "N")


write_tsv(data_fasta_tab, 
          paste(input_dir, "/snpAlignment_subset.fasta", sep = ""),
          col_names = F)

system(paste("awk '{print \">\"$1; $1=\"\"; print $0}' OFS= ", input_dir, "/snpAlignment_subset.fasta | gzip > ", input_dir, "/snpAlignment_subset.fasta.gz", sep = ""))

system(paste("rm ", input_dir, "/snpAlignment_subset.fasta", sep = ""))
           