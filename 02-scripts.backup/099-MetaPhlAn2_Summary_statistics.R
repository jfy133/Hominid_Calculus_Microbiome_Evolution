library(tidyverse)

## Input estimated species assigned reads and clean up
raw_mp2 <- read_tsv("../04-analysis/screening/metaphlan2/output_readcounts/metaphlan2_merged_estimatedreadcount_table_20190401.txt", comment = "") %>%
  filter(grepl("s__", ID)) %>% 
  filter(!grepl("t__", ID)) %>%
  rowwise() %>%
  mutate(ID = str_split_fixed(ID, "\\|s__", n = 2)[2])

## Calculate 0.01 percent of total number of species hits per sample 
thresholds <- colSums(raw_mp2[2:ncol(raw_mp2)]) * 0.0001

tibble(Min_Support_Reads_Threshold_MP2 = round(thresholds, 0), 
       Individual = names(thresholds)) %>% 
  write_tsv(paste("../00-documentation.backup/05-MetaPhlAn2_taxonomic_binning_summary_statistics_", format(Sys.Date(), "%Y%m%d"), ".tsv"))


