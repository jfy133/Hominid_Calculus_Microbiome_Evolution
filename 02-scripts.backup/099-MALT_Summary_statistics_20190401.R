## In BASH
# DB="<SELECT_HERE>" && grep -e "Loading MEGAN File:" -e "Total reads:" -e "With hits:" -e "Alignments:" -e "Assig. Taxonomy" -e "Min-supp. changes" -e "Numb. Tax. classes:" -e "Class. Taxonomy:" -e "Num. of queries:" -e "Aligned queries:" -e "Num. alignments:" -e "MinSupport set to:" /projects1/microbiome_calculus/evolution/04-analysis/screening/malt/nt/*log | cut -d":" -f 2-99 > /projects1/microbiome_calculus/evolution/00-documentation.backup/99-maltAlignedReadsSummary_raw_"$DB"_$(date "+%Y%m%d").txt

library(tidyverse)

## or nt
db <- "refseq_bacarchhomo_gcs_20181122" ## nt or refseq_bacarchhomo_gcs_20181122

if (db == "nt") {
  raw_data <- read_lines(paste("/home/fellows/projects1/microbiome_calculus/evolution/00-documentation.backup/99-maltAlignedReadsSummary_raw_", db,"_", format(Sys.time(), "%Y%m%d"), ".txt", sep = "")) %>% 
    as_tibble()
} else if (db == "refseq_bacarchhomo_gcs_20181122") {
  raw_data <- read_lines(paste("/home/fellows/projects1/microbiome_calculus/evolution/00-documentation.backup/99-maltAlignedReadsSummary_raw_", db,"_", format(Sys.time(), "%Y%m%d"), ".txt", sep = "")) %>% 
    as_tibble()
}


clean_data <- raw_data %>% 
  separate(value, into = c("variable", "value"), sep = ":") %>% 
  mutate(variable = str_trim(variable), value = str_trim(value), value = gsub(",", "", value)) %>% 
  mutate(objnum = cumsum(variable == "Loading MEGAN File")) %>%
  spread(variable, value) %>%
  select(`Loading MEGAN File`, 
         `Num. of queries`,
         `With hits`,
         `Alignments`,
         `Assig. Taxonomy`,
         `MinSupport set to`) %>%
  rename(File = `Loading MEGAN File`,
         Input_Reads = `Num. of queries`,
         Input_Reads_With_Hits = `With hits`,
         Total_Alignments = `Alignments`,
         Input_Reads_With_Taxonomy = `Assig. Taxonomy`,
         Min_Support_Reads_Threshold_MALT = `MinSupport set to`) %>%
  mutate(Input_Reads = as.numeric(Input_Reads),
         Input_Reads_With_Hits = as.numeric(Input_Reads_With_Hits),
         Total_Alignments = as.numeric(Total_Alignments),
         Input_Reads_With_Taxonomy = as.numeric(Input_Reads_With_Taxonomy),
         Min_Support_Reads_Threshold_MALT = as.numeric(Min_Support_Reads_Threshold_MALT)) %>%
  mutate(Percentage_Hits = Input_Reads_With_Hits / Input_Reads * 100,
         Percentage_Hits_with_Taxonomy = Input_Reads_With_Taxonomy / Input_Reads * 100) %>%
         filter(!is.na(Input_Reads)) %>%
  distinct()

write_tsv(clean_data, paste("/home/fellows/projects1/microbiome_calculus/evolution/00-documentation.backup/05-MALT_taxonomic_binning_summary_statistics_", db, "_", format(Sys.time(), "%Y%m%d"),".tsv", sep = ""))

          