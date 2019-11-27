library(tidyverse)
library(ggalluvial)
library(gridExtra)

################
## Load files ##
################

main_dir <- "../04-analysis/screening/presenceabsence_intersection/"

individual_thresholds <- read_tsv("../04-analysis/screening/presenceabsence_intersection/04-coremicrobiome_samplenumberthresholdcalculation_hostpopulations_alldbs_allthresholdlevels_alltaxalevels_20190204.tsv")

files_taxapassedthresholdstats <- list.files(main_dir, pattern = '05-coremicrobiome_samplenumberthresholdpassed_hostpopulation_*', full.names = T, recursive = F)
files_taxapassedthresholdslist <- list.files(main_dir, pattern = '06-coremicrobiome_microbialtaxapassingthreshold_hostgenus_*', full.names = T, recursive = F)
files_intersectiontable <- list.files(main_dir, pattern = "08a-coremicrobiome_presenceabsence_upsettable_allsoftware_*", full.names = T, recursive = F) %>% .[grepl("hostgenus",.)]

collate_data <- function(file_list){
  result <- data_frame(filename = file_list) %>%
    mutate(file_contents = purrr::map(filename, ~ suppressWarnings(suppressMessages(as_tibble(read_tsv(.)))))) %>%
    rowwise()  %>%
    mutate(Grouping = str_split(filename, "_")[[1]][[5]],
           SoftwareIncluded = str_split(filename, "_")[[1]][[6]],
           MALT_DB = str_split(filename, "_")[[1]][[7]],
           MetaPhlAn2_DB = str_split(filename, "_")[[1]][[8]],
           Threshold = str_split(filename, "_")[[1]][[9]]) %>%
    select(-filename) %>%
    unnest()
  return(result)
}


#################################################
## List of each Taxa and population prevalence ##
#################################################

summary_taxapopulationprevalence <- collate_data(files_taxapassedthresholdslist) %>%
  select(Software, MALT_DB, MetaPhlAn2_DB, Tax_Level, Threshold, Host_Genus, Taxon, No_Pops_Present, No_Pops_Group, Proportion_Pop_Groups_Presence) %>%
  mutate(MALT_DB = gsub("maltdb", "", MALT_DB),
           MetaPhlAn2_DB = gsub("metaphlan2db", "", MetaPhlAn2_DB),
           MetaPhlAn2_DB = if_else(Software == "MALT", NA_character_, MetaPhlAn2_DB),
           MALT_DB = if_else(Software == "MetaPhlAn2", NA_character_, MALT_DB)
           ) %>%
  rename(Populations_Taxa_Present_In = No_Pops_Present,
         Populations_in_Study = No_Pops_Group,
         Proportion_Taxa_Present_In = Proportion_Pop_Groups_Presence) %>%
  replace_na(list(MALT_DB = "", MetaPhlAn2_DB = "")) %>%
  unite(Database, MALT_DB, MetaPhlAn2_DB, sep = "") 

# write_tsv(paste("../00-documentation.backup/24-intersection_proktaxapassingthresholdstaxalist_", format(Sys.Date(), "%Y%m%d"), ".tsv", sep = ""), summary_taxapopulationprevalence ) %>% 
#   print()

###############################################################
## Summary of No. of Taxa at each Min Support per Population ##
###############################################################

data_taxapassedthresholdstats <-  data_frame(filename = files_taxapassedthresholdstats) %>%
  mutate(file_contents = purrr::map(filename, ~ suppressWarnings(suppressMessages(as_tibble(read_tsv(.)))))) %>%
  rowwise()  %>%
  mutate(Grouping = str_split(filename, "_")[[1]][[5]],
         SoftwareIncluded = str_split(filename, "_")[[1]][[6]],
         MALT_DB = str_split(filename, "_")[[1]][[7]],
         MetaPhlAn2_DB = str_split(filename, "_")[[1]][[8]],
         Min_Support = str_split(filename, "_")[[1]][[9]]) %>%
  select(-filename) %>%
  unnest()

summary_taxapassedthresholdstats <- data_taxapassedthresholdstats %>% 
  spread(Min_Support, No_Taxa_Per_Pop) %>% 
  filter(Software == "MALT") %>% 
  mutate(MALT_DB = gsub("maltdb", "", MALT_DB)) %>%
  select(-Grouping, -SoftwareIncluded, -MetaPhlAn2_DB) %>%
  select(Software, MALT_DB, Tax_Level, Host_Genus, Host_Common, Env, `0.01`, `0.02`, `0.05`)

# write_tsv(summary_taxapassedthresholdstats, paste("../00-documentation.backup/23-intersection_proktaxapassingthresholdsstats_", format(Sys.Date(), "%Y%m%d"), ".tsv", sep = ""))

summary_taxapassedthresholdstats2 <- summary_taxapassedthresholdstats %>%
  mutate(Host_Genus = as_factor(Host_Genus)) %>% 
  mutate(Host_Genus = fct_relevel(Host_Genus, c("Control", "Alouatta", "Gorilla", "Pan", "Homo"))) %>% 
  gather(Threshold, Count, c(`0.05`, `0.02`, `0.01`))

summary_taxapassedthresholdstats_plot <- ggplot(summary_taxapassedthresholdstats2, aes(Host_Genus, Count, fill = Threshold)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Dark2") + 
  theme_minimal(base_family = "Roboto", base_size = 7) + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(MALT_DB ~ Tax_Level)

 ggsave(paste("11-coremicrobiome_presenceabsence_filterthresholdbarpot_nt_alldb_alltaxlevel_", format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), 
        summary_taxapassedthresholdstats_plot, 
        path = "../04-analysis/screening/presenceabsence_intersection/", 
        device = cairo_pdf, 
        width = 3.5, 
        height = 3.5, 
        units = "in",
        dpi = 600)

#################################################################################
## Summary of effect of core microbiome calculation each min support threshold ##
#################################################################################


data_intersectiontables <- data_frame(filename = files_intersectiontable) %>%
  mutate(file_contents = purrr::map(filename, ~ suppressWarnings(suppressMessages(as_tibble(read_tsv(.)))))) %>%
  rowwise() %>%
  mutate(MALT_DB = str_split(filename, "_")[[1]][[7]],
         MetaPhlAn2_DB = str_split(filename, "_")[[1]][[8]],
         Min_Support = str_split(filename, "_")[[1]][[9]],
         Grouping = str_split(filename, "_")[[1]][[10]],
         Taxon_Level = str_split(filename, "_")[[1]][[11]]) %>%
  select(-filename) %>%
  unnest()


## Get which taxa are in which combination depending on threshold

summary_intersectiontables <- data_intersectiontables %>% 
  filter(Combination %in% c("Alouatta:Gorilla:Pan:Homo", "Gorilla:Pan:Homo", "Pan:Homo")) %>% 
  mutate(Combination = case_when(Combination == "Pan:Homo" ~ "Hominini",
                                 Combination == "Gorilla:Pan:Homo" ~ "Homininae",
                                 Combination == "Alouatta:Gorilla:Pan:Homo" ~ "Anthropoid",
                                 TRUE ~ "Other Comb.")) %>%
  mutate(Present = T) %>%
  filter(Software == "MALT", Tax_Level == "genus") %>%
  spread(Min_Support, Combination) %>%
  replace_na(list(`0.01` = "Other Comb.",`0.02` = "Other Comb.",`0.05` = "Other Comb.")) %>%
  filter(MALT_DB == "maltdbnt")


group_order <- c("Hominini", "Homininae", "Anthropoid", "Other Comb.")  
summary_intersectiontables$`0.01` <- factor(summary_intersectiontables$`0.01`, levels = group_order)
summary_intersectiontables$`0.02` <- factor(summary_intersectiontables$`0.02`, levels = group_order)
summary_intersectiontables$`0.05` <- factor(summary_intersectiontables$`0.05`, levels = group_order)

## Load and clean OTU table to get average number of alignments for each taxon

data_cleaner <- function(x, meta_dat, supp_col){
  colnames(x) <- gsub("_S.*_L.*_R1_.*.fastq.combined.fq.prefixed.extractunmapped.bam", "", colnames(x)) 
  colnames(x) <- gsub("_S.*_L.*_R1_.*.fastq.merged.prefixed.hg19unmapped", "", colnames(x))
  colnames(x) <- gsub("_S0_L001_R1_001.fastq.extractunmapped.bam", "", colnames(x))
  colnames(x) <- gsub(".mp2profile", "", colnames(x))
  colnames(x)[1] <- "Taxon"
  x <- x %>% 
    gather(Individual, Value, 2:ncol(x)) %>% 
    mutate(Value = as.numeric(Value)) %>%
    left_join(select(meta_dat, Individual, one_of(c(supp_col)), Age)) %>%
    rename_("Threshold" = supp_col) %>%
    mutate(Threshold = as.numeric(Threshold)) %>%
    mutate(Filter_Passed = if_else(Value == 0, 
                                   1,
                                   if_else(Value >= Threshold, 1, 0))) %>% 
    filter(Filter_Passed == 1) %>% #,
    #Age != "ModernDay") %>% 
    select(Taxon, Individual, Value)
  return(x)
}


raw_malt_genus <- read_tsv("../04-analysis/screening/megan.backup/Evolution-Comparison_MEGAN_20180817-ex_absolute_genus_prokaryotes_summarised.txt")

metadata <- read_tsv("../00-documentation.backup/02-calculus_microbiome-deep_evolution-individualscontrolssources_metadata_20190228.tsv.tsv") %>%
  rename(Individual = `#SampleID`)

bad_samples <- read_tsv("../00-documentation.backup/06-sourcetracker_samples_to_filter_out_181112.tsv") 

minsupp_malt <- "Min_Support_Reads_Threshold_MALT"

full_data <- data_cleaner(raw_malt_genus, metadata, minsupp_malt)

full_data_meta <- left_join(full_data, 
                            metadata %>% 
                              select(Individual, 
                                     Env, 
                                     Host_Genus, 
                                     Host_Common, 
                                     SourceSink, 
                                     Sample_or_Control,
                                     Individual_Seq_Depth)) %>% 
  filter(#Sample_or_Control == "Sample", 
    SourceSink == "sink")

readcount_summary <- full_data_meta %>% filter(!Individual %in% pull(bad_samples)) %>%
  filter(Host_Genus %in% c("Alouatta", "Gorilla", "Pan", "Homo")) %>%
  group_by(Taxon) %>%
  summarise(Average_Reads = mean(Value)) %>%
  mutate()

## Get threshold range and see how to categorize data

readcount_summary %>% 
  select(Average_Reads) %>% 
  ggplot(., aes(x = Average_Reads)) + 
  geom_histogram(bins = 60) + 
  scale_x_continuous(limits = c(0,7000), breaks = seq(0,7000,200)) + 
  theme(axis.text.x = element_text(angle = 90))

## Assign cateogories 

summary_intersectiontables <- summary_intersectiontables %>% 
  left_join(readcount_summary)

summary_intersectiontables <- summary_intersectiontables %>% mutate(
  Avg_Read_Count = case_when(Average_Reads < 100 ~ "0-99",
                             Average_Reads >= 100 & Average_Reads <= 199 ~ "100-199",
                             Average_Reads >= 200 & Average_Reads <= 399 ~ "200-399",
                             Average_Reads >= 400 & Average_Reads <= 799 ~ "400-799",
                             Average_Reads >= 800 & Average_Reads <= 1599 ~ "800-1599",
                             Average_Reads >= 1600 & Average_Reads <= 3199 ~ "1600-3199",
                             Average_Reads >= 3200 & Average_Reads <= 6399 ~ "3200-6399",
                             Average_Reads >= 6400 ~ "> 6400")
  )

summary_intersectiontables$Avg_Read_Count <- factor(summary_intersectiontables$Avg_Read_Count, 
                                                    levels = c("0-99", "100-199", "200-399", "400-799", "800-1599", "1600-3199", "3200-6399", "> 6400"))
  

## Assign Source Categories

out <- read_tsv("../00-documentation.backup/99-results_bacdive_ohne_references.tsv.gz")

filter(out, field == "sample_type", search_term == "")

for (i in summary_intersectiontables$Taxon) {print(i)}

filter(out, field == "sample_type", search_term %in% list)


## BacDive by Eye
isolation_cats <- c(Acidipropionibacterium = "environmental", 
  Acidivorax = "environmental",
  Actinomyces = "oral",
  Arthrobacer = "environmental",
  Bacillus = "environmental",
  Bifidobacterium = "faeces",
  Bordetella = "environmental/mammalian",
  Brachybacterium = "environmental",
  Campylobacter = "mammalian",
  Capnocytophaga = "oral",
  Cellulomonas = "environmental",
  Clostridioides = "environmental",
  Clostridium = "environmental",
  Corynebacterium = "environmental/mammalian/oral",
  Desulfomicrobium = "environmental",
  Desulfovibrio = "environmental",
  Fretibacterium = "oral",
  Fusobacterium = "oral",
  Kocuria = "environmental",
  Microbacterium = "environmental",
  Micromonospora = "environmental",
  Mogibacterium = NA,
  Mycobacterium = "environmental",
  Neisseria = "oral",
  Nocardia = "environmental/mammalian",
  Olsenella = "oral",
  Ottowia = "environmental",
  Paenibacillus = "environmental",
  Parvimonas = NA,
  Porphyromonas = "oral",
  Prevotella = "oral",
  Propionibacterium = "environmental",
  Pseudopropionibacterium = "oral",
  Rhodococcus = "environmental",
  Selenomonas = "oral",
  Streptococcus = "oral",
  Streptomyces = "environmental",
  Tannerella = "oral",
  Tessaracoccus = "environmental",
  Treponema = "oral"
) %>% bind_rows() %>% gather(Taxon, Source)

## HOMD
isolation_cats2 <- c(Acidipropionibacterium = "oral", 
                    Acidivorax = "oral",
                    Actinomyces = "oral",
                    Arthrobacter = "environmental",
                    Bacillus = "oral",
                    Bifidobacterium = "oral",
                    Bordetella = "oral",
                    Brachybacterium = "environmental",
                    Campylobacter = "oral",
                    Capnocytophaga = "oral",
                    Cellulomonas = "environmental",
                    Clostridioides = "environmental",
                    Clostridium = "environmental",
                    Corynebacterium = "oral",
                    Desulfomicrobium = "oral",
                    Desulfovibrio = "oral",
                    Fretibacterium = "oral",
                    Fusobacterium = "oral",
                    Kocuria = "oral",
                    Microbacterium = "oral",
                    Micromonospora = "oral",
                    Mogibacterium = "oral",
                    Mycobacterium = "oral",
                    Neisseria = "oral",
                    Nocardia = "environmental",
                    Olsenella = "oral",
                    Ottowia = "oral",
                    Paenibacillus = "oral",
                    Parvimonas = "oral",
                    Porphyromonas = "oral",
                    Prevotella = "oral",
                    Propionibacterium = "environmental",
                    Pseudopropionibacterium = "oral",
                    Rhodococcus = "oral",
                    Selenomonas = "oral",
                    Streptococcus = "oral",
                    Streptomyces = "environmental",
                    Tannerella = "oral",
                    Tessaracoccus = "environmental",
                    Treponema = "oral") %>% 
  bind_rows() %>% 
  gather(Taxon, Source2)

## BacDive by eye condensed and google for NAs
isolation_cats3 <- c(Acidipropionibacterium = "environmental", 
                    Acidovorax = "environmental",
                    Actinomyces = "oral",
                    Arthrobacter = "environmental",
                    Bacillus = "environmental",
                    Bifidobacterium = "faeces",
                    Bordetella = "diverse",
                    Brachybacterium = "environmental",
                    Campylobacter = "faeces",
                    Capnocytophaga = "oral",
                    Cellulomonas = "environmental",
                    Clostridioides = "environmental",
                    Clostridium = "environmental",
                    Corynebacterium = "diverse",
                    Desulfomicrobium = "environmental",
                    Desulfovibrio = "environmental",
                    Fretibacterium = "oral",
                    Fusobacterium = "oral",
                    Kocuria = "environmental",
                    Microbacterium = "environmental",
                    Micromonospora = "environmental",
                    Mogibacterium = "oral",
                    Mycobacterium = "environmental",
                    Neisseria = "oral",
                    Nocardia = "diverse",
                    Olsenella = "oral",
                    Ottowia = "environmental",
                    Paenibacillus = "environmental",
                    Parvimonas = "oral",
                    Porphyromonas = "oral",
                    Prevotella = "oral",
                    Propionibacterium = "environmental",
                    Pseudopropionibacterium = "oral",
                    Rhodococcus = "environmental",
                    Selenomonas = "oral",
                    Streptococcus = "oral",
                    Streptomyces = "environmental",
                    Tannerella = "oral",
                    Tessaracoccus = "environmental",
                    Treponema = "oral"
) %>% bind_rows() %>% gather(Taxon, Source3)

summary_intersectiontables <- left_join(summary_intersectiontables, isolation_cats)
summary_intersectiontables2 <- left_join(summary_intersectiontables, isolation_cats2)
summary_intersectiontables3 <- left_join(summary_intersectiontables, isolation_cats3)


## Plot

### Read count
summary_intersectiontables %>% 
  ggplot(aes(axis1 = `0.01`, axis2 = `0.02`, axis3 = `0.05`), y = "Combination") +
  geom_alluvium(aes(fill = Avg_Read_Count)) +
  geom_stratum() +
  geom_text(stat = "stratum", label.strata = TRUE, size = 3) +
  ylab("Number of Taxa") +
  xlab("Minimum Support Threshold") +
  scale_x_continuous(breaks = 1:3, 
                     labels = c("0.01", "0.02", "0.05"), 
                     position = "top") +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  theme_minimal(base_family = "Roboto", base_size = 7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        )

### By Taxon
summary_intersectiontables %>% 
  ggplot(aes(axis1 = `0.01`, axis2 = `0.02`, axis3 = `0.05`), y = "Combination") +
  geom_alluvium(aes(fill = Taxon)) +
  geom_stratum() +
  geom_text(stat = "stratum", label.strata = TRUE, size = 3) +
  ylab("Number of Taxa") +
  xlab("Minimum Support Threshold") +
  scale_x_continuous(breaks = 1:3, 
                     labels = c("0.01", "0.02", "0.05"), 
                     position = "top") +
  theme_minimal(base_family = "Roboto", base_size = 7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )


### Bacdive and Eye
summary_intersectiontables %>% 
  ggplot(aes(axis1 = `0.01`, axis2 = `0.02`, axis3 = `0.05`), y = "Combination") +
  geom_alluvium(aes(fill = Source)) +
  geom_stratum() +
  geom_text(stat = "stratum", label.strata = TRUE, size = 3) +
  ylab("Number of Taxa") +
  xlab("Minimum Support Threshold") +
  scale_x_continuous(breaks = 1:3, 
                     labels = c("0.01", "0.02", "0.05"), 
                     position = "top") +
  theme_minimal(base_family = "Roboto", base_size = 7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

## HOMD
summary_intersectiontables2 %>% 
  ggplot(aes(axis1 = `0.01`, axis2 = `0.02`, axis3 = `0.05`), y = "Combination") +
  geom_alluvium(aes(fill = Source2)) +
  geom_stratum() +
  geom_text(stat = "stratum", label.strata = TRUE, size = 3) +
  ylab("Number of Taxa") +
  xlab("Minimum Support Threshold") +
  scale_x_continuous(breaks = 1:3, 
                     labels = c("0.01", "0.02", "0.05"), 
                     position = "top") +
  theme_minimal(base_family = "Roboto", base_size = 7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

## BacDive and Eye subset and Google
alluvial_plot <- summary_intersectiontables3 %>% 
  ggplot(aes(axis1 = `0.01`, axis2 = `0.02`, axis3 = `0.05`), y = "Combination") +
  geom_alluvium(aes(fill = Source3)) +
  geom_stratum() +
  geom_text(stat = "stratum", label.strata = TRUE, size = 2, angle = -45) +
  ylab("Number of Taxa") +
  xlab("Minimum Support Threshold") +
  scale_x_continuous(breaks = 1:3, 
                     labels = c("0.01", "0.02", "0.05"), 
                     position = "top") +
  theme_minimal(base_family = "Roboto", base_size = 7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

alluvial_plot
  
  
# ggsave(paste("12-coremicrobiome_presenceabsence_filterthresholdalluvial_malt_nt_genus_", format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), alluvial_plot, path = "../04-analysis/screening/presenceabsence_intersection/", device = cairo_pdf, width = 3.5, height = 3.5, units = "in", dpi = 600)


#####################################################################################
## Corrobration of Core Microbiome between Databases and Software at 0.05% support ##
#####################################################################################

## Compare number of taxa found in a combination by both MALT and MetaPhlAn2 
compare_software <- function(dat, min_supp, tax_level){
  dat %>% 
    filter(Min_Support == min_supp, Tax_Level == tax_level) %>% 
    mutate(Present = T) %>%
    spread(Software, Present, fill = F) %>%
    mutate(Where_Found = case_when(MALT & MetaPhlAn2 ~ "both",
                                   MALT & !MetaPhlAn2 ~ "MALT_unique",
                                   !MALT & MetaPhlAn2~ "MetaPhlAn2_unique")
    ) %>%
    group_by(Combination, MALT_DB) %>%
    count(Where_Found) %>%
    spread(Where_Found, n) %>%
    replace_na(list(both = 0, MALT_unique = 0, MetaPhlAn2_unique = 0)) %>%
    ungroup() %>%
    rename(MALT_Database = MALT_DB) %>%
    mutate(MALT_Database = gsub("maltdb", "", MALT_Database)) %>%
    arrange(MALT_Database)
}

## Compare number of taxa found in a combination when classified in both nt and 
## refseq databases
compare_database <- function(dat, db_col, soft, min_supp, tax_level){
  dat %>% 
    filter(Software == soft, Min_Support == min_supp, Tax_Level == tax_level) %>%
    mutate(Present = T) %>%
    spread(db_col, Present, fill = F) %>%
    mutate(Where_Found = case_when(maltdbnt & maltdbrefseq ~ "both",
                                   maltdbnt & !maltdbrefseq ~ "nt_unique",
                                   !maltdbnt & maltdbrefseq ~ "refseq_unique")
    ) %>%
    group_by(Combination) %>%
    count(Where_Found) %>%
    spread(Where_Found, n) %>%
    replace_na(list(both = 0, nt_unique = 0, refseq_unique = 0)) %>%
    print
}

## Raw taxa counts per software/taxa level/database/combination
per_combination_taxa_counts <- data_intersectiontables %>% 
  rename(MALT_Database = MALT_DB) %>%
  mutate(MALT_Database = gsub("maltdb", "", MALT_Database)) %>%
  group_by(Software, MALT_Database, Taxon_Level, Combination, Min_Support) %>% 
  filter(Min_Support == 0.05) %>%
  summarise(No_Taxa = n()) %>% 
  print(n = 99)

## Corroboration functions
software_comparison_genus <- compare_software(data_intersectiontables, 0.05, "genus")
software_comparison_species <- compare_software(data_intersectiontables, 0.05, "species")
database_comparison_genus <- compare_database(data_intersectiontables, "MALT_DB", "MALT", 0.05, "genus")
database_comparison_species <- compare_database(data_intersectiontables, "MALT_DB","MALT", 0.05, "species")

## Proportion corroborated for MALT NT 0.05 in MetaPhlan2, and remove taxa in MP2 not in MALT
per_combination_taxa_counts %>% 
  filter(Software == "MALT", MALT_Database == "nt", Taxon_Level == "genus") %>% 
  select(Combination, No_Taxa) %>% 
  rename(Overall_Taxa_MALT_Nt = No_Taxa) %>% 
  right_join(software_comparison_species %>% filter(MALT_Database == "nt")) %>% 
  mutate(Percent_MetaPhlAn2_also_in_MALT = Both / Overall_Taxa_MALT_Nt * 100) %>% 
  select(Combination, Both, MALT_only, MetaPhlAn2_only, Overall_Taxa_MALT_Nt, Percent_MetaPhlAn2_also_in_MALT) %>% 
  filter(!is.na(Overall_Taxa_MALT_Nt)) %>%
  print(n = 99)

per_combination_taxa_counts %>% 
  filter(Software == "MALT", MALT_Database == "nt", Taxon_Level == "species") %>% 
  select(Combination, No_Taxa) %>% 
  rename(Overall_Taxa_MALT_Nt = No_Taxa) %>% 
  right_join(software_comparison_species %>% filter(MALT_Database == "nt")) %>% 
  mutate(Percent_MetaPhlAn2_also_in_MALT = Both / Overall_Taxa_MALT_Nt * 100) %>% 
  select(Combination, Both, MALT_only, MetaPhlAn2_only, Overall_Taxa_MALT_Nt, Percent_MetaPhlAn2_also_in_MALT) %>% 
  filter(!is.na(Overall_Taxa_MALT_Nt)) %>%
  print(n = 99)

## As there is overlap when comparing either databases and/or software, lets
## investigate where they are going with an alluvial plot again

core_combinations <- c("Pan:Homo", "Gorilla:Pan:Homo", "Alouatta:Gorilla:Pan:Homo", "Not_Core")


data_intersectiontables %>% 
  filter(MALT_DB == "maltdbnt",
         Min_Support == "0.05",
         Tax_Level == "genus") %>%
  filter(Combination %in% core_combinations) %>%
  spread(Software, Combination) %>% 
  replace_na(list(MALT = "Not_Core", MetaPhlAn2 = "Not_Core")) %>%
  mutate(MALT = factor(MALT, levels = core_combinations),
         MetaPhlAn2 = factor(MetaPhlAn2, levels = core_combinations)) %>%
  ggplot(aes(axis1 = MALT, axis2 = MetaPhlAn2), y = "Combination") +
  geom_alluvium(color = "darkgray") +
  geom_stratum() +
  geom_text(stat = "stratum", label.strata = TRUE, size = 2) +
  ylab("Number of Taxa") +
  xlab("Software") +
  scale_x_continuous(breaks = 1:2, 
                     labels = c("MALT", "MetaPhlAn2"), 
                     position = "top") +
  theme_minimal(base_family = "Roboto", base_size = 7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

## What about overlaps between MALT databasE?




compare_maltdb_core <- function(min_tax, tax_level){
  my_table <- data_intersectiontables %>% 
    filter(Software == "MALT",
           Min_Support == min_tax,
           Tax_Level == tax_level) %>% 
    filter(Combination %in% core_combinations) %>%
    mutate(MALT_DB = gsub("maltdb", "", MALT_DB)) %>%
    spread(MALT_DB, Combination) %>% 
    replace_na(list(nt = "Not_Core", refseq = "Not_Core")) %>%
    mutate(nt = factor(nt, levels = core_combinations),
           refseq = factor(refseq, levels = core_combinations)) %>%
    arrange(nt, refseq)
    
    
  my_plot <- ggplot(my_table, aes(axis1 = nt, axis2 = refseq), y = "Combination") +
    geom_alluvium(color = "darkgray") +
    geom_stratum() +
    geom_text(stat = "stratum", label.strata = TRUE, size = 2) +
    ylab("Number of Taxa") +
    xlab("Software") +
    scale_x_continuous(breaks = 1:2, 
                       labels = c("nt", "refseq"), 
                       position = "top") +
    theme_minimal(base_family = "Roboto", base_size = 7) +
    theme(axis.ticks.y = element_line(),
          panel.grid.major = element_blank(),
          legend.position = "right"
    )
  
  return(list("plot" = my_plot, "table" = my_table))
}

## Note - species doesn't work as MALT for refseq doesn't get any of the 
## combinations.
dat <- compare_maltdb_core("0.05", "genus")

final_dbcompare_plot <- grid.arrange(dat$plot, tableGrob(dat$table %>% select(Taxon, nt, refseq), rows = NULL), nrow = 1)
  
ggsave(paste("13-coremicrobiome_presenceabsence_databasecomparisonalluvial_malt_0.05_genus_", format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), 
       final_dbcompare_plot, 
       path = "../04-analysis/screening/presenceabsence_intersection/", 
       device = cairo_pdf, 
       width = 7, 
       height = 3.5, 
       units = "in", 
       dpi = 600)
