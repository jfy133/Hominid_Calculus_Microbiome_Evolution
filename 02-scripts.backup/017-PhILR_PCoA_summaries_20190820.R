#!/user/bin/env Rscript

## Infrastructure - libs
### Libraries

library(tidyverse)

## Functions

main_dir <- "../04-analysis/screening/philr.backup/"

collect_files <- function(string) {
  list.files(main_dir, pattern = string, full.names = T, recursive = F) %>% 
    enframe(name = NULL, value = "files") %>%
    filter(grepl("badsamplesout", files) & 
             grepl("noSources", files) & 
             grepl("noControls", files))
}

collect_data <- function(file_list){
  as_tibble(file_list) %>%
    mutate(file_contents = purrr::map(files, ~ suppressWarnings(suppressMessages(as_tibble(read_tsv(.)))))) %>%
    rowwise()  %>%
    mutate(Database = str_split(files, "_")[[1]][[5]],
           Taxonomic_Level = str_split(files, "_")[[1]][[6]],
           Sources_Included = str_split(files, "_")[[1]][[7]],
           Controls_Included = str_split(files, "_")[[1]][[8]],
           LowPreservedSamples_Included = str_split(files, "_")[[1]][[9]],
           Min_Support = str_split(files, "_")[[1]][[11]]) %>%
    select(-files) %>%
    unnest()
}

## Load data

files_betadisperanova <- collect_files("07-betadisper_homogenity_anova_*")
files_betadisperpermutest <- collect_files("08-betadisper_homogenity_permutest_*")
files_permanovainitial <- collect_files("10a-permanova_initial_hostgenuscentroids_*")
files_permanovabootstrap <- collect_files("11a-permanova_bootstrap_hostgenuscentroids_*")

## Data wrangling

data_betadisperanova <- collect_data(files_betadisperanova) %>%
  filter(Taxonomic_Level == "genus" & Min_Support == "0.07" | 
           Taxonomic_Level == "species" & Min_Support == "0.04") %>%
  select(-Sources_Included, -Controls_Included, -LowPreservedSamples_Included) %>%
  mutate(p.value = if_else(!is.na(p.value), 
                           format(p.value, scientific = F), 
                           NA_character_)) %>% 
  mutate(p.value = as.double(p.value), 
         test = "anova") %>%
  rename(Test = test,
         Term = term, 
         `Degrees of Freedom` = df, 
         `Sum of Squares` = sumsq,
         `Mean of Squares` = meansq,
         `Test Statistic` = statistic, 
         `p value` = p.value
    ) %>%
  select(Database, Taxonomic_Level, Min_Support, Test, everything()) %>%
  print(n = 99) %>%
  split(.$Test)

data_betadisperpermutest <- collect_data(files_betadisperpermutest) %>%
  filter(Taxonomic_Level == "genus" & Min_Support == "0.07" | 
           Taxonomic_Level == "species" & Min_Support == "0.04") %>%
  rename(comb_dat  = `Permutation test for homogeneity of multivariate dispersions`) %>% 
  filter(grepl("Groups|Residuals", comb_dat)) %>% 
  mutate(comb_dat = gsub("\\s+", " ", comb_dat)) %>% 
  separate(comb_dat, 
           c("term", "df", "sumsq", "meansq", "statistic", "No. Permutations", "p.value", "sig"), 
           sep = " ",
           fill = "right") %>% 
  mutate(df = as.numeric(df),
         sumsq = as.numeric(sumsq),
         meansq = as.numeric(meansq),
         statistic = as.numeric(statistic),
         `No. Permutations` = as.numeric(`No. Permutations`),
         p.value = as.numeric(p.value), 
         test = "permutest"
         ) %>%
  select(-Sources_Included, -Controls_Included, -LowPreservedSamples_Included, -sig) %>%
  rename(Test = test,
         Term = term, 
         `No. Permutations (permutest)` = `No. Permutations`,
         `Degrees of Freedom` = df, 
         `Sum of Squares` = sumsq,
         `Mean of Squares` = meansq,
         `Test Statistic` = statistic, 
         `p value` = p.value
  ) %>%
  select(Database, Taxonomic_Level, Min_Support, Test, everything()) %>%
  print(n = 99) %>%
  split(.$Test)

## single permanova (not tabular input data so extra cleanup)
data_permanovasingle <- collect_data(files_permanovainitial) %>% 
  rename(comb_dat  = `Call:`) %>% 
  filter(Taxonomic_Level == "genus" & Min_Support == "0.07" | 
           Taxonomic_Level == "species" & Min_Support == "0.04") %>%
  filter(grepl("groups|Residuals|Total", comb_dat), !grepl("adonis", comb_dat)) %>% 
  mutate(comb_dat = gsub("\\s+", " ", comb_dat)) %>% 
  separate(comb_dat, 
           c("term", "df", "sumsq", "meansq", "F.Model", "R2", "Pr(>F)", "sig"), 
           sep = " ",
           fill = "right") %>% 
  mutate(df = as.numeric(df),
         sumsq = as.numeric(sumsq),
         meansq = as.numeric(meansq),
         `Pr(>F)` = as.numeric(`Pr(>F)`),
         R2 = as.numeric(R2),
         test = "adonis"
  ) %>% 
  rename(Test = test,
         Term = term, 
         `Degrees of Freedom` = df, 
         `Sum of Squares` = sumsq,
         `Mean of Squares` = meansq,
         `Test Statistic` = F.Model, 
         `R Squared` = R2,
         `p value` = `Pr(>F)`,
         ) %>%
  select(-Sources_Included, -Controls_Included, 
         -LowPreservedSamples_Included, -sig) %>%
  select(Database, Taxonomic_Level, Min_Support, Test, everything()) %>%
  print(n = 999) %>%
  split(.$Test)

summary_permanovabootstrap <- collect_data(files_permanovabootstrap) %>% 
  group_by(Database, 
           Taxonomic_Level, 
           Sources_Included, 
           Controls_Included, 
           LowPreservedSamples_Included, 
           Min_Support) %>%
  summarise_all(list(Mean = mean, SD = sd), na.rm = TRUE) %>%
  ungroup() %>%
  select(-contains("Included"), -contains("Run")) %>% 
  gather(Metric, Value, 4:ncol(.)) %>% 
  separate(Metric, c("Test", "Statistic", "Blub", "Summary_Metric"), 
           sep = "_") %>% 
  select(-Blub) %>%
  spread(Summary_Metric, Value) %>%
  print(n = 999) %>%
  split(.$Test)



write_tsv(data_betadisperanova$anova, paste("../00-documentation.backup/19a-philr_hostgenus_noSources_noControls_badsamplesOut_withinvariation_betadispersion_singlerun_summary_anova_", format(Sys.Date(), "%Y%m%d"), ".tsv", sep = ""))

write_tsv(data_betadisperpermutest$permutest, paste("../00-documentation.backup/19b-philr_hostgenus_noSources_noControls_badsamplesOut_withinvariation_betadispersion_singlerun_summary_permutest_", format(Sys.Date(), "%Y%m%d"), ".tsv", sep = ""))

write_tsv(data_permanovasingle$adonis, paste("../00-documentation.backup/21-philr_hostgenus_noSources_noControls_badsamplesOut_withinvariation_permanova_singlerun_summary_adonis_", format(Sys.Date(), "%Y%m%d"), ".tsv", sep = ""))

write_tsv(summary_permanovabootstrap$Adonis, paste("../00-documentation.backup/22a-philr_hostgenus_noSources_noControls_badsamplesOut_withinvariation_permanova_bootstrap_summary_adonis_", format(Sys.Date(), "%Y%m%d"), ".tsv", sep = ""))

write_tsv(summary_permanovabootstrap$BetaDispersion, paste("../00-documentation.backup/22a-philr_hostgenus_noSources_noControls_badsamplesOut_withinvariation_permanova_bootstrap_summary_betadisper_", format(Sys.Date(), "%Y%m%d"), ".tsv", sep = ""))

