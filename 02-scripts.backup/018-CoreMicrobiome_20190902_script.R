#' ---
#' title: "Core Microbiome Estimation"
#' output: html_notebook
#' editor_options: 
#'  chunk_output_type: inline
#' ---
#' 
#' # Preparation
#' 
#' ## For script conversion
#' 
#' As this notebook became more and more parameterised, I decided to make it
#' easy to convert into a script. But for the script, we also need to define
#' input arguments for the options in chunk 1.
#' 
#' The options in this chunk are our default selection. 
#' 

#' 
#' Reasons for these settings:
#' 
#' * Min support is quite strict 
#'  * to ensure only taxon with good read abundance selected and remove spurious hits
#' * Only half of an individuals in a population needs to have a taxon to be considered core
#'  * We need to consider that different individuals will have different biofilm grow stages, and this accounts for this
#' * 2/3s of the populations of a host genus needs to have a taxon as core to be a host genus taxon
#'  * This allows us to be sure it is found across a reasonable number of populations that is unlikely to be related from local contamination.
#'  * For Gorillas and Chimps this means it has to be found across at least differnt 2 populations at minimum
#'  * For humans this means 4/6 of the populations which means it either has to be in both ancient/modern samples or _all_ ancient populations that passed preservation thresholds.
#' 
#' Patterns we saw on each threshold for why the above
#'  * Increasing Min. Support to 0.05 meant highly likely environmental contaminants ( _Pseudomonas_ ) were not considered as core (often dropped out suggesting low level contamination of suprious hits)
#'  * Increasing minimum number of individuals from 0.5 to 0.6 saw very different number of changes in number of genera, and saw Mycobacterium being lost from the 'control' core, despite being a likely environmental contaminant. 0.8 saw quite a few clear oral genera being lost (e.g. _Porphyromonas_) 
#'  * Changing the fraction of populations per host genus saw changing from 0.5 - 0.66 single genus loss of 'Tessaracoccus' as anthropoid core (but most isolates from that genus is from activated sludge), increasing to 0.75 saw a loss of a lot of clear highly abundant oral taxa (e.g. Capnocytophaga, Ottowia, Prevotella, Filifactor etc.)
#' 
#' Script settings.
#' 
## ------------------------------------------------------------------------

 args = commandArgs(trailingOnly = TRUE)
 
 if (length(args) == 0 | args[1] == "" | args[1] == "-h" | args[1] == "--help") {
  cat("Usage: 018-CoreMicrobiome_20XXXXXX_script.R <db> <minsupp_multiplier> <fraction_individuals> <fraction_populations> <drop_single> \n")
  cat("db: nt, refseq in quotes!\n")
  cat("minsupp_multiplier: an integer to multiple 0.01% by \n")
  cat("fraction_individuals: e.g. 0.8 for 80% \n")
  cat("fraction_populations: e.g. 0.66 for 66% \n")
  cat("drop_single: T or F to drop populations with single individuals")
  cat("save_upsetr: T or F to save upsetR plot? At some point the package broke and doesn't work with ggsave and pdf devices do wierd stuff. Option to not print but still generate all other otutput")
  stop()
 } else if (length(args) == 6) {
  db <- args[1]
  minsupp_multiplier <- as.numeric(args[2])
  fraction_individuals <- as.numeric(args[3])
  fraction_populations = as.numeric(args[4])
  drop_single = args[5]
  save_upsetr = args[6]
  script <- T
 }

print(args)

minsupp_threshold <- 0.01 * minsupp_multiplier

#' 
#' 
#' ## Data Loading
#' 
#' Load data. Note: removed all MetaPhlan2 as no 'bad sample' removal method.
#' 
## ------------------------------------------------------------------------
 cat("loading data\n")

library(tidyverse)
library(VennDiagram)
library(UpSetR)
library(vegan) ## for some non base distance metrics
library(viridis) ## for colour 
library(scales) ## for colour
library(amap) ## alternative for correlation distance matrix
library(patchwork)

metadata <- read_tsv("../00-documentation.backup/02-calculus_microbiome-deep_evolution-individualscontrolssources_metadata_20190523.tsv") %>%
 rename(Individual = `#SampleID`)


## Bad Samples
if (db == "nt") {
  
  data_out_sampfil <- read_tsv("../04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninwithinfluctuationSDvariation_nt_fractionOralThreshold_50_20190509.tsv")
  
  bad_samples <- data_out_sampfil %>% 
      filter(!withinfluctuationvariation_pass) %>%
      filter(!grepl("ARS|SRR|ERR|LIB|EXB", sample)) %>%
      pull(sample)
} else if (db == "refseq") {
  
  data_out_sampfil <- read_tsv("../04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninwithinfluctuationSDvariation_refseq_fractionOralThreshold_65_20190509.tsv")
  
  bad_samples <- data_out_sampfil %>% 
      filter(!withinfluctuationvariation_pass) %>%
      filter(!grepl("ARS|SRR|ERR|LIB|EXB", sample)) %>%
      pull(sample)
}

## Contaminants
if (db == "nt") {
 contaminants_species <- read_tsv("../04-analysis/screening/decontam.backup/decontam_taxa_to_remove_megan_nt_species_combined_0.99_190411.tsv")
 contaminants_genus <- read_tsv("../04-analysis/screening/decontam.backup/decontam_taxa_to_remove_megan_nt_genus_combined_0.99_190411.tsv")
} else if (db == "refseq") {
 contaminants_species <- read_tsv("../04-analysis/screening/decontam.backup/decontam_taxa_to_remove_megan_refseq_species_combined_0.99_190411.tsv")
 contaminants_genus <- read_tsv("../04-analysis/screening/decontam.backup/decontam_taxa_to_remove_megan_refseq_genus_combined_0.99_190411.tsv")
}

## Load OTU tables

if (db == "nt") {
 raw_malt_genus <- read_tsv("../04-analysis/screening/megan.backup/Evolution-Comparison_MEGAN_20190401-ex_absolute_genus_prokaryotes_summarised_nt.txt")

## At Species Level
 raw_malt_species <- read_tsv("../04-analysis/screening/megan.backup/Evolution-Comparison_MEGAN_20190401-ex_absolute_species_prokaryotes_summarised_nt.txt")

 
} else if (db == "refseq") {
 raw_malt_genus <- read_tsv("../04-analysis/screening/megan.backup/Evolution-Comparison_MEGAN_20190410-ex_absolute_genus_prokaryotes_summarised_refseq.txt")
 
 ## At Species Level
 raw_malt_species <- read_tsv("../04-analysis/screening/megan.backup/Evolution-Comparison_MEGAN_20190410-ex_absolute_species_prokaryotes_summarised_refseq.txt")

}

## Set foundation minsupport columns
if (db == "nt") {
  minsupp_malt <- "Min_Support_Reads_Threshold_MALT"
} else if (db == "refseq") {
  minsupp_malt <- "Min_Support_Reads_Threshold_MALT_refseq"
}



#' 
#' ## Data Clean-Up
#' 
#' Clean up to remove column name cruft, and convert to long format.
#' 
## ------------------------------------------------------------------------
cat("cleaning data\n")

## Cleaner function
data_cleaner_new <- function(x, meta_dat, supp_col, ms_mltpr) {
  supp_col <- enquo(supp_col)
  ms_mltpr <- enquo(ms_mltpr)
  
  colnames(x) <- gsub("_S.*_L.*_R1_.*.fastq.combined.fq.prefixed.extractunmapped.bam", "", colnames(x))
  colnames(x) <- gsub("_S.*_L.*_R1_.*.fastq.merged.prefixed.hg19unmapped", "", colnames(x))
  colnames(x) <- gsub("_S0_L001_R1_001.fastq.extractunmapped.bam", "", colnames(x))
  colnames(x)[1] <- "Taxon"
  x <- x %>%
    gather(Individual, Value, 2:ncol(x)) %>%
    mutate(Value = as.numeric(Value)) %>%
    left_join(select(meta_dat, Individual, one_of(c(!! supp_col)), Age)) %>%
    rename(Orig_Threshold = !! supp_col) %>%
    mutate(Multiplier = as.numeric(!! ms_mltpr)) %>%
    mutate(Threshold = Orig_Threshold * Multiplier) %>%
    mutate(Threshold = as.numeric(Threshold)) %>%
    mutate(Filter_Passed = if_else(Value == 0,
      1,
      if_else(Value >= Threshold, 1, 0)
    )) %>%
    filter(Filter_Passed == 1) %>% # ,
    # Age != "ModernDay") %>%
    select(Taxon, Individual, Value)
  return(x)
}

## Add extra information and clean up the OTU table
clean_malt_genus <- data_cleaner_new(
  raw_malt_genus, metadata,
  minsupp_malt, minsupp_multiplier
) %>%
  mutate(Software = "MALT", Tax_Level = "genus")

clean_malt_species <- data_cleaner_new(
  raw_malt_species, metadata,
  minsupp_malt, minsupp_multiplier
) %>%
  mutate(Software = "MALT", Tax_Level = "species")

full_data <- bind_rows(
  clean_malt_genus, clean_malt_species
)


## Extra formatting mMetadata
metadata$Env <- factor(metadata$Env, levels = c(
  "Howler_Monkey",
  "Gorilla_1",
  "Gorilla_2",
  "Gorilla_3",
  "Chimp_1",
  "Chimp_2",
  "Chimp_3",
  "Chimp_4",
  "Neanderthal",
  "PreagriculturalHuman_1",
  "PreagriculturalHuman_2",
  "PreantibioticHuman_1",
  "PreantibioticHuman_2",
  "ModernDayHuman_1",
  "ModernDayHuman_2",
  "ExtractionControl",
  "LibraryControl",
  "ruralGut",
  "sediment",
  "skin",
  "subPlaque",
  "supPlaque",
  "urbanGut",
  "EnvironmentalControl"
))


metadata$Host_General <- factor(metadata$Host_General, levels = c(
  "Howler",
  "Gorilla",
  "Chimp",
  "Neanderthal",
  "PreagriculturalHuman",
  "PreantibioticHuman",
  "ModernDayHuman",
  "ExtractionControl",
  "LibraryControl",
  "ruralGut",
  "sediment",
  "skin",
  "subPlaque",
  "supPlaque",
  "urbanGut",
  "EnvironmentalControl"
))

metadata$Host_Common <- factor(metadata$Host_Common, levels = c(
  "Alouatta",
  "Gorilla",
  "Pan",
  "Homo (Neanderthal)",
  "Homo (Modern Human)",
  "ExtractionControl",
  "LibraryControl",
  "Plaque",
  "Gut",
  "Skin",
  "EnvironmentalControl",
  "Sediment"
))

metadata$Host_Genus <- factor(metadata$Host_Genus, levels = c(
  "Alouatta",
  "Gorilla",
  "Pan",
  "Homo",
  "Control",
  NA
))


#' 
#' 
#' ## Data Processing
#' 
#' Now we can do a presence/absence conversion. For the moment we can assume any
#' value that isn't 0 counts.
#' 
## ------------------------------------------------------------------------
 cat("processing data\n")

full_data <- full_data %>% mutate(Presence = if_else(Value > 0, 1, 0))

#' 
#' Next we need to add our metadata, so we can ask if a particular taxa is
#' present in all individuals of a particular host group. We can also remove
#' the controls and sources.
#' 
#' If we want to remove particular samples, we can also do this at this step.
#' 
## ------------------------------------------------------------------------
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

## Remove bad samples
full_data_meta <- full_data_meta %>% 
  filter(!Individual %in% bad_samples)

## Drop populations with only single individuals
envs_to_drop <- c()

if (drop_single) {
  envs_to_drop <- full_data_meta %>% 
    select(Individual, Env) %>% 
    distinct() %>% 
    group_by(Env) %>% 
    summarise(Individuals_N = n()) %>% 
    filter(Individuals_N <= 1) %>%
    pull(Env)
  
  full_data_meta <- full_data_meta %>% 
    filter(!Env %in% as.vector(envs_to_drop))
}





#' 
#' # Shannon Index
#' 
#' Alpha diversity considering number of species and evenness on the MALT/genus 
#' data
#' 
## ------------------------------------------------------------------------
 cat("shannon\n")
calculate_shannon_index <- function(data, met_data, soft, tax) {
 wide <- data %>% filter(Software == soft, Tax_Level == tax) %>%
  select(Taxon, Individual, Value) %>% 
  spread(Taxon, Value, fill = 0)
 matr <- as.matrix(wide[,2:ncol(wide)])
 rownames(matr) <- wide$Individual
 matr_tb <- tibble(diversity(matr), rownames = names(diversity(matr)))
 matr_tb <- matr_tb %>% rename(Individual = rownames, 
                Shannon_Index = `diversity(matr)`)
 final <- left_join(matr_tb, met_data)
 final <- final %>% rename(Population = Env)
 return(final)
}


my_colours <- c("#1f78b4", 
                 "#6a3d9a", 
                 "#6a3d9a", 
                 "#6a3d9a", 
                 "#33a02c", 
                 "#33a02c", 
                 "#33a02c", 
                 "#33a02c", 
                 "#e31a1c", 
                 "#ff7f00", 
                 "#ff7f00", 
                 "#ff7f00", 
                 "#ff7f00", 
                 "#ff7f00", 
                 "#ff7f00", 
                 "#d9d9d9", 
                 "#d9d9d9",
                 "#d9d9d9", 
                 "#d9d9d9", 
                 "#d9d9d9", 
                 "#d9d9d9", 
                 "#d9d9d9", 
                 "#d9d9d9",
                 "#d9d9d9")

names(my_colours) <- c("Howler_Monkey", 
           "Gorilla_1", 
           "Gorilla_2", 
           "Gorilla_3", 
           "Chimp_1", 
           "Chimp_2", 
           "Chimp_3", 
           "Chimp_4",
           "Neanderthal", 
           "PreagriculturalHuman_1", 
           "PreagriculturalHuman_2", 
           "PreantibioticHuman_1", 
           "PreantibioticHuman_2", 
           "ModernDayHuman_1", 
           "ModernDayHuman_2", 
           "ExtractionControl", 
           "LibraryControl", 
           "ruralGut", 
           "urbanGut",
           "sediment", 
           "skin", 
           "subPlaque",
           "supPlaque",
           "EnvironmentalControl"
           )

 my_shapes <- c(8,
                0,
                1,
                2,
                0,
                1,
                2,
                5,
                11,
                0,
                12,
                1,
                10,
                2,
                6,
                10,
                13,
                7,
                8,
                14,
                3,
                4,
                12,
                0)

names(my_shapes) <- c("Howler_Monkey", 
           "Gorilla_1", 
           "Gorilla_2", 
           "Gorilla_3", 
           "Chimp_1", 
           "Chimp_2", 
           "Chimp_3", 
           "Chimp_4",
           "Neanderthal", 
           "PreagriculturalHuman_1", 
           "PreagriculturalHuman_2", 
           "PreantibioticHuman_1", 
           "PreantibioticHuman_2", 
           "ModernDayHuman_1", 
           "ModernDayHuman_2", 
           "ExtractionControl", 
           "LibraryControl", 
           "ruralGut", 
           "sediment", 
           "skin", 
           "subPlaque",
           "supPlaque", 
           "urbanGut",
           "EnvironmentalControl")

shannon_data_malt_genus <- calculate_shannon_index(full_data_meta, metadata, 
                                                   "MALT", "genus")
shannon_data_malt_species <- calculate_shannon_index(full_data_meta, metadata, 
                                                     "MALT", "species")

plot_shannon <- function(x) {
 plot_one <- x %>% filter(Sample_or_Control == "Sample") %>% 
    ggplot(aes(Host_Genus, Shannon_Index)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(width = 0.1, size = 3, aes(colour = Population, 
                                         shape = Population)) +
  scale_colour_manual(values = my_colours) +
  scale_shape_manual(values = my_shapes) +
  xlab("Host Group") +
  ylab("Shannon's Diversity Index") +
  ggtitle("With Modern Day Humans") +
  theme_minimal(base_size = 7, base_family = "Roboto") +
  ggplot(x %>% filter(Age != "ModernDay"), aes(Host_Genus, 
                                               Shannon_Index)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(width = 0.1, size = 3, aes(colour = Population, 
                                         shape = Population)) +
  xlab("Host Group") +
  ylab("") +
  ggtitle("Without Modern Day Humans") +
  scale_colour_manual(values = my_colours) +
  scale_shape_manual(values = my_shapes) +
  theme_minimal(base_size = 7, base_family = "Roboto")
}

shannon_plot_malt_genus <- plot_shannon(shannon_data_malt_genus)
shannon_plot_malt_species <- plot_shannon(shannon_data_malt_species)

if (script == F) {
 shannon_plot_malt_genus
 shannon_plot_malt_species
}


ggsave(paste("01a-coremicrobiome_presenceabsence_shannondiversity_analysis_malt_",
  db,
  "_",
  minsupp_threshold,
  "_fracinds",
  fraction_individuals,
  "_fracpops",
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
  "_genus_",
  format(Sys.Date(), "%Y%m%d"),
  ".pdf",
  sep = ""
),
plot = shannon_plot_malt_genus,
"../04-analysis/screening/presenceabsence_intersection.backup/",
device = cairo_pdf,
width = 7,
height = 3.5,
units = "in",
dpi = 600
)

ggsave(paste("01b-coremicrobiome_presenceabsence_shannondiversity_analysis_malt_", 
 db,
 "_",
             minsupp_threshold,
             "_fracinds",
             fraction_individuals,
             "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
             "_species_", 
             format(Sys.Date(), "%Y%m%d"),
".pdf", 
sep = ""),
plot = shannon_plot_malt_species,
       "../04-analysis/screening/presenceabsence_intersection.backup/", 
device = cairo_pdf, 
       width = 7, 
       height = 3.5, 
       units = "in", 
       dpi = 600)




#' 
#' "4. The Shannon index
#' increases as both the richness and the evenness of the community increase. The
#' fact that the index incorporates both components of biodiversity can be seen as
#' both a strength and a weakness. It is a strength because it provides a simple,
#' synthetic summary, but it is a weakness because it makes it difficult to compare
#' communities that differ greatly in richness"
#' 
#' http://biology.kenyon.edu/courses/biol229/diversity.pdf.
#' 
#' So in this case the humans are possible less rich and less even, suggesting
#' they have more dominant handful of taxa than Gorillas/Pan/Howlers. 
#' Interestingly, this is not driven by just the Modern individuals, as ancient
#' samples are also falling around the mean.
#' 
#' # Equitability
#' 
#' To test for Equibility
#' 
## ------------------------------------------------------------------------
cat("equitability\n")

## From https://cran.r-project.org/web/packages/vegan/vignettes/diversity-vegan.pdf

calculate_equitability <- function(data, met_data, soft, tax) {
 wide <- data %>% filter(Software == soft, Tax_Level == tax) %>%
  select(Taxon, Individual, Value) %>% 
  spread(Taxon, Value, fill = 0)
 matr <- as.matrix(wide[,2:ncol(wide)])
 rownames(matr) <- wide$Individual
 equitability <- diversity(matr)/log(specnumber(matr))
 temp_tib <- tibble(equitability, rownames = names(equitability))
 temp_tib <- temp_tib %>% rename(Individual = rownames, Equitability = `equitability`)
 temp_tib <- left_join(temp_tib, met_data)
 result <- rename(temp_tib, Population = Env)
 return(result)
}

equitability_data_malt_genus <- calculate_equitability(full_data_meta, 
                                                       metadata, 
                                                       "MALT", 
                                                       "genus")
equitability_data_malt_species <- calculate_equitability(full_data_meta, 
                                                         metadata, 
                                                         "MALT", 
                                                         "species")


plot_equitability <- function(x) {
 ggplot(x %>% filter(Sample_or_Control == "Sample"), aes(Host_Genus, Equitability)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(width = 0.1, size = 2, aes(colour = Population, 
                                         shape = Population)) +
  scale_colour_manual(values = my_colours) +
  scale_shape_manual(values = my_shapes) +
  ggtitle("Without Modern Day Humans") +
  xlab("Host Genus") +
  theme_minimal(base_size = 7, base_family = "Roboto") +
  ggplot(x %>% filter(Age != "ModernDay"), aes(Host_Genus, Equitability)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(width = 0.1, size = 2, aes(colour = Population, shape = 
                                           Population)) +
  ylab("") +
  xlab("Host Genus") +
  ggtitle("Without Modern Day Humans") +
  scale_colour_manual(values = my_colours) +
  scale_shape_manual(values = my_shapes) +
  theme_minimal(base_size = 7, base_family = "Roboto")
}

equitability_plot_malt_genus <- plot_equitability(equitability_data_malt_genus)
equitability_plot_malt_species <- plot_equitability(equitability_data_malt_species)

if (script == F) {
 equitability_plot_malt_genus
 equitability_plot_malt_species
}

ggsave(paste("02a-coremicrobiome_presenceabsence_equitability_analysis_malt", 
 db,
 "_", 
  minsupp_threshold,
 "_fracinds",
 fraction_individuals,
 "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
 "_genus_", format(Sys.Date(), "%Y%m%d"),
".pdf", 
sep = ""), 
plot = equitability_plot_malt_genus, "../04-analysis/screening/presenceabsence_intersection.backup/", 
device = cairo_pdf, width = 7, 
height = 3.5, 
units = "in", 
dpi = 600)

ggsave(paste("02b-coremicrobiome_presenceabsence_equitability_analysis_malt_",
 db,
 "_", 
minsupp_threshold,
 "_fracinds",
 fraction_individuals,
 "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
 "_species_", format(Sys.Date(), "%Y%m%d"),
".pdf", 
sep = ""), 
plot = equitability_plot_malt_species, "../04-analysis/screening/presenceabsence_intersection.backup/", 
device = cairo_pdf, width = 7, 
height = 3.5, 
units = "in", 
dpi = 600)


#' 
#' Where the higher evenness (here implemented as Pielou's evennes index) 
#' is closer to 1. The less evenness, the closer to 0.
#' 
#' However, this suggests that the evenness is actually slightly higher in modern
#' humans, albeit similar to Gorillas (Pan is reduced) - so richness is main
#' driver of the Shannon scores for humans, while Pan could be more due to having
#' low evenness.
#' 
#' It is worth noting it seems that ModernDay humans seem to have more evenness
#' than of ancient humans who have distinctly less.
#' 
#' http://www.tiem.utk.edu/~gross/bioed/bealsmodules/shannonDI.html
#' 
#' # Individual Taxon Abundance Distributions
#' 
#' The purpose of this plot is akin to the equitibility plot above, as it is
#' asking for each taxa across the first 50 most abundant taxa, what is the 
#' proportion of alignments to each taxa. If the distribution is more uneven, you 
#' should see a higher proportion of the most abundant taxa (so taller peak) then
#' very low abundant tails, whereas more even groups should see smaller peak at
#' the most abundant and higher proportion tails. Note this odes
#' 
## ------------------------------------------------------------------------
cat("taxon abundance distribution\n")

my_colours_abundist <- c("#F7756C", 
        "#A2A500", 
        "#00BF7D", 
        "#00B0F6", 
        "#E76BF2", 
        "#6a3d9a", 
        "#ff7f00", 
        "#a6cee3", 
        "#b2df8a", 
        "#fb9a99", 
        "#fdbf6f",
        "#66c2a5")

names(my_colours_abundist) <- c("Alouatta",
                 "Gorilla",
                "Pan",
                "Homo (Neanderthal)",
                "Homo (Modern Human)",
                "ExtractionControl",
                "LibraryControl",
                "Plaque",
                "Gut",
                "Skin", 
                "Sediment",
                "EnvironmentalControl")

sorted_data_meta <- full_data_meta %>% 
 filter(Value != 0) %>%
 arrange(Env, desc(Value))

sample_totals <- sorted_data_meta %>%
 group_by(Individual, Software, Tax_Level) %>%
 summarise(Total_Alignments = sum(Value))

sorted_data_meta <- left_join(sorted_data_meta, sample_totals) %>%
 mutate(Proportion_Reads = Value / Total_Alignments * 100)

distribution_plotter <- function(in_dat, soft, tax_level) {
 trans_dat <- filter(in_dat, Software == soft, Tax_Level == tax_level, Value != 0) %>%
  group_by(Individual, Value) %>%
  arrange(Host_Common, Individual, desc(Value)) %>%
  left_join(summarise(group_by(., Individual), 
            Total_Alignments = sum(Value))) %>%
  ungroup() %>%
  mutate(Proportion_Alignments = Value / Total_Alignments * 100) %>%
  select(Host_Common, Individual, Proportion_Alignments) %>%
  group_by(Individual) %>%
  mutate(rowid = row_number()) %>%
  filter(rowid <= 50)

 ggplot(trans_dat, aes(rowid, Proportion_Alignments, group = Individual, colour = Host_Common)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ Host_Common) +
  scale_colour_manual(values = my_colours_abundist) +
  theme_minimal(base_family = "Roboto", base_size = 7) +
  ylab("Most to Least Abundant Taxon") +
  theme(axis.text.x = element_blank())
}

distribution_plot_malt_genus <- distribution_plotter(full_data_meta, "MALT", "genus")
distribution_plot_malt_species <- distribution_plotter(full_data_meta, "MALT", "species")


if (script == F) {
 distribution_plot_malt_genus
 distribution_plot_malt_species
}

ggsave(paste("03a-coremicrobiome_presenceabsence_taxonabundancedistribution_analysis_malt_",
   db,
   "_", 
  minsupp_threshold,
   "_fracinds",
   fraction_individuals,
   "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
   "_genus_", format(Sys.Date(), "%Y%m%d"),
  ".pdf", 
  sep = ""), 
plot = distribution_plot_malt_genus, "../04-analysis/screening/presenceabsence_intersection.backup/", 
device = cairo_pdf, 
width = 3.5, 
height = 3.5, 
units = "in", 
dpi = 600)

ggsave(paste("03b-coremicrobiome_presenceabsence_taxonabundancedistribution_analysis_malt_",
   db,
   "_", 
  minsupp_threshold,
   "_fracinds",
   fraction_individuals,
   "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
   "_species_", format(Sys.Date(), "%Y%m%d"),
  ".pdf", 
  sep = ""), 
plot = distribution_plot_malt_species, "../04-analysis/screening/presenceabsence_intersection.backup/", 
device = cairo_pdf, 
width = 3.5, 
height = 3.5, 
units = "in", 
dpi = 600)



#' 
#' # Presence Absence Calculations
#' 
#' To ensure we are finding taxa that are indeed core to each group, yet 
#' accounting for low peservation or CoDa-related drop out, we will require a 
#' taxon to be present in a majority of individuals across each of the populations 
#' of each host.
#' 
#' To identify this threshold, we first need to see how many individuals we have
#' in each host population,
#' 
## ------------------------------------------------------------------------
cat("presence absence inds threshold calculus\n")

full_data_meta %>% 
 select(Individual, Env) %>% 
 distinct() %>% 
 group_by(Env) %>% 
 summarise(No_Individuals = n()) %>%
 mutate(Threshold = fraction_individuals,
 Individuals_Passing_Threshold = No_Individuals * fraction_individuals) %>%
 write_tsv(paste("../04-analysis/screening/presenceabsence_intersection.backup/04-coremicrobiome_samplenumberthresholdcalculation_hostpopulations_maltdb", db,"_allthresholdlevels_alltaxalevels_fracinds", fraction_individuals,
 "_fracpops", 
 fraction_populations,
 "_singleindpopsdropped",
  drop_single,
 "_", 
 format(Sys.Date(), "%Y%m%d"), ".tsv", 
sep = ""))

if (script == F) {
 full_data_meta %>% 
  select(Individual, Env) %>% 
  distinct() %>% 
  group_by(Env) %>% 
  summarise(No_Individuals = n()) %>%
  mutate(Threshold = fraction_individuals,
 Individuals_Passing_Threshold = No_Individuals * fraction_individuals)
 
 full_data_meta %>% 
  select(Individual, Host_Common) %>% 
  distinct() %>% 
  group_by(Host_Common) %>% 
  summarise(No_Individuals = n()) %>%
  mutate(Passing_Proportion = No_Individuals * fraction_individuals)
 
 full_data_meta %>% 
  select(Individual, Host_Genus) %>% 
  distinct() %>% 
  group_by(Host_Genus) %>% 
  summarise(No_Individuals = n()) %>%
  mutate(Passing_Proportion = No_Individuals * fraction_individuals)
}


#' 
#' Now we can count how many times a particular taxa was found in a particular 
#' population and make a proportion.
#' 
## ------------------------------------------------------------------------
cat("presence absence taxa passing threshold\n")

## Number individuals a taxa is in per population
summary_data <- full_data_meta %>% 
 filter(Presence == 1) %>%
 group_by(Software, Tax_Level, Env, Host_Common, Host_Genus, Taxon) %>%
 summarise(No_Individuals_Present = n())

## Number of individuals in each population overall
summary_meta <- full_data_meta %>% select(Individual, Env, Host_Genus) %>% 
 distinct() %>%
 group_by(Env) %>%
 summarise(No_Individuals_Pop = n())

summary_data <- left_join(summary_data, summary_meta)

## Proportion of individuals in each population a taxa is present
summary_data <- summary_data %>% 
 mutate(Presence_Proportion_Per_Pop = round(No_Individuals_Present / No_Individuals_Pop * 100, digits = 1))

#' 
#' If we set our threshold as a taxon needs to be in 80% of individuals of a 
#' particular population, we can ask how many taxa are 'core' to each population.
#' 
## ------------------------------------------------------------------------
cat("Summarise taxa core to population\n")

if (script == F) {
 summary_data %>% 
  filter(Presence_Proportion_Per_Pop >= (fraction_individuals * 100)) %>%
  group_by(Software, Tax_Level, Env, Host_Common, Host_Genus) %>%
  summarise(No_Taxa_Per_Pop = n())
}

summary_data %>% 
 filter(Presence_Proportion_Per_Pop >= (fraction_individuals * 100)) %>%
 group_by(Software, Tax_Level, Env, Host_Common, Host_Genus) %>%
 summarise(No_Taxa_Per_Pop = n()) %>%
 write_tsv(paste("../04-analysis/screening/presenceabsence_intersection.backup/05-coremicrobiome_samplenumberthresholdpassed_hostpopulation_maltdb",
 db,
 "_",
 minsupp_threshold,
 "_alltaxalevels_fracinds", fraction_individuals,
 "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
 "_", 
 format(Sys.Date(), "%Y%m%d"), ".tsv", 
sep = ""))

majority_population_taxa <- summary_data %>% 
 filter(Presence_Proportion_Per_Pop >= (fraction_individuals * 100))

#' 
#' Finally, we can calculate the proportion of _populations_ of a host group, that 
#' haveparticular taxa in 80% of the groups individuals. 
#' 
## ------------------------------------------------------------------------
cat("Populations passing host genus threshold\n")

host_genus_count <- full_data_meta %>% select(Env, Host_Genus) %>% 
 distinct() %>%
 group_by(Host_Genus) %>%
 summarise(No_Pops_Group = n())

group_pop_prop_data_hostgenus <- summary_data %>% 
 filter(Presence_Proportion_Per_Pop >= (fraction_individuals * 100)) %>%
 group_by(Taxon, Software, Tax_Level, Host_Genus) %>%
 summarise(No_Pops_Present = n()) %>%
 left_join(host_genus_count) %>%
 mutate(Proportion_Pop_Groups_Presence = round(No_Pops_Present / No_Pops_Group * 100, digits = 1))


#' 
#' However,to look for taxa that are common to across multiple populations of each 
#' Host genus, i.e. not just to specific populations, we can filter any taxa that 
#' is not in 66% (i.e. 2/3) of the populations of a host group, and optionally 
#' remove any taxa that are considered lab contaminants. Two thirds are here 
#' selected as most of the Host Genera have a minimum of three subpopulations. Thus
#' 2/3rds will allow us to select taxa which are indeed found in more than one
#' population and not unique to that specific population.
#' 
## ------------------------------------------------------------------------
cat("Final core in each host genus calculation\n")



final_taxa_hostgenus <- filter(group_pop_prop_data_hostgenus, Proportion_Pop_Groups_Presence >= fraction_populations * 100)

cat("gap\n")

if (script == F) {
 final_taxa_hostgenus %>% 
  filter(!Taxon %in% pull(contaminants_genus)) %>% 
  filter(!Taxon %in% pull(contaminants_species))
}

cat("Contaminant removal and save\n")

final_taxa_hostgenus <- final_taxa_hostgenus %>% 
 filter(!Taxon %in% pull(contaminants_genus)) %>% 
 filter(!Taxon %in% pull(contaminants_species)) %>%
 write_tsv(paste("../04-analysis/screening/presenceabsence_intersection.backup/06-coremicrobiome_microbialtaxapassingthreshold_hostgenus_allsoftware_maltdb", db ,"_", minsupp_threshold,
 "_fracinds",
 fraction_individuals,
 "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
 "_", 
 format(Sys.Date(), "%Y%m%d"), ".tsv", 
sep = ""))


#' 
#' With this we can make the data ready for intersection analysis via UpSetR. 
#' To do this, for each particular software and taxonomy level, we pull the taxa 
#' defined for each group that matches the thresholds defined above of a taxa
#' being present in at least 80% of indivduals of a population, and 2/3s of 
#' the populations of a host genus.
#' 
#' In addition, for finer scale analysis, we will look for taxa found at just 80%
#' of each population only and remove unique population taxa. This allows us to 
#' look at which taxa are shared from which population - but also allows us
#' to find taxa in Neanderthals - which are a particular interest of this project.
#' 
#' 
## ------------------------------------------------------------------------
cat("Intersections between each host genus\n")

## For Host Genus UpSetR with 80% Individuals per Prop, and 66% of populations
venn_prep_hostgenus <- function(dat, soft, tax) {
 venn_data <- list(
  Alouatta = dat %>%
   filter(Software == soft, Tax_Level == tax) %>%
   filter(Host_Genus == "Alouatta") %>%
   pull(Taxon),
  Gorilla = dat %>% 
   filter(Software == soft, Tax_Level == tax) %>%
   filter(Host_Genus == "Gorilla") %>% 
   pull(Taxon),
  Pan = dat %>% 
   filter(Software == soft, Tax_Level == tax) %>%
   filter(Host_Genus == "Pan") %>% 
   pull(Taxon),
  Homo = dat %>% 
   filter(Software == soft, Tax_Level == tax) %>%
   filter(Host_Genus == "Homo") %>% 
   pull(Taxon),
  Control = dat %>% 
   filter(Software == soft, Tax_Level == tax) %>%
   filter(Host_Genus == "Control") %>% 
   pull(Taxon)
 )
 return(venn_data)
}

venn_malt_species_hostgenus <- venn_prep_hostgenus(final_taxa_hostgenus, "MALT", "species")
venn_malt_genus_hostgenus <- venn_prep_hostgenus(final_taxa_hostgenus, "MALT", "genus")


## For population only
venn_prep_env <- function(data, soft, tax) {
 result <- list()
 for (pop in data$Env %>% unique) {
   res <- list(filter(data, Env == pop, Software == soft, Tax_Level == tax)$Taxon)
   result[paste(pop)] <- res
 }
 return(result)
}

venn_malt_species_hostenv <- venn_prep_env(majority_population_taxa, "MALT", "species")
venn_malt_genus_hostenv <- venn_prep_env(majority_population_taxa, "MALT", "genus")
 


#' 
#' # Plotting
#' 
#' ## UpSet Plot
#' 
#' Alternatively from a VennDiagram we could generate an UpSet plot.
#' 
## ------------------------------------------------------------------------
cat("Upset\n")

#upset(fromList(venn_data))
 malt_genus_hostgenus <- upset(fromList(venn_malt_genus_hostgenus), nsets = length(venn_malt_genus_hostgenus), nintersects = NA, keep.order = T)
 malt_species_hostgenus <- upset(fromList(venn_malt_species_hostgenus), nsets = length(venn_malt_species_hostgenus), nintersects = NA, keep.order = T)
 
 malt_genus_hostenv <- upset(fromList(venn_malt_genus_hostenv), nsets = length(venn_malt_genus_hostenv), nintersects = NA, keep.order = T)
 malt_species_hostenv <- upset(fromList(venn_malt_species_hostenv), nsets = length(venn_malt_species_hostenv), nintersects = NA, keep.order = T)



cairo_pdf(filename = paste("../04-analysis/screening/presenceabsence_intersection.backup/07a-coremicrobiome_presenceabsence_upsetplot_malt_",
     db,
     "_", 
    minsupp_threshold,
     "_fracinds",
     fraction_individuals,
     "_fracpops", 
    fraction_populations,
    "_singleindpopsdropped",
    drop_single,
     "_hostgenus_genus_", format(Sys.Date(), "%Y%m%d"),
    ".pdf", 
    sep = ""),
  width = 3.5, 
  height = 3.5)
malt_genus_hostgenus
dev.off()

cairo_pdf(filename = paste("../04-analysis/screening/presenceabsence_intersection.backup/07a-coremicrobiome_presenceabsence_upsetplot_malt_",
     db,
     "_", 
    minsupp_threshold,
     "_fracinds",
     fraction_individuals,
     "_fracpops", 
    fraction_populations,
    "_singleindpopsdropped",
    drop_single,
     "_hostgenus_species_", format(Sys.Date(), "%Y%m%d"),
    ".pdf", 
    sep = ""),
  width = 3.5, 
  height = 3.5)
malt_species_hostgenus
dev.off()

cairo_pdf(filename = paste("../04-analysis/screening/presenceabsence_intersection.backup/07a-coremicrobiome_presenceabsence_upsetplot_malt_",
     db,
     "_", 
    minsupp_threshold,
     "_fracinds",
     fraction_individuals,
     "_fracpops", 
    fraction_populations,
    "_singleindpopsdropped",
    drop_single,
     "_hostenv_genus_", format(Sys.Date(), "%Y%m%d"),
    ".pdf", 
    sep = ""),
  width = 14, 
  height = 3.5)
malt_genus_hostenv
dev.off()

cairo_pdf(filename = paste("../04-analysis/screening/presenceabsence_intersection.backup/07a-coremicrobiome_presenceabsence_upsetplot_malt_",
     db,
     "_", 
    minsupp_threshold,
     "_fracinds",
     fraction_individuals,
     "_fracpops", 
    fraction_populations,
    "_singleindpopsdropped",
    drop_single,
     "_hostenv_species_", format(Sys.Date(), "%Y%m%d"),
    ".pdf", 
    sep = ""),
  width = 14, 
  height = 3.5)
malt_species_hostenv
dev.off()
 
# ## To save hostgenus
# if (isTRUE(save_upsetr)) {
#   ggsave(paste("07a-coremicrobiome_presenceabsence_upsetplot_malt_",
#      db,
#      "_", 
#     minsupp_threshold,
#      "_fracinds",
#      fraction_individuals,
#      "_fracpops", 
#     fraction_populations,
#     "_singleindpopsdropped",
#     drop_single,
#      "_hostgenus_genus_", format(Sys.Date(), "%Y%m%d"),
#     ".pdf", 
#     sep = ""), 
#   plot = upset(fromList(venn_malt_genus_hostgenus), nsets = length(venn_malt_genus_hostgenus), nintersects = NA, keep.order = T)
#   , "../04-analysis/screening/presenceabsence_intersection.backup/", 
#   device = cairo_pdf, 
#   width = 3.5, 
#   height = 3.5, 
#   units = "in", 
#   dpi = 600)
#   
#   ggsave(paste("07b-coremicrobiome_presenceabsence_upsetplot_malt",
#      db,
#      "_", 
#     minsupp_threshold,
#      "_fracinds",
#      fraction_individuals,
#      "_fracpops", 
#     fraction_populations,
#     "_singleindpopsdropped",
#     drop_single,
#      "_hostgenus_species_", format(Sys.Date(), "%Y%m%d"),
#     ".pdf", 
#     sep = ""), 
#   plot = upset(fromList(venn_malt_species_hostgenus), nsets = length(venn_malt_species_hostgenus), nintersects = NA, keep.order = T)
#   , "../04-analysis/screening/presenceabsence_intersection.backup/", 
#   device = cairo_pdf, 
#   width = 3.5, 
#   height = 3.5, 
#   units = "in", 
#   dpi = 600)
#   
#   
#   ## To save hostcommon
#   ggsave(paste("07c-coremicrobiome_presenceabsence_upsetplot_malt_",
#    db,
#    "_", 
#   minsupp_threshold,
#      "_fracinds",
#      fraction_individuals,
#      "_fracpops", 
#     fraction_populations,
#     "_singleindpopsdropped",
#     drop_single,
#      "_hostenv_genus_", format(Sys.Date(), "%Y%m%d"),
#     ".pdf", 
#     sep = ""), 
#   plot = upset(fromList(venn_malt_genus_hostenv), nsets = length(venn_malt_genus_hostenv), nintersects = NA, keep.order = T)
#   , "../04-analysis/screening/presenceabsence_intersection.backup/", 
#   device = cairo_pdf, 
#   width = 3.5, 
#   height = 3.5, 
#   units = "in", 
#   dpi = 600)
#   
#   ggsave(paste("07d-coremicrobiome_presenceabsence_upsetplot_malt",
#      db,
#      "_", 
#     minsupp_threshold,
#      "_fracinds",
#      fraction_individuals,
#      "_fracpops", 
#     fraction_populations,
#     "_singleindpopsdropped",
#     drop_single,
#      "_hostenv_species_", format(Sys.Date(), "%Y%m%d"),
#     ".pdf", 
#     sep = ""), 
#   plot = upset(fromList(venn_malt_species_hostenv), nsets = length(venn_malt_species_hostenv), nintersects = NA, keep.order = T)
#   , "../04-analysis/screening/presenceabsence_intersection.backup/", 
#   device = cairo_pdf, 
#   width = 3.5, 
#   height = 3.5, 
#   units = "in", 
#   dpi = 600)
# }

#' 
#' # Generation Per Grouping Taxa List 
#' 
#' Although we have the numbers of taxa that overlap per group, we still don't 
#' know what taxa these are.
#' 
#' We here extract the these by calling the indicies of the overlapped taxa
#' in each list to pull the names from the input list itself.
#' 
## ------------------------------------------------------------------------
 cat("intersection\n")

## Modified from https://github.com/hms-dbmi/UpSetR/issues/85#issuecomment-415480954
overlapGroups <- function(listInput, sort = TRUE) {
 listInputmat  <- fromList(listInput) == 1
 listInputunique <- unique(listInputmat)
 grouplist <- list()
 for (i in 1:nrow(listInputunique)) {
  currentRow <- listInputunique[i,]
  myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
  attr(myelements, "groups") <- currentRow
  grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
  myelements
 }
 if (sort) {
  grouplist <- grouplist[order(sapply(grouplist, 
                                      function(x) length(x)), decreasing = TRUE)]
 }
 attr(grouplist, "elements") <- unique(unlist(listInput))

 core_per_groups <- c()
 
 for (i in names(grouplist)) {
  tmp_list <- list()
  tmp_list[1] <- list(attributes(grouplist)$elements[grouplist[[i]]])
  names(tmp_list) <- paste(i)
  core_per_groups <- append(core_per_groups, tmp_list[1])
 }
 core_per_groups <- as_tibble(stack(core_per_groups)) %>%
  select(ind, values) %>%
  rename(Combination = ind, Taxon = values)
 return(core_per_groups)
}

overlaptaxa_malt_genus_hostgenus <- overlapGroups(venn_malt_genus_hostgenus) %>% 
  mutate(Software = "MALT", Tax_Level = "genus")
overlaptaxa_malt_species_hostgenus <- overlapGroups(venn_malt_species_hostgenus) %>% 
  mutate(Software = "MALT", Tax_Level = "species")


overlaptaxa_malt_genus_hostenv <- overlapGroups(venn_malt_genus_hostenv) %>% 
  mutate(Software = "MALT", Tax_Level = "genus")
overlaptaxa_malt_species_hostenv <- overlapGroups(venn_malt_species_hostenv) %>% 
  mutate(Software = "MALT", Tax_Level = "species")

core_taxa_list_hostgenus <- bind_rows(overlaptaxa_malt_genus_hostgenus, 
     overlaptaxa_malt_species_hostgenus) %>% 
 mutate(Taxon = gsub(" ", "_", 
 Taxon))

core_taxa_list_hostenv <- bind_rows(overlaptaxa_malt_genus_hostenv, 
     overlaptaxa_malt_species_hostenv) %>% 
 mutate(Taxon = gsub(" ", "_", 
 Taxon))


## To save for host genus
bind_rows(overlaptaxa_malt_genus_hostgenus) %>% 
 mutate(Taxon = gsub(" ", "_", 
 Taxon)) %>%
 write_tsv(paste("../04-analysis/screening/presenceabsence_intersection.backup/08a-coremicrobiome_presenceabsence_upsettable_allsoftware_maltdb",
 db,
 "_", minsupp_threshold,
 "_fracinds",
 fraction_individuals,
 "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
 "_hostgenus_genus_", format(Sys.Date(), "%Y%m%d"), ".tsv", 
sep = ""))

bind_rows(overlaptaxa_malt_species_hostgenus) %>% 
 mutate(Taxon = gsub(" ", "_", 
 Taxon)) %>%
  write_tsv(paste("../04-analysis/screening/presenceabsence_intersection.backup/08b-coremicrobiome_presenceabsence_upsettable_allsoftware_maltdb",
 db,
"_", minsupp_threshold,
 "_fracinds",
 fraction_individuals,
 "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
 "_hostgenus_species_", format(Sys.Date(), "%Y%m%d"), ".tsv", 
sep = ""))

## To save for host species
bind_rows(overlaptaxa_malt_genus_hostenv) %>% 
 mutate(Taxon = gsub(" ", "_", 
 Taxon)) %>%
 write_tsv(paste("../04-analysis/screening/presenceabsence_intersection.backup/08c-coremicrobiome_presenceabsence_upsettable_allsoftware_maltdb",
 db,
 "_", minsupp_threshold,
 "_fracinds",
 fraction_individuals,
 "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
 "_hostenv_genus_", format(Sys.Date(), "%Y%m%d"), ".tsv", 
sep = ""))

bind_rows(overlaptaxa_malt_species_hostenv) %>% 
 mutate(Taxon = gsub(" ", "_", 
 Taxon)) %>%
  write_tsv(paste("../04-analysis/screening/presenceabsence_intersection.backup/08d-coremicrobiome_presenceabsence_upsettable_allsoftware_maltdb",
 db,
"_", minsupp_threshold,
 "_fracinds",
 fraction_individuals,
 "_fracpops", 
  fraction_populations,
  "_singleindpopsdropped",
  drop_single,
 "_hostenv_species_", format(Sys.Date(), "%Y%m%d"), ".tsv", 
sep = ""))

#' 
#' And save this list
#' 
#' # Grouping Taxa Relative Abundances
#' 
#' Define pseudo count and CLR transform function for calculating ratio abundances
#' 
## ------------------------------------------------------------------------
cat("lCLR\n")

pse_clr_pca <- function(in_data, soft, tax) {
 filt_data <- filter(in_data, Software == soft, Tax_Level == tax)
 ## firstly convert to a matrix
 wide_data <- filt_data %>% 
  dplyr::select(Taxon, Individual, Value) %>% 
  spread(Taxon, Value, fill = 0)
 matrix <- as.matrix(wide_data[2:ncol(wide_data)])
 rownames(matrix) <- wide_data$Individual
 matrix.pse <- matrix[,2:ncol(matrix)] <- matrix[,2:ncol(matrix)] + 1
 matrix.clr <- t(apply(matrix.pse, 1, function(x) {log(x) - mean(log(x))}))
 return(as_tibble(matrix.clr, rownames = "Individual") %>% 
      gather(Taxon, Value, 2:ncol(.)))
}

#' 
#' Apply CLR transform to each software and taxonomic level data
#' 
## ------------------------------------------------------------------------
clr_malt_species <- pse_clr_pca(full_data_meta, "MALT", "species")
clr_malt_genus <- pse_clr_pca(full_data_meta, "MALT", "genus")


#' 
#' Some plotting parameters
#' 
## ------------------------------------------------------------------------
my_colours <- c("#1f78b4", 
        "#6a3d9a", 
        "#33a02c",
        "#ff7f00", 
        "#d9d9d9" 
        )

names(my_colours) <- c("Alouatta", 
           "Gorilla", 
           "Pan", 
           "Homo", 
           "Control"
           )

#' 
#' Violin plots for showing variation in relative ratio of taxa in samples
#' 
#' 
## ----eval = FALSE--------------------------------------------------------
## 
## clr_ratio_plotter <- function(in_dat, met_dat, tax_list) {
##  plot_list <- list()
##  for(i in names(tax_list)) {
##   plot_dat <- suppressMessages(left_join(in_dat, met_dat) %>%
##    filter(Taxon %in% unlist(tax_list[i]))) #%>%
##    #filter(Age != "ModernDay")
##   plot_dat$Taxon <- factor(unlist(tax_list[i]), levels = unlist(tax_list[i]))
##   my_plot <- ggplot(plot_dat, aes(Host_Genus, Value, fill = Host_Genus)) +
##    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
##    geom_jitter(size = 1, width = 0.2) +
##    scale_fill_manual(values = my_colours) +
##    xlab("Group") +
##    ylab("Centered Log Ratio") +
##    facet_wrap(~ Taxon) +
##    theme_minimal(base_family = "Roboto", base_size = 7) +
##    theme(axis.text.x = element_text(angle = 90, hjust = 1))
##   #plot_list[[i]] <- my_plot
##  }
##  return(my_plot)
## }
## 
## #clr_violins_malt_genus <- clr_ratio_plotter(clr_malt_genus, metadata, overlaptaxa_malt_genus)
## #clr_violins_malt_species <- clr_ratio_plotter(clr_malt_genus, metadata, overlaptaxa_malt_genus)
## 
## 
## ## Potential main figure example
## paper_list <- list(Example = factor((core_taxa_list %>%
##                     filter(Combination == "Alouatta:Gorilla:Pan:Homo",
##                        Software == "MALT",
##                        Tax_Level == "genus") %>%
##                     pull(Taxon))),
##                   levels = ((core_taxa_list %>%
##                          filter(Combination == "Alouatta:Gorilla:Pan:Homo",
##                             Software == "MALT",
##                             Tax_Level == "genus") %>%
##                          pull(Taxon)))
##                   )
## )
## 
## clr_violins_malt_genus_paper <- clr_ratio_plotter(clr_malt_genus, metadata, paper_list)
## 
## if (script == F) {
##  clr_violins_malt_genus_paper
## }

#' 
#' Misc Test to see if modern humans are artificially inflating the core
#' taxa relative abundances consistantly across all taxa, or this is indeed
#' a homo specific thing
#' 
## ----eval = F------------------------------------------------------------
## ## Get mean CLR of with today
## with_today <- left_join(clr_malt_genus, metadata) %>%
##    filter(Taxon %in% paper_list$Example) %>%
##  group_by(Host_Genus, Taxon) %>%
##  summarise(Mean_CLR = mean(Value)) %>%
##  mutate(With_ModernDay = "Y")
## 
## ## Get mean CLR of without today
## without_today <- left_join(clr_malt_genus, metadata) %>%
##    filter(Taxon %in% paper_list$Example, Age != "ModernDay") %>%
##  group_by(Host_Genus, Taxon) %>%
##  summarise(Mean_CLR = mean(Value)) %>%
##  mutate(With_ModernDay = "N")
## 
## ## Join Together
## 
## if (script == F) {
##  bind_rows(with_today, without_today) %>% spread(Host_Genus, Mean_CLR) %>%
##   arrange(Taxon)
## }

#' 
#' 
#' # Script generation
#' 
#' Use the following to generate a new fast script.
#' 
## ----eval = F------------------------------------------------------------
## cat("script gen\n")
## 
## knitr::purl("../02-scripts.backup/018-CoreMicrobiome_20190902.Rmd", "../02-scripts.backup/018-CoreMicrobiome_20190902_script.R", documentation = 2)

#' 
