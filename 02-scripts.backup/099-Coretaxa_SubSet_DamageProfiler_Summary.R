#! /usr/bin/env Rscript

######LIBRARIES########

library(tidyverse) ## for general data cleaning
library(data.table) ## for large table loading
library(taxize) ## for NCBI taxonomic info collection
library(ape) ## for tree manipulation
library(phyloseq) ## for data format as input for PhILR
library(philr) ## for data transform
library(vegan) ## for statistical testing
library(patchwork) ## for further visualisation assistance
library(cowplot) ## further visualisation assistance

######MAP########


data_meta_raw <- read_csv("00-documentation.backup/01-calculus_microbiome_deep_evolution_samplescontrols_metadata_20191112.csv")

data_meta_filtered <- filter(data_meta_raw, !grepl("Blank", Group)) %>%
  select(Group, Host_Common, Latitude, Longitude) %>% 
  unique

data_meta_filtered$Group <- factor(data_meta_filtered$Group, levels = c("Howler_Monkey", 
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
                                                                        "ExtractionBlank", 
                                                                        "LibraryBlank", 
                                                                        "ruralGut", 
                                                                        "sediment", 
                                                                        "skin", 
                                                                        "subPlaque",
                                                                        "supPlaque", 
                                                                        "urbanGut",
                                                                        "EnvironmentalControl")
)

common_colours <- c(Alouatta = "#1f78b4", Gorilla = "#6a3d9a", Pan = "#33a02c", 
                    `Homo (Neanderthal)` = "#fdbf6f", 
                    `Homo (Modern Human)` = "#ff7f00", ExtractionControl = "darkgrey", 
                    LibraryControl = "darkgrey", Plaque = "darkgrey", Gut = "darkgrey", 
                    Skin = "darkgrey", Sediment = "darkgrey", EnvironmentalControl = "darkgrey")

common_shapes <- c(Alouatta = 8, Gorilla = 0, Pan = 1, 
                   `Homo (Neanderthal)` = 2, `Homo (Modern Human)` = 6,
                   ExtractionControl = 3, LibraryControl = 7, Plaque = 9, 
                   Gut = 7, Skin = 9, Sediment = 10, EnvironmentalControl = 12)

## Metadata
raw_metadata <- read_tsv("00-documentation.backup/02-calculus_microbiome-deep_evolution-individualscontrolssources_metadata_20190523.tsv") %>%
  rename(Individual = `#SampleID`)

######DAMAGE########

data_longdamage <- read_tsv("00-documentation.backup/14-damageprofiler_screening_3p_5p_summaries_20190503.tsv")

selected_taxa_damage <- c("streptococcus_gordonii", "fusobacterium_nucleatum", 
                          "pseudopropionibacterium_propionicum", 
                          "tannerella_forysthia", "treponema_denticola",
                          "porphyromonas_gingivalis")

selected_groups <- c("Neanderthal", "PreagriculturalHuman", "ModernDayHuman")

selected_taxa_damage_subset <- c("streptococcus_gordonii", 
                                 "pseudopropionibacterium_propionicum", 
                                 "tannerella_forysthia",
                                 "treponema_denticola")


data_longdamage$Strand <- factor(data_longdamage$Strand, 
                                 levels = c("5pC>T", "3pG>A"))

data_longdamage$Position <- factor(data_longdamage$Position, 
                                   levels = c(1:25, -25:-1))

data_longdamage_plot <- data_longdamage %>% 
  filter(Taxon %in% selected_taxa_damage,
         Group %in% selected_groups) %>% 
  mutate(Group = as_factor(Group),
         Group = fct_relevel(Group, selected_groups),
         Taxon = as_factor(Taxon),
         Taxon = fct_relevel(Taxon, selected_taxa_damage)) %>%
  arrange(Group) %>%
  mutate(Sample = as_factor(Sample)) %>%
  filter(Strand == "5pC>T")

group_colours <- c(Neanderthal = "#e41a1c",
                   PreagriculturalHuman = "#ff7f00",
                   ModernDayHuman = "#ff7f00")

group_linetype <- c(Neanderthal = "solid",
                    PreagriculturalHuman = "solid",
                    ModernDayHuman = "dashed")

selected_samples <- c("PES001.B0101", "EMN001.A0101", "JAE006.A0101")

# From https://stackoverflow.com/a/16249622
capitalise_first <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

data_longdamage_final <- data_longdamage_plot %>%
  filter(Strand == "5pC>T",
         Group %in% selected_groups,
         Sample %in% selected_samples,
         Taxon %in% selected_taxa_damage_subset) %>%
  mutate(Taxon = str_replace(Taxon, "_", " ("),
         Taxon = str_replace(Taxon, "$", ")"),
         Taxon = capitalise_first(Taxon),
         Sample = case_when(grepl("PES", Sample) ~ "Neanderthal (PES001)",
                            grepl("EMN", Sample) ~ "Ancient Human (EMN001)",
                            grepl("JAE", Sample) ~ "Modern Day Human (JAE006)"),
         Sample = factor(Sample, levels = c("Neanderthal (PES001)",
                                            "Ancient Human (EMN001)",
                                            "Modern Day Human (JAE006)"))
         )


######damage plot########

damage_plot <- ggplot(data_longdamage_final, aes(Position, Frequency, 
                                  group = Taxon, 
                                  colour = Taxon, 
                                  linetype = Taxon)) +
  geom_line(size = 0.7) +
  xlab("Position on 5' strand") +
  ylab("Frequency C>T Subtitutions") +
  labs(tag = "c") +
  theme_minimal(base_size = 7) +
  scale_colour_manual(values = c("#7570b3", "#1b9e77","#d95f02", "#e7298a")) +
  scale_linetype_manual(values = c("solid", "solid","dashed", "dashed")) +
  facet_wrap(~ Sample, scales = "free_x") +
  theme(text  = element_text(family = "Roboto"),
        legend.position = "right") +
  guides(colour = guide_legend("Reference Genome"),
         linetype = guide_legend("Reference Genome"))

damage_plot

ggplot2::ggsave("Damage_Only.pdf", plot = damage_plot, device = cairo_pdf, width = 7, height = 3.5, units = "in")
ggplot2::ggsave("Damage_Only.png", plot = damage_plot, width = 7, height = 3.5, units = "in", dpi = 600)

