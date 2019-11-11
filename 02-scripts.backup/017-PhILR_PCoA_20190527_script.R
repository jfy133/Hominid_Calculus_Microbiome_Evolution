#' ---
#' title: "PhILR PCoA and PERMANOVA"
#' ---
#' 
#' # Preamble
#' 
#' Following the introducory vignette: http://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.html
#' 
#' Purpose of this notebook is to statistically show that indeed our groups
#' are distinguishable at a compositional level (i.e. all individuals of all
#' group could be derived from the same population, rather than individuals
#' from each group coming from group-specific populations).
#' 
#' We will run at Genus level as we do not yet know the divergence level
#' of different species which might be mapped to human-specific taxa due to
#' lack of non-human microbiota references.
#' 
#' # Preparation
#' 
#' ## For script conversion
#' 
#' As this notebook became more and more parameterised, I decided to make it
#' easy to convert into a script. But for the script, we also need to define
#' input arguments for the options in chunk 1.
#' 

#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
args = commandArgs(trailingOnly = TRUE)

if (args[1] == "" | args[1] == "-h" | args[1] == "--help") {
 cat("Usage: 017-PhILR_PCoA_20XXXXXX_script.R <db> <tax_level> <sources> <controls> <bad_samples> <sample_filter> <minsupp_multiplier>\n")
 cat("db: nt or refseq \n")
 cat("tax_level: genus or species \n") 
 cat("sources: noSources or withSources \n")
 cat("controls: noControls or withControls \n")
 cat("bad_samples: in or out \n")
 cat("sample_filter: sourcetracker, onepcburnin, twocburnin, fivepcburnin, tenpcburnin, withinvariation or none \n")
 cat("minsupp_multiplier: number to multiply foundation 0.01% minimum support")
 stop()
} else if (length(args) == 7) {
 db <- args[1]
 tax_level <- args[2]
 sources <- args[3]
 controls <- args[4]
 bad_samples <- args[5]
 sample_filter <- args[6]
 minsupp_multiplier <- as.numeric(args[7])
 script <- T
}

cat(args)

minsupp_multiplier <- as.numeric(minsupp_multiplier)
minsupp_threshold <- 0.01 * minsupp_multiplier
 

#' 
#' 
#' ## Data Loading and Cleaning
#' 
#' Load libraries
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
library(tidyverse) ## for general data cleaning
library(taxize) ## for NCBI taxonomic info collection
library(ape) ## for tree manipulation
library(phyloseq) ## for data format as input for PhILR
library(philr) ## for data transform
library(vegan) ## for statistical testing
library(ggtree) ## for tree visualisation
library(patchwork) ## for further visualisation assistance
library(broom) ## for saving stats summary tables
library(plotly) ## for interactive plots
library(pairwiseAdonis) ## for PERMANOVA bonferoni post-hoc

#' 
#' Load already generated data from MEGAN and metadata. We also need to export the 
#' same data with as a tree from MEGAN with the option: file > export > tree.
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------

## The tree related to the OTU table
if (tax_level == "genus" & db == "nt") {
 otu_tree <- read.tree("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/megan.backup/Evolution-Comparison_20190401_nt_prokaryotes_genus.nwk")
} else if (tax_level == "species" & db == "nt") {
 otu_tree <- read.tree("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/megan.backup/Evolution-Comparison_20190401_nt_prokaryotes_species.nwk")
} else if (tax_level == "genus" & db == "refseq") {
 otu_tree <- read.tree("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/megan.backup/Evolution-Comparison_20190410_refseq_prokaryotes_genus.nwk")
} else if (tax_level == "species" & db == "refseq") {
  otu_tree <- read.tree("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/megan.backup/Evolution-Comparison_20190410_refseq_prokaryotes_species.nwk")
}


## OTU tables
if (tax_level == "genus" & db == "nt") {
 otu_table <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/megan.backup/Evolution-Comparison_MEGAN_20190401-ex_absolute_genus_prokaryotes_summarised_nt.txt")
} else if (tax_level == "species" & db == "nt") {
 otu_table <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/megan.backup/Evolution-Comparison_MEGAN_20190401-ex_absolute_species_prokaryotes_summarised_nt.txt")
} else if (tax_level == "genus" & db == "refseq") {
 otu_table <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/megan.backup/Evolution-Comparison_MEGAN_20190410-ex_absolute_genus_prokaryotes_summarised_refseq.txt")
} else if (tax_level == "species" & db == "refseq") {
 otu_table <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/megan.backup/Evolution-Comparison_MEGAN_20190410-ex_absolute_species_all_summarised_refseq.txt")
}


## Predicted contaminant taxa to remove
if (tax_level == "genus" & db == "nt") {
 taxa_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/decontam.backup/decontam_taxa_to_remove_megan_nt_genus_combined_0.99_190411.tsv")
} else if (tax_level == "species" & db == "nt") {
 taxa_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/decontam.backup/decontam_taxa_to_remove_megan_nt_species_combined_0.99_190411.tsv")
} else if (tax_level == "genus" & db == "refseq") {
 taxa_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/decontam.backup/decontam_taxa_to_remove_megan_refseq_genus_combined_0.99_190411.tsv")
} else if (tax_level == "species" & db == "refseq") {
 taxa_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/decontam.backup/decontam_taxa_to_remove_megan_refseq_species_combined_0.99_190411.tsv")
}

## Metadata
raw_metadata <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/00-documentation.backup/02-calculus_microbiome-deep_evolution-individualscontrolssources_metadata_20190523.tsv")

## Bad samples to remove

if (sample_filter == "sourcetracker") {
 samples_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/sourcetracker.backup/sourcetracker_filtering_results_190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "nt" && sample_filter == "onepcburnin") {
 samples_to_remove <- read_tsv("home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninfilter1pc_nt_fractionOralThreshold_50_20190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "nt" && sample_filter == "twopcburnin") {
 samples_to_remove <- read_tsv("home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninfilter2pc_nt_fractionOralThreshold_50_20190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "nt" && sample_filter == "fivepcburnin") {
 samples_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninfilter5pc_nt_fractionOralThreshold_50_20190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "nt" && sample_filter == "tenpcburnin") {
 samples_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninfilter10pc_nt_fractionOralThreshold_50_20190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "nt" && sample_filter == "withinvariation") {
 samples_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninwithinfluctuationSDvariation_nt_fractionOralThreshold_50_20190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "refseq" && sample_filter == "onepcburnin") {
 samples_to_remove <- read_tsv("home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninfilter1pc_refseq_fractionOralThreshold_65_20190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "refseq" && sample_filter == "twopcburnin") {
 samples_to_remove <- read_tsv("home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninfilter2pc_refseq_fractionOralThreshold_65_20190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "refseq" && sample_filter == "fivepcburnin") {
 samples_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninfilter5pc_refseq_fractionOralThreshold_65_20190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "refseq" && sample_filter == "tenpcburnin") {
 samples_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninfilter10pc_refseq_fractionOralThreshold_65_20190509.tsv") %>% 
  filter(more_env == T)
} else if (db == "refseq" && sample_filter == "withinvariation") {
 samples_to_remove <- read_tsv("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/cumulative_decay.backup/cumulativeproportiondecay_burninwithinfluctuationSDvariation_refseq_fractionOralThreshold_65_20190509.tsv") %>% 
  filter(more_env == T)
}




#' 
#' Clean up to remove samples not required and then remove any OTUs that
#' now have no counts. Also remove OTUs that are likely lab contaminants
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
data_cleaner <- function(x) {
 colnames(x) <- gsub("_S.*_L.*_R1_.*.fastq.combined.fq.prefixed.extractunmapped.bam","", colnames(x))
 colnames(x) <- gsub("_S.*_L00.*_R1_.*.fastq.merged.prefixed.hg19unmapped", "", colnames(x))
 colnames(x) <- gsub("_S.*_L00.*_R1_.*.fastq.extractunmapped.bam", "", colnames(x))
 colnames(x) <- gsub("_S.*_L.*_R1_.*.fastq.merged", "", colnames(x))
 colnames(x) <- gsub(".prefixed.hg19unmapped", "", colnames(x))
 colnames(x) <- gsub("_S0_L003_R1_001.sorted.bam.unmapped", "", colnames(x))
 colnames(x)[1] <- "Taxon"
 return(x)
}

## Remove col cruft of Metadata
raw_metadata <- rename(raw_metadata, Individual = `#SampleID`)

## Remove bad sources from OTU table
if (bad_samples == "in") {
 otu_table <- otu_table %>% 
 data_cleaner
} else if (bad_samples == "out") {
 otu_table <- otu_table %>% 
  data_cleaner %>% 
  dplyr::select(-one_of(samples_to_remove %>% 
                         left_join(raw_metadata, 
                                   by = c("sample" = "Individual")) %>%
                         select(sample, SourceSink, Sample_or_Control) %>% 
                         filter(SourceSink == "sink", 
                                Sample_or_Control == "Sample") %>% 
                          pull(sample)) 
 )
}
 


## Conditional filtering out of sources and/or Controls
if (sources == "withSources") {
 NA
} else if (sources == "noSources") {
 otu_table <- otu_table %>% 
 dplyr::select(Taxon, one_of(filter(raw_metadata, SourceSink == "sink") %>% 
             pull(Individual)))
}

if (controls == "withControls") {
 NA
} else if (controls == "noControls") {
 otu_table <- otu_table %>% 
 dplyr::select(Taxon, one_of(filter(raw_metadata, 
                                    Sample_or_Control == "Sample") %>% 
             pull(Individual)), 
             contains("ARS"))
}


## Filter taxa not passing min support threshold 
if (db == "nt") {
 otu_table <- otu_table %>% 
  gather(Individual, Value, 2:ncol(.)) %>% 
  left_join(select(raw_metadata, Individual, Min_Support_Reads_Threshold_MALT)) %>%
  mutate(Threshold = Min_Support_Reads_Threshold_MALT * minsupp_multiplier) %>%
  mutate(Threshold = as.numeric(Threshold)) %>%
  mutate(Filter_Passed = if_else(Value >= Threshold, 1, 0)) %>% 
  filter(Filter_Passed == 1) %>%
  select(Taxon, Individual, Value) %>%
  spread(Individual, Value, fill = 0)
} else if (db == "refseq") {
 otu_table <- otu_table %>% 
  gather(Individual, Value, 2:ncol(.)) %>% 
  left_join(select(raw_metadata, Individual, Min_Support_Reads_Threshold_MALT_refseq)) %>%
  mutate(Threshold = Min_Support_Reads_Threshold_MALT_refseq * minsupp_multiplier) %>%
  mutate(Threshold = as.numeric(Threshold)) %>%
  mutate(Filter_Passed = if_else(Value >= Threshold, 1, 0)) %>% 
  filter(Filter_Passed == 1) %>%
  select(Taxon, Individual, Value) %>%
  spread(Individual, Value, fill = 0)
}
 
## Convert to matrix
otu_matrix <- as.matrix(dplyr::select(otu_table, -Taxon))
rownames(otu_matrix) <- otu_table$Taxon

## Remove any taxa that were unique to the bad samples
pos_otus <- rowSums(otu_matrix)
pos_otus <- pos_otus[pos_otus != 0] 

## Remove lab contaminants
otu_matrix <- subset(otu_matrix, !rownames(otu_matrix) %in% (taxa_to_remove %>% pull))
otu_matrix_final <- subset(otu_matrix, rownames(otu_matrix) %in% names(pos_otus))

rownames(otu_matrix_final) <- gsub(" ", "_", rownames(otu_matrix_final))

#' 
#' ## Make Taxonomy Table
#' 
#' We also need to make our taxonomy table at species level. This takes a long 
#' to run, so I've generated ones in the past and re-load them when needed.
#' 
## ----eval = FALSE--------------------------------------------------------
## taxonomy_to_tibble <- function(x) {
##  y <- gsub("\\[", "", x)
##  y <- gsub("\\]", "", y)
##  y <- gsub("\\(.*\\)", "", y)
## 
##  ttt_uid <- taxize::get_uid(y, db = "ncbi") ## inspired by myTAI::taxonomy()
##  ttt_result <- taxize::classification(ttt_uid)[[1]] ## inspired by myTAI::taxonomy()
## 
##  if ( is.na(ttt_result)[1] ) {
##   ttt_out <- x
##  } else {
##   ttt_out <- as_tibble(ttt_result %>%
##    mutate(taxon = paste(x)) %>%
##    filter(rank %in% c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
##    dplyr::select(name, rank, taxon) %>%
##    spread(rank, name))
##  }
##  return(ttt_out)
##   Sys.sleep(2)
## }
## 
## # taxonomy_summary <- rownames(otu_matrix_final) %>% map(., .f = taxonomy_to_tibble)
## 
## taxonomy_summary <- tibble(taxon = character(),
##               superkingdom = character(),
##               kingdom = character(),
##               phylum = character(),
##               class = character(),
##               order = character(),
##               family = character(),
##               genus = character(),
##               species = character()
##               )
## 
## n <- 0
## tot <- nrow(otu_matrix_final)
## fail_taxa <- c()
## 
## ## Temp due to curl timeout, change number beginning øf range depending on
## ## where n got to at time out
## otu_matrix_final_sub <- otu_matrix_final[903:nrow(otu_matrix_final),]
## 
## ## If time out, replace otu_matrix_final with otu_matrix_final_sub and re-run
## for(i in rownames(otu_matrix_final_sub)) {
##  print(paste(i," - ", format(round((n / tot) * 100), nsmall = 2), "%", sep = ""))
##  n <- n + 1
##  taxonomy_temp <-taxonomy_to_tibble(i)
## 
##  if (length(taxonomy_temp) = 1) {
##   fail_taxa <-append(fail_taxa, i)
##  } else {
##    taxonomy_summary <- bind_rows(taxonomy_summary, taxonomy_temp)
##  }
## }
## 
## taxonomy_summary <- taxonomy_summary %>% distinct()
## 

#' 
#' And save
#' 
## ----eval = FALSE--------------------------------------------------------
## save(taxonomy_summary, file = paste("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/0-evolution-philr_taxonomytable", "_", db,"_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".robj", sep = ""))

#' 
#' For the failed taxa, manually look up paths and add taxonomic paths and add them
#' and save.
#' 
## ----eval = FALSE--------------------------------------------------------
## ## To check for the ones still missing
## missing <- setdiff(rownames(otu_matrix_final), rownames(taxonomy_table_final))
## 
## out <- tibble(taxon = character(),
##               superkingdom = character(),
##               kingdom = character(),
##               phylum = character(),
##               class = character(),
##               order = character(),
##               family = character(),
##               genus = character(),
##               species = character()
##               )
## 
## n <- 0
## tot <- length(missing)
## fail_taxa <- c()
## 
## for(i in missing[1:length(missing)]) {
##  print(paste(i," - ", format(round((n / tot) * 100), nsmall = 2), "%", sep = ""))
##  n <- n + 1
##  taxonomy_temp <- taxonomy_to_tibble(i)
## 
##  if (length(taxonomy_temp) == 1) {
##   fail_taxa <- append(fail_taxa, i)
##  } else {
##    out <- bind_rows(out, taxonomy_temp)
##  }
## }
## 
## taxonomy_summary <- bind_rows(taxonomy_summary, out) %>% distinct()
## 
## 
## ## Some manual additions
## temp <- read_tsv("~/Downloads/Genus_extra.csv")
## taxonomy_summary <- bind_rows(taxonomy_summary, temp) %>% distinct()
## 
## save(taxonomy_summary, file = "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/0-evolution-philr_taxonomytable_withsources_withControls_all_20190404.robj")

#' 
#' Load a previously made taxonomy table
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
if (db == "nt") {
 load("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/0-evolution-philr_taxonomytable_withsources_withblanks_all_20190404.robj")
} else if (db == "refseq") {
 load("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/0-evolution-philr_taxonomytable_refseq_species_withSources_withControls_20190115.robj")
}

taxonomy_table <- as.data.frame(taxonomy_summary[2:ncol(taxonomy_summary)])
rownames(taxonomy_table) <- taxonomy_summary$taxon 
taxonomy_table <- as.matrix.data.frame(taxonomy_table)


## To make sure OTU names match tree, which in import replaced spaces with 
## undercores
#rownames(taxonomy_table) <- gsub(" ", "_", rownames(taxonomy_table))

## Filter for genus level
if (tax_level == "genus") {
 taxonomy_table_final <- taxonomy_table[,c("superkingdom", 
                                           "kingdom", 
                                           "phylum", 
                                           "class", 
                                           "order", 
                                           "family", 
                                           "genus")] %>% 
  as.data.frame.matrix
 rownames(taxonomy_table_final) <- c()
 taxonomy_table_final <- subset(taxonomy_table_final, 
                                genus %in% gsub("_", 
                                                " ", 
                                                rownames(otu_matrix_final))) %>% 
   unique()
 taxonomy_table_final <- unique(taxonomy_table_final)
 row.names(taxonomy_table_final) <- taxonomy_table_final[,"genus"]
 taxonomy_table_final <- as.matrix(taxonomy_table_final)
 taxonomy_table_final <- unique(taxonomy_table_final)

 } else if (tax_level == "species") {
 taxonomy_table_final <- taxonomy_table
}

#' 
#' 
#' ## Check Tree Format Valid
#' 
#' As need a purely bifurcating tree (thus do not allow polytomy), and need a
#' rooted tree we should check this with:
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
is.rooted(otu_tree)
is.binary.tree(otu_tree)

#' 
#' If we don't have a binary tree, solve polytomies with 
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
otu_tree <- multi2di(otu_tree)
is.binary.tree(otu_tree)

#' 
#' Or if this doesn't work we can remove singletons (as in single-intermediate 
#' nodes between a leaf and a parent node) with
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
otu_tree <- collapse.singles(otu_tree) # from http://blog.phytools.org/2018/05/when-phylogeny-fails-isbinary-but-is.html
is.binary.tree(otu_tree)


#' 
#' ## Prepare Metdata
#' 
#' Prepare sample metadata
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
meta_data <- as.data.frame(raw_metadata[2:ncol(raw_metadata)], 
              row.names = raw_metadata$Individual, 
              colnames = col.names(raw_metadata[2:col(raw_metadata),]))
rownames(meta_data) <- raw_metadata$Individual


#' 
#' 
#' # Clustering Analysis
#' 
#' ## PhILR Transform
#' 
#' ### Make PhyloSeq object
#' 
#' Make our Phyloseq object.
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
ps_data <- phyloseq(otu_table(otu_matrix_final, taxa_are_rows = TRUE), tax_table(taxonomy_table_final), phy_tree(otu_tree), sample_data(meta_data))

#' 
#' ### Zero Replacement
#' 
#' We also need to a run zero removal procedure as the *LR methods require
#' positive numbers. Here just doing pseudocount (note might be flawed - see Fodor 
#' paper), but here just following PhILR tutorial. Have used zCompositions
#' in the past but that was with random guessing of parameters (or using
#' recommended defaults)
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
ps_data <- transform_sample_counts(ps_data, function(x) x + 1)

#' 
#' ### Final Checks
#' 
#' We fix the tree nodes with and something else from the PhILR tutorial I don't 
#' fully understand with name balancing and transpose to sample row/taxa column
#' composition convention.
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
phy_tree(ps_data) <- makeNodeLabel(phy_tree(ps_data), method = "number", prefix = 'n')
name.balance(phy_tree(ps_data), tax_table(ps_data), 'n1')

otu.table <- t(otu_table(ps_data))
tree <- phy_tree(ps_data)


#' 
#' Some last checks
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
otu.table <- t(otu_table(ps_data))
tree <- phy_tree(ps_data)
metadata <- sample_data(ps_data)
tax <- tax_table(ps_data)

otu.table[1:2,1:2] # OTU Table
tree
head(metadata,2)

#' 
#' ### Run PhILR
#' 
#' Now run PhILR (hopefully)... here goes...
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
ps_data.philr <- philr(otu.table, 
            tree,
            part.weights = "enorm.x.gm.counts",
            ilr.weights = "uniform")

ps_data.philr[1:5,1:5]

#' 
#' And we can now do ordination within the PhILR space with 
#' Euclidean distances
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
ps_data.dist <- dist(ps_data.philr, method = "euclidean")
ps_data.pcoa <- ordinate(ps_data, 'PCoA', distance = ps_data.dist)

#' 
#' ### Plot PhILR
#' 
#' Pre-plotting formatting
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
env_colours <- c("#1f78b4", 
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

names(env_colours) <- c("Howler_Monkey", 
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

 env_shapes <- c(8,
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

 
 names(env_shapes) <- c("Howler_Monkey", 
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
          )
 
common_colours <- c(Alouatta = "#1f78b4", Gorilla = "#6a3d9a", Pan = "#33a02c", 
          `Homo (Neanderthal)` = "#ff7f00", 
          `Homo (Modern Human)` = "#ff7f00", ExtractionControl = "#d9d9d9", 
          LibraryControl = "#d9d9d9", Plaque = "#d9d9d9", Gut = "#d9d9d9", 
          Skin = "#d9d9d9", Sediment = "#d9d9d9", EnvironmentalControl = "#d9d9d9")

common_shapes <- c(Alouatta = 8, Gorilla = 0, Pan = 1, 
          `Homo (Neanderthal)` = 2, `Homo (Modern Human)` = 6,
          ExtractionControl = 10, LibraryControl = 13, Plaque = 9, 
          Gut = 4, Skin = 14, Sediment = 7, EnvironmentalControl = 12)


 
 

#' 
#' 
#' And plot with host populations
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------

## Convert to tidy data and add metadata
ps_data.pcoa_tib <- as_tibble(ps_data.pcoa$vectors, rownames = "Individual") %>%
 left_join(as_tibble(metadata, rownames = "Individual"))

ps_data.pcoa_tib$Env <- factor(ps_data.pcoa_tib$Env, levels = c("Howler_Monkey", 
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
           )

## Extract percentage variation
percentage_data <- round(100 * ps_data.pcoa$values[1] / sum(ps_data.pcoa$values[1]), 2)
rownames(percentage_data) <- rownames(ps_data.pcoa$value)

## Colour/Shapes by sample type and source/sink
pcoa_plot_hostpops <- ggplot(ps_data.pcoa_tib, aes(Axis.1, Axis.2, colour = Env, shape = Env, text = Individual)) +
 geom_point(size = 2) +
 scale_colour_manual(values = env_colours) +
 scale_shape_manual(values = env_shapes) +
 xlab(paste("Axis 1 (", percentage_data[1,], "%)", sep = "")) +
 ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
 theme_minimal(base_family = "Roboto", base_size = 7) +
 labs(shape = "Host Population", colour = "Host Population") +
 guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))


if (script == F) {
 pcoa_plot_hostpops
}

ggsave(paste("01a-philr_pcoa_malt_euclidean_axis1axis2_hostpopulations_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = pcoa_plot_hostpops, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 3.5, units = "in", dpi = 600)


## Axis 3 and 2
pcoa_plot_hostpops <- ggplot(ps_data.pcoa_tib, aes(Axis.3, Axis.2, colour = Env, shape = Env)) +
 geom_point(size = 2) +
 scale_colour_manual(values = env_colours) +
 scale_shape_manual(values = env_shapes) +
 xlab(paste("Axis 3 (", percentage_data[3,], "%)", sep = "")) +
 ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
 theme_minimal(base_family = "Roboto", base_size = 7) +
 labs(shape = "Host Population", colour = "Host Population") + 
 guides(shape = guide_legend(ncol = 2), colour = guide_legend(ncol = 2))



if (script == F) {
pcoa_plot_hostpops
}

ggsave(paste("01b-philr_pcoa_malt_euclidean_axis3axis2_hostpopulations_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = pcoa_plot_hostpops, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 3.5, units = "in", dpi = 600)




#' 
#' And plot without host populations (but Neanderthals indicated)
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------

## Set factor
ps_data.pcoa_tib$Env <- factor(ps_data.pcoa_tib$Env, levels = c("Howler_Monkey", 
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
           )

## Column name change for proper plotting
ps_data.pcoa_tib$Host_Common <- factor(ps_data.pcoa_tib$Host_Common, 
                    levels = c("Alouatta", "Gorilla", "Pan", "Homo (Modern Human)", "Homo (Neanderthal)",
                         "ExtractionControl", "LibraryControl", "Plaque", "Gut", "Skin", "Sediment", "EnvironmentalControl"))

## Colour/Shapes by sample type and source/sink
pcoa_plot_hostgenus <- ggplot(ps_data.pcoa_tib, aes(Axis.1, Axis.2, colour = Host_Common, shape = Host_Common, fill = Host_Common)) +
 geom_point(size = 2) +
 scale_colour_manual(values = common_colours) +
 scale_fill_manual(values = common_colours) +
 scale_shape_manual(values = common_shapes) +
 xlab(paste("Axis 1 (", percentage_data[1,], "%)", sep = "")) +
 ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
 theme_minimal(base_family = "Roboto", base_size = 7) +
 labs(shape = "Host Genus", colour = "Host Genus", fill = "Host Genus") +
 guides(shape = guide_legend(ncol = 2), 
        colour = guide_legend(ncol = 2),
        fill = guide_legend(ncol = 2))


if (script == F) {
pcoa_plot_hostgenus
}

ggsave(paste("02a-philr_pcoa_malt_euclidean_axis1axis2_hostgenus_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = pcoa_plot_hostgenus, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 3.5, units = "in", dpi = 600)


## Axis 3 and 2

pcoa_plot_hostgenus <- ggplot(ps_data.pcoa_tib, aes(Axis.3, Axis.2, colour = Host_Common, shape = Host_Common, fill = Host_Common)) +
 geom_point(size = 2) +
 scale_colour_manual(values = common_colours) +
 scale_fill_manual(values = common_colours) +
 scale_shape_manual(values = common_shapes) +
 xlab(paste("Axis 3 (", percentage_data[3,], "%)", sep = "")) +
 ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
 theme_minimal(base_family = "Roboto", base_size = 7) +
 labs(shape = "Host Genus", colour = "Host Genus", fill = "Host Genus") +
 guides(shape = guide_legend(ncol = 2), 
        colour = guide_legend(ncol = 2),
        fill = guide_legend(ncol = 2))

if (script == F) {
pcoa_plot_hostgenus
}

ggsave(paste("02b-philr_pcoa_malt_euclidean_axis3axis2_hostgenus_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = pcoa_plot_hostgenus, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 3.5, units = "in", dpi = 600)


#' 
#' Host genera again, but removing less important sources
#' 
## ------------------------------------------------------------------------
## get names of less important sources

sourcecommon_to_filterout <- c("ExtractionControl", "LibraryControl", 
                               "Skin", "Sediment")
samplecommon_to_filterout <-  raw_metadata %>% 
  filter(Host_Common %in% sourcecommon_to_filterout) %>%
  pull(Individual)

## Colour/Shapes by sample type and source/sink
pcoa_plot_hostgenus_reduced <- ggplot(ps_data.pcoa_tib %>%
                                        filter(!Individual %in% 
                                                 samplecommon_to_filterout), 
                                      aes(Axis.1, Axis.2, 
                                          colour = Host_Common, 
                                          shape = Host_Common, 
                                          fill = Host_Common)) +
 geom_point(size = 2) +
 scale_colour_manual(values = common_colours) +
 scale_fill_manual(values = common_colours) +
 scale_shape_manual(values = common_shapes) +
 xlab(paste("Axis 1 (", percentage_data[1,], "%)", sep = "")) +
 ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
 theme_minimal(base_family = "Roboto", base_size = 7) +
 labs(shape = "Host Genus", colour = "Host Genus", fill = "Host Genus")

if (script == F) {
pcoa_plot_hostgenus_reduced
}

ggsave(paste("02a-philr_pcoa_malt_euclidean_axis1axis2_hostgenus_", db , "_", tax_level, "_", sources, "reduced_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = pcoa_plot_hostgenus_reduced, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 3.5, units = "in", dpi = 600)


## Axis 3 and 2

pcoa_plot_hostgenus_reduced <- ggplot(ps_data.pcoa_tib %>%
                                        filter(!Individual %in% 
                                                 samplecommon_to_filterout), 
                                      aes(Axis.3, Axis.2, 
                                          colour = Host_Common, 
                                          shape = Host_Common, 
                                          fill = Host_Common)) +
 geom_point(size = 2) +
 scale_colour_manual(values = common_colours) +
 scale_fill_manual(values = common_colours) +
 scale_shape_manual(values = common_shapes) +
 xlab(paste("Axis 3 (", percentage_data[3,], "%)", sep = "")) +
 ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
 theme_minimal(base_family = "Roboto", base_size = 7) +
 labs(shape = "Host Genus", colour = "Host Genus", fill = "Host Genus")

if (script == F) {
pcoa_plot_hostgenus_reduced
}

ggsave(paste("02b-philr_pcoa_malt_euclidean_axis3axis2_hostgenus_", db , "_", tax_level, "_", sources, "reduced_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = pcoa_plot_hostgenus_reduced, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 3.5, units = "in", dpi = 600)

#' 
#' 
#' Optional plot of whether sample from oral cavity or not (this will break
#' depending on filtering options set in Data Loading/Cleaning)
#' 
## ----fig.height=5, fig.width=7-------------------------------------------
## Colour/Shapes by sample type and source/sink
pcoa_plot_isoral_axis1axis2 <- ggplot(ps_data.pcoa_tib, aes(Axis.1, Axis.2, colour = Is_Oral, shape = Env)) +
 geom_point(size = 2) +
 scale_shape_manual(values = env_shapes) +
 xlab(paste("Axis 1 (", percentage_data[1,], "%)", sep = "")) +
 ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
 theme_minimal(base_family = "Roboto", base_size = 7) +
 labs(shape = "Typical Isolation Source", colour = "Host Population") +
 guides(shape = guide_legend(ncol = 2), 
        colour = guide_legend(ncol = 2),
        fill = guide_legend(ncol = 2))

if (script == F) {
pcoa_plot_isoral_axis1axis2
}

pcoa_plot_isoral_axis3axis2 <- ggplot(ps_data.pcoa_tib, aes(Axis.3, Axis.2, colour = Is_Oral, shape = Env)) +
 geom_point(size = 2) +
 scale_shape_manual(values = env_shapes) +
 xlab(paste("Axis 3 (", percentage_data[3,], "%)", sep = "")) +
 ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
 theme_minimal(base_size = 7, base_family = "Roboto") +
 labs(colour = "Host Population", shape = "Typical Isolation Source") +
 guides(colour = guide_legend(ncol = 2),
        shape = guide_legend(ncol = 2), 
        fill = guide_legend(ncol = 2))

if (script == F) {
pcoa_plot_isoral_axis3axis2
}

ggsave(paste("03a-philr_pcoa_malt_euclidean_axis1axis2_isoral_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = pcoa_plot_isoral_axis1axis2, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 5, units = "in", dpi = 600)

ggsave(paste("03b-philr_pcoa_malt_euclidean_axis3axis2_isoral_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = pcoa_plot_isoral_axis3axis2, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 5, units = "in", dpi = 600)


#' 
#' Additional plot of each Host Genus separately ordinated
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
final_plots_a1a2 <- list()
final_plots_a3a2 <- list()

for (i in filter(raw_metadata, 
                 Sample_or_Control == "Sample", 
                 SourceSink == "sink") %>% 
     select(Host_Genus) %>% 
     distinct() %>% 
     pull()
     ) {
 print(i)
 sample_list <- filter(raw_metadata, Host_Genus == i) %>% 
  pull(Individual)
 temp_otu.table <- subset(otu.table, rownames(otu.table) %in% sample_list)
 temp_philr <- philr(temp_otu.table, 
    tree,
    part.weights = "enorm.x.gm.counts",
    ilr.weights = "uniform")
 temp_philr_dist <- dist(temp_philr, method = "euclidean")
 temp_philr_pcoa <- ordinate(temp_otu.table, 'PCoA', distance = temp_philr_dist)
 temp_philr_pcoa_tib <- as_tibble(temp_philr_pcoa$vectors, rownames = "Individual") %>%
  left_join(as_tibble(metadata, rownames = "Individual"))
 
 temp_philr_pcoa_tib$Env <- factor(temp_philr_pcoa_tib$Env, levels = c("Howler_Monkey", 
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
 )
 
 ## Extract percentage variation
 percentage_data <- round(100 * temp_philr_pcoa$values[1] / sum(temp_philr_pcoa$values[1]), 2)
 rownames(percentage_data) <- rownames(temp_philr_pcoa$value)
 

 
 ## Colour/Shapes by sample type and source/sink
 temp_plot_a1a2 <- ggplot(temp_philr_pcoa_tib, aes(Axis.1, Axis.2, colour = Env, shape = Env)) +
  geom_point(size = 2) +
  scale_colour_manual(values = env_colours) +
  scale_shape_manual(values = env_shapes) +
  xlab(paste("Axis 1 (", percentage_data[1,], "%)", sep = "")) +
  ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
  theme_minimal(base_size = 7, base_family = "Roboto") +
  theme_minimal(base_family = "Roboto", base_size = 7) +
   labs(shape = "Host Population", colour = "Host Population")
  
 final_plots_a1a2[[i]] <- temp_plot_a1a2
 
 ## Colour/Shapes by sample type and source/sink
 temp_plot_a3a2 <- ggplot(temp_philr_pcoa_tib, aes(Axis.3, Axis.2, colour = Env, shape = Env)) +
  geom_point(size = 2) +
  scale_colour_manual(values = env_colours) +
  scale_shape_manual(values = env_shapes) +
  #ylim(-1500,1500) +
  xlab(paste("Axis 3 (", percentage_data[3,], "%)", sep = "")) +
  ylab(paste("Axis 2 (", percentage_data[2,], "%)", sep = "")) +
  theme_minimal(base_size = 7, base_family = "Roboto") +
  theme_minimal(base_family = "Roboto", base_size = 7) +
  labs(shape = "Host Population", colour = "Host Population")

   
 final_plots_a3a2[[i]] <- temp_plot_a3a2
}

combined_hostgenusseparated_a1a2_plot <- final_plots_a1a2$Alouatta + final_plots_a1a2$Gorilla + final_plots_a1a2$Pan + final_plots_a1a2$Homo
combined_hostgenusseparated_a3a2_plot <- final_plots_a3a2$Alouatta + final_plots_a3a2$Gorilla + final_plots_a3a2$Pan + final_plots_a3a2$Homo



if (script == F) {
combined_hostgenusseparated_a1a2_plot
combined_hostgenusseparated_a3a2_plot
}
 
ggsave(paste("04a-philr_pcoa_malt_euclidean_axis1axis2_hostgenusseparated_", db , "_",tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = combined_hostgenusseparated_a1a2_plot, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 5, units = "in", dpi = 600)

ggsave(paste("04b-philr_pcoa_malt_euclidean_axis3axis2_hostgenusseparated_", db , "_",tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = combined_hostgenusseparated_a3a2_plot, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 5, units = "in", dpi = 600)


#' 
#' 
#' ## Hierarchical Clustering
#' 
#' Moving away from the PhILR Tutorial, what happens if we run hierarchical 
#' cluster analysis instead, using the distance matrix made in the previous
#' method.
#' 
## ----eval = F------------------------------------------------------------
## ## Cluster
## ps_data.philr_hclust <- hclust(ps_data.dist, method = "ward.D2")
## 
## ## Convert hclust to phylo object for plotting
## ps_data.philr_clust_phylo <- as.phylo(ps_data.philr_hclust)
## 
## ## Plot
## ps_data.philr_clust_dendro <- ggtree(as.phylo(ps_data.philr_clust_phylo), layout = "rectangular", right = TRUE, ladderize = TRUE) %<+%
##  raw_metadata +
##  geom_tiplab(aes(colour = Env), size = 3, hjust = -0.1) +
##  geom_tippoint(aes(shape = Env, colour = Env), stroke = 1.1, size = 2) +
##  scale_shape_manual(values = env_shapes) +
##  scale_color_manual(values = env_colours) +
##  ggtitle("ward.D2 Hierarchical clustering of Euclidean Distances") +
##  theme_tree2() +
##  theme(text = element_text(family = "Roboto", size = 7))
## 
## if (script == F) {
## ps_data.philr_clust_dendro
## }
## 
## ggsave(paste("05-philr_pcoa_malt_euclidean_hclustwardD2_hostpopulation_", db , "_",tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = ps_data.philr_clust_dendro, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 7, height = 7, units = "in", dpi = 600)

#' 
#' ## Summary
#' 
#' From both the PCoA and the clustering of the distances between as represented
#' in the hierarchical clustering that there does seem to be clustering by Host
#' genus.
#' 
#' Interesting to note, however, is that this doesn't follow host relationship 
#' as Gorillas and Humans are more similar to each other than Humans and Chimps.
#' 
#' So likely there is some other driver rather than host phylogeny. However,
#' we should properly test whether indeed the groups are significantly different.
#' 
#' # Multivariate Analysis
#' 
#' Next we want to calculate if the groups by host are truly distinct. For
#' this people typically use PERMANOVA based on a distance matrix used to visualise
#' the groups. It is nice as it is semi-parametric (i.e. doesn't need normal
#' distributions), yet still allows downstream typical downstream statistical
#' tests.
#' 
#' Howver, we first need to check some assumptions related to it as described in
#' Anderson 2017 [(Wiley StatsRef: Statistics Reference Online) ](https://onlinelibrary.wiley.com/doi/full/10.1002/9781118445112.stat07841).
#' We can assume we have a good dissimilarity measure, as we have followed current
#' best practises for CoDa data and also are using the 'simplest' distance measure
#' which is Euclidean.
#' 
#' ## Dispersion Homogenity Test
#' 
#' However, the second assumption is that there is homogenity of dispersions among
#' groups (or beta-dispersion). While the method is robust when this is applied
#' on balanced designs (i.e. same sample sizes per group), it is not robust
#' when you have different numbers of sample sizes (for example: 
#' http://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html, 
#' or originally published Anderson and Walsh (2013 in Ecological Monographs 83(4))). 
#' 
#' Anderson thus recommends using the 'PERMDISP' test which is implemented as 
#' `betadisper()` in R's `vegan` library (originally published in Anderson et al. 
#' (2006) Biometrics 62(1) 245-253)
#' 
#' First we need to make a vector for the grouping of the sample, with the order 
#' matching that in the distance matrix.
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
metadata_subset <- raw_metadata %>% 
 filter(Individual %in% rownames(otu.table)) %>% 
 dplyr::select(Individual, Host_Genus)

ind_groups <- metadata_subset %>% pull(Host_Genus)

#' 
#' Now we can run the `betadisper()` function to get the average distances from
#' the centroid
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
ps_data.betadisper <- betadisper(ps_data.dist, ind_groups, type = "centroid")
ps_data.betadisper

## to save
sink(paste("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/06-betadisper_homogenity_", db , "_",tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".txt", sep = ""))
ps_data.betadisper
sink(file = NULL)


#' 
#' To check that we are violating the assumption of homogenity of 
#' variances (i.e. the amount beta-dispersion does not vary between each group), 
#' we can check this with an ANOVA.
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
## to save

ps_data.betadisper_anova <- anova(ps_data.betadisper) %>% tidy

ps_data.betadisper_anova

write_tsv(ps_data.betadisper_anova, paste("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/07-betadisper_homogenity_anova_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".txt", sep = ""))


#' 
#' Here, see we do have a significant difference between each group.
#' 
#' We can double-double check this with an alternative method suggested by the
#' `vegan` library developers, were we can use a permutation test.
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
ps_data.betadisper_permutest <- permutest(ps_data.betadisper , pairwise = TRUE)
ps_data.betadisper_permutest

# to save
sink(paste("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/08-betadisper_homogenity_permutest_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".txt", sep = ""))

ps_data.betadisper_permutest
sink(file = NULL)

#' 
#' For species and genus level, the F and P value suggest that the variances in the 
#' dispersions are indeed significantly different (i.e. unlikely to have been 
#' picked from the same distribution by chance).
#' 
#' According to the `vegan` help page on the `betadisper` function we can plot 
#' this with
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------

betadisper_scatterplot <- plot(ps_data.betadisper, main = 'PCoA Beta Dispersion of Groups', xlab = "Axis 1", ylab = "Axis 2", cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7, cex.sub = 0.7)

if (script == F) {
 betadisper_scatterplot
}

## Note: can't #ggsave plot()

genus_colours <- c(Alouatta = "#F7756C", Gorilla = "#A2A500", Pan = "#00BF7D", 
          Homo = "#E76BF2", Control = "#6a3d9a", Plaque = "#fdbf6f", 
          Gut = "#a6cee3", Skin = "#fb9a99", Sediment = "#b2df8a")


betadisper_boxplot <- tibble(ps_data.betadisper$distances, rownames = names(ps_data.betadisper$distances)) %>% 
 rename(Individual = rownames, Distance_to_centroid = `ps_data.betadisper$distances`) %>% 
 left_join(raw_metadata) %>% 
 rename(Host = Host_Genus) %>%
 ggplot(aes(Host, Distance_to_centroid, fill = Host)) + 
 geom_boxplot() + 
 xlab("Host Genus") +
 ylab("Distance to Centroid") +
 scale_fill_manual(values = genus_colours) +
 theme_minimal(base_size = 7, base_family = "Roboto") +
 labs(fill = "Host Genus") +
 theme(axis.text.x = element_text(angle = 30))

if (script == F) {
 betadisper_boxplot
}
ggsave(paste("09-philr_betadisper_malt_boxplot_hostgenus_", db , "_",tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = betadisper_boxplot, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 3.5, height = 3.5, units = "in", dpi = 600)


#' 
#' Where we can explicitly see large variations in dispersion.
#' 
#' ## PERMANOVA
#' 
#' The PERMANOVA null hypothesis is that
#' under the assumption of 'exchangeability of the sample units among the groups'
#' (as in if you randomly re-assign a sample to a different group), the centroids 
#' (X~Y mean) among the groups is equvalient for all groups. Thus, if the null 
#' hypothesis is true then any observered differences would be the same as if you 
#' randomly allocated sample units to the groups. In contrast, if it is false, 
#' then the observed differences is actually from different centroids between 
#' the groups.
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
permanova_factors <- factor(ind_groups)
original_result <- adonis(ps_data.dist ~ ind_groups) %>% print

## to save
sink(paste("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/10a-permanova_initial_hostgenuscentroids_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".txt", sep = ""))
 original_result
sink(file = NULL)

## Store values for later
original_f <- original_result$aov.tab$F.Model[1]
original_p <- original_result$aov.tab$`Pr(>F)`[1]


#' 
#' The results here would suggest that we do a significant difference between 
#' the centroid of each group, suggesting that they are distinct. However, we 
#' violate of assumption of homogenity between the beta-dispersion of each group. 
#' This is possibly due to us having different sample sizes (see next section).
#' 
#' ## PERMANOVA with Multiple Correction Testing 
#' 
#' We also need to perform multiple-testing correction as we are using the same
#' group over and over in this test, if I understand correctly.
#' 
## ------------------------------------------------------------------------
original_result_posthoc <- pairwise.adonis(ps_data.dist, ind_groups)

original_result_posthoc <- original_result_posthoc %>% as_tibble()
  
write_tsv(original_result_posthoc, 
          paste("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/10b-permanova_initial_hostgenuscentroids_multiplecorrectiontesting_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".txt", sep = ""))




#' 
#' 
#' 
#' ## PERMANOVA with bootstrapping
#' 
#' While a method to account for this issue has been published (Anderson et al. 
#' 2017 Australian and New Zealand Journal of Statistics), it is very recent and 
#' I could not find an R implementation fo it.
#' 
#' We could instead try bootstrapping each group down to the same number of 
#' individuals, run multiple times, then compare the differences in the result.
#' 
#' https://www.researchgate.net/post/How_to_run_PERMANOVA_in_R_to_compare_community_comp_between_sites_when_data_is_heteroscedastic_with_unbalanced_designs_RE_Anderson_et_al_2017
#' 
#' To get the number of individuals to bootstrap down to
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
as_tibble(ind_groups) %>% group_by(value) %>% summarise(No_Individuals = n())

#' 
#' Lets say 10 per group, and drop Alouatta as our outgroup. Now lets run the 
#' PERMANOVA multiple times with randomly selected individuals and store the 
#' output of each into vectors.
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
## Give an OTU table, a table with samples and gorups
boot_permanova <- function(tab, meta, group, ind, tr, n_boot) {
 
 ## Assign empty vectors for results
 betadisp_f_vector <- c()
 betadisp_p_vector <- c()
 #betadisp_r_vector <- c()
 adonis_f_vector <- c()
 adonis_p_vector <- c()
 adonis_r2_vector <- c()
 adonis_rres_vector <- c()
 adonis_mct_vector <- list()


 prog <- txtProgressBar(min = 1, max = n_boot, style = 3)

 
 ## Loop with number of bootstraps to perform
 for (i in 1:n_boot) {
  setTxtProgressBar(prog, i)
  ## get list of 10 individuals per group
  boot_table <- meta %>% 
   group_by_(group) %>% 
   sample_n(10, replace = TRUE)
  ## extract individuals and group
  boot_inds <- suppressMessages(boot_table %>% dplyr::select(ind) %>% pull())
  boot_groups <- boot_table %>% dplyr::select(group) %>% pull()
  ## subset otu table to just individuals above
  sub_tab <- tab[boot_inds,]
  ## run PhILR transform on subset otu table
  philr_data <- suppressMessages(philr(sub_tab, 
     tr,
     part.weights = "enorm.x.gm.counts",
     ilr.weights = "uniform"))
  ## calculate euclidiean distances from transformed data
  philr_dist <- dist(philr_data, method = "euclidean")
  ## perform beta-dispersion test, for legacy
  betadisp_anova <- anova(betadisper(philr_dist, boot_groups, type = "centroid"))
  ## run PERMANOVA on this subset data
  result <- adonis(philr_dist ~ boot_groups)
  ## run pairwise PERMANOVA 
  result_bonferroni <- pairwise.adonis(philr_dist, boot_groups) %>% as_tibble()
  ## add F and P values to results vectors
  betadisp_f_vector <- append(betadisp_f_vector, betadisp_anova$`F value`[1])
  betadisp_p_vector <- append(betadisp_p_vector, betadisp_anova$`Pr(>F)`[1])
  adonis_f_vector <- append(adonis_f_vector, result$aov.tab$F.Model[1])
  adonis_p_vector <- append(adonis_p_vector, result$aov.tab$`Pr(>F)`[1])
  adonis_r2_vector <- append(adonis_r2_vector, result$aov.tab$R2[1])
  adonis_rres_vector <- append(adonis_rres_vector, result$aov.tab$R2[2])
  adonis_mct_vector[[i]] <- result_bonferroni

 }
 ## return vectors
 return(list("betadisp_f_values" = betadisp_f_vector, 
       "betadisp_p_values" = betadisp_p_vector,
       "adonis_f_vector" = adonis_f_vector,
       "adonis_p_vector" = adonis_p_vector, 
       "adonis_r2_vector" = adonis_r2_vector,
       "adonis_rres_vector" = adonis_rres_vector,
       "adonis_mct_vector" = adonis_mct_vector))
}

bootstraps <- 100

balanced_permanova_results <- boot_permanova(otu.table, 
                       metadata_subset %>% 
                        filter(Host_Genus != "Alouatta"), 
                       "Host_Genus", 
                       "Individual", 
                       tree, 
                       bootstraps)

## to save
balanced_permanova_results_tib <- balanced_permanova_results
balanced_permanova_results_tib$adonis_mct_vector <- NULL

balanced_permanova_results_tib <- balanced_permanova_results_tib %>%
  as_tibble() %>% 
  mutate(run = row_number()) %>% 
  select(run, everything()) %>%
  rename(Run = run,
     BetaDispersion_F_Value = betadisp_f_values,
     BetaDispersion_p_Value = betadisp_p_values,
     Adonis_pseudoF_Value = adonis_f_vector,
     Adonis_p_Value = adonis_p_vector,
     Adonis_RsquaredGroups_Value = adonis_r2_vector,
     Adonis_RsquaredResiduals_Value = adonis_rres_vector)

balanced_permanova_results_mct_tib <- balanced_permanova_results$adonis_mct_vector %>% 
  enframe() %>% 
  unnest %>%
  select(-sig, -name) %>%
  mutate(pairs = as.character(pairs)) %>%
  group_by(pairs) %>%
  summarise(mean_SumsOfSqs = mean(SumsOfSqs),
            sd_SumsOfSqs = sd(SumsOfSqs),
            mean_F.Model = mean(F.Model),
            sd_F.Model = sd(F.Model),
            mean_R2 = mean(R2),
            sd_R2 = sd(R2),
            mean_p.value = mean(p.value),
            sd_p.value = sd(p.value),
            mean_p.adjusted = mean(p.adjusted),
            sd_p.adjusted = sd(p.adjusted)
            ) %>%
  mutate(bootstraps = bootstraps)

write_tsv(balanced_permanova_results_tib, paste("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/11a-permanova_bootstrap_hostgenuscentroids_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"), ".tsv", sep = ""))

write_tsv(balanced_permanova_results_mct_tib, 
          paste("/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/11b-permanova_bootstrap_hostgenuscentroids_multiplecorrectiontesting_", db , "_", tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"), ".tsv", sep = ""))
 

#' 
#' ### R2 Value
#' 
#' R-squared is a statistical measure of how close the data are to the fitted 
#' regression line. It is also known as the coefficient of determination, or the 
#' coefficient of multiple determination for multiple regression.
#' 
#' The definition of R-squared is fairly straight-forward; it is the percentage 
#' of the response variable variation that is explained by a linear model. Or:
#' 
#' R-squared = Explained variation / Total variation
#' 
#' R-squared is always between 0 and 100%:
#' 
#' 0% indicates that the model explains none of the variability of the response 
#' data around its mean.
#' 100% indicates that the model explains all the variability of the response data 
#' around its mean.
#' 
#' From http://blog.minitab.com/blog/adventures-in-statistics-2/regression-analysis-how-do-i-interpret-r-squared-and-assess-the-goodness-of-fit
#' 
#' R2 "[...] is usually interpreted as the proportion of variance in the dependent
#' variable "explained" by the independent variables in the regression equation"
#' From: http://www.polisci.msu.edu/jacoby/icpsr/reg2015/handouts/multreg1/R-sqd%20and%20ANOVA,%202015%20ICPSR%20Handout.pdf
#' 
#' To plot a summary of the values generated after running the 'bootstrapped' 
#' PERMANOVA
#' 
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
adonis_r2_plot <- ggplot(as_tibble(balanced_permanova_results$adonis_r2_vector), aes(value)) +
 geom_histogram(alpha = 0.5) +
 theme_minimal(base_size = 7, base_family = "roboto") +
 ggtitle("Bootstrapped PERMANOVA R squared value/ residuals summary")

adonis_rres_plot <- ggplot(as_tibble(balanced_permanova_results$adonis_rres_vector), aes(value)) +
 geom_histogram(alpha = 0.5) +
 theme_minimal(base_size = 7, base_family = "roboto")

if (script == F) {
 adonis_r2_plot
}

ggsave(paste("12-philr_bootstrappermanova_malt_r2distribution_hostgenus_", db , "_",tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = adonis_r2_plot, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 3.5, height = 3.5, units = "in", dpi = 600)



#' 
#' 
#' ### F Value
#' 
#' "The PERMANOVA pseudo F statistic for testing the number of differences in the 
#' positions of the group centroids in the space of the chosen dissimilarity 
#' measure is given by"
#' 
#' F = (SumSquaresAmonggroups/SumSequencesWithinGroup)·[(Nrows−gowersMeanMatrix)/(gowersMeanMatrix−1)]
#' 
#' However, for a more simple explanation for stupid non-sciencey people like me:
#' 
#' F = variation between sample means / variation within the samples
#' 
#' Note: In general, an F-statistic is a ratio of two quantities that are expected 
#' to be roughly equal under the null hypothesis, which produces an F-statistic of 
#' approximately 1. In order to reject the null hypothesis that the group means 
#' are equal, we need a high F-value. 
#' 
#' From: http://blog.minitab.com/blog/adventures-in-statistics-2/understanding-analysis-of-variance-anova-and-the-f-test
#' 
#' To plot a summary of the values generated after running the 'bootstrapped' 
#' PERMANOVA
#' 
## ----fig.height=3.5, fig.width=7-----------------------------------------
f_mean <- mean(balanced_permanova_results$adonis_f_vector)
f_median <- median(balanced_permanova_results$adonis_f_vector)
p_mean <- mean(balanced_permanova_results$adonis_f_vector)

adonis_f_plot <- ggplot(as_tibble(balanced_permanova_results$adonis_f_vector), aes(value)) +
 geom_histogram(alpha = 0.5) +
 geom_vline(xintercept = original_f, colour = "red", linetype = "dashed", size = 1, alpha = 0.5) +
 geom_vline(xintercept = f_mean, colour = "blue", linetype = "dashed", size = 1, alpha = 0.5) +
 annotate("text", x = original_f, y = 11, label = "Original", family = "Roboto", size = 3, colour = "red", vjust = 0.5) +
 annotate("text", x = f_mean, y = 11, label = "Mean", family = "Roboto", size = 3, colour = "blue") +
 xlab(expression(paste("Pseudo ", italic(F), " value"))) +
 ylab("Count") +
 theme_minimal(base_size = 7, base_family = "roboto") +
 theme()

if (script == F) {
 adonis_f_plot
}

ggsave(paste("13-philr_bootstrappermanova_malt_pseudofdistribution_hostgenus_", db , "_",tax_level, "_", sources, "_", controls, "_badsamples", bad_samples, "_", sample_filter, "_", minsupp_threshold, "_" , format(Sys.Date(), "%Y%m%d"),".pdf", sep = ""), plot = adonis_f_plot, "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/screening/philr.backup/", device = cairo_pdf, width = 3.5, height = 3.5, units = "in", dpi = 600)


#' 
#' 
#' ## Multivariate Analysis Conclusion
#' 
#' While initially we were possibly confounding our analysis due to the combination 
#' of variation in the beta-dispersion of the samples in each group, and also an
#' imbalanced design (i.e. sample size) PERMANOVA analysis, bootstrapping on
#' random combination of 10 individuals per group (note: only for Gorillas, Chimps 
#' and Humans) to make a balanced design (which PERMANOVA is robust against),
#' we still had significant (p = 0.001) pseudo-F values indicating that the 
#' differences in the centroids of the data occuring despite the samples deriving
#' from the same population would only occur <5% of the time. So, we have a low
#' probability that the null hypothesis of no differences in the centroids 
#' of the groups is true.
#' 
#' However, it should be noted the R2 (or the variation explained by the variable),
#' is still relatively low. So there is probably better partitioning of the data
#' that could descibe this.
#' 
#' ## Script generation
#' 
#' Use the following to generate a new fast script.
#' 

#' 
