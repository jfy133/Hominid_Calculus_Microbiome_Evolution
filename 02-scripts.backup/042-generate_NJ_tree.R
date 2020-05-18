#! /usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

if (args[1] %in% c("-h", "--help", "-?", "-help")) {
  cat("
Usage: Rscript generate_NJ_tree.R <in_fasta> <minimum_no_snps> <substitution_model> <number_of_bootstraps> <samples_to_exclude>

      Generates a bootstrapped neighbour-joining tree, with prior sample 
      removal based on minimum number of SNP threshold. 
      
      Input file is a fasta or fasta.gz alignment. Valid substitution models 
      are those valid with ape's dist.dna() function. Samples to exclude should 
      be either 'none' to keep all, or sample names in form of comma separated
      list e.g. Sample1,Sample2,Sample3
      
      SNPs are point mutations only. Output is a newick file in the same folder 
      as your input file.

      Required R packages: ape, ade4, adegenet
      
      Note: if a .fasta.gz file is supplied, please ignore 'wrong file 
      extension' errors.
")
  quit()
  
} else if (length(args) < 4 | length(args) > 6) {
  cat("
Usage: Rscript generate_NJ_tree.R <in_fasta> <minimum_no_snps> <substitution_model> <number_of_bootstraps> <samples_to_exclude>

      Generates a bootstrapped neighbour-joining tree, with prior sample 
      removal based on minimum number of SNP threshold. 
      
      Input file is a fasta or fasta.gz alignment. Valid substitution models 
      are those valid with ape's dist.dna() function. Samples to exclude should 
      be either 'none' to keep all, or sample names in form of comma separated
      list e.g. Sample1,Sample2,Sample3
      
      SNPs are point mutations only. Output is a newick file in the same folder 
      as your input file.

      Required R packages: ape, ade4, adegenet
      
      Note: if a .fasta.gz file is supplied, please ignore 'wrong file 
      extension' errors.

")
  stop("Too few or too many input files given")
  
}

## Set variables

in_file <- as.character(args[1])
min_snps <- as.numeric(args[2])
model = as.character(args[3])
n_boots <- as.numeric(args[4])
exclude <- as.character(args[5])

## Testing variables
 # in_file <- "/home/fellows/projects1/microbiome_calculus/evolution/04-analysis/deep/multivcfanalyzer/initial_single_genome/output/initial_single_genome_2X_0.7_0.7/Porphyromonas_gingivalis_ATCC_33277/snpAlignment.fasta.gz"
 # min_snps <- 1000
 # model <- "JC69"
 # n_boots <- 100
 # exclude <- "none"

## Load libraries
cat("\nLOADING LIBRARIES\n")
library(ape)
library(ade4)
library(adegenet)

## Set functions
cat("SETTING FUNCTIONS\n")

count_nonN_nucleotides <- function(dnabin_obj){
  out_list <- list()
  
  for (i in 1:nrow(dnabin_obj)) {
    total <- sum(base.freq(dnabin_obj[i,], freq = T))
    out_list[row.names(dnabin_obj[i,])] <- total
    
  }
  return(out_list)
}

## Load data
cat("LOADING DATA (IF GZIP IGNORE WRONG EXTENSION)\n")
raw_fasta <- fasta2DNAbin(file = in_file, quiet = T)

cat("CALCULATING SAMPLE FILTERING STATS\n")

## Filter alignment to remove samples requested to be excluded
if (exclude != "none") { 
  vec_exclude <- unlist(strsplit(exclude, ",")) 
  subset_fasta <- raw_fasta[!rownames(raw_fasta) %in% vec_exclude,]
  samplesexcluded <- "samplesexcludedT" 
} else {
  subset_fasta <- raw_fasta
  samplesexcluded <- "samplesexcludedF" 
  }
  

## Calculate number of non-N nucleotides per sample
alignment_stats <- count_nonN_nucleotides(subset_fasta)

## Filter alignments to remove samples with less than min_snps
subset_fasta <- subset_fasta[alignment_stats >= min_snps,]

## Report samples kept/lost
sample_report <- data.frame(unlist(alignment_stats), 
                            names(unlist(alignment_stats)))

names(sample_report) <- c("Number_of_Positions", "Sample_Name")
row.names(sample_report) <- NULL
sample_report <- sample_report[order(-sample_report$Number_of_Positions),]

sample_report$Passed_Filter <- ifelse(sample_report$Number_of_Positions >= min_snps, 
                                              T, 
                                              F)

write.csv(sample_report,
          file = paste0(in_file,
                        "_minSamplePositions",
                        min_snps,
                        "_methodNJ_model",
                        model,
                        "_bootstraps",
                        n_boots, "_",
                        samplesexcluded,
                        "_pairwiseDel_sampleFilteringReport",
                        format(Sys.Date(), "%Y%m%d"),
                        ".csv")
)

## Report base frequencies for all (this allows removal of samples that violate
## model assumptions)
cat("GENERATING PER-SAMPLE BASE FREQUENCY TABLE\n")
basefreq_table <- list()

for (i in 1:nrow(subset_fasta)) {
  basefreq <- base.freq(subset_fasta[i,])
  basefreq_table[[rownames(subset_fasta[i,])]] <- basefreq  
}

basefreq_report <- t(as.data.frame(basefreq_table))

write.csv(basefreq_report,
          file = paste0(in_file,
                        "_minSamplePositions",
                        min_snps,
                        "_methodNJ_model",
                        model,
                        "_bootstraps",
                        n_boots, "_",
                        samplesexcluded,
                        "_pairwiseDel_baseFreqReport",
                        format(Sys.Date(), "%Y%m%d"),
                        ".csv")
)

## Report number of pairwise overlapping bases (i.e. how many nucleotides
## shared between two samples, of these which are same and different bases)
cat("CALCULATING SHARED NON-N POSITION STATISTICS\n")
comb_pairs <- combn(row.names(subset_fasta), 2, function(x) paste(x))

overlap_table <- data.frame()

## For each combination, run Ftab for non-N mutation count matrix, convert to
## long, summarise to total of shared and different bases, total all.  
for (i in 1:ncol(comb_pairs)) {
  overlap <- as.data.frame(Ftab(subset_fasta[comb_pairs[,i],]))
  sample_x <- rownames(subset_fasta[comb_pairs[,i],])[1]
  sample_y <- rownames(subset_fasta[comb_pairs[,i],])[2]
  overlap$sample_2 <- rownames(overlap)
  rownames(overlap) <- NULL
  overlap <- reshape(overlap, 
                     direction = "long", 
                     varying = list(names(overlap)[1:4]), 
                     v.names = "count", 
                     idvar = "sample_2", 
                     timevar = "sample_1", 
                     times = c("a", "c", "g", "t"))
  rownames(overlap) <- NULL
  overlap$snp <- overlap$sample_2 == overlap$sample_1
  result <- aggregate(overlap$count, by = list(snp = overlap$snp), FUN = sum)
  result$combination <- paste0(sample_x, "_", sample_y)
  overlap_table <- rbind(overlap_table, as.data.frame(result))
}

overlap_table <- reshape(overlap_table, 
                         idvar = "combination", 
                         timevar = "snp", 
                         direction = "wide")

colnames(overlap_table) <- c("Combination", 
                             "total_different_bases", 
                             "total_shared_bases")

overlap_table$total_overlapping_bases <- overlap_table$total_different_bases + overlap_table$total_shared_bases

overlap_table <- overlap_table[order(-overlap_table$total_overlapping_bases),]

write.csv(overlap_table,
          file = paste0(in_file,
                        "_minSamplePositions",
                        min_snps,
                        "_methodNJ_model",
                        model,
                        "_bootstraps",
                        n_boots, "_",
                        samplesexcluded,
                        "_pairwiseDel_overlappingNucleotidesReport",
                        format(Sys.Date(), "%Y%m%d"),
                        ".csv")
)


cat("BUILDING NJ DISTANCE TREE\n")

#subset_fasta <- subset_fasta[!grepl("OME003", rownames(subset_fasta)),]

## Build NJ distance tree based on Jukes-Cantor and boostrap
D <- dist.dna(subset_fasta, model = "JC69", pairwise.deletion = T)


tre <- nj(D)

cat("BOOTSTRAPPING NJ DISTANCE TREE\n")

boots <- boot.phylo(tre, 
                    subset_fasta, 
                    function(e) nj(dist.dna(e, 
                                            model = "JC69", 
                                            pairwise.deletion = T)),
                    B = n_boots)

## Add bootstraps to tree
cat("TREE ANNOTATION AND CLEANUP\n")

tre$node.label <- boots

## Fix negative branch lengths
## Reason why negative: https://www.sequentix.de/gelquest/help/neighbor_joining_method.htm
## Fix: http://boopsboops.blogspot.com/2010/10/negative-branch-lengths-in-neighbour.html
tre$edge.length[tre$edge.length < 0] <- 0

# Save tree in newick format
write.tree(tre,
           file = paste0(in_file,
                         "_minSamplePositions",
                         min_snps,
                         "_methodNJ_model",
                         model,
                         "_bootstraps",
                         n_boots, "_",
                         samplesexcluded,
                         "_pairwiseDel_",
                         format(Sys.Date(), "%Y%m%d"),
                         ".nwk")
           )

cat("FINISHED. CHECK FOR EXPECTED OUTPUT!")

