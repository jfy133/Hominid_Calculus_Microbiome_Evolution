library(tidyverse)
library(readxl)

## Load two tables, one from Slon 2017 SI and other from ENA project
ena_data <- read_tsv("../00-documentation.backup/99-Slon2017_PRJEB18629.txt")
paper_data <- read_xlsx("../00-documentation.backup/99-aam9695_DataFileS1.xlsx")

## Combine two tables and find files with >= 10 million reads
combined_data <- left_join(paper_data, ena_data, by = c("Indexed library ID" = "library_name"))
filtered_data <- filter(combined_data, read_count >= 10000000)

## Save file
write.table(filtered_data$fastq_ftp, 
            file = "../00-documentation.backup/99-Slon2017_AccessionsToDownload_2.tsv", 
            row.names = FALSE, 
            quote = FALSE, 
            col.names = FALSE)
