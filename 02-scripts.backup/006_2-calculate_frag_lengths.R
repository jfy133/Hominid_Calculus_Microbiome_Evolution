#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

input_file = paste(sub("/$", "", args[1]),"/damageprofiler/lgdistribution.txt", sep="")

raw_data <- read.table(input_file, header=TRUE, sep="\t")

median_calculator <- function(x, y){
  long_data <- rep(x, y)
  my_median <- median(long_data)
  return(my_median)
}

print(median_calculator(raw_data$Length,raw_data$Occurrences))