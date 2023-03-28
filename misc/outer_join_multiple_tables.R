args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(tidyverse)

outpath <- args[1]
cat(length(args), "\n", args, "\n")

table_list = list()
for (fname in args[2:length(args)]){
  f <- read.delim(fname,  header = T, as.is = T, check.names = F, sep = "\t")
  cat("Read", fname, "\n")
  colnames(f)[1] <- "id"
  table_list[[fname]] <- f
}
cat("Read", length(table_list), "files\n")
merged_table <- table_list %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="id"), .)

write.table(merged_table, file = outpath, sep = "\t", quote = F, row.names = F) 
cat("Finished!\nDimensions of the merged table:", dim(merged_table), "\n")
