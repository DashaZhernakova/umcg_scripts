setwd("/Users/Dasha/work/UMCG/data/NEXT/metabolites/")

d <- read.delim("20230112_residuals_mtb_time_corrected_for_genetics.txt", header = T, sep = "\t", as.is = T)
d$NEXT_ID_geno <- gsub("_2", "", d$NEXT_ID)

samples <- read.delim("genotype_ids.txt", header = F, as.is = T)

d <- d[d$NEXT_ID_geno %in% samples[,1],]

save_dedup <- function(d, who, when){
  d_timepoint <- d[d$Plasma_Sample_type == who & d$Plasma_Timepoint == when,]
  
  counts <- table(d_timepoint$NEXT_ID_geno)
  duplicates <- names(counts[counts > 1])
  
  rows_to_remove <- which(d_timepoint$NEXT_ID_geno %in% duplicates, )[-1]
  cat("Removing ", length(rows_to_remove), "rows\n")
  dedup <- d_timepoint[-rows_to_remove,]
  
  dedup$NEXT_ID <- dedup$NEXT_ID_geno
  dedup <- subset(dedup, select = -c(PARENT_SAMPLE_NAME,NEXT_Mother_ID, Family_ID, NEXT_ID_geno ,Plasma_Sample_type, Plasma_Timepoint))
  write.table(t(dedup), file = paste0("metabo_", who, "_", when, ".filtered.txt"), sep = "\t", quote = F, col.names = F)
  
}

table(d[,c("Plasma_Sample_type", "Plasma_Timepoint")])
save_dedup(d, "Mother", "B")
save_dedup(d, "Mother", "P12")
save_dedup(d, "Mother", "P28")
save_dedup(d, "Baby", "B")

