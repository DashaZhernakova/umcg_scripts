setwd("/Users/Dasha/work/UMCG/data/NEXT/metabolites/")

d <- read.delim("20230112_residuals_mtb_time_corrected_for_genetics.txt", header = T, sep = "\t", as.is = T, check.names = FALSE)
d$NEXT_ID_geno <- gsub("_2", "", d$NEXT_ID)

samples <- read.delim("genotype_ids.txt", header = F, as.is = T)

d <- d[d$NEXT_ID_geno %in% samples[,1],]

save_dedup <- function(d, who, when){
  d_timepoint <- d[d$Plasma_Sample_type == who & d$Plasma_Timepoint == when,]
  
  counts <- table(d_timepoint$NEXT_ID_geno)
  duplicates <- names(counts[counts > 1])
  
  rows_to_remove <- which(d_timepoint$NEXT_ID_geno %in% duplicates, )[-1]
  cat("Removing ", length(rows_to_remove), "rows\n")
  if (length(rows_to_remove) > 0) {
    dedup <- d_timepoint[-rows_to_remove,]
  } else {
    dedup <- d_timepoint
  }
  dedup$NEXT_ID <- dedup$NEXT_ID_geno
  dedup <- subset(dedup, select = -c(PARENT_SAMPLE_NAME,NEXT_Mother_ID, Family_ID, NEXT_ID_geno ,Plasma_Sample_type, Plasma_Timepoint))
  write.table(t(dedup), file = paste0("metabo_", who, "_", when, ".filtered.txt"), sep = "\t", quote = F, col.names = F)
  
  row.names(dedup) <- dedup$NEXT_ID
  dedup <- subset(dedup, select=-NEXT_ID)
  return (dedup)
}

table(d[,c("Plasma_Sample_type", "Plasma_Timepoint")])
dedup = NULL
dedup <- save_dedup(d, "Mother", "B")
adjust_for_PCs_sex(dedup, "Mother", "B")

dedup = NULL
dedup <- save_dedup(d, "Mother", "P12")
adjust_for_PCs_sex(dedup, "Mother", "P12")

dedup = NULL
dedup <- save_dedup(d, "Mother", "P28")
adjust_for_PCs_sex(dedup,"Mother", "P28")

dedup = NULL
dedup <- save_dedup(d, "Baby", "B")
adjust_for_PCs_sex(dedup, "Baby", "B")



adjust_for_PCs_sex <- function(d, who, when){
  if (who == "Mother"){
    covar <- read.table("../metabolites/mothers.PC1-2.txt", header = F,  as.is = T, check.names = F)
    row.names(covar) <- covar[,2]
    
    ids <- intersect(row.names(d), row.names(covar))
    length(ids)
    
    m <- cbind(d[ids,], covar[ids,c(3,4)])
    predictors <- c("PC1", "PC2")
    colnames(m) <- c(colnames(d), predictors)
    
    lm_formula <- "metabo ~ PC1 + PC2"
  } else if (who == "Baby"){
    covar <- read.table("../metabolites/babies_PC+sex.txt", header = T,  as.is = T, row.names = 1, check.names = F)
    
    ids <- intersect(row.names(d), row.names(covar))
    length(ids)
    
    m <- cbind(d[ids,], covar[ids,])
    predictors <- c("PC1", "PC2", "sex")
    colnames(m) <- c(colnames(d), predictors)
    
    lm_formula <- "metabo ~ PC1 + PC2 + sex"
  }
  
  d_adj <- data.frame(matrix(nrow = length(ids), ncol = ncol(d)))
  row.names(d_adj) <- ids
  colnames(d_adj) <- colnames(d)
  
  for (metabo in colnames(d)){
    subs <- m[,c(metabo, predictors)]
    colnames(subs)[1] <- "metabo"
    
    d_adj[,metabo] <- residuals(lm(as.formula(lm_formula), data = subs))
  }
  d_adj_t <- as.data.frame(t(d_adj))
  write.table(d_adj_t, file = paste0("metabo_", who, "_", when, ".filtered.adjPCs.t.txt"), sep = "\t", quote = F, col.names = NA)
  
}
