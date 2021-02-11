library(feather)
args <- commandArgs(trailingOnly = TRUE)
setwd('/groups/umcg-lld/scr01/dasha/methylation/')

traits_path <- "data/ILLUMINA450K_All2Blood_mValue.Selection.Dasen.QuantileNormalized.LLDsubset.LLD_ids.scaled.feather"
gene_table_path <- "data/meth_probe2gene.txt"

pheno_path <- args[1]
out_path <- args[2]
permute <- args[3]

print(paste("methylation path:", traits_path))
print(paste("phenotype path:", pheno_path))
print(paste("output path:", out_path))
print(paste("permute labels?:", permute))
#################

traits <- as.data.frame(read_feather(traits_path))
row.names(traits) <- traits$sample
traits$sample <- NULL

pheno <- read.delim(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))
nrow(traits_m)

gene_table <- read.delim(gene_table_path, header = T, sep = "\t", as.is = T, check.names = F)

phenos <- colnames(traits_m)
phenos2 <- colnames(pheno_m)[10:ncol(pheno_m)]
res_table <- matrix(NA, nrow = length(phenos) * length(phenos2), ncol = 8)
colnames(res_table) <- c("pheno1", "gene", "pheno2", "beta", "se", "tval", "pval", "N")
print("Using the following phenotypes:")
print(phenos2)

pheno_m$gender <- as.factor(pheno_m$gender)
if ("SMK13" %in% colnames(pheno_m)) pheno_m$SMK13 <- as.factor(pheno_m$SMK13)
if ("SMK14" %in% colnames(pheno_m)) pheno_m$SMK14 <- as.factor(pheno_m$SMK14)
if ("parents_diab" %in% colnames(pheno_m)) pheno_m$parents_diab <- as.factor(pheno_m$parents_diab)
if ("father_diab" %in% colnames(pheno_m)) pheno_m$father_diab <- as.factor(pheno_m$father_diab)
if ("mother_diab" %in% colnames(pheno_m)) pheno_m$mother_diab <- as.factor(pheno_m$mother_diab)

if (permute){
  print("NB! Running a permutation!!!")
  pheno_m <- pheno_m[sample(nrow(pheno_m)),]
}
##########

cnt <- 1
n <- length(phenos) * length(phenos2)
for (g in seq_along(colnames(traits_m))) {
  probe <- traits_m[,g]
  g_id <- colnames(traits_m)[g]
  g_name <- gene_table[gene_table[, 1] == g_id, 2]
  if (length(g_name) == 0) g_name <- g_id

  if (cnt %% 10000 == 0){
          print(paste0("Processed ", cnt, " pairs (", format(100 * cnt / n, digits = 4), "%)"))
  }
  
  for (p_name in phenos2) {
        merged_tab <- cbind(probe, pheno_m[, c(p_name, "age", "gender", "ba", "eo", "er", "gr", "ly", "mo", "tr")])
        colnames(merged_tab) <- c("pheno1", "pheno2", "age", "gender", "ba", "eo", "er", "gr", "ly", "mo", "tr")
        merged_tab$pheno1 <- as.numeric(as.character(merged_tab$pheno1))
        merged_tab <- na.omit(merged_tab)
        lm_res <- summary(lm(pheno1 ~ age + gender + pheno2 + ba + eo + er + gr + ly + mo + tr, data = merged_tab))$coefficients	
        res_table[cnt, ] <- c(g_id, g_name, p_name, as.numeric(lm_res[4, ]), nrow(merged_tab))
        cnt <- cnt + 1
  }
}


write.table(res_table, file = out_path, sep = "\t", quote = F, col.names = NA)

