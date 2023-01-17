args <- commandArgs(trailingOnly = TRUE)

fname <- args[1]
d <- read.delim(fname, header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)
#d <- read.delim("/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/HMO_M1.INT.txt", header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)
covar <- read.table("/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/qc/mothers.PC1-2.txt", header = F,  as.is = T, check.names = F)
row.names(covar) <- covar[,2]

ids <- intersect(row.names(d), row.names(covar))

m <- cbind(d[ids,], covar[ids,c(3,4)])
colnames(m) <- c(colnames(d), "PC1", "PC2")

d_adj <- data.frame(matrix(nrow = length(ids), ncol = ncol(d)))
row.names(d_adj) <- ids
colnames(d_adj) <- colnames(d)

for (hmo in colnames(d)){
  subs <- m[,c(hmo, "PC1", "PC2")]
  colnames(subs)[1] <- "HMO"
  
  d_adj[,hmo] <- residuals(lm(HMO ~ PC1 + PC2, data = subs))
}
d_adj_t <- as.data.frame(t(d_adj))
write.table(d_adj_t, paste0(fname, ".adjPCs.t.txt"), sep = "\t", quote = F, col.names = NA)




