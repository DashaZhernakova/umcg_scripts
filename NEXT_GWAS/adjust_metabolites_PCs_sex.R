
fname <- "../metabolites/metabo_Mother_P28.filtered.txt"
d <- as.data.frame(t(read.delim(fname, header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)))
covar <- read.table("../metabolites/mothers.PC1-2.txt", header = F,  as.is = T, check.names = F)
row.names(covar) <- covar[,2]

ids <- intersect(row.names(d), row.names(covar))
length(ids)

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

# Babies



fname <- "../metabolites/metabo_Baby_B.filtered.txt"
d <- as.data.frame(t(read.delim(fname, header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)))
covar <- read.table("../metabolites/babies_PC+sex.txt", header = T,  as.is = T, row.names = 1, check.names = F)

ids <- intersect(row.names(d), row.names(covar))
length(ids)

m <- cbind(d[ids,], covar[ids,])
colnames(m) <- c(colnames(d), "PC1", "PC2", "sex")

d_adj <- data.frame(matrix(nrow = length(ids), ncol = ncol(d)))
row.names(d_adj) <- ids
colnames(d_adj) <- colnames(d)

for (hmo in colnames(d)){
  subs <- m[,c(hmo, "PC1", "PC2", "sex")]
  colnames(subs)[1] <- "HMO"
  subs$sex <- as.factor(subs$sex)
  
  d_adj[,hmo] <- residuals(lm(HMO ~ PC1 + PC2 + sex, data = subs))
}
d_adj_t <- as.data.frame(t(d_adj))
write.table(d_adj_t, paste0(fname, ".adjPCs.t.txt"), sep = "\t", quote = F, col.names = NA)



