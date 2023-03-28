args <- commandArgs(trailingOnly = TRUE)

fname <- args[1]
varexp_cutoff <- 0.9

d <- read.delim(fname, check.names = F, header = T, row.names = 1)

pca <- prcomp(d, center = TRUE,scale = TRUE)
exp_var <- summary(pca)$importance[3,]
write.table(pca$x[,which(exp_var < varexp_cutoff)], file = paste0(fname, ".PCs_", varexp_cutoff, ".txt"), sep = "\t", quote = F, col.names = NA)


