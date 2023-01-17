args <- commandArgs(trailingOnly = TRUE)

cat("Current directory: ", args[1], "\n")
setwd(args[1])

# Mothers
eigenvec <- read.table("all_chr.mothers.no_rel.pruned_r2_0.2.eigenvec",  header = F)
colnames(eigenvec) = c(c("FID","IID"), paste0("PC", seq(5)))
eigenvec <- eigenvec[,seq(1,7)]
pdf("mothers_genotype_PCA.pdf", width = 10, height = 10, useDingbats = F)
pairs(eigenvec[,3:ncol(eigenvec)], pch = 16, col = "darkblue")
dev.off()

plot(eigenvec$PC1, eigenvec$PC3, pch = 16, col = "darkblue")
pc1_outliers <- get_outlier_threshold(eigenvec$PC1, cutoff = 10)
pc2_outliers <- get_outlier_threshold(eigenvec$PC2, cutoff = 10)
pc3_outliers <- get_outlier_threshold(eigenvec$PC3, cutoff = 10)
abline(v = pc1_outliers[[1]], col = "red")
abline(v = pc1_outliers[[2]], col = "red")

abline(h = pc3_outliers[[1]], col = "red")
abline(h = pc3_outliers[[2]], col = "red")

outliers <- eigenvec[eigenvec$PC3 > 0.4 | eigenvec$PC1 > 0.4,]
outliers <- eigenvec[eigenvec$PC3 > 0.4 ,]
outliers 




# Babies
eigenvec <- read.table("all_chr.babies.no_rel.pruned_r2_0.2.eigenvec",  header = F)
colnames(eigenvec) = c(c("FID","IID"), paste0("PC", seq(5)))
eigenvec <- eigenvec[,seq(1,7)]
pdf("babies_genotype_PCA.pdf", width = 10, height = 10, useDingbats = F)
pairs(eigenvec[,3:ncol(eigenvec)], pch = 16, col = "darkblue")
dev.off()

plot(eigenvec$PC1, eigenvec$PC3, pch = 16, col = "darkblue")
pc1_outliers <- get_outlier_threshold(eigenvec$PC1, cutoff = 10)
pc2_outliers <- get_outlier_threshold(eigenvec$PC2, cutoff = 10)
pc3_outliers <- get_outlier_threshold(eigenvec$PC3, cutoff = 10)
abline(v = pc1_outliers[[1]], col = "red")
abline(v = pc1_outliers[[2]], col = "red")

abline(h = pc3_outliers[[1]], col = "red")
abline(h = pc3_outliers[[2]], col = "red")

outliers <- eigenvec[eigenvec$PC3 > 0.4 | eigenvec$PC1 > 0.4,]
outliers <- eigenvec[eigenvec$PC3 > 0.4 ,]
outliers 



get_outlier_threshold <- function(d, cutoff = 3){
  q1 <- quantile(d, probs = 0.25)
  q3 <- quantile(d, probs = 0.75)
  iqr <- q3 - q1
  return (list(q1 - cutoff*iqr, q3 + cutoff*iqr))
}
